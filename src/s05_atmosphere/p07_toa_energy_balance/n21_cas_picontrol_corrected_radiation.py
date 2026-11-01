# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import netCDF4 as nc
import numpy as np
import pandas as pd


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
forced_source = root / 'output/s05_atmosphere/p07_toa_energy_balance/n20_cas_toa_radiation'
destination = root / 'output/s05_atmosphere/p07_toa_energy_balance/n21_cas_picontrol_corrected_radiation'
assert not destination.exists(), '正式输出存在，须先归档。'
model, member = 'CAS-ESM2-0', 'r1i1p1f1'
variables = ['rsdt', 'rsut', 'rlut']
quantities = variables + ['asr', 'negative_olr', 'eei']
chunk_months = 24
# 主算按升高分支parent/child均为0的相对计数解释；日历标签冲突另作敏感性。
fit_first, fit_last = 1, 280
fit_cases = [('relative_start_1_280', 1, 280, 0),
             ('calendar_match_200_479', 200, 479, 199),
             ('all550_relative_start', 1, 550, 0),
             ('all550_calendar_match', 1, 550, 199)]
monthly_parts, inventory, coefficient_rows, sensitivity_rows, case_parts = [], [], [], [], []
with np.load(forced_source / 'native_weights.npz') as saved:
    latitude, longitude, weights = saved['latitude'], saved['longitude'], saved['weight']
forced = pd.read_csv(forced_source / 'annual_global_energy.csv', float_precision='round_trip')
forced_metadata = json.loads((forced_source / 'metadata.json').read_text(encoding='utf-8'))
for entry in forced_metadata['source_files']:
    if entry['experiment'] == '1pctCO2':
        assert entry['branch_time_in_parent'] == entry['branch_time_in_child'] == 0


# %% 3. 读取550年原生全球控制通量

for variable in variables:
    paths = sorted(Path('F:/sea/cmip6' + variable + '-picontrol').rglob(
        f'{variable}_Amon_{model}_piControl_{member}_*.nc'))
    assert len(paths) == 1
    for path in paths:
        with nc.Dataset(path) as dataset:
            assert dataset.source_id == model and dataset.experiment_id == 'piControl'
            assert dataset.variant_label == member and dataset.variable_id == variable
            field, time = dataset[variable], dataset['time']
            assert field.dimensions == ('time', 'lat', 'lon') and field.units == 'W m-2'
            np.testing.assert_array_equal(dataset['lat'][:], latitude)
            np.testing.assert_array_equal(dataset['lon'][:], longitude)
            dates = nc.num2date(time[:], time.units, calendar=time.calendar)
            years = np.array([date.year for date in dates])
            months = np.array([date.month for date in dates])
            np.testing.assert_array_equal(years * 12 + months - 1, np.arange(12, 551 * 12))
            bounds = nc.num2date(dataset[time.bounds][:].reshape(-1), time.units,
                                calendar=time.calendar).reshape(-1, 2)
            days = np.array([(pair[1] - pair[0]).total_seconds() / 86400 for pair in bounds])
            assert np.all(days > 0)
            global_mean = np.empty(len(dates))
            for offset in range(0, len(dates), chunk_months):
                values = np.ma.asarray(field[offset:offset + chunk_months], dtype=np.float64)
                assert not np.ma.getmaskarray(values).any()
                values = np.asarray(values)
                assert np.isfinite(values).all()
                global_mean[offset:offset + len(values)] = np.einsum('tyx,yx->t', values, weights, optimize=True)
            monthly_parts.append(pd.DataFrame({'model': model, 'variable': variable,
                                               'source_year': years, 'month': months,
                                               'days': days, 'global_Wm2': global_mean}))
            inventory.append({'path': str(path), 'variable': variable, 'experiment': 'piControl',
                              'member': member, 'dimensions': list(field.dimensions),
                              'shape': list(field.shape), 'units': field.units,
                              'calendar': time.calendar, 'time_units': time.units,
                              'first_year': int(years[0]), 'last_year': int(years[-1]),
                              'months': len(dates), 'tracking_id': dataset.tracking_id,
                              'bytes': path.stat().st_size, 'mtime_ns': path.stat().st_mtime_ns})
        print(variable, len(dates), 'control months complete', flush=True)


# %% 4. 完整控制年平均与辐射记账

monthly = pd.concat(monthly_parts, ignore_index=True)
assert not monthly.duplicated(['variable', 'source_year', 'month']).any()
annual_parts = []
for variable in variables:
    selected = monthly.loc[monthly.variable.eq(variable)].copy()
    selected['weighted_flux'] = selected.global_Wm2 * selected.days
    annual = selected.groupby('source_year').agg(months=('month', 'count'), days=('days', 'sum'),
                                                weighted=('weighted_flux', 'sum')).reset_index()
    assert annual.months.eq(12).all()
    annual[variable] = annual.weighted / annual.days
    annual_parts.append(annual[['source_year', variable]])
control = annual_parts[0].merge(annual_parts[1], on='source_year').merge(annual_parts[2], on='source_year')
control['asr'] = control.rsdt - control.rsut
control['negative_olr'] = -control.rlut
control['eei'] = control.asr + control.negative_olr
np.testing.assert_array_equal(control.source_year, np.arange(1, 551))


# %% 5. 对应控制背景OLS拟合及窗口敏感性

for case, first, last, year_offset in fit_cases:
    selected = control.loc[control.source_year.between(first, last)]
    assert len(selected) == last - first + 1
    origin = (first + last) / 2
    design = np.column_stack([np.ones(len(selected)), selected.source_year.to_numpy() - origin])
    corrected = forced[['model', 'phase', 'source_year', 'simulation_year'] + quantities].copy()
    corrected['control_year'] = corrected.simulation_year + year_offset
    corrected['case'] = case
    for quantity in quantities:
        intercept, slope = np.linalg.lstsq(design, selected[quantity].to_numpy(), rcond=None)[0]
        corrected[quantity + '_control_fit'] = intercept + slope * (corrected.control_year - origin)
        corrected[quantity + '_relative_control'] = corrected[quantity] - corrected[quantity + '_control_fit']
        for phase in ['increase', 'removal']:
            mask = corrected.phase.eq(phase)
            corrected.loc[mask, quantity + '_relative_control_smooth11'] = corrected.loc[
                mask, quantity + '_relative_control'].rolling(11, center=True, min_periods=6).mean()
        coefficient_rows.append({'case': case, 'quantity': quantity, 'fit_first': first, 'fit_last': last,
                                 'origin': origin, 'intercept_Wm2': intercept,
                                 'slope_Wm2_per_year': slope, 'evaluation_year_offset': year_offset})
        for window, start, stop in [('initial', 1, 10), ('peak_year140', 140, 140), ('late', 271, 280)]:
            sensitivity_rows.append({'case': case, 'quantity': quantity, 'window': window,
                                     'mean_relative_control_Wm2': corrected.loc[
                                         corrected.simulation_year.between(start, stop), quantity + '_relative_control'].mean()})
    for suffix in ['_relative_control', '_relative_control_smooth11']:
        np.testing.assert_allclose(corrected['eei' + suffix],
                                   corrected['asr' + suffix] + corrected['negative_olr' + suffix],
                                   atol=1e-10, rtol=0)
    case_parts.append(corrected)
all_cases = pd.concat(case_parts, ignore_index=True)
primary = all_cases.loc[all_cases.case.eq('relative_start_1_280')].copy()
differences = []
for case, first, last, offset in fit_cases[1:]:
    alternate = all_cases.loc[all_cases.case.eq(case)].set_index('simulation_year')
    for quantity in ['eei', 'asr', 'negative_olr']:
        difference = alternate[quantity + '_relative_control'].to_numpy() - primary[quantity + '_relative_control'].to_numpy()
        differences.append({'case': case, 'quantity': quantity, 'max_abs_CAS_difference_Wm2': float(np.max(np.abs(difference))),
                            'max_abs_recovery_mean_difference_Wm2': float(np.max(np.abs(difference))) / 3})


# %% 6. 保存主算、敏感性及无法按字面解码的分支属性

destination.mkdir(parents=True)
monthly.to_csv(destination / 'monthly_control.csv', index=False)
control.to_csv(destination / 'annual_control.csv', index=False)
primary.to_csv(destination / 'corrected_annual.csv', index=False)
all_cases.to_csv(destination / 'alignment_sensitivity_annual.csv', index=False)
pd.DataFrame(coefficient_rows).to_csv(destination / 'control_fit_coefficients.csv', index=False)
pd.DataFrame(sensitivity_rows).to_csv(destination / 'alignment_sensitivity_endpoints.csv', index=False)
pd.DataFrame(differences).to_csv(destination / 'alignment_sensitivity_maximum.csv', index=False)
pd.DataFrame(inventory).to_csv(destination / 'source_inventory.csv', index=False)
metadata = {'status': 'passed_with_branch_origin_limitation', 'model': model, 'member': member,
            'control_years': [1, 550], 'primary_fit_years': [fit_first, fit_last],
            'primary_mapping': 'simulation_year maps to control source_year1-280; interpretation of zero-valued parent/child elapsed counters.',
            'alignment_limitation': 'Not an independently validated branch mapping. Forced parent_time_units declares1850 origin while actual control uses0001; CDR branch dates also conflict with actual coordinates. The relative-start mapping is provisional and calendar200-479 plus all550 alternatives are retained.',
            'fit': 'Equal-year OLS intercept+slope; subtract fitted control background, not random control fluctuations; no forced initial centering.',
            'raw_dimensions': 'time xlat128 xlon256; control550 complete years; forced280 complete years',
            'spatial_weights': 'Exactly reuse forced native latitude-bound spherical area weights after coordinate checks.',
            'annual_weights': 'CF month-bound durations; all global cells and12months required.',
            'smoothing': 'Saved centered11year separately within each phase, min_periods6.',
            'sensitivity_cases': [row[0] for row in fit_cases],
            'source_files': inventory,
            'small_inputs': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
                             for path in [forced_source / 'annual_global_energy.csv', forced_source / 'metadata.json', forced_source / 'native_weights.npz']]}
(destination / 'metadata.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
print(pd.DataFrame(differences).to_string(index=False), flush=True)
print(destination.resolve(), flush=True)
