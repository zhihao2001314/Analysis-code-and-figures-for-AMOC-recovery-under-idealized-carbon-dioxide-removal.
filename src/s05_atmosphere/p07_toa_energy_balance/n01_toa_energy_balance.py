# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import netCDF4 as nc
import numpy as np
import pandas as pd


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
destination = root / 'output/s05_atmosphere/p07_toa_energy_balance/n01_toa_energy_balance'
assert not destination.exists(), '正式产品已存在，须先归档。'
destination.mkdir(parents=True)
registry_path = root / 'output/s05_atmosphere/p05_amoc_recovery_composites/n49_group_naming/group_labels.json'
registry = json.loads(registry_path.read_text(encoding='utf-8'))
models = [model for members in registry['members'].values() for model in members]
starts = {'CESM2': [1, 1], 'CNRM-ESM2-1': [1850, 1990],
          'NorESM2-LM': [1, 141], 'GFDL-ESM4': [6, 141],
          'ACCESS-ESM1-5': [101, 241], 'CanESM5': [1850, 1991],
          'MIROC-ES2L': [1850, 1989], 'UKESM1-0-LL': [1850, 1990]}
members = {'CESM2': 'r1i1p1f1', 'CNRM-ESM2-1': 'r1i1p1f2',
           'NorESM2-LM': 'r1i1p1f1', 'GFDL-ESM4': 'r1i1p1f1',
           'ACCESS-ESM1-5': 'r1i1p1f1', 'CanESM5': 'r1i1p2f1',
           'MIROC-ES2L': 'r1i1p1f2', 'UKESM1-0-LL': 'r1i1p1f2'}
variables = ['rsdt', 'rsut', 'rlut']
phases = ['increase', 'removal']
experiments = ['1pctCO2', '1pctCO2-cdr']
inventory, monthly_parts, grid_report = [], [], []
reference_grids = {}


# %% 3. 按文件内部标识清点来源，不按所在模式目录判断

for variable in variables:
    for suffix in ['', 'ad']:
        for path in sorted(Path('F:/sea/cmip6' + variable + suffix).rglob('*.nc')):
            with nc.Dataset(path) as dataset:
                model = dataset.source_id
                experiment = dataset.experiment_id
                assert dataset.variable_id == variable
                assert model in models and experiment in experiments
                assert dataset.variant_label == members[model]
                field = dataset.variables[variable]
                time = dataset.variables['time']
                dates = nc.num2date(time[:], time.units, calendar=time.calendar)
                source_year = np.array([date.year for date in dates])
                source_month = np.array([date.month for date in dates])
                assert field.dimensions == ('time', 'lat', 'lon')
                assert field.units == 'W m-2'
                phase_index = experiments.index(experiment)
                first_year = starts[model][phase_index]
                # 多保留GFDL上升原始1–5年，只用于1–140与6–145窗口敏感性。
                read_first = 1 if model == 'GFDL-ESM4' and phase_index == 0 else first_year
                take = np.flatnonzero((source_year >= read_first) & (source_year <= first_year + 139))
                inventory.append({'model': model, 'experiment': experiment, 'variable': variable,
                                  'member': dataset.variant_label, 'grid': dataset.grid_label,
                                  'path': str(path), 'bytes': path.stat().st_size,
                                  'dimensions': str(field.dimensions), 'shape': str(field.shape),
                                  'calendar': time.calendar, 'units': field.units,
                                  'first_date': str(dates[0]), 'last_date': str(dates[-1]),
                                  'selected_months': len(take)})
                if not len(take):
                    continue
                latitude = np.asarray(dataset.variables['lat'][:], dtype=float)
                longitude = np.asarray(dataset.variables['lon'][:], dtype=float)
                np.testing.assert_allclose(np.diff(longitude), 360 / len(longitude), atol=1e-4, rtol=0)
                if model == 'CNRM-ESM2-1':
                    nodes, latitude_weight = np.polynomial.legendre.leggauss(len(latitude))
                    np.testing.assert_allclose(latitude, np.rad2deg(np.arcsin(nodes)), atol=1e-4, rtol=0)
                    area_method = 'Gaussian quadrature; verified Gaussian nodes and uniform longitudes'
                else:
                    latitude_bounds = np.asarray(dataset.variables[dataset.variables['lat'].bounds][:], dtype=float)
                    assert latitude_bounds.shape == (len(latitude), 2)
                    np.testing.assert_allclose(latitude_bounds[1:, 0], latitude_bounds[:-1, 1], atol=1e-6)
                    np.testing.assert_allclose(latitude_bounds[[0, -1], [0, 1]], [-90, 90], atol=1e-6)
                    latitude_weight = np.diff(np.sin(np.deg2rad(latitude_bounds)), axis=1)[:, 0]
                    longitude_bounds = np.asarray(dataset.variables[dataset.variables['lon'].bounds][:], dtype=float)
                    if model == 'NorESM2-LM' and phase_index == 0:
                        # 上升文件把周期首界裁到0、末界拉到360；中心仍为0..357.5的均匀周期网格。
                        expected_width = np.full(len(longitude), 2.5)
                        expected_width[[0, -1]] = [1.25, 3.75]
                        np.testing.assert_allclose(np.diff(longitude_bounds, axis=1)[:, 0], expected_width, atol=1e-4)
                    else:
                        np.testing.assert_allclose(np.diff(longitude_bounds, axis=1)[:, 0], 360 / len(longitude), atol=1e-4)
                    area_method = 'Native latitude-bound spherical areas; uniform cyclic longitude-center widths (Nor upward bounds have clipped seam)'
                assert np.all(latitude_weight > 0)
                np.testing.assert_allclose(latitude_weight.sum(), 2, atol=1e-6, rtol=0)
                weights = latitude_weight[:, None] * np.ones((1, len(longitude)))
                weights /= weights.sum()
                if model not in reference_grids:
                    reference_grids[model] = (latitude.copy(), longitude.copy(), weights.copy())
                    np.savez(destination / ('weights_' + model + '.npz'), lat=latitude, lon=longitude, weight=weights)
                    grid_report.append({'model': model, 'nlat': len(latitude), 'nlon': len(longitude),
                                        'method': area_method, 'normalized_area': float(weights.sum())})
                for current, reference in zip([latitude, longitude, weights], reference_grids[model]):
                    np.testing.assert_allclose(current, reference, atol=1e-12, rtol=0)
                time_bounds = np.asarray(dataset.variables[time.bounds][:], dtype=float)
                assert time_bounds.shape == (len(time), 2)
                bounds_dates = nc.num2date(time_bounds[take].reshape(-1), time.units, calendar=time.calendar).reshape(-1, 2)
                days = np.array([(pair[1] - pair[0]).total_seconds() / 86400 for pair in bounds_dates])
                assert np.all(days > 0)
                for local, index in enumerate(take):
                    assert bounds_dates[local, 0].year == source_year[index]
                    assert bounds_dates[local, 0].month == source_month[index] and bounds_dates[local, 0].day == 1
                    assert (bounds_dates[local, 1].year * 12 + bounds_dates[local, 1].month
                            == source_year[index] * 12 + source_month[index] + 1)
                global_mean = np.empty(len(take))
                for offset in range(0, len(take), 120):
                    indices = take[offset:offset + 120]
                    values = np.ma.asarray(field[indices, :, :], dtype=np.float64)
                    assert not np.ma.getmaskarray(values).any(), '全球TOA场缺测，不能缩小分母或填零。'
                    values = np.asarray(values)
                    assert np.isfinite(values).all()
                    global_mean[offset:offset + len(indices)] = np.einsum('tyx,yx->t', values, weights, optimize=True)
                monthly_parts.append(pd.DataFrame({'model': model, 'phase': phases[phase_index], 'variable': variable,
                                                   'source_year': source_year[take], 'month': source_month[take],
                                                   'days': days, 'global_Wm2': global_mean,
                                                   'path': str(path), 'file_index': take}))
            print(model, experiment, variable, len(take), 'months', flush=True)

pd.DataFrame(inventory).to_csv(destination / 'source_inventory.csv', index=False)
pd.DataFrame(grid_report).to_csv(destination / 'grid_weights.csv', index=False)
monthly_long = pd.concat(monthly_parts, ignore_index=True)
assert not monthly_long.duplicated(['model', 'phase', 'variable', 'source_year', 'month']).any()
monthly_long.to_csv(destination / 'monthly_global_components_long.csv', index=False)


# %% 4. 按实际日历日数计算完整年均，记录真实源年与实验年

annual_parts, coverage_rows, sensitivity_rows = [], [], []
for model in models:
    for phase_index, phase in enumerate(phases):
        first_year = starts[model][phase_index]
        selected = monthly_long[(monthly_long.model == model) & (monthly_long.phase == phase)]
        for variable in variables:
            part = selected[selected.variable == variable]
            selected_years = part[(part.source_year >= first_year) & (part.source_year < first_year + 140)]
            ids = selected_years.source_year.to_numpy() * 12 + selected_years.month.to_numpy()
            ids.sort()
            expected = np.arange(first_year * 12 + 1, (first_year + 140) * 12 + 1)
            missing = np.setdiff1d(expected, ids)
            assert len(ids) == len(np.unique(ids))
            coverage_rows.append({'model': model, 'phase': phase, 'variable': variable,
                                  'source_first': first_year, 'source_last': first_year + 139,
                                  'months': len(ids), 'missing_months': len(missing),
                                  'missing_source_ym': ';'.join([str(int(item // 12)) + '-' + str(int(item % 12 + 1)) for item in missing - 1])})
            # 当前只有UK移除阶段的rsdt完全缺失；不以常数或上升气候态代替。
            if len(missing):
                assert model == 'UKESM1-0-LL' and phase == 'removal' and variable == 'rsdt' and len(ids) == 0
        wide = selected.pivot(index=['source_year', 'month'], columns='variable', values='global_Wm2').reindex(columns=variables)
        duration = selected.groupby(['source_year', 'month']).days.agg(['min', 'max'])
        np.testing.assert_allclose(duration['min'], duration['max'], atol=1e-12)
        wide['days'] = duration['min']
        wide['asr'] = wide.rsdt - wide.rsut
        wide['eei'] = wide.asr - wide.rlut
        wide['negative_olr'] = -wide.rlut
        for source_year, year_part in wide.groupby(level='source_year'):
            np.testing.assert_array_equal(year_part.index.get_level_values('month'), np.arange(1, 13))
            record = {'model': model, 'phase': phase, 'source_year': source_year,
                      'simulation_year': int(source_year - first_year + 1 + 140 * phase_index),
                      'days': year_part.days.sum(), 'months': 12}
            for variable in variables + ['asr', 'negative_olr', 'eei']:
                # 固定12月分母；任何缺月/缺变量保留NaN，禁止sum跳过缺测。
                record[variable] = float(np.average(year_part[variable], weights=year_part.days))
            if model == 'GFDL-ESM4' and phase == 'increase':
                sensitivity_rows.append(record.copy())
            if source_year >= first_year:
                annual_parts.append(record)

annual = pd.DataFrame(annual_parts).sort_values(['model', 'simulation_year']).reset_index(drop=True)
assert len(annual) == 2240
pd.DataFrame(coverage_rows).to_csv(destination / 'coverage.csv', index=False)
pd.DataFrame(sensitivity_rows).to_csv(destination / 'gfdl_source_year_sensitivity.csv', index=False)


# %% 5. 初始参照异常、分阶段平滑与固定成员组平均

for model in models:
    indices = annual.index[annual.model == model]
    initial = annual.loc[indices].query('simulation_year <= 10')
    for variable in variables + ['asr', 'negative_olr', 'eei']:
        baseline = float(initial[variable].mean())
        annual.loc[indices, variable + '_anomaly'] = annual.loc[indices, variable] - baseline
        # 不跨CO2转折拼接平滑，尤其不能混合GFDL的不同源窗口。
        for phase in phases:
            phase_indices = annual.index[(annual.model == model) & (annual.phase == phase)]
            annual.loc[phase_indices, variable + '_smooth11'] = annual.loc[phase_indices, variable].rolling(11, center=True, min_periods=6).mean()
            annual.loc[phase_indices, variable + '_anomaly_smooth11'] = annual.loc[phase_indices, variable + '_anomaly'].rolling(11, center=True, min_periods=6).mean()
    group = next(group for group, group_members in registry['members'].items() if model in group_members)
    annual.loc[indices, 'group'] = group

group_rows = []
for group, group_members in registry['members'].items():
    for year in range(1, 281):
        part = annual[(annual.group == group) & (annual.simulation_year == year)]
        record = {'group': group, 'simulation_year': year, 'expected_members': len(group_members),
                  'available_members': int(part.eei.notna().sum())}
        for variable in ['eei', 'asr_anomaly', 'negative_olr_anomaly', 'eei_anomaly', 'eei_smooth11']:
            record[variable] = float(part[variable].mean()) if part[variable].notna().sum() == len(group_members) else np.nan
        group_rows.append(record)
pd.DataFrame(group_rows).to_csv(destination / 'group_annual.csv', index=False)
annual.to_csv(destination / 'annual_global_energy.csv', index=False)


# %% 6. 未平滑端点与过零诊断，不自动归因于AMOC

endpoint_rows, event_rows = [], []
windows = {'initial': [1, 10], 'peak_co2': [131, 140], 'early_removal': [141, 150], 'late': [271, 280]}
for model in models:
    model_part = annual[annual.model == model]
    for label, (first, last) in windows.items():
        part = model_part[model_part.simulation_year.between(first, last)]
        record = {'model': model, 'window': label, 'first_year': first, 'last_year': last,
                  'valid_eei_years': int(part.eei.notna().sum())}
        for variable in variables + ['asr', 'negative_olr', 'eei', 'asr_anomaly', 'negative_olr_anomaly', 'eei_anomaly']:
            record[variable] = float(part[variable].mean()) if part[variable].notna().sum() == 10 else np.nan
        endpoint_rows.append(record)
    removal = model_part[model_part.phase == 'removal']
    negative = removal[removal.eei_smooth11 < 0]
    peak = model_part.loc[model_part.eei_smooth11.idxmax()]
    minimum = model_part.loc[model_part.eei_smooth11.idxmin()]
    event_rows.append({'model': model, 'first_negative_smoothed_removal_year': None if negative.empty else int(negative.iloc[0].simulation_year),
                       'max_smooth11_Wm2': peak.eei_smooth11, 'max_year': int(peak.simulation_year),
                       'min_smooth11_Wm2': minimum.eei_smooth11, 'min_year': int(minimum.simulation_year)})
pd.DataFrame(endpoint_rows).to_csv(destination / 'endpoint_summary.csv', index=False)
pd.DataFrame(event_rows).to_csv(destination / 'event_summary.csv', index=False)
closure = annual.eei - annual.rsdt + annual.rsut + annual.rlut
assert closure.abs().max() < 1e-10
metadata = {'status': 'passed_with_explicit_missing_UKESM_removal_rsdt',
            'formula': 'EEI = global(rsdt - rsut - rlut); positive downward', 'units': 'W m-2',
            'models': models, 'groups': registry['members'], 'group_registry': str(registry_path),
            'starts_increase_removal': starts, 'no_drift_correction': True, 'interpolation': False,
            'annual_weighting': 'CF time-bound duration; 12 complete months required',
            'spatial_weighting': 'Global land plus ocean, native latitude-bound areas; CNRM verified Gaussian quadrature',
            'smoothing': '11-year centered mean, at least 6 years, separately within each phase',
            'baseline': 'Each model increase simulation years 1-10; subtraction is not piControl correction',
            'gfdl_note': 'User convention: upward source6-145->sim1-140; CDR141-280->sim141-280. Not an uninterrupted absolute calendar.',
            'ukesm_note': 'Local rsdt for CDR absent; no removal EEI or ASR, no fixed-four-member group mean after140.',
            'annual_rows': len(annual), 'complete_eei_years': int(annual.eei.notna().sum()),
            'max_identity_error_Wm2': float(closure.abs().max()),
            'code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(metadata, ensure_ascii=False, indent=2), flush=True)
