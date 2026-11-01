# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import netCDF4 as nc
import numpy as np
import pandas as pd


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
destination = root / 'output/s05_atmosphere/p07_toa_energy_balance/n20_cas_toa_radiation'
assert not destination.exists(), '正式输出存在，须先归档。'
model, member = 'CAS-ESM2-0', 'r1i1p1f1'
variables = ['rsdt', 'rsut', 'rlut']
experiments = ['1pctCO2', '1pctCO2-cdr']
first_years = [200, 340]
chunk_months = 24
monthly_parts, inventory, annual_parts = [], [], []


# %% 3. 核对实际维度与原生全球面积权重

reference_path = next(Path('F:/sea/cmip6rsdtad').rglob('rsdt_Amon_CAS-ESM2-0_1pctCO2_r1i1p1f1_*.nc'))
with nc.Dataset(reference_path) as dataset:
    latitude = np.asarray(dataset['lat'][:], dtype=float)
    longitude = np.asarray(dataset['lon'][:], dtype=float)
    latitude_bounds = np.asarray(dataset[dataset['lat'].bounds][:], dtype=float)
    longitude_bounds = np.asarray(dataset[dataset['lon'].bounds][:], dtype=float)
np.testing.assert_allclose(latitude_bounds[1:, 0], latitude_bounds[:-1, 1], atol=1e-5)
np.testing.assert_allclose(latitude_bounds[[0, -1], [0, 1]], [-90, 90], atol=1e-5)
np.testing.assert_allclose(np.diff(longitude), 360 / len(longitude), atol=1e-4)
np.testing.assert_allclose(np.diff(longitude_bounds, axis=1)[:, 0], 360 / len(longitude), atol=1e-4)
latitude_weights = np.sin(np.deg2rad(latitude_bounds[:, 1])) - np.sin(np.deg2rad(latitude_bounds[:, 0]))
weights = latitude_weights[:, None] * np.ones((1, len(longitude)))
weights /= weights.sum()
assert weights.shape == (128, 256) and np.all(weights > 0)
np.testing.assert_allclose(weights.sum(), 1., atol=1e-12)


# %% 4. 分块读取两个140年强迫阶段，CF月界时长加权

for variable in variables:
    for phase_index, (experiment, suffix) in enumerate(zip(experiments, ['ad', ''])):
        paths = sorted(Path('F:/sea/cmip6' + variable + suffix).rglob(
            f'{variable}_Amon_{model}_{experiment}_{member}_*.nc'))
        first_year = first_years[phase_index]
        for path in paths:
            with nc.Dataset(path) as dataset:
                assert dataset.source_id == model and dataset.experiment_id == experiment
                assert dataset.variant_label == member and dataset.variable_id == variable
                field, time = dataset[variable], dataset['time']
                assert field.dimensions == ('time', 'lat', 'lon') and field.units == 'W m-2'
                np.testing.assert_array_equal(dataset['lat'][:], latitude)
                np.testing.assert_array_equal(dataset['lon'][:], longitude)
                dates = nc.num2date(time[:], time.units, calendar=time.calendar)
                years = np.array([date.year for date in dates])
                months = np.array([date.month for date in dates])
                selected = np.flatnonzero((years >= first_year) & (years < first_year + 140))
                assert len(selected) == 1680
                np.testing.assert_array_equal(years[selected] * 12 + months[selected] - 1,
                                               np.arange(first_year * 12, (first_year + 140) * 12))
                bounds_dates = nc.num2date(dataset[time.bounds][selected].reshape(-1), time.units,
                                          calendar=time.calendar).reshape(-1, 2)
                days = np.array([(pair[1] - pair[0]).total_seconds() / 86400 for pair in bounds_dates])
                assert np.all(days > 0)
                global_mean = np.empty(len(selected))
                for offset in range(0, len(selected), chunk_months):
                    indices = selected[offset:offset + chunk_months]
                    values = np.ma.asarray(field[indices], dtype=np.float64)
                    assert not np.ma.getmaskarray(values).any()
                    values = np.asarray(values)
                    assert np.isfinite(values).all()
                    global_mean[offset:offset + len(indices)] = np.einsum('tyx,yx->t', values, weights, optimize=True)
                monthly_parts.append(pd.DataFrame({'model': model, 'variable': variable,
                                                   'phase': 'increase' if phase_index == 0 else 'removal',
                                                   'source_year': years[selected], 'month': months[selected],
                                                   'simulation_year': years[selected] - first_year + 1 + phase_index * 140,
                                                   'days': days, 'global_Wm2': global_mean}))
                inventory.append({'path': str(path), 'variable': variable, 'experiment': experiment,
                                  'member': member, 'dimensions': list(field.dimensions),
                                  'shape': list(field.shape), 'units': field.units, 'calendar': time.calendar,
                                  'time_units': time.units, 'first_selected': str(dates[selected[0]]),
                                  'last_selected': str(dates[selected[-1]]), 'selected_months': len(selected),
                                  'parent_experiment_id': dataset.parent_experiment_id,
                                  'parent_time_units': dataset.parent_time_units,
                                  'branch_time_in_parent': float(dataset.branch_time_in_parent),
                                  'branch_time_in_child': float(dataset.branch_time_in_child),
                                  'bytes': path.stat().st_size, 'mtime_ns': path.stat().st_mtime_ns})
            print(variable, experiment, len(selected), 'months complete', flush=True)


# %% 5. 年平均、长短波分解与初始异常，尚不作piControl订正

monthly = pd.concat(monthly_parts, ignore_index=True)
assert not monthly.duplicated(['variable', 'simulation_year', 'month']).any()
for variable in variables:
    selected = monthly.loc[monthly.variable.eq(variable)].copy()
    selected['weighted_flux'] = selected.global_Wm2 * selected.days
    annual_variable = selected.groupby('simulation_year').agg(
        source_year=('source_year', 'first'), months=('month', 'count'),
        days=('days', 'sum'), weighted=('weighted_flux', 'sum')).reset_index()
    assert annual_variable.months.eq(12).all()
    annual_variable[variable] = annual_variable.weighted / annual_variable.days
    annual_parts.append(annual_variable[['simulation_year', 'source_year', variable]])
annual = annual_parts[0].merge(annual_parts[1], on=['simulation_year', 'source_year']).merge(
    annual_parts[2], on=['simulation_year', 'source_year'])
annual.insert(0, 'model', model)
annual['phase'] = np.where(annual.simulation_year <= 140, 'increase', 'removal')
annual['asr'] = annual.rsdt - annual.rsut
annual['negative_olr'] = -annual.rlut
annual['eei'] = annual.asr + annual.negative_olr
for quantity in variables + ['asr', 'negative_olr', 'eei']:
    annual[quantity + '_initial_anomaly'] = annual[quantity] - annual.loc[annual.simulation_year.between(1, 10), quantity].mean()
    for phase in ['increase', 'removal']:
        mask = annual.phase.eq(phase)
        for suffix in ['', '_initial_anomaly']:
            annual.loc[mask, quantity + suffix + '_smooth11'] = annual.loc[mask, quantity + suffix].rolling(
                11, center=True, min_periods=6).mean()
np.testing.assert_allclose(annual.eei, annual.rsdt - annual.rsut - annual.rlut, atol=1e-11)
np.testing.assert_array_equal(annual.simulation_year, np.arange(1, 281))


# %% 6. 保存原值和完整处理口径

destination.mkdir(parents=True)
monthly.to_csv(destination / 'monthly_global_long.csv', index=False)
annual.to_csv(destination / 'annual_global_energy.csv', index=False)
pd.DataFrame(inventory).to_csv(destination / 'source_inventory.csv', index=False)
np.savez_compressed(destination / 'native_weights.npz', latitude=latitude, longitude=longitude, weight=weights)
metadata = {'status': 'raw_forced_complete_control_pending', 'model': model,
            'member': member, 'units': 'W m-2', 'formula': 'EEI=rsdt-rsut-rlut; ASR=rsdt-rsut; LW contribution=-rlut',
            'dimensions': 'raw:time xlat128 xlon256; annual:year280',
            'phase_years': {'increase_source': [200, 339], 'removal_source': [340, 479]},
            'spatial_weights': 'Native latitude-bound spherical cell areas, uniform cyclic longitude.',
            'annual_weights': 'CF month-bound durations; twelve complete months and all global cells required.',
            'drift_correction': 'none; no CAS piControl radiation files located yet',
            'control_alignment': 'CAS forced parent attributes retained; their mixed calendar origins require checking against actual control files, not a guessed window.',
            'smoothing': 'Centered11year separately within phases, min_periods6; no initial re-zeroing for raw series.',
            'source_files': inventory,
            'code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(destination / 'metadata.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
print(annual.loc[annual.simulation_year.isin([1, 140, 280]), ['simulation_year', 'eei', 'asr', 'negative_olr']].to_string(index=False))
print(destination.resolve(), flush=True)
