# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import xarray as xr
from scipy.integrate import trapezoid
from scipy.stats import false_discovery_control
from threadpoolctl import threadpool_limits


# %% 2. 路径与固定六模式

root = Path('F:/sea/Antic sea ice constraint')
base = root / 'output/s05_atmosphere/p05_amoc_recovery_composites'
sample_source = base / 'n72_six_model_annual_circulation_precipitation_samples'
statistics_source = base / 'n73_six_model_precipitation_significance'
destination = base / 'n74_six_model_circulation_composites'
assert not destination.exists(), '正式结果存在，须先归档。'
models = ['GFDL-ESM4', 'NorESM2-LM', 'CAS-ESM2-0', 'CESM2', 'CNRM-ESM2-1', 'UKESM1-0-LL']
groups = {'recovered': models[:3], 'not_recovered': models[3:]}
model, member = 'CAS-ESM2-0', 'r1i1p1f1'
periods, source_starts = ['early', 'later'], [200, 470]
pressure_levels = np.array([1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000,
                            30000, 40000, 50000, 60000, 70000, 85000, 92500, 100000], dtype=float)
latitude = np.arange(-88., 90., 2.)
radius, gravity = 6371000., 9.80665
replicates, block_lengths, alpha = 9999, [2, 3, 4], .05
inventory, integration_rows, independent_rows, phases = [], [], [], []
with xr.open_dataset(sample_source / 'available_circulation_samples.nc') as dataset:
    previous_psi = dataset.psi.load()
with xr.open_dataset(sample_source / 'annual_precipitation_samples.nc') as dataset:
    days = dataset.day_weight.sel(model=models).load()
with xr.open_dataarray(sample_source / 'group_precipitation_changes.nc') as dataset:
    precipitation_changes = dataset.load()


# %% 3. CAS完整初末十年va和ps积分，不重复计算已有五模式

for period_index, period in enumerate(periods):
    experiment = ['1pctCO2', '1pctCO2-cdr'][period_index]
    va_root = Path('F:/sea/cmip6vaad/CAS' if period_index == 0 else 'F:/sea/cmip6va/CAS')
    va_paths = sorted(va_root.glob(f'va_Amon_{model}_{experiment}_{member}_*.nc'))
    # 新下载的升高段ps也在cmip6ps/CAS，按明确实验名筛选，不按文件夹名推断实验。
    ps_paths = sorted(Path('F:/sea/cmip6ps/CAS').glob(f'ps_Amon_{model}_{experiment}_{member}_*.nc'))
    annual_fields = []
    for year_index, year in enumerate(range(source_starts[period_index], source_starts[period_index] + 10)):
        fields = {}
        for variable, paths in [('va', va_paths), ('ps', ps_paths)]:
            pieces = []
            for path in paths:
                with xr.open_dataset(path, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as dataset:
                    assert dataset.attrs['source_id'] == model and dataset.attrs['experiment_id'] == experiment
                    assert dataset.attrs['variant_label'] == member
                    expected_dimensions = ('time', 'plev', 'lat', 'lon') if variable == 'va' else ('time', 'lat', 'lon')
                    assert dataset[variable].dims == expected_dimensions
                    assert dataset[variable].attrs['units'] == ('m s-1' if variable == 'va' else 'Pa')
                    take = dataset.time.dt.year == year
                    if not bool(take.any()):
                        continue
                    selected = dataset[variable].isel(time=take)
                    if variable == 'va':
                        assert dataset.plev.attrs['units'] == 'Pa'
                        selected = selected.sel(plev=pressure_levels)
                    pieces.append(selected.load().astype(float))
                    inventory.append({'model': model, 'member': member, 'period': period, 'year': year,
                                      'variable': variable, 'path': str(path), 'dimensions': str(dataset[variable].dims),
                                      'shape': str(dataset[variable].shape), 'units': dataset[variable].attrs['units'],
                                      'months': int(take.sum()), 'calendar': dataset.time.dt.calendar,
                                      'file_bytes': path.stat().st_size, 'mtime_ns': path.stat().st_mtime_ns})
            fields[variable] = xr.concat(pieces, 'time').sortby('time')
            np.testing.assert_array_equal(fields[variable].time.dt.month, np.arange(1, 13))
            assert fields[variable].time.dt.year.values.tolist() == [year] * 12
        wind, surface = fields['va'], fields['ps']
        np.testing.assert_array_equal(wind.lat, surface.lat)
        np.testing.assert_array_equal(wind.lon, surface.lon)
        np.testing.assert_array_equal(wind.time, surface.time)
        assert wind.shape == (12, 17, 128, 256) and np.isfinite(surface).all()
        ps = surface.values
        velocity = wind.where(wind.plev <= surface).transpose('time', 'plev', 'lat', 'lon').values
        valid = np.isfinite(velocity)
        assert valid[:, 0].all() and not np.any((~valid[:, :-1]) & valid[:, 1:])
        filled = velocity.copy()
        for level_index in range(1, len(pressure_levels)):
            filled[:, level_index] = np.where(valid[:, level_index], filled[:, level_index], filled[:, level_index - 1])
        integral = np.zeros((12, 128, 256))
        profiles = [integral.mean(axis=2).copy()]
        for level_index in range(len(pressure_levels) - 1):
            width = pressure_levels[level_index + 1] - pressure_levels[level_index]
            thickness = np.clip(ps - pressure_levels[level_index], 0., width)
            bottom_velocity = filled[:, level_index] + (filled[:, level_index + 1] - filled[:, level_index]) * thickness / width
            integral += .5 * (filled[:, level_index] + bottom_velocity) * thickness
            profiles.append(integral.mean(axis=2).copy())
        psi_monthly = np.stack(profiles, axis=1) * (2 * np.pi * radius * np.cos(np.deg2rad(wind.lat.values)) / gravity)[None, None] / 1e9
        month_days = wind.time.dt.days_in_month.values.astype(float)
        assert month_days.sum() == float(days.sel(model=model, period=period, year_in_window=year_index + 1))
        field = xr.DataArray(np.average(psi_monthly, axis=0, weights=month_days), dims=['plev', 'lat'],
                             coords={'plev': pressure_levels, 'lat': wind.lat}).interp(lat=latitude)
        assert np.isfinite(field).all() and np.max(np.abs(psi_monthly[:, 0])) == 0
        annual_fields.append(field.expand_dims(year_in_window=[year_index + 1]))
        integration_rows.append({'model': model, 'period': period, 'year': year,
                                  'top_boundary_absmax': float(abs(psi_monthly[:, 0]).max()),
                                  'monthly_absmax_1e9kgs': float(abs(psi_monthly).max())})
        if year_index == 0:
            tropical = np.abs(wind.lat.values) <= 14
            upper = pressure_levels <= 30000
            upper_velocity = velocity[:, upper][:, :, tropical]
            assert np.isfinite(upper_velocity).all()
            independent = trapezoid(upper_velocity.mean(axis=3), x=pressure_levels[upper], axis=1)
            independent *= (2 * np.pi * radius / gravity * np.cos(np.deg2rad(wind.lat.values[tropical])))[None] / 1e9
            independent_annual = np.average(independent, axis=0, weights=month_days)
            for target in [-10., 10.]:
                expected = np.interp(target, wind.lat.values[tropical], independent_annual)
                actual = float(field.sel(plev=30000, lat=target))
                assert abs(actual - expected) < 1e-8
                independent_rows.append({'model': model, 'period': period, 'year': year, 'lat': target,
                                         'pressure_Pa': 30000, 'independent_psi': expected, 'saved_psi': actual,
                                         'error_1e9kgs': actual - expected})
        print(model, period, 'circulation year', year_index + 1, 'of10 complete', flush=True)
    phases.append(xr.concat(annual_fields, 'year_in_window').expand_dims(period=[period]))
cas_psi = xr.concat(phases, 'period').expand_dims(model=[model])
psi = xr.concat([previous_psi, cas_psi], 'model').sel(model=models).transpose('model', 'period', 'year_in_window', 'plev', 'lat').rename('psi').drop_encoding()
assert psi.shape == (6, 2, 10, 17, 89)


# %% 4. 六模式初末气候态与三对三组均值，核对初始强度量级

climate = psi.weighted(days).mean('year_in_window', skipna=False)
member_changes = climate.sel(period='later', drop=True) - climate.sel(period='early', drop=True)
strength_rows = []
for selected_model in models:
    for period in periods:
        field = climate.sel(model=selected_model, period=period)
        record = {'model': selected_model, 'period': period,
                  'NH500_max_1e9kgs': float(field.sel(plev=50000, lat=slice(0, 30)).max()),
                  'SH500_min_1e9kgs': float(field.sel(plev=50000, lat=slice(-30, 0)).min()),
                  'tropical_max_1e9kgs': float(field.sel(plev=slice(20000, 70000), lat=slice(-30, 30)).max()),
                  'tropical_min_1e9kgs': float(field.sel(plev=slice(20000, 70000), lat=slice(-30, 30)).min()),
                  'bottom_residual_absmax_1e9kgs': float(abs(field.sel(plev=100000)).max())}
        assert record['NH500_max_1e9kgs'] > 0 and record['SH500_min_1e9kgs'] < 0
        if period == 'early':
            assert 30 < record['tropical_max_1e9kgs'] < 300 and 30 < -record['tropical_min_1e9kgs'] < 300
        strength_rows.append(record)
group_climate = xr.concat([climate.sel(model=group_members).mean('model', skipna=False).expand_dims(product=[group])
                           for group, group_members in groups.items()], 'product').rename('psi')
group_changes = group_climate.sel(period='later', drop=True) - group_climate.sel(period='early', drop=True)
contrast = (group_changes.sel(product='not_recovered', drop=True) - group_changes.sel(product='recovered', drop=True)).expand_dims(product=['contrast'])
group_changes = xr.concat([group_changes, contrast], 'product').rename('psi')
destination.mkdir(parents=True)
xr.Dataset({'psi': psi, 'day_weight': days}).to_netcdf(destination / 'annual_circulation_samples.nc')
climate.rename('psi').to_netcdf(destination / 'member_circulation_climatology.nc')
member_changes.rename('psi').to_netcdf(destination / 'member_circulation_changes.nc')
group_climate.to_netcdf(destination / 'group_circulation_climatology.nc')
group_changes.to_netcdf(destination / 'group_circulation_changes.nc')
pd.DataFrame(inventory).to_csv(destination / 'cas_source_inventory.csv', index=False)
pd.DataFrame(integration_rows).to_csv(destination / 'cas_integration_qc.csv', index=False)
pd.DataFrame(independent_rows).to_csv(destination / 'cas_independent_300hPa_check.csv', index=False)
pd.DataFrame(strength_rows).to_csv(destination / 'six_model_strength_audit.csv', index=False)


# %% 5. 复用六模式降水抽样索引，对环流组差做同口径检验

with np.load(statistics_source / 'bootstrap_year_indices.npz') as saved:
    indices_archive = {name: saved[name].copy() for name in saved.files}
display = psi.sel(plev=slice(10000, 100000))
values = display.values.reshape(6, 2, 10, -1)
day_array = days.values
observed = group_changes.sel(product='contrast', plev=display.plev).values.ravel()
coefficients = np.array([-1/3] * 3 + [1/3] * 3)
probabilities = np.empty((3, values.shape[-1]))
for block_index, length in enumerate(block_lengths):
    indices = indices_archive['L' + str(length)]
    weights = np.empty((6, 2, replicates, 10))
    for model_index in range(6):
        for period_index in range(2):
            counts = np.eye(10)[indices[model_index, period_index]].sum(axis=1) * day_array[model_index, period_index][None]
            weights[model_index, period_index] = counts / counts.sum(axis=1, keepdims=True)
    for offset in range(0, values.shape[-1], 256):
        last = min(offset + 256, values.shape[-1])
        with threadpool_limits(limits=4):
            members_bootstrap = np.stack([weights[index, 1] @ values[index, 1, :, offset:last]
                                          - weights[index, 0] @ values[index, 0, :, offset:last] for index in range(6)])
            bootstrap = np.einsum('m,mbc->bc', coefficients, members_bootstrap, optimize=True)
            errors = bootstrap - observed[offset:last][None]
        probabilities[block_index, offset:last] = (1 + np.sum(abs(errors) >= abs(observed[offset:last])[None], axis=0)) / (replicates + 1)
    print('Six-model psi contrast bootstrap block', length, 'complete', flush=True)
p_max = probabilities.max(axis=0)
q_values = false_discovery_control(p_max, method='bh')
shape = (display.sizes['plev'], display.sizes['lat'])
coordinates = {'plev': display.plev, 'lat': display.lat}
results = xr.Dataset()
for name, values_flat in [('delta', observed), ('p_max', p_max), ('q', q_values),
                          ('tested', np.ones_like(q_values, dtype=bool)), ('significant', q_values < alpha),
                          ('nonsignificant', q_values >= alpha)]:
    results[name] = xr.DataArray(values_flat.reshape(shape), dims=['plev', 'lat'], coords=coordinates)
results['p_by_block'] = xr.DataArray(probabilities.reshape((3,) + shape), dims=['block_length', 'plev', 'lat'],
                                     coords={'block_length': block_lengths, **coordinates})
results = results.expand_dims(product=['contrast'])
results.attrs = {'units': '1e9 kg s-1', 'test': 'Centered two-sided circular-year-block bootstrap; max p across L2/3/4, then BH over displayed100-1000hPa x88S-88N cells.',
                 'inference': 'Temporal sampling conditional on six fixed model identities.'}
results.to_netcdf(destination / 'psi_significance.nc')
metadata = {'status': 'passed', 'models': models, 'groups': groups, 'group_sizes': [3, 3],
            'CAS_starts': source_starts, 'raw_CAS_dimensions': {'va': 'time,plev19,lat128,lon256;m s-1', 'ps': 'time,lat128,lon256;Pa'},
            'output_dimensions': 'psi(model6,period2,year_in_window10,plev17,lat89)',
            'psi_method': 'Inherited10hPa-top, terrain-clipped trapezoidal pressure integral of zonal mean meridional wind;2pi*a*cos(lat)/g; units1e9kg/s; no mass adjustment.',
            'strength_gate': 'All six initial tropical extrema30-300 in1e9kg/s with expected signs. Unit/order screen only, not observational validation.',
            'bootstrap': '9999 centered two-sided circular year blocks; same independent model/period year draws as n73 precipitation; max p over2/3/4years then per-panel global BH q<0.05.',
            'psi_tested_panel': 'Only non-recovery minus recovery contrast; no new significance marks on the two individual circulation group panels.',
            'precipitation_source': str(statistics_source), 'drift_correction': 'none',
            'contrast': 'Non-recovery minus recovery residual responses, each group three fixed equal-weight models.',
            'limitations': ['Each endpoint only10annual samples; temporal inference is conditional on fixed model identities.',
                           '10hPa finite top and nonzero near-surface mass residual; no column mass correction.',
                           'CAS branch metadata anomalies remain an input provenance caveat, not repaired by offsets.'],
            'source_files': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in
                             [sample_source / 'available_circulation_samples.nc', sample_source / 'annual_precipitation_samples.nc',
                              statistics_source / 'bootstrap_year_indices.npz', statistics_source / 'pr_significance.nc']]}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print(pd.DataFrame(strength_rows).to_string(index=False), flush=True)
print('Psi contrast significant:', int(np.sum(q_values < alpha)), 'of', len(q_values), flush=True)
print(destination.resolve(), flush=True)
