# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import netCDF4 as nc
import numpy as np
import pandas as pd
import xarray as xr
from scipy.integrate import trapezoid


# %% 2. 路径、固定六模式及原始维序

root = Path('F:/sea/Antic sea ice constraint')
base = root / 'output/s05_atmosphere/p05_amoc_recovery_composites'
destination = base / 'n72_six_model_annual_circulation_precipitation_samples'
assert not destination.exists(), '正式结果存在，须先归档。'
models = ['GFDL-ESM4', 'NorESM2-LM', 'CAS-ESM2-0', 'CESM2', 'CNRM-ESM2-1', 'UKESM1-0-LL']
groups = {'recovered': models[:3], 'not_recovered': models[3:]}
periods = ['early', 'later']
starts = {'GFDL-ESM4': [1, 271], 'NorESM2-LM': [1, 271], 'CAS-ESM2-0': [200, 470],
          'CESM2': [1, 131], 'CNRM-ESM2-1': [1850, 2120], 'UKESM1-0-LL': [1850, 2120]}
members = {'CAS-ESM2-0': 'r1i1p1f1', 'UKESM1-0-LL': 'r1i1p1f2'}
pressure_levels = np.array([1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000,
                            30000, 40000, 50000, 60000, 70000, 85000, 92500, 100000], dtype=float)
latitude, longitude = np.arange(-88., 90., 2.), np.arange(-180., 180., 2.)
radius, gravity = 6371000., 9.80665
pr_parts, day_parts, psi_parts, inventory, strength_rows, independent_rows = [], [], [], [], [], []


# %% 3. 复用四模式原逐年降水与流函数，不改其窗口和计算

with xr.open_dataset(base / 'n71_circulation_precipitation_significance/annual_precipitation_samples.nc') as dataset:
    original_pr = dataset.pr.load()
    original_days = dataset.day_weight.load()
for model in ['GFDL-ESM4', 'NorESM2-LM', 'CESM2', 'CNRM-ESM2-1']:
    pr_parts.append(original_pr.sel(model=[model]))
    day_parts.append(original_days.sel(model=[model]))
    with xr.open_dataset(base / 'n35_global_meridional_samples' / f'{model}_samples.nc') as dataset:
        psi_parts.append(dataset.psi.sel(model=[model], season='Annual', drop=True).load())


# %% 4. CAS和UK降水重建完整初末十年，原生年平均后双线性到2度

for model in ['CAS-ESM2-0', 'UKESM1-0-LL']:
    model_phases, model_day_phases = [], []
    for period_index, period in enumerate(periods):
        folder = Path('F:/sea/cmip6prad' if period_index == 0 else 'F:/sea/cmip6pr')
        paths = sorted(folder.rglob(f'pr_Amon_{model}_{["1pctCO2", "1pctCO2-cdr"][period_index]}_{members[model]}_*.nc'))
        annual_maps, annual_days = [], []
        for year_index, year in enumerate(range(starts[model][period_index], starts[model][period_index] + 10)):
            pieces, dates_all = [], []
            for path in paths:
                with xr.open_dataset(path, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as dataset:
                    assert dataset.attrs['source_id'] == model and dataset.attrs['variant_label'] == members[model]
                    assert dataset.pr.dims == ('time', 'lat', 'lon') and dataset.pr.attrs['units'] == 'kg m-2 s-1'
                    take = dataset.time.dt.year == year
                    if not bool(take.any()):
                        continue
                    pieces.append(dataset.pr.isel(time=take).load().astype(float) * 86400.)
                    inventory.append({'model': model, 'variable': 'pr', 'period': period, 'year': year,
                                      'path': str(path), 'dimensions': str(dataset.pr.dims),
                                      'shape': str(dataset.pr.shape), 'months': int(take.sum()),
                                      'units': dataset.pr.attrs['units'], 'calendar': dataset.time.dt.calendar,
                                      'bytes': path.stat().st_size, 'mtime_ns': path.stat().st_mtime_ns})
            monthly = xr.concat(pieces, 'time').sortby('time')
            np.testing.assert_array_equal(monthly.time.dt.month, np.arange(1, 13))
            assert monthly.time.dt.year.values.tolist() == [year] * 12 and np.isfinite(monthly).all()
            days = monthly.time.dt.days_in_month.astype(float)
            native = monthly.weighted(days).mean('time', skipna=False)
            native = native.assign_coords(lon=(native.lon + 180) % 360 - 180).sortby('lon').sortby('lat')
            periodic = xr.concat([native.isel(lon=[-1]).assign_coords(lon=native.lon.values[[-1]] - 360),
                                  native, native.isel(lon=[0]).assign_coords(lon=native.lon.values[[0]] + 360)], 'lon')
            common = periodic.interp(lat=latitude, lon=longitude).transpose('lat', 'lon')
            assert np.isfinite(common).all()
            annual_maps.append(common.expand_dims(year_in_window=[year_index + 1]))
            annual_days.append(float(days.sum()))
        model_phases.append(xr.concat(annual_maps, 'year_in_window').expand_dims(period=[period]))
        model_day_phases.append(xr.DataArray(annual_days, dims='year_in_window',
                                            coords={'year_in_window': np.arange(1, 11)}).expand_dims(period=[period]))
        print(model, period, '10 annual precipitation fields complete', flush=True)
    pr_parts.append(xr.concat(model_phases, 'period').expand_dims(model=[model]))
    day_parts.append(xr.concat(model_day_phases, 'period').expand_dims(model=[model]))


# %% 5. UK经向风地形截断积分，风纬度145与气压纬度144先对齐

model = 'UKESM1-0-LL'
model_phases = []
for period_index, period in enumerate(periods):
    experiment = ['1pctCO2', '1pctCO2-cdr'][period_index]
    va_paths = sorted(Path('F:/sea/cmip6vaad' if period_index == 0 else 'F:/sea/cmip6va').rglob(
        f'va_Amon_{model}_{experiment}_{members[model]}_*.nc'))
    ps_paths = sorted(Path('F:/sea/cmip6psad' if period_index == 0 else 'F:/sea/cmip6ps').rglob(
        f'ps_Amon_{model}_{experiment}_{members[model]}_*.nc'))
    annual_samples = []
    for year_index, year in enumerate(range(starts[model][period_index], starts[model][period_index] + 10)):
        fields = {}
        for variable, paths in [('va', va_paths), ('ps', ps_paths)]:
            pieces = []
            for path in paths:
                with xr.open_dataset(path, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as dataset:
                    assert dataset.attrs['source_id'] == model and dataset.attrs['variant_label'] == members[model]
                    expected_dims = ('time', 'plev', 'lat', 'lon') if variable == 'va' else ('time', 'lat', 'lon')
                    assert dataset[variable].dims == expected_dims
                    assert dataset[variable].attrs['units'] == ('m s-1' if variable == 'va' else 'Pa')
                    take = dataset.time.dt.year == year
                    if not bool(take.any()):
                        continue
                    selected = dataset[variable].isel(time=take)
                    if variable == 'va':
                        # 丢弃精确南北极的交错风点；目标88S–88N仍由内部点插值，不外推ps到极点。
                        selected = selected.sel(plev=pressure_levels, lat=slice(-88.75, 88.75))
                    pieces.append(selected.load().astype(float))
                    inventory.append({'model': model, 'variable': variable, 'period': period, 'year': year,
                                      'path': str(path), 'dimensions': str(dataset[variable].dims),
                                      'shape': str(dataset[variable].shape), 'months': int(take.sum()),
                                      'units': dataset[variable].attrs['units'], 'calendar': dataset.time.dt.calendar,
                                      'bytes': path.stat().st_size, 'mtime_ns': path.stat().st_mtime_ns})
            fields[variable] = xr.concat(pieces, 'time').sortby('time')
            np.testing.assert_array_equal(fields[variable].time.dt.month, np.arange(1, 13))
            assert fields[variable].time.dt.year.values.tolist() == [year] * 12
        wind = fields['va']
        # UK地面气压与va在经纬向均交错；周期扩展经度，避免0°风点被当作越界。
        native_ps = fields['ps']
        periodic_ps = xr.concat([native_ps.isel(lon=[-1]).assign_coords(lon=native_ps.lon.values[[-1]] - 360),
                                 native_ps, native_ps.isel(lon=[0]).assign_coords(lon=native_ps.lon.values[[0]] + 360)], 'lon')
        surface = periodic_ps.interp(lat=wind.lat, lon=wind.lon)
        assert np.isfinite(surface).all()
        ps = surface.values
        velocity = wind.where(wind.plev <= surface).transpose('time', 'plev', 'lat', 'lon').values
        valid = np.isfinite(velocity)
        assert valid[:, 0].all() and not np.any((~valid[:, :-1]) & valid[:, 1:])
        filled = velocity.copy()
        for level_index in range(1, len(pressure_levels)):
            filled[:, level_index] = np.where(valid[:, level_index], filled[:, level_index], filled[:, level_index - 1])
        integral = np.zeros((12, wind.sizes['lat'], wind.sizes['lon']))
        profiles = [integral.mean(axis=2).copy()]
        for level_index in range(len(pressure_levels) - 1):
            width = pressure_levels[level_index + 1] - pressure_levels[level_index]
            thickness = np.clip(ps - pressure_levels[level_index], 0., width)
            bottom_velocity = filled[:, level_index] + (filled[:, level_index + 1] - filled[:, level_index]) * thickness / width
            integral += .5 * (filled[:, level_index] + bottom_velocity) * thickness
            profiles.append(integral.mean(axis=2).copy())
        psi_monthly = np.stack(profiles, axis=1) * (2 * np.pi * radius * np.cos(np.deg2rad(wind.lat.values)) / gravity)[None, None] / 1e9
        days = wind.time.dt.days_in_month.values.astype(float)
        field = xr.DataArray(np.average(psi_monthly, axis=0, weights=days), dims=['plev', 'lat'],
                             coords={'plev': pressure_levels, 'lat': wind.lat}).interp(lat=latitude)
        assert np.isfinite(field).all() and np.max(np.abs(psi_monthly[:, 0])) == 0
        annual_samples.append(field.expand_dims(year_in_window=[year_index + 1]))
        if year_index == 0:
            tropical = np.abs(wind.lat.values) <= 14
            upper = pressure_levels <= 30000
            upper_velocity = velocity[:, upper][:, :, tropical]
            assert np.isfinite(upper_velocity).all()
            independent = trapezoid(upper_velocity.mean(axis=3), x=pressure_levels[upper], axis=1)
            independent *= (2 * np.pi * radius / gravity * np.cos(np.deg2rad(wind.lat.values[tropical])))[None] / 1e9
            independent_annual = np.average(independent, axis=0, weights=days)
            for target in [-10., 10.]:
                expected = np.interp(target, wind.lat.values[tropical], independent_annual)
                actual = float(field.sel(plev=30000, lat=target))
                assert abs(actual - expected) < 1e-8
                independent_rows.append({'model': model, 'period': period, 'year': year, 'lat': target,
                                         'pressure_Pa': 30000, 'independent_psi': expected, 'saved_psi': actual,
                                         'error_1e9kgs': actual - expected})
        print(model, period, 'circulation year', year_index + 1, 'of10 complete', flush=True)
    model_phases.append(xr.concat(annual_samples, 'year_in_window').expand_dims(period=[period]))
psi_parts.append(xr.concat(model_phases, 'period').expand_dims(model=[model]))


# %% 6. 保存六模式降水与已具备的五模式环流，不把恢复组缺失成员排除后平均

precipitation = xr.concat(pr_parts, 'model').sel(model=models).transpose('model', 'period', 'year_in_window', 'lat', 'lon').rename('pr')
day_weight = xr.concat(day_parts, 'model').sel(model=models).transpose('model', 'period', 'year_in_window').rename('day_weight')
circulation = xr.concat(psi_parts, 'model').transpose('model', 'period', 'year_in_window', 'plev', 'lat').rename('psi').drop_encoding()
uk = circulation.sel(model=model)
uk_days = day_weight.sel(model=model)
uk_climate = uk.weighted(uk_days).mean('year_in_window', skipna=False)
for period in periods:
    field = uk_climate.sel(period=period)
    record = {'model': model, 'period': period,
              'NH500_max_1e9kgs': float(field.sel(plev=50000, lat=slice(0, 30)).max()),
              'SH500_min_1e9kgs': float(field.sel(plev=50000, lat=slice(-30, 0)).min()),
              'tropical_max_1e9kgs': float(field.sel(plev=slice(20000, 70000), lat=slice(-30, 30)).max()),
              'tropical_min_1e9kgs': float(field.sel(plev=slice(20000, 70000), lat=slice(-30, 30)).min()),
              'bottom_residual_absmax_1e9kgs': float(abs(field.sel(plev=100000)).max())}
    assert record['NH500_max_1e9kgs'] > 0 and record['SH500_min_1e9kgs'] < 0
    if period == 'early':
        assert 30 < record['tropical_max_1e9kgs'] < 300 and 30 < -record['tropical_min_1e9kgs'] < 300
    strength_rows.append(record)
climate = precipitation.weighted(day_weight).mean('year_in_window', skipna=False)
group_climate = xr.concat([climate.sel(model=members_in_group).mean('model', skipna=False).expand_dims(product=[group])
                           for group, members_in_group in groups.items()], 'product')
changes = group_climate.sel(period='later', drop=True) - group_climate.sel(period='early', drop=True)
contrast = (changes.sel(product='not_recovered', drop=True) - changes.sel(product='recovered', drop=True)).expand_dims(product=['contrast'])
changes = xr.concat([changes, contrast], 'product')
peak_rows = []
for group in groups:
    for period in periods:
        field = group_climate.sel(product=group, period=period, lat=slice(-20, 20))
        indices = field.values.argmax(axis=0)
        for index, lon in enumerate(longitude):
            peak_rows.append({'product': group, 'period': period, 'season': 'Annual', 'lon': lon,
                              'peak_lat': float(field.lat.values[indices[index]]),
                              'peak_pr_mm_day': float(field.values[indices[index], index])})
destination.mkdir(parents=True)
xr.Dataset({'pr': precipitation, 'day_weight': day_weight}).to_netcdf(destination / 'annual_precipitation_samples.nc')
xr.Dataset({'psi': circulation, 'day_weight': day_weight.sel(model=circulation.model)}).drop_encoding().to_netcdf(destination / 'available_circulation_samples.nc')
group_climate.rename('pr').to_netcdf(destination / 'group_precipitation_climatology.nc')
changes.rename('pr').to_netcdf(destination / 'group_precipitation_changes.nc')
pd.DataFrame(peak_rows).to_csv(destination / 'rainband_peak_gridpoints.csv', index=False)
pd.DataFrame(inventory).to_csv(destination / 'new_model_source_inventory.csv', index=False)
pd.DataFrame(strength_rows).to_csv(destination / 'uk_strength_audit.csv', index=False)
pd.DataFrame(independent_rows).to_csv(destination / 'uk_independent_300hPa_check.csv', index=False)
metadata = {'status': 'precipitation_six_models_complete_circulation_CAS_pending', 'models': models,
            'groups': groups, 'starts': starts, 'periods': periods, 'season': 'Annual',
            'missing_circulation_models': ['CAS-ESM2-0'],
            'CAS_required_inputs': {'va': ['1pctCO2 years200-209', '1pctCO2-cdr years470-479'],
                                    'ps': ['1pctCO2 years200-209', '1pctCO2-cdr years470-479']},
            'raw_dimensions': {'new_pr': 'time,lat,lon;kg m-2 s-1', 'UK_va': 'time,plev,lat145,lon192;m s-1',
                               'UK_ps': 'time,lat144,lon192;Pa'},
            'UK_staggered_grid': 'Trim exact-pole v rows; periodic bilinear interpolation of ps from half-offset lat/lon to internal v grid. Common88S-88N lies strictly inside retained native wind support; no pole extrapolation.',
            'psi_method': 'Same n35:2pi*a*cos(phi)/g integral of v from10hPa to p, terrain-clipped trapezoid, lowest valid v held to surface; divide1e9; no column mass correction.',
            'pr_method': 'CF monthly day weighting; annual native fields; periodic bilinear2degree mapping; no drift correction.',
            'grouping': 'Three fixed models each; no partial-membership circulation group means computed while CAS is absent.',
            'source_outputs': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in
                               [base / 'n71_circulation_precipitation_significance/annual_precipitation_samples.nc'] +
                               [base / 'n35_global_meridional_samples' / f'{model}_samples.nc' for model in
                                ['GFDL-ESM4', 'NorESM2-LM', 'CESM2', 'CNRM-ESM2-1']]]}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print(pd.DataFrame(strength_rows).to_string(index=False), flush=True)
print(destination.resolve(), flush=True)
