# %% 1. 导入依赖

from pathlib import Path
import json
import numpy as np
import pandas as pd
import xarray as xr


# %% 2. 参数与已有来源

root = Path('F:/sea/Antic sea ice constraint')
base = root / 'output/s05_atmosphere/p05_amoc_recovery_composites'
output = base / 'n32_nh_meridional_circulation'
assert not output.exists(), '正式结果存在，先归档。'
models = ['GFDL-ESM4', 'NorESM2-LM', 'CESM2', 'CNRM-ESM2-1']
members = ['r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f2']
starts = [[1, 271], [1, 271], [1, 131], [1850, 2120]]
periods = ['early', 'later']
pressure_levels = np.array([1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000,
                            40000, 50000, 60000, 70000, 85000, 92500, 100000], dtype=float)
latitude = np.arange(-30., 90., 2.)
radius = 6371000.
gravity = 9.80665
records = json.loads((base / 'n13_endpoint10_moisture_budget/metadata.json').read_text(encoding='utf-8'))['files']
with xr.open_dataset(base / 'n28_nh_endpoint_samples/GFDL-ESM4.nc') as reference:
    map_latitude = reference.lat.values.copy()
    map_longitude = reference.lon.values.copy()
audit = []
qc = []
model_samples = []
output.mkdir(parents=True)


# %% 3. 逐模式、逐端点重读全经度ua与va，恢复真实压力坐标

for model_index, model in enumerate(models):
    phase_fields = []
    for period_index, period in enumerate(periods):
        start = starts[model_index][period_index]
        expected_ids = np.array([year * 12 + month for year in range(start, start + 10) for month in [1, 2, 12]])
        fields = {}
        if model != 'NorESM2-LM':
            ps_parts = []
            for filename in sorted(set(r['file'] for r in records if r['model'] == model and r['period'] == period and r['variable'] == 'ps')):
                with xr.open_dataset(filename, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as source:
                    take = source.time.dt.year.isin(range(start, start + 10)) & source.time.dt.month.isin([1, 2, 12])
                    if not bool(take.any()):
                        continue
                    assert source.ps.dims == ('time', 'lat', 'lon') and source.ps.attrs['units'] == 'Pa'
                    native_ps = source.ps.isel(time=take).sel(lat=slice(-35, 90)).load()
                    ids = native_ps.time.dt.year.values * 12 + native_ps.time.dt.month.values
                    ps_parts.append(native_ps.assign_coords(time=ids).reset_coords(drop=True))
                    audit.append({'model': model, 'period': period, 'variable': 'ps', 'file': filename, 'dimensions': list(source.ps.dims), 'months': int(take.sum())})
            surface_pressure = xr.concat(ps_parts, 'time').sortby('time')
            assert np.array_equal(surface_pressure.time.values, expected_ids)

        for variable in ['ua', 'va']:
            parts = []
            embedded_ps_parts = []
            for filename in sorted(set(r['file'] for r in records if r['model'] == model and r['period'] == period and r['variable'] == variable)):
                with xr.open_dataset(filename, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as source:
                    take = source.time.dt.year.isin(range(start, start + 10)) & source.time.dt.month.isin([1, 2, 12])
                    if not bool(take.any()):
                        continue
                    assert source.attrs['source_id'] == model and source.attrs['variant_label'] == members[model_index]
                    assert source[variable].attrs['units'] == 'm s-1'
                    vertical = 'lev' if model == 'NorESM2-LM' else 'plev'
                    assert source[variable].dims == ('time', vertical, 'lat', 'lon')
                    native = source[variable].isel(time=take).sel(lat=slice(-35, 90)).load().astype(float)
                    if model == 'NorESM2-LM':
                        native_ps = source.ps.isel(time=take).sel(lat=native.lat).load()
                        offset = source.a.values * float(source.p0)
                        fraction = source.b.values
                        converted = np.full((native.sizes['time'], len(pressure_levels), native.sizes['lat'], native.sizes['lon']), np.nan)
                        for month_index in range(native.sizes['time']):
                            pressure = offset[:, None, None] + fraction[:, None, None] * native_ps.values[month_index][None]
                            wind = native.values[month_index]
                            if np.nanmean(pressure[0]) > np.nanmean(pressure[-1]):
                                pressure = pressure[::-1]
                                wind = wind[::-1]
                            assert np.all(np.diff(pressure, axis=0) > 0)
                            for level_index, target in enumerate(pressure_levels):
                                lower = np.clip(np.sum(pressure <= target, axis=0) - 1, 0, len(offset) - 2)
                                p0 = np.take_along_axis(pressure, lower[None], axis=0)[0]
                                p1 = np.take_along_axis(pressure, (lower + 1)[None], axis=0)[0]
                                w0 = np.take_along_axis(wind, lower[None], axis=0)[0]
                                w1 = np.take_along_axis(wind, (lower + 1)[None], axis=0)[0]
                                valid = (target >= pressure[0]) & (target <= pressure[-1]) & (target <= native_ps.values[month_index])
                                converted[month_index, level_index] = np.where(valid, w0 + (w1 - w0) * (np.log(target) - np.log(p0)) / (np.log(p1) - np.log(p0)), np.nan)
                        field = xr.DataArray(converted, dims=('time', 'plev', 'lat', 'lon'), coords={'time': native.time, 'plev': pressure_levels, 'lat': native.lat, 'lon': native.lon})
                        ids = native.time.dt.year.values * 12 + native.time.dt.month.values
                        embedded_ps_parts.append(native_ps.assign_coords(time=ids).reset_coords(drop=True))
                    else:
                        field = native.assign_coords(plev=np.rint(native.plev)).sel(plev=pressure_levels)
                    parts.append(field.reset_coords(drop=True))
                    audit.append({'model': model, 'period': period, 'variable': variable, 'file': filename, 'dimensions': list(source[variable].dims), 'months': int(take.sum())})
            monthly = xr.concat(parts, 'time').sortby('time')
            ids = monthly.time.dt.year.values * 12 + monthly.time.dt.month.values
            assert np.array_equal(ids, expected_ids)
            month_days = monthly.time.dt.days_in_month.values.astype(float)
            if model == 'NorESM2-LM':
                surface_pressure = xr.concat(embedded_ps_parts, 'time').sortby('time')
            ps = surface_pressure.sel(time=ids).interp(lat=monthly.lat, lon=monthly.lon).values
            assert np.isfinite(ps).all()
            monthly = monthly.where(xr.DataArray(pressure_levels, dims='plev', coords={'plev': pressure_levels}) <= xr.DataArray(ps, dims=('time', 'lat', 'lon'), coords={'time': monthly.time, 'lat': monthly.lat, 'lon': monthly.lon}))
            wind = monthly.values
            valid = np.isfinite(wind)
            assert valid[:, 0].all()
            assert not np.any((~valid[:, :-1]) & valid[:, 1:]), '压力坐标风场出现内部缺层。'
            zonal_coverage = valid.mean(axis=3)
            zonal_wind = monthly.mean('lon', skipna=True).where(xr.DataArray(zonal_coverage >= .8, dims=('time','plev','lat'), coords={'time': monthly.time, 'plev': pressure_levels, 'lat': monthly.lat}))

            # 全球纬向平均才用于标准经圈流函数；不把NA区域平均冒充Hadley环流。
            if variable == 'va':
                filled = wind.copy()
                for level_index in range(1, len(pressure_levels)):
                    filled[:, level_index] = np.where(valid[:, level_index], filled[:, level_index], filled[:, level_index - 1])
                integral = np.zeros((wind.shape[0], wind.shape[2], wind.shape[3]))
                profiles = [integral.mean(axis=2).copy()]
                for level_index in range(len(pressure_levels) - 1):
                    width = pressure_levels[level_index + 1] - pressure_levels[level_index]
                    thickness = np.clip(ps - pressure_levels[level_index], 0., width)
                    lower_wind = filled[:, level_index] + (filled[:, level_index + 1] - filled[:, level_index]) * thickness / width
                    integral += .5 * (filled[:, level_index] + lower_wind) * thickness
                    profiles.append(integral.mean(axis=2).copy())
                streamfunction = np.stack(profiles, axis=1) * (2 * np.pi * radius * np.cos(np.deg2rad(monthly.lat.values)) / gravity)[None, None] / 1.e9
                field = xr.DataArray(streamfunction, dims=('time','plev','lat'), coords={'time': monthly.time, 'plev': pressure_levels, 'lat': monthly.lat}).interp(lat=latitude)
                fields['psi'] = field
                qc.append({'model': model, 'period': period, 'psi_10hPa_max_abs': float(abs(streamfunction[:, 0]).max()),
                           'psi_max_abs_1e9kgs': float(abs(streamfunction).max()), 'coverage_min': float(zonal_coverage.min())})
            else:
                fields['u'] = zonal_wind.interp(lat=latitude)
                sector_longitude = (monthly.lon + 180) % 360 - 180
                fields['u_atlantic'] = monthly.where((sector_longitude >= -60) & (sector_longitude <= 0), drop=True).mean('lon', skipna=True).interp(lat=latitude)
                map_field = monthly.sel(plev=25000, drop=True).assign_coords(lon=monthly.lon % 360).sortby('lon')
                cyclic = xr.concat([map_field.isel(lon=[-1]).assign_coords(lon=[float(map_field.lon[-1])-360]), map_field,
                                    map_field.isel(lon=[0]).assign_coords(lon=[float(map_field.lon[0])+360])], 'lon')
                fields['u250'] = cyclic.interp(lat=map_latitude, lon=map_longitude).rename(lat='map_lat')
            print(model, period, variable, 'global longitudes read', flush=True)

        monthly_result = xr.Dataset(fields)
        annual = []
        annual_days = []
        for year in range(start, start + 10):
            take = monthly_result.time.dt.year == year
            subset = monthly_result.sel(time=take)
            assert subset.sizes['time'] == 3
            annual.append(subset.weighted(subset.time.dt.days_in_month).mean('time', skipna=False))
            annual_days.append(float(subset.time.dt.days_in_month.sum()))
        phase = xr.concat(annual, pd.Index(range(1, 11), name='year_in_window'))
        phase['day_weight'] = ('year_in_window', annual_days)
        phase_fields.append(phase.expand_dims(period=[period]))
    result = xr.concat(phase_fields, 'period').expand_dims(model=[model])
    result.u.attrs['units'] = 'm s-1'
    result.psi.attrs['units'] = '1e9 kg s-1'
    result.attrs = {'source': 'Full-longitude monthly ua/va; no piControl correction', 'psi': '2pi*a*cos(lat)/g times integral of zonal v dp from10hPa down; surface-clipped trapezoid, lowest valid v held to surface, no column mass correction.',
                    'boundary': 'psi(10hPa)=0; no measured circulation above10hPa included.', 'zonal_u': 'Valid-atmosphere mean shown only with at least80percent longitude support.',
                    'atlantic_u': '60W–0 sector mean, not a meridional mass streamfunction.'}
    result.to_netcdf(output / (model + '_samples.nc'))
    model_samples.append(result)


# %% 4. 保存成员气候态、变化与固定组差

samples = xr.concat(model_samples, 'model')
climate = samples.drop_vars('day_weight').weighted(samples.day_weight).mean('year_in_window', skipna=False)
change = climate.sel(period='later', drop=True) - climate.sel(period='early', drop=True)
groups = xr.concat([change.sel(model=models[:2]).mean('model', skipna=False).expand_dims(product=['recovered']),
                    change.sel(model=models[2:]).mean('model', skipna=False).expand_dims(product=['not_recovered'])], 'product')
contrast = (groups.sel(product='not_recovered', drop=True) - groups.sel(product='recovered', drop=True)).expand_dims(product=['contrast'])
group_climate = xr.concat([climate.sel(model=models[:2]).mean('model', skipna=False).expand_dims(product=['recovered']),
                           climate.sel(model=models[2:]).mean('model', skipna=False).expand_dims(product=['not_recovered'])], 'product')
climate.to_netcdf(output / 'member_climatology.nc')
change.to_netcdf(output / 'member_changes.nc')
xr.concat([groups, contrast], 'product').to_netcdf(output / 'group_changes.nc')
group_climate.to_netcdf(output / 'group_climatology.nc')
pd.DataFrame(audit).to_csv(output / 'input_audit.csv', index=False)
pd.DataFrame(qc).to_csv(output / 'integration_qc.csv', index=False)
(output / 'metadata.json').write_text(json.dumps({'models': models, 'members': members, 'starts': starts, 'period': '10 calendar-year DJFs, Jan/Feb/Dec, dayweighted',
    'pressure_levels_Pa': pressure_levels.tolist(), 'source_dimensions': {'pressure_models': ['time','plev','lat','lon'], 'NorESM': ['time','lev','lat','lon'], 'ps': ['time','lat','lon']},
    'output_dimensions': {'sections': ['model','period','year_in_window','plev','lat'], 'maps': ['model','period','year_in_window','map_lat','lon']},
    'streamfunction': result.attrs['psi'], 'boundary': result.attrs['boundary'], 'interpretation': 'Zonal mean circulation diagnostic, not proof that Hadley changes cause Atlantic eddy-driven jet or continental precipitation.'}, ensure_ascii=False, indent=2), encoding='utf-8')
print(output, flush=True)
