# %% 1. 导入依赖

from pathlib import Path
import json
import numpy as np
import pandas as pd
import xarray as xr


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
output = root / 'output/s05_atmosphere/p02_ua/n23_eight_model_pressure_thermal'
assert not output.exists()
models = ['ACCESS-ESM1-5', 'CanESM5', 'CESM2', 'CNRM-ESM2-1', 'GFDL-ESM4', 'MIROC-ES2L', 'NorESM2-LM', 'UKESM1-0-LL']
new_models = ['MIROC-ES2L', 'NorESM2-LM']
starts = {'MIROC-ES2L': (1850, 2119), 'NorESM2-LM': (1, 271)}
members = {'MIROC-ES2L': 'r1i1p1f2', 'NorESM2-LM': 'r1i1p1f1'}
periods = ['early', 'later']
pressure = [85000, 70000, 50000]
radius = 6371000.
gas_constant = 287.05
rotation_rate = 7.292115e-5
with xr.open_dataset(root / 'output/s05_atmosphere/p02_ua/n14_matched_temperature_wind/matched_temperature_wind.nc') as source:
    existing = source.load()
with xr.open_dataset(root / 'output/s05_atmosphere/p02_ua/n12_endpoint10_temperature_gradient_jet/temperature_gradient_jet.nc') as source:
    main = source.load()
latitude = main.lat.values
longitude = main.lon.values
qc = []
products = {}


# %% 3. 读取新补两模式相同窗口温度与风

for variable in ['ta', 'ua']:
    model_fields = []
    for model in new_models:
        period_fields = []
        for period_index, period in enumerate(periods):
            experiment = '1pctCO2' if period == 'early' else '1pctCO2-cdr'
            folder = Path(f'F:/sea/cmip6{variable}' + ('ad' if period == 'early' else ''))
            table = 'AERmon' if variable == 'ua' and model == 'NorESM2-LM' else 'Amon'
            paths = sorted(folder.rglob(f'{variable}_{table}_{model}_{experiment}_{members[model]}_*.nc'))
            paths = [path for path in paths if '(' not in path.name or not path.with_name(path.name.split(' (')[0] + '.nc').exists()]
            start = starts[model][period_index]
            parts = []
            for path in paths:
                with xr.open_dataset(path, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as source:
                    selected = (source.time.dt.year >= start) & (source.time.dt.year < start + 10) & source.time.dt.month.isin([1, 2, 12])
                    if int(selected.sum()) == 0:
                        continue
                    assert source.attrs['source_id'] == model
                    assert source.attrs['experiment_id'] == experiment
                    assert source.attrs['variant_label'] == members[model]
                    if table == 'Amon':
                        assert source[variable].dims == ('time', 'plev', 'lat', 'lon')
                        assert source.plev.attrs['units'] == 'Pa'
                        assert source[variable].attrs['units'] == ('K' if variable == 'ta' else 'm s-1')
                        field = source[variable].isel(time=selected).assign_coords(plev=np.rint(source.plev)).sel(plev=pressure, lat=slice(15, 90)).load().reset_coords(drop=True)
                    else:
                        # Nor原生混合坐标，逐月由配套ps计算真实压力，不向地下外推。
                        assert source.ua.dims == ('time', 'lev', 'lat', 'lon')
                        assert source.ps.dims == ('time', 'lat', 'lon') and source.ps.attrs['units'] == 'Pa'
                        native = source.ua.isel(time=selected).sel(lat=slice(15, 90)).load()
                        ps = source.ps.isel(time=selected).sel(lat=slice(15, 90)).load()
                        offset = source.a.values * float(source.p0)
                        fraction = source.b.values
                        values = np.full((native.sizes['time'], 3, native.sizes['lat'], native.sizes['lon']), np.nan)
                        for month_index in range(native.sizes['time']):
                            monthly_u = native.isel(time=month_index).values.astype(float)
                            monthly_p = offset[:, None, None] + fraction[:, None, None] * ps.isel(time=month_index).values[None]
                            if np.nanmean(monthly_p[0]) > np.nanmean(monthly_p[-1]):
                                monthly_p = monthly_p[::-1]
                                monthly_u = monthly_u[::-1]
                            assert np.all(np.diff(monthly_p, axis=0) > 0)
                            for level_index, target in enumerate(pressure):
                                lower = np.clip(np.sum(monthly_p <= target, axis=0) - 1, 0, len(offset) - 2)
                                low_p = np.take_along_axis(monthly_p, lower[None], axis=0)[0]
                                high_p = np.take_along_axis(monthly_p, (lower + 1)[None], axis=0)[0]
                                low_u = np.take_along_axis(monthly_u, lower[None], axis=0)[0]
                                high_u = np.take_along_axis(monthly_u, (lower + 1)[None], axis=0)[0]
                                valid = (target >= monthly_p[0]) & (target <= monthly_p[-1]) & np.isfinite(low_u) & np.isfinite(high_u)
                                weight = (np.log(target) - np.log(low_p)) / (np.log(high_p) - np.log(low_p))
                                values[month_index, level_index] = np.where(valid, low_u + weight * (high_u - low_u), np.nan)
                        field = xr.DataArray(values, dims=('time', 'plev', 'lat', 'lon'), coords={'time': native.time, 'plev': pressure, 'lat': native.lat, 'lon': native.lon})
                    parts.append(field)
                    qc.append({'file': str(path), 'model': model, 'variable': variable, 'period': period, 'dimensions': list(source[variable].dims), 'member': source.attrs['variant_label'], 'first_year': start, 'last_year': start + 9, 'selected_months': int(selected.sum())})
            monthly = xr.concat(parts, dim='time').sortby('time')
            ids = monthly.time.dt.year.values * 12 + monthly.time.dt.month.values
            expected = np.array([year * 12 + month for year in range(start, start + 10) for month in [1, 2, 12]])
            assert np.array_equal(np.sort(ids), np.sort(expected)) and len(np.unique(ids)) == 30
            climate = monthly.weighted(monthly.time.dt.days_in_month).mean('time').where(monthly.notnull().all('time'))
            climate = climate.assign_coords(lon=climate.lon % 360).sortby('lon')
            left = climate.isel(lon=[-1]).assign_coords(lon=[float(climate.lon[-1]) - 360])
            right = climate.isel(lon=[0]).assign_coords(lon=[float(climate.lon[0]) + 360])
            climate = xr.concat([left, climate, right], dim='lon').interp(lat=latitude, lon=longitude)
            period_fields.append(climate.expand_dims(period=[period]))
            print(variable, model, period, len(ids), flush=True)
        model_fields.append(xr.concat(period_fields, dim='period').expand_dims(model=[model]))
    added = xr.concat(model_fields, dim='model').transpose('model', 'period', 'plev', 'lat', 'lon')
    added = added.where(added.notnull().all('period'))
    products[variable] = xr.concat([existing[variable], added], dim='model').sel(model=models).transpose('model', 'period', 'plev', 'lat', 'lon')


# %% 4. 同层梯度与垂直切变计算

temperature = products['ta']
wind = products['ua']
gradient = -temperature.differentiate('lat') / (np.pi * radius / 180.) * 1.e6
joint850 = gradient.sel(plev=85000).notnull().all('period') & main.u850.notnull().all('period')
gradient850 = gradient.sel(plev=85000, drop=True).where(joint850)
integral = gradient.sel(plev=[50000, 70000, 85000]).assign_coords(plev=np.log([50000., 70000., 85000.])).integrate('plev')
coriolis = 2. * rotation_rate * np.sin(np.deg2rad(temperature.lat))
thermal = (gas_constant / coriolis * integral / 1.e6).transpose('model', 'period', 'lat', 'lon')
actual = wind.sel(plev=50000, drop=True) - wind.sel(plev=85000, drop=True)
joint = thermal.notnull().all('period') & actual.notnull().all('period')
delta_thermal = (thermal.sel(period='later') - thermal.sel(period='early')).where(joint)
delta_actual = (actual.sel(period='later') - actual.sel(period='early')).where(joint)
statistics = []
for region in ['North Atlantic-Europe', 'North Pacific']:
    mask = ((wind.lon >= 300) | (wind.lon <= 60)) if region == 'North Atlantic-Europe' else ((wind.lon >= 140) & (wind.lon <= 240))
    for model in models:
        predicted = delta_thermal.sel(model=model).where(mask, drop=True).sel(lat=slice(30, 70))
        observed = delta_actual.sel(model=model).where(mask, drop=True).sel(lat=slice(30, 70))
        valid = np.isfinite(predicted.values) & np.isfinite(observed.values)
        weights = np.broadcast_to(np.cos(np.deg2rad(predicted.lat.values))[:, None], predicted.shape)[valid]
        weights = weights / weights.sum()
        x = predicted.values[valid]
        y = observed.values[valid]
        centered_x = x - np.sum(weights * x)
        centered_y = y - np.sum(weights * y)
        correlation = np.sum(weights * centered_x * centered_y) / np.sqrt(np.sum(weights * centered_x**2) * np.sum(weights * centered_y**2))
        statistics.append({'model': model, 'region': region, 'area_weighted_spatial_r': correlation, 'rms_prediction_error_m_s': np.sqrt(np.sum(weights * (x - y)**2)), 'valid_cells': int(valid.sum())})


# %% 5. 独立数组检查与保存

manual_gradient = -np.gradient(temperature.values, latitude, axis=3) / (np.pi * radius / 180.) * 1.e6
gradient_error = float(np.nanmax(abs(manual_gradient - gradient.values)))
numpy_integral = np.trapezoid(gradient.sel(plev=[50000, 70000, 85000]).values, x=np.log([50000., 70000., 85000.]), axis=2)
numpy_thermal = gas_constant / coriolis.values[None, None, :, None] * numpy_integral / 1.e6
thermal_error = float(np.nanmax(abs(numpy_thermal[:, 1] - numpy_thermal[:, 0] - delta_thermal.values)))
wind_error = float(abs(wind.sel(plev=85000, drop=True) - main.u850).max())
assert gradient_error < 1.e-10 and thermal_error < 1.e-10 and wind_error < 1.e-8
assert bool(gradient850.notnull().any(['period', 'lat', 'lon']).all('model'))
output.mkdir(parents=True)
result = xr.Dataset({'ta': temperature, 'ua': wind, 'gradient': gradient, 'ta850': temperature.sel(plev=85000, drop=True), 'u850': main.u850, 'delta_gradient850': gradient850.sel(period='later') - gradient850.sel(period='early'), 'delta_u850': main.delta_u850, 'delta_thermal_shear': delta_thermal, 'delta_actual_shear': delta_actual})
result.attrs = {'periods': 'early1pctCO2 years1–10; laterCDR131–140 (continuous271–280); DJF Jan/Feb/Dec30months each,day-weighted', 'regions': 'NA60W–60E; NP140E–120W; jet squared-wind centroid30–70N unchanged', 'gradient_units': 'K per1000km; -dT/dy', 'thermal_wind': 'Dry approximation Rd/f integral(-dT/dy)dlnp from500 to850hPa,500/700/850hPa three-level trapezoid; compares vertical shear not jet latitude', 'mask': 'Two-period common valid values; no belowground filling; thermal diagnostic additionally common temperature/shear mask', 'missing_temperature_models': 'none', 'limitations': 'Spatial shear correlation is not cross-model jet-shift correlation or causal attribution; no significance inferred', 'drift_correction': 'none'}
result.to_netcdf(output / 'pressure_matched_patterns.nc')
pd.DataFrame(statistics).to_csv(output / 'thermal_shear_diagnostics.csv', index=False)
pd.read_csv(root / 'output/s05_atmosphere/p02_ua/n12_endpoint10_temperature_gradient_jet/jet_later10_early10.csv').to_csv(output / 'jet_later10_early10.csv', index=False)
(output / 'qc.json').write_text(json.dumps({'new_files': qc, 'six_model_source': 'output/s05_atmosphere/p02_ua/n14_matched_temperature_wind', 'dimensions': {'ta': list(temperature.dims), 'ua': list(wind.dims), 'delta_thermal_shear': list(delta_thermal.dims)}, 'numpy_gradient_max_error': gradient_error, 'numpy_thermal_max_error_m_s': thermal_error, 'main850_wind_max_error_m_s': wind_error, 'method': result.attrs}, indent=2), encoding='utf-8')
print(pd.DataFrame(statistics).to_string(index=False), flush=True)
