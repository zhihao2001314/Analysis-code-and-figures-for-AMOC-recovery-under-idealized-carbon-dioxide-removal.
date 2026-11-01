# %% 1. 导入依赖

from pathlib import Path
import json
import numpy as np
import pandas as pd
import xarray as xr


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
output = root / 'output/s05_atmosphere/p02_ua/n27_endpoint10_jet_sections'
assert not output.exists()
base = root / 'output/s05_atmosphere/p02_ua/n12_endpoint10_temperature_gradient_jet'
inputs = json.loads((base / 'input_qc.json').read_text(encoding='utf-8'))['files']
jets = pd.read_csv(base / 'jet_later10_early10.csv')
with xr.open_dataset(base / 'temperature_gradient_jet.nc') as source:
    main = source.load()
with xr.open_dataset(root / 'output/s05_atmosphere/p02_ua/n23_eight_model_pressure_thermal/pressure_matched_patterns.nc') as source:
    three_levels = source.load()
models = list(main.model.values)
periods = ['early', 'later']
regions = ['North Atlantic-Europe', 'North Pacific']
pressure_levels = [100000, 92500, 85000, 70000, 60000, 50000, 40000, 30000, 25000, 20000, 15000, 10000]
latitude = main.lat.values
longitude = main.lon.values
radius = 6371000.
coverage_threshold = 0.8
products = {}
qc = []


# %% 3. 读取同窗口12层温度和风

for variable in ['ta', 'ua']:
    model_fields = []
    for model in models:
        period_fields = []
        for period in periods:
            reference = [record for record in inputs if record['variable'] == 'ua' and record['model'] == model and record['period'] == period][0]
            start = int(reference['first_year'])
            member = reference['variant']
            experiment = '1pctCO2' if period == 'early' else '1pctCO2-cdr'
            folder = Path(f'F:/sea/cmip6{variable}' + ('ad' if period == 'early' else ''))
            table = 'AERmon' if model == 'NorESM2-LM' and variable == 'ua' else 'Amon'
            paths = sorted(folder.rglob(f'{variable}_{table}_{model}_{experiment}_{member}_*.nc'))
            paths = [path for path in paths if 'yearmean' not in path.name.lower() and ('(' not in path.name or not path.with_name(path.name.split(' (')[0] + '.nc').exists())]
            parts = []
            for path in paths:
                with xr.open_dataset(path, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as source:
                    selected = (source.time.dt.year >= start) & (source.time.dt.year < start + 10) & source.time.dt.month.isin([1, 2, 12])
                    if int(selected.sum()) == 0:
                        continue
                    assert source.attrs['source_id'] == model and source.attrs['experiment_id'] == experiment and source.attrs['variant_label'] == member
                    if table == 'Amon':
                        assert source[variable].dims == ('time', 'plev', 'lat', 'lon') and source.plev.attrs['units'] == 'Pa'
                        field = source[variable].isel(time=selected).assign_coords(plev=np.rint(source.plev)).sel(plev=pressure_levels, lat=slice(15, 90)).load().reset_coords(drop=True)
                    else:
                        assert source.ua.dims == ('time', 'lev', 'lat', 'lon') and source.ps.dims == ('time', 'lat', 'lon')
                        native = source.ua.isel(time=selected).sel(lat=slice(15, 90)).load()
                        ps = source.ps.isel(time=selected).sel(lat=slice(15, 90)).load()
                        offset = source.a.values * float(source.p0)
                        fraction = source.b.values
                        values = np.full((native.sizes['time'], len(pressure_levels), native.sizes['lat'], native.sizes['lon']), np.nan)
                        for month_index in range(native.sizes['time']):
                            native_u = native.isel(time=month_index).values.astype(float)
                            native_p = offset[:, None, None] + fraction[:, None, None] * ps.isel(time=month_index).values[None]
                            if np.nanmean(native_p[0]) > np.nanmean(native_p[-1]):
                                native_p = native_p[::-1]
                                native_u = native_u[::-1]
                            assert np.all(np.diff(native_p, axis=0) > 0)
                            for level_index, target in enumerate(pressure_levels):
                                lower = np.clip(np.sum(native_p <= target, axis=0) - 1, 0, len(offset) - 2)
                                low_p = np.take_along_axis(native_p, lower[None], axis=0)[0]
                                high_p = np.take_along_axis(native_p, (lower + 1)[None], axis=0)[0]
                                low_u = np.take_along_axis(native_u, lower[None], axis=0)[0]
                                high_u = np.take_along_axis(native_u, (lower + 1)[None], axis=0)[0]
                                valid = (target >= native_p[0]) & (target <= native_p[-1]) & np.isfinite(low_u) & np.isfinite(high_u)
                                weight = (np.log(target) - np.log(low_p)) / (np.log(high_p) - np.log(low_p))
                                values[month_index, level_index] = np.where(valid, low_u + weight * (high_u - low_u), np.nan)
                        field = xr.DataArray(values, dims=('time', 'plev', 'lat', 'lon'), coords={'time': native.time, 'plev': pressure_levels, 'lat': native.lat, 'lon': native.lon})
                    parts.append(field)
                    qc.append({'model': model, 'variable': variable, 'period': period, 'file': str(path), 'member': member, 'native_dimensions': list(source[variable].dims), 'first_year': start, 'last_year': start + 9, 'selected_months': int(selected.sum())})
            monthly = xr.concat(parts, dim='time').sortby('time')
            ids = monthly.time.dt.year.values * 12 + monthly.time.dt.month.values
            expected = np.array([year * 12 + month for year in range(start, start + 10) for month in [1, 2, 12]])
            assert len(np.unique(ids)) == 30 and np.array_equal(np.sort(ids), np.sort(expected))
            # 剖面沿用n23全月有效口径，低层缺测不补值；与n25区域指数的可用月份平均不同。
            climate = monthly.weighted(monthly.time.dt.days_in_month).mean('time').where(monthly.notnull().all('time'))
            climate = climate.assign_coords(lon=climate.lon % 360).sortby('lon')
            left = climate.isel(lon=[-1]).assign_coords(lon=[float(climate.lon[-1]) - 360])
            right = climate.isel(lon=[0]).assign_coords(lon=[float(climate.lon[0]) + 360])
            climate = xr.concat([left, climate, right], dim='lon').interp(lat=latitude, lon=longitude)
            period_fields.append(climate.expand_dims(period=[period]))
            print(variable, model, period, flush=True)
        model_fields.append(xr.concat(period_fields, dim='period').expand_dims(model=[model]))
    merged = xr.concat(model_fields, dim='model').transpose('model', 'period', 'plev', 'lat', 'lon')
    products[variable] = merged.where(merged.notnull().all('period'))


# %% 4. 同经度支持的温度梯度剖面与剪切廓线

temperature = products['ta']
wind = products['ua']
gradient = -temperature.differentiate('lat') / (np.pi * radius / 180.) * 1.e6
sections = []
shears = []
for region in regions:
    mask = ((temperature.lon >= 300) | (temperature.lon <= 60)) if region == regions[0] else ((temperature.lon >= 140) & (temperature.lon <= 240))
    region_gradient = gradient.where(mask, drop=True)
    region_wind = wind.where(mask, drop=True)
    common = region_gradient.notnull().all('period') & region_wind.notnull().all('period')
    coverage = common.mean('lon')
    gradient_profile = region_gradient.where(common).mean('lon').where(coverage >= coverage_threshold)
    wind_profile = region_wind.where(common).mean('lon').where(coverage >= coverage_threshold)
    sections.append(xr.Dataset({'gradient_profile': gradient_profile, 'wind_profile': wind_profile, 'delta_gradient_profile': gradient_profile.sel(period='later') - gradient_profile.sel(period='early'), 'delta_wind_profile': wind_profile.sel(period='later') - wind_profile.sel(period='early'), 'valid_longitude_fraction': coverage}).expand_dims(region=[region]))
    predicted = three_levels.delta_thermal_shear.where(mask, drop=True)
    actual = three_levels.delta_actual_shear.where(mask, drop=True)
    upper = (three_levels.ua.sel(period='later', plev=50000, drop=True) - three_levels.ua.sel(period='early', plev=50000, drop=True)).where(mask, drop=True)
    lower = three_levels.delta_u850.where(mask, drop=True)
    common_shear = predicted.notnull() & actual.notnull() & upper.notnull() & lower.notnull()
    shear_coverage = common_shear.mean('lon')
    fields = {}
    for name, field in [('thermal_shear', predicted), ('actual_shear', actual), ('upper_wind', upper), ('lower_wind', lower)]:
        fields[name] = field.where(common_shear).mean('lon').where(shear_coverage >= coverage_threshold)
    fields['valid_shear_longitude_fraction'] = shear_coverage
    shears.append(xr.Dataset(fields).expand_dims(region=[region]))
section_data = xr.concat(sections, dim='region')
shear_data = xr.concat(shears, dim='region')


# %% 5. 保留主急流并描述风带形状变化

shape_rows = []
for region in regions:
    for model in models:
        profile = main.regional_u850.sel(region=region, model=model, lat=slice(30, 70))
        early = profile.sel(period='early').values
        later = profile.sel(period='later').values
        lat = profile.lat.values
        change = later - early
        matrix = np.column_stack([early, -np.gradient(early, lat), np.ones(len(lat))])
        coefficients = np.linalg.lstsq(matrix, change, rcond=None)[0]
        row = jets[(jets.region == region) & (jets.model == model)].iloc[0]
        record = {'model': model, 'region': region, 'main_jet_shift_degree': row.jet_shift_degree, 'early_centroid_degree': row.early_jet_latitude, 'later_centroid_degree': row.later_jet_latitude, 'linear_translation_degree': coefficients[1], 'amplitude_fraction': coefficients[0], 'background_change_m_s': coefficients[2], 'relative_fit_residual_energy': float(np.sum((change - matrix @ coefficients)**2) / np.sum(change**2))}
        for period, values in [('early', early), ('later', later)]:
            centroid = np.sum(values**2 * lat) / np.sum(values**2)
            record[f'{period}_peak_grid_latitude'] = float(lat[np.argmax(values)])
            record[f'{period}_peak_speed_m_s'] = float(np.max(values))
            record[f'{period}_squared_wind_width_degree'] = float(np.sqrt(np.sum(values**2 * (lat - centroid)**2) / np.sum(values**2)))
        shape_rows.append(record)


# %% 6. QC与保存

wind_error = float(abs(wind.sel(plev=[85000, 70000, 50000]) - three_levels.ua).max())
temperature_error = float(abs(temperature.sel(plev=[85000, 70000, 50000]) - three_levels.ta).max())
numpy_gradient = -np.gradient(temperature.values, latitude, axis=3) / (np.pi * radius / 180.) * 1.e6
gradient_error = float(np.nanmax(abs(gradient.values - numpy_gradient)))
identity_error = float(abs(shear_data.upper_wind - shear_data.lower_wind - shear_data.actual_shear).max())
assert wind_error < 1.e-8 and temperature_error < 1.e-8 and gradient_error < 1.e-10 and identity_error < 1.e-10
output.mkdir(parents=True)
result = xr.merge([section_data, shear_data])
result['main_u850_profile'] = main.regional_u850
result['main_delta_u850_profile'] = main.regional_u850.sel(period='later') - main.regional_u850.sel(period='early')
result.attrs = {'periods': 'CDR131–140 minus1pctCO2 1–10; DJF Jan/Feb/Dec30months each, day-weighted; no stabilization', 'regions': 'NA60W–60E; NP140E–120W', 'sections': 'Cellwise -dT/dy then regional common ta-gradient/ua longitude mean per pressure, >=80% coverage; two-period and all-month valid; no smoothing/filling', 'shear': 'Saved n23 dry three-level850–500hPa thermal-wind estimate and actual shear, same longitude support and >=80% coverage; upper_wind is actual500hPa change', 'main_jet': 'n12 regional u850 profile and squared-wind centroid unchanged; section/shear masks do not replace main masks', 'shape_fit': 'First-order Δu=a*u_early−shift*du_early/dlat+background; residual retained, not exact translation or causal allocation', 'limits': 'Discrete peak grid spacing about2.79degrees; no significance/CI inferred; thermal-wind constrains vertical shear not absolute jet latitude', 'drift_correction': 'none'}
result.to_netcdf(output / 'jet_temperature_sections.nc')
pd.DataFrame(shape_rows).to_csv(output / 'wind_shape_summary.csv', index=False)
jets.to_csv(output / 'main_jet_later10_early10.csv', index=False)
(output / 'qc.json').write_text(json.dumps({'files': qc, 'source_dimensions': {'ta': list(temperature.dims), 'ua': list(wind.dims)}, 'output_dimensions': {name: list(result[name].dims) for name in result.data_vars}, 'n23_wind_max_error_m_s': wind_error, 'n23_temperature_max_error_K': temperature_error, 'numpy_gradient_max_error': gradient_error, 'upper_minus_lower_shear_identity_error_m_s': identity_error, 'method': result.attrs}, indent=2), encoding='utf-8')
print(pd.DataFrame(shape_rows).to_string(index=False), flush=True)
