# %% 1. 导入依赖

from pathlib import Path
import json
import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import pearsonr


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
source_directory = root / 'output/s05_atmosphere/p02_ua/n12_endpoint10_temperature_gradient_jet'
output = root / 'output/s05_atmosphere/p02_ua/n14_matched_temperature_wind'
assert not output.exists()
models = ['ACCESS-ESM1-5', 'CanESM5', 'CESM2', 'CNRM-ESM2-1', 'GFDL-ESM4', 'UKESM1-0-LL']
regions = ['North Atlantic-Europe', 'North Pacific']
periods = ['early', 'later']
starts = {'ACCESS-ESM1-5': (101, 371), 'CanESM5': (1850, 2121), 'CESM2': (1, 131), 'CNRM-ESM2-1': (1850, 2120), 'GFDL-ESM4': (1, 271), 'UKESM1-0-LL': (1850, 2120)}
pressure = [85000, 70000, 50000]
radius = 6371000.
gas_constant = 287.05
rotation_rate = 7.292115e-5
with xr.open_dataset(source_directory / 'temperature_gradient_jet.nc') as dataset:
    previous = dataset.load()
latitude = previous.lat.values
longitude = previous.lon.values
jets = pd.read_csv(source_directory / 'jet_later10_early10.csv')
qc = []
products = {}


# %% 3. 读取同窗口等压面温度和风

for variable in ['ta', 'ua']:
    model_fields = []
    for model in models:
        period_fields = []
        for period_index, period in enumerate(periods):
            folder = Path(f'F:/sea/cmip6{variable}' + ('ad' if period == 'early' else ''))
            experiment = '1pctCO2' if period == 'early' else '1pctCO2-cdr'
            start = starts[model][period_index]
            paths = sorted(folder.rglob(f'{variable}_Amon_{model}_{experiment}_*.nc'))
            paths = [path for path in paths if '(' not in path.name or not path.with_name(path.name.split(' (')[0] + '.nc').exists()]
            parts = []
            for path in paths:
                with xr.open_dataset(path, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as dataset:
                    assert dataset[variable].dims == ('time', 'plev', 'lat', 'lon')
                    assert dataset.plev.attrs['units'] == 'Pa'
                    selected = (dataset.time.dt.year >= start) & (dataset.time.dt.year < start + 10) & dataset.time.dt.month.isin([1, 2, 12])
                    if int(selected.sum()) == 0:
                        continue
                    field = dataset[variable].isel(time=selected).assign_coords(plev=np.rint(dataset.plev)).sel(plev=pressure, lat=slice(15, 90)).load().reset_coords(drop=True)
                    parts.append(field)
                    qc.append({'file': str(path), 'model': model, 'variable': variable, 'period': period, 'dimensions': list(dataset[variable].dims), 'months': int(selected.sum()), 'first_year': start, 'last_year': start + 9})
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
            print(variable, model, period, flush=True)
        model_fields.append(xr.concat(period_fields, dim='period').expand_dims(model=[model]))
    data = xr.concat(model_fields, dim='model').transpose('model', 'period', 'plev', 'lat', 'lon')
    products[variable] = data.where(data.notnull().all('period'))


# %% 4. 检查同层梯度和热成风近似

temperature = products['ta']
wind = products['ua']
gradient = -temperature.differentiate('lat') * 1000. / (np.pi * radius / 180.) * 1000.
delta_gradient = gradient.sel(period='later') - gradient.sel(period='early')
delta_wind = wind.sel(period='later') - wind.sel(period='early')
# 热成风约束垂直风切变，不直接决定850hPa急流纬度；用实际T而非虚温作干空气近似。
layer_gradient = gradient.sel(plev=[50000, 70000, 85000]).assign_coords(plev=np.log([50000., 70000., 85000.])).integrate('plev') / np.log(85000. / 50000.)
coriolis = 2. * rotation_rate * np.sin(np.deg2rad(gradient.lat))
thermal_shear = gas_constant / coriolis * layer_gradient / 1.e6 * np.log(85000. / 50000.)
actual_shear = wind.sel(plev=50000, drop=True) - wind.sel(plev=85000, drop=True)
joint = thermal_shear.notnull().all('period') & actual_shear.notnull().all('period')
thermal_shear = thermal_shear.where(joint)
actual_shear = actual_shear.where(joint)
delta_thermal = thermal_shear.sel(period='later') - thermal_shear.sel(period='early')
delta_actual = actual_shear.sel(period='later') - actual_shear.sel(period='early')
flanks = []
thermal_statistics = []
for region in regions:
    longitude_mask = ((gradient.lon >= 300) | (gradient.lon <= 60)) if region == regions[0] else ((gradient.lon >= 140) & (gradient.lon <= 240))
    for model in models:
        jet_row = jets[(jets.model == model) & (jets.region == region)].iloc[0]
        for label, field in [('surface', previous.delta_surface_gradient), ('850hPa', delta_gradient.sel(plev=85000, drop=True))]:
            profile = field.sel(model=model).where(longitude_mask, drop=True).mean('lon')
            for width in [5, 10, 15]:
                south = float(profile.sel(lat=slice(jet_row.early_jet_latitude - width, jet_row.early_jet_latitude)).mean('lat'))
                north = float(profile.sel(lat=slice(jet_row.early_jet_latitude, jet_row.early_jet_latitude + width)).mean('lat'))
                flanks.append({'model': model, 'region': region, 'temperature_level': label, 'width_degree': width, 'contrast': north - south, 'jet_shift_degree': jet_row.jet_shift_degree})
        predicted = delta_thermal.sel(model=model).where(longitude_mask, drop=True).sel(lat=slice(30, 70))
        actual = delta_actual.sel(model=model).where(longitude_mask, drop=True).sel(lat=slice(30, 70))
        # 面积权重空间相关，使用实际剪切和预测剪切的共同格点；不是跨模式急流位移相关。
        valid = np.isfinite(predicted.values) & np.isfinite(actual.values)
        weight = np.broadcast_to(np.cos(np.deg2rad(actual.lat.values))[:, None], actual.shape)[valid]
        x = predicted.values[valid]
        y = actual.values[valid]
        weight = weight / weight.sum()
        x_mean = np.sum(weight * x)
        y_mean = np.sum(weight * y)
        spatial_r = np.sum(weight * (x - x_mean) * (y - y_mean)) / np.sqrt(np.sum(weight * (x - x_mean)**2) * np.sum(weight * (y - y_mean)**2))
        thermal_statistics.append({'model': model, 'region': region, 'area_weighted_spatial_r': spatial_r, 'rms_prediction_error_m_s': np.sqrt(np.sum(weight * (x - y)**2)), 'actual_shear_rms_m_s': np.sqrt(np.sum(weight * y**2)), 'valid_cells': int(valid.sum())})
flank_table = pd.DataFrame(flanks)
correlations = []
for region in regions:
    for label in ['surface', '850hPa']:
        for width in [5, 10, 15]:
            selected = flank_table[(flank_table.region == region) & (flank_table.temperature_level == label) & (flank_table.width_degree == width)]
            correlation = pearsonr(selected.contrast, selected.jet_shift_degree)
            correlations.append({'region': region, 'temperature_level': label, 'width_degree': width, 'n_models': 6, 'pearson_r': float(correlation.statistic), 'p_descriptive': float(correlation.pvalue)})


# %% 5. 分离风廓线的近似平移与强度变化

# 对全部八模式一视同仁，固定30–70N；线性一阶形状拟合不是另选主急流指标。
fit_rows = []
fit_fields = []
for region in regions:
    for model in previous.model.values:
        profile = previous.regional_u850.sel(region=region, model=model, lat=slice(30, 70))
        early = profile.sel(period='early').values
        change = profile.sel(period='later').values - early
        derivative = np.gradient(early, profile.lat.values)
        matrix = np.column_stack([early, -derivative, np.ones_like(early)])
        coefficients, residuals, rank, singular = np.linalg.lstsq(matrix, change, rcond=None)
        fitted = matrix @ coefficients
        residual = change - fitted
        fit_rows.append({'model': str(model), 'region': region, 'amplitude_fraction': coefficients[0], 'linear_translation_degree': coefficients[1], 'background_change_m_s': coefficients[2], 'relative_residual_energy': float(np.sum(residual**2) / np.sum(change**2)), 'design_condition_number': float(np.linalg.cond(matrix)), 'linear_translation_note': 'First-order fit, not exact displacement; poorly represented changes retain residual'})
        fit_fields.append(xr.Dataset({'actual_change': ('lat', change), 'fitted_change': ('lat', fitted), 'residual': ('lat', residual)}, coords={'lat': profile.lat}).expand_dims(case=[f'{region}|{model}']))


# %% 6. 保存正式结果和核对

assert float(abs(wind.sel(plev=85000, drop=True) - previous.u850.sel(model=models)).max()) < 1.e-8
assert float(np.nanmax(abs(gradient.values + np.gradient(temperature.values, latitude, axis=3) / (np.pi * radius / 180.) * 1.e6))) < 1.e-8
output.mkdir(parents=True)
result = xr.Dataset({'ta': temperature, 'ua': wind, 'gradient': gradient, 'delta_gradient': delta_gradient, 'delta_u': delta_wind, 'delta_thermal_shear': delta_thermal, 'delta_actual_shear': delta_actual})
result.attrs = {'periods': 'early1pctCO2 years1–10; laterCDR131–140; DJF calendar months1/2/12; 30months each day-weighted', 'regions': 'NA60W–60E; NP140E–120W; fixed definitions, not tuned by correlation', 'gradient_units': 'K per1000km; -dT/dy', 'thermal_wind': 'Rd/f integral[-dT/dy] dln(p) from500 to850hPa; three levels500/700/850; dry-temperature approximation; compares actual zonal wind shear, not geostrophic wind', 'limitations': 'No humidity correction, ageostrophic or numerical closure; cannot infer jet-latitude control from shear pattern correlation', 'missing_models': 'MIROC-ES2L and NorESM2-LM absent from local paired ta inventory, not excluded by result', 'drift_correction': 'none'}
result.to_netcdf(output / 'matched_temperature_wind.nc')
xr.concat(fit_fields, dim='case').to_netcdf(output / 'wind_shape_fits.nc')
flank_table.to_csv(output / 'flank_metrics.csv', index=False)
pd.DataFrame(correlations).to_csv(output / 'same_six_model_correlations.csv', index=False)
pd.DataFrame(thermal_statistics).to_csv(output / 'thermal_shear_diagnostics.csv', index=False)
pd.DataFrame(fit_rows).to_csv(output / 'wind_shape_diagnostics.csv', index=False)
(output / 'input_qc.json').write_text(json.dumps({'files': qc, 'method': result.attrs}, indent=2), encoding='utf-8')
print(pd.DataFrame(correlations).to_string(index=False), flush=True)
print(pd.DataFrame(thermal_statistics).to_string(index=False), flush=True)
print(pd.DataFrame(fit_rows).to_string(index=False), flush=True)
