# %% 1. 导入依赖

from pathlib import Path
import json
import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import pearsonr, spearmanr


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
output = root / 'output/s05_atmosphere/p02_ua/n12_endpoint10_temperature_gradient_jet'
assert not output.exists()
models = ['ACCESS-ESM1-5', 'CanESM5', 'CESM2', 'CNRM-ESM2-1', 'GFDL-ESM4', 'MIROC-ES2L', 'NorESM2-LM', 'UKESM1-0-LL']
regions = ['North Atlantic-Europe', 'North Pacific']
periods = ['early', 'later']
start_years = {'ACCESS-ESM1-5': (101, 371), 'CanESM5': (1850, 2121), 'CESM2': (1, 131), 'CNRM-ESM2-1': (1850, 2120), 'GFDL-ESM4': (1, 271), 'MIROC-ES2L': (1850, 2119), 'NorESM2-LM': (1, 271), 'UKESM1-0-LL': (1850, 2120)}
with xr.open_dataset(root / 'output/s05_atmosphere/p02_ua/n06_article_jet_later_early/jet_climatology.nc') as source:
    grid_latitude = source.lat.values
    grid_longitude = source.lon.values
km_per_degree = np.pi * 6371. / 180.
qc = []
products = {}


# %% 3. 读取两段各10年气候态并统一网格

for variable in ['tas', 'ua']:
    model_fields = []
    for model in models:
        period_fields = []
        for period_index, period in enumerate(periods):
            folder = Path(f'F:/sea/cmip6{variable}' + ('ad' if period == 'early' else ''))
            experiment = '1pctCO2' if period == 'early' else '1pctCO2-cdr'
            table_name = 'AERmon' if variable == 'ua' and model == 'NorESM2-LM' else 'Amon'
            start = start_years[model][period_index]
            paths = sorted(folder.rglob(f'{variable}_{table_name}_{model}_{experiment}_*.nc'))
            paths = [path for path in paths if 'yearmean' not in path.name.lower() and ('(' not in path.name or not path.with_name(path.name.split(' (')[0] + '.nc').exists())]
            monthly_parts = []
            for path in paths:
                with xr.open_dataset(path, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as source:
                    selected = (source.time.dt.year >= start) & (source.time.dt.year < start + 10) & source.time.dt.month.isin([1, 2, 12])
                    if int(selected.sum()) == 0:
                        continue
                    if variable == 'tas':
                        assert source.tas.dims == ('time', 'lat', 'lon') and source.tas.attrs['units'] == 'K'
                        field = source.tas.isel(time=selected).sel(lat=slice(15, 90)).load()
                    elif table_name == 'Amon':
                        assert source.ua.dims == ('time', 'plev', 'lat', 'lon') and source.plev.attrs['units'] == 'Pa'
                        field = source.ua.isel(time=selected).assign_coords(plev=np.rint(source.plev)).sel(plev=85000, drop=True).sel(lat=slice(15, 90)).load()
                    else:
                        assert source.ua.dims == ('time', 'lev', 'lat', 'lon') and source.ps.dims == ('time', 'lat', 'lon')
                        assert source.ps.attrs['units'] == 'Pa'
                        monthly_wind = source.ua.isel(time=selected).sel(lat=slice(15, 90)).load()
                        monthly_ps = source.ps.isel(time=selected).sel(lat=slice(15, 90)).load()
                        pressure_offset = source.a.values * float(source.p0)
                        pressure_fraction = source.b.values
                        values = np.full((monthly_wind.sizes['time'], monthly_wind.sizes['lat'], monthly_wind.sizes['lon']), np.nan)
                        for month_index in range(monthly_wind.sizes['time']):
                            native_wind = monthly_wind.isel(time=month_index).values.astype(float)
                            native_pressure = pressure_offset[:, None, None] + pressure_fraction[:, None, None] * monthly_ps.isel(time=month_index).values[None]
                            if np.nanmean(native_pressure[0]) > np.nanmean(native_pressure[-1]):
                                native_pressure = native_pressure[::-1]
                                native_wind = native_wind[::-1]
                            assert np.all(np.diff(native_pressure, axis=0) > 0)
                            lower = np.clip(np.sum(native_pressure <= 85000, axis=0) - 1, 0, len(pressure_offset) - 2)
                            low_p = np.take_along_axis(native_pressure, lower[None], axis=0)[0]
                            high_p = np.take_along_axis(native_pressure, (lower + 1)[None], axis=0)[0]
                            low_u = np.take_along_axis(native_wind, lower[None], axis=0)[0]
                            high_u = np.take_along_axis(native_wind, (lower + 1)[None], axis=0)[0]
                            valid = (85000 >= native_pressure[0]) & (85000 <= native_pressure[-1]) & np.isfinite(low_u) & np.isfinite(high_u)
                            fraction = (np.log(85000) - np.log(low_p)) / (np.log(high_p) - np.log(low_p))
                            values[month_index] = np.where(valid, low_u + fraction * (high_u - low_u), np.nan)
                        field = xr.DataArray(values, dims=('time', 'lat', 'lon'), coords={'time': monthly_wind.time, 'lat': monthly_wind.lat, 'lon': monthly_wind.lon})
                    monthly_parts.append(field.reset_coords(drop=True))
                    qc.append({'model': model, 'variable': variable, 'period': period, 'file': str(path), 'dimensions': list(source[variable].dims), 'months': int(selected.sum()), 'first_year': start, 'last_year': start + 9, 'variant': source.attrs.get('variant_label')})
            monthly = xr.concat(monthly_parts, dim='time').sortby('time')
            ids = monthly.time.dt.year.values * 12 + monthly.time.dt.month.values
            expected = np.array([year * 12 + month for year in range(start, start + 10) for month in [1, 2, 12]])
            assert np.array_equal(np.sort(ids), np.sort(expected)) and len(np.unique(ids)) == 30
            climate = monthly.weighted(monthly.time.dt.days_in_month).mean('time').where(monthly.notnull().all('time'))
            climate = climate.assign_coords(lon=climate.lon % 360).sortby('lon')
            left = climate.isel(lon=[-1]).assign_coords(lon=[float(climate.lon[-1]) - 360])
            right = climate.isel(lon=[0]).assign_coords(lon=[float(climate.lon[0]) + 360])
            climate = xr.concat([left, climate, right], dim='lon').interp(lat=grid_latitude, lon=grid_longitude)
            period_fields.append(climate.expand_dims(period=[period]))
            print(variable, model, period, len(ids), flush=True)
        model_fields.append(xr.concat(period_fields, dim='period').expand_dims(model=[model]))
    data = xr.concat(model_fields, dim='model').transpose('model', 'period', 'lat', 'lon')
    products[variable] = data.where(data.notnull().all('period'))


# %% 4. 计算急流、温度异常与梯度空间结构

wind = products['ua']
temperature = products['tas']
gradient = -temperature.differentiate('lat') * 1000. / km_per_degree
jets = []
flanks = []
profiles = []
wind_profiles = []
for region in regions:
    longitude_mask = ((wind.lon >= 300) | (wind.lon <= 60)) if region == regions[0] else ((wind.lon >= 140) & (wind.lon <= 240))
    wind_profile = wind.where(longitude_mask, drop=True).mean('lon')
    search = wind_profile.sel(lat=slice(30, 70))
    jet_latitude = (search**2 * search.lat).sum('lat') / (search**2).sum('lat')
    gradient_profile = gradient.where(longitude_mask, drop=True).mean('lon')
    profiles.append(gradient_profile.expand_dims(region=[region]))
    wind_profiles.append(wind_profile.expand_dims(region=[region]))
    for model in models:
        early_latitude = float(jet_latitude.sel(model=model, period='early'))
        later_latitude = float(jet_latitude.sel(model=model, period='later'))
        jets.append({'model': model, 'region': region, 'early_jet_latitude': early_latitude, 'later_jet_latitude': later_latitude, 'jet_shift_degree': later_latitude - early_latitude})
        delta = gradient_profile.sel(model=model, period='later') - gradient_profile.sel(model=model, period='early')
        for width in [5, 10, 15]:
            south = float(delta.sel(lat=slice(early_latitude - width, early_latitude)).mean('lat'))
            north = float(delta.sel(lat=slice(early_latitude, early_latitude + width)).mean('lat'))
            flanks.append({'model': model, 'region': region, 'flank_width_degree': width, 'jet_shift_degree': later_latitude - early_latitude, 'south_gradient_change': south, 'north_gradient_change': north, 'north_minus_south_gradient_change': north - south})
jet_table = pd.DataFrame(jets)
flank_table = pd.DataFrame(flanks)
statistics = []
for region in regions:
    for width in [5, 10, 15]:
        selected = flank_table[(flank_table.region == region) & (flank_table.flank_width_degree == width)]
        x = selected.north_minus_south_gradient_change.to_numpy()
        y = selected.jet_shift_degree.to_numpy()
        assert np.isfinite(x).all() and np.isfinite(y).all()
        leave_one = []
        for index in range(8):
            mask = np.arange(8) != index
            leave_one.append(float(pearsonr(x[mask], y[mask]).statistic))
        correlation = pearsonr(x, y)
        statistics.append({'region': region, 'width_degree': width, 'pearson_r': float(correlation.statistic), 'p_descriptive': float(correlation.pvalue), 'spearman_rho': float(spearmanr(x, y).statistic), 'leave_one_min': min(leave_one), 'leave_one_max': max(leave_one), 'same_sign_count': int(np.sum(x * y > 0))})


# %% 5. 独立公式核对与保存结果

derivative_error = float(np.nanmax(abs(-np.gradient(temperature.values, grid_latitude, axis=2) * 1000. / km_per_degree - gradient.values)))
assert derivative_error < 1.e-10
manual_errors = []
for region in regions:
    region_profile = xr.concat(wind_profiles, dim='region').sel(region=region, lat=slice(30, 70))
    manual_jet = np.nansum(region_profile.values**2 * region_profile.lat.values[None, None], axis=-1) / np.nansum(region_profile.values**2, axis=-1)
    selected = jet_table[jet_table.region == region].set_index('model').loc[models]
    manual_errors.append(float(np.max(abs(manual_jet - selected[['early_jet_latitude', 'later_jet_latitude']].values))))
assert max(manual_errors) < 1.e-10
previous = pd.read_csv(root / 'output/s05_atmosphere/p02_ua/n06_article_jet_later_early/jet_later_early.csv')
comparison = jet_table.merge(previous[['model', 'region', 'cdr_end_minus_early_shift']], on=['model', 'region'])
previous_max_difference = float(abs(comparison.jet_shift_degree - comparison.cdr_end_minus_early_shift).max())
result = xr.Dataset({'tas': temperature, 'u850': wind, 'surface_gradient': gradient, 'delta_tas': temperature.sel(period='later') - temperature.sel(period='early'), 'delta_u850': wind.sel(period='later') - wind.sel(period='early'), 'delta_surface_gradient': gradient.sel(period='later') - gradient.sel(period='early'), 'regional_surface_gradient': xr.concat(profiles, dim='region'), 'regional_u850': xr.concat(wind_profiles, dim='region')})
result.surface_gradient.attrs['units'] = 'K per 1000 km'
result.delta_surface_gradient.attrs['units'] = 'K per 1000 km'
result.attrs = {'early': '1pctCO2 years1–10', 'later': 'CDR years131–140; continuous years271–280; no stabilization months', 'season': 'DJF calendar months1/2/12 in each ten-year window, 30months each, day-weighted', 'jet_method': '850hPa squared regional-zonal-mean wind centroid over30–70N; no peak search or smoothing; period climatology first', 'regions': 'NA-Europe60W–60E; NP140E–120W', 'grid': 'CanESM5 native Gaussian latitude/longitude; cyclic bilinear interpolation', 'mask': 'Variable-specific common valid mask for two periods only; no below-ground filling', 'gradient': '-d(tas)/dy with actual nonuniform latitude spacing; K/1000km; positive change means stronger poleward cooling', 'flank_metric': 'Fixed5/10/15degree bands about early regional jet latitude; north mean minus south mean, equal latitude-grid-point weighting', 'interpretation': 'Regional horizontal jet lines are not local jet axes; gradient association is not thermal-wind closure or causal proof', 'drift_correction': 'none'}
output.mkdir(parents=True)
result.to_netcdf(output / 'temperature_gradient_jet.nc')
jet_table.to_csv(output / 'jet_later10_early10.csv', index=False)
flank_table.to_csv(output / 'jet_flank_gradient.csv', index=False)
pd.DataFrame(statistics).to_csv(output / 'flank_width_sensitivity.csv', index=False)
(output / 'input_qc.json').write_text(json.dumps({'files': qc, 'method': result.attrs, 'numpy_gradient_error': derivative_error, 'numpy_jet_error': max(manual_errors), 'previous_cdr_endpoint_shift_max_difference': previous_max_difference, 'all_temperature_finite': bool(np.isfinite(temperature).all())}, indent=2), encoding='utf-8')
print(jet_table.to_string(index=False), flush=True)
print(flank_table[flank_table.flank_width_degree == 10].to_string(index=False), flush=True)
print(pd.DataFrame(statistics).to_string(index=False), flush=True)
