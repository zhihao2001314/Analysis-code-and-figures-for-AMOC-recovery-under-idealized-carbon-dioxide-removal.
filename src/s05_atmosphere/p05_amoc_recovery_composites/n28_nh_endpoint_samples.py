# %% 1. 导入依赖

from pathlib import Path
import json
import numpy as np
import pandas as pd
import xarray as xr


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
base = root / 'output/s05_atmosphere/p05_amoc_recovery_composites'
jet_base = root / 'output/s05_atmosphere/p02_ua'
output = base / 'n28_nh_endpoint_samples'
assert not output.exists(), '正式结果存在，先归档。'
models = ['GFDL-ESM4', 'NorESM2-LM', 'CESM2', 'CNRM-ESM2-1']
members = {'GFDL-ESM4': 'r1i1p1f1', 'NorESM2-LM': 'r1i1p1f1', 'CESM2': 'r1i1p1f1', 'CNRM-ESM2-1': 'r1i1p1f2'}
starts = {'GFDL-ESM4': [1, 271], 'NorESM2-LM': [1, 271], 'CESM2': [1, 131], 'CNRM-ESM2-1': [1850, 2120]}
periods = ['early', 'later']
provenance = json.loads((base / 'n01_endpoint_groups/qc.json').read_text(encoding='utf-8'))['inputs']
wind_provenance = json.loads((jet_base / 'n12_endpoint10_temperature_gradient_jet/input_qc.json').read_text(encoding='utf-8'))['files']
audit = []
errors = []
with xr.open_dataset(base / 'n01_endpoint_groups/endpoint_climatology.nc') as source:
    reference = source.sel(model=models).load()
grid_record = next(row for row in provenance if row['model'] == 'CNRM-ESM2-1' and row['variable'] == 'pr' and row['period'] == 'early')
with xr.open_dataset(grid_record['file'], decode_times=False) as grid_source:
    lower_latitudes = grid_source.lat.values[(grid_source.lat.values >= 0) & (grid_source.lat.values < float(reference.lat.min()))]
latitude = np.concatenate([lower_latitudes, reference.lat.values])
assert np.isin(reference.lat.values, latitude).all()
longitude = reference.lon.values
target_lon, target_lat = np.meshgrid(longitude, latitude)
target_xyz = np.column_stack([(np.cos(np.deg2rad(target_lat)) * np.cos(np.deg2rad(target_lon))).ravel(),
                             (np.cos(np.deg2rad(target_lat)) * np.sin(np.deg2rad(target_lon))).ravel(), np.sin(np.deg2rad(target_lat)).ravel()])
products = {}


# %% 3. 从已核验来源重读每窗30月，保存十个日历年冬季样本

for model in models:
    variables = {}
    for variable in ['pr', 'evspsbl', 'ua']:
        phases = []
        weights_by_phase = []
        monthly_phases = []
        coordinates = []
        records = wind_provenance if variable in ['tas', 'ua'] else provenance
        for period_index, period in enumerate(periods):
            start = starts[model][period_index]
            paths = sorted(set(row['file'] for row in records if row['model'] == model and row['variable'] == variable and row['period'] == period))
            assert paths
            parts = []
            for filename in paths:
                with xr.open_dataset(filename, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as source:
                    assert source.attrs['source_id'] == model and source.attrs['variant_label'] == members[model]
                    selected = source.time.dt.year.isin(range(start, start + 10)) & source.time.dt.month.isin([12, 1, 2])
                    assert bool(selected.any())
                    if variable == 'ua' and model == 'NorESM2-LM':
                        assert source.ua.dims == ('time', 'lev', 'lat', 'lon')
                        native = source.ua.isel(time=selected).sel(lat=slice(-5, 90)).load()
                        ps = source.ps.isel(time=selected).sel(lat=native.lat).load()
                        offset = source.a.values * float(source.p0)
                        fraction = source.b.values
                        values = np.full((native.sizes['time'], native.sizes['lat'], native.sizes['lon']), np.nan)
                        for month in range(native.sizes['time']):
                            pressure = offset[:, None, None] + fraction[:, None, None] * ps.isel(time=month).values[None]
                            wind = native.isel(time=month).values.astype(float)
                            if np.nanmean(pressure[0]) > np.nanmean(pressure[-1]):
                                pressure = pressure[::-1]
                                wind = wind[::-1]
                            assert np.all(np.diff(pressure, axis=0) > 0)
                            lower = np.clip(np.sum(pressure <= 85000, axis=0) - 1, 0, len(offset) - 2)
                            low_p = np.take_along_axis(pressure, lower[None], axis=0)[0]
                            high_p = np.take_along_axis(pressure, (lower + 1)[None], axis=0)[0]
                            low_u = np.take_along_axis(wind, lower[None], axis=0)[0]
                            high_u = np.take_along_axis(wind, (lower + 1)[None], axis=0)[0]
                            valid = (85000 >= pressure[0]) & (85000 <= pressure[-1]) & np.isfinite(low_u) & np.isfinite(high_u)
                            amount = (np.log(85000) - np.log(low_p)) / (np.log(high_p) - np.log(low_p))
                            values[month] = np.where(valid, low_u + amount * (high_u - low_u), np.nan)
                        field = xr.DataArray(values, dims=('time', 'lat', 'lon'), coords={'time': native.time, 'lat': native.lat, 'lon': native.lon})
                    elif variable == 'ua':
                        assert source.ua.dims == ('time', 'plev', 'lat', 'lon')
                        field = source.ua.isel(time=selected).assign_coords(plev=np.rint(source.plev)).sel(plev=85000, drop=True).sel(lat=slice(-5, 90)).load()
                    else:
                        field = source[variable].isel(time=selected).load()
                        if variable == 'tas':
                            field = field.sel(lat=slice(-5, 90))
                    parts.append(field.reset_coords(drop=True).astype(float))
                    if variable == 'siconc':
                        lat_name, lon_name = ('latitude', 'longitude') if model == 'NorESM2-LM' else ('lat', 'lon')
                        native_lat = source[lat_name].values
                        native_lon = source[lon_name].values
                    audit.append({'model': model, 'period': period, 'variable': variable, 'file': filename,
                                  'dimensions': list(source[variable].dims), 'selected_months': int(selected.sum())})
            monthly = xr.concat(parts, 'time', join='exact').sortby('time')
            ids = monthly.time.dt.year.values * 12 + monthly.time.dt.month.values
            expected = np.array([year * 12 + month for year in range(start, start + 10) for month in [1, 2, 12]])
            assert np.array_equal(ids, expected)
            if variable in ['pr', 'evspsbl']:
                monthly = monthly * 86400.
            monthly = monthly.assign_coords(lon=monthly.lon % 360).sortby('lon')
            left = monthly.isel(lon=[-1]).assign_coords(lon=[float(monthly.lon[-1]) - 360])
            right = monthly.isel(lon=[0]).assign_coords(lon=[float(monthly.lon[0]) + 360])
            monthly = xr.concat([left, monthly, right], 'lon').interp(lat=latitude, lon=longitude)
            monthly_phases.append(monthly)
        name = 'u850' if variable == 'ua' else variable
        support = monthly_phases[0].notnull().all('time').values & monthly_phases[1].notnull().all('time').values
        support = xr.DataArray(support, dims=('lat', 'lon'), coords={'lat': latitude, 'lon': longitude})
        for phase_index, monthly in enumerate(monthly_phases):
            yearly_fields = []
            yearly_days = []
            start = starts[model][phase_index]
            for year in range(start, start + 10):
                selected = monthly.sel(time=monthly.time.dt.year == year)
                assert selected.sizes['time'] == 3
                yearly_fields.append(selected.weighted(selected.time.dt.days_in_month).mean('time').where(support))
                yearly_days.append(float(selected.time.dt.days_in_month.sum()))
            phases.append(xr.concat(yearly_fields, pd.Index(range(1, 11), name='year_in_window')))
            weights_by_phase.append(yearly_days)
        variables[name] = xr.concat(phases, pd.Index(periods, name='period'))
        day_weights = xr.DataArray(weights_by_phase, dims=('period', 'year_in_window'), coords={'period': periods, 'year_in_window': range(1, 11)})
        recovered_mean = variables[name].weighted(day_weights).mean('year_in_window')
        error = float(abs(recovered_mean.sel(lat=reference.lat, lon=reference.lon) - reference[name].sel(model=model, drop=True)).max())
        assert error < 1.e-6, (model, name, error)
        errors.append({'model': model, 'variable': name, 'original_climatology_max_error': error})
        print(model, name, '10+10 samples; reference error', error, flush=True)
    dataset = xr.Dataset(variables)
    dataset['pe'] = dataset.pr - dataset.evspsbl
    dataset['day_weight'] = day_weights
    dataset['calendar_year'] = xr.DataArray([np.arange(start, start + 10) for start in starts[model]], dims=('period', 'year_in_window'))
    dataset.attrs = {'model': model, 'member': members[model], 'sample': 'Each calendar-year Jan/Feb/Dec dayweighted;10samples per endpoint, not30independent winters',
                     'units': 'pr/E/PE mmday-1,u850 m/s', 'statistical_warning': 'Short transient windows; within-window stationarity must not be assumed automatically'}
    products[model] = dataset


# %% 4. 保存样本、来源和与旧结果的一致性

output.mkdir(parents=True)
for model, dataset in products.items():
    dataset.to_netcdf(output / (model + '.nc'))
pd.DataFrame(audit).to_csv(output / 'input_audit.csv', index=False)
pd.DataFrame(errors).to_csv(output / 'climatology_reproduction.csv', index=False)
(output / 'qc.json').write_text(json.dumps({'models': models, 'dimensions': {'samples': ['period', 'year_in_window', 'lat', 'lon'], 'days': ['period', 'year_in_window']},
    'max_error': max(row['original_climatology_max_error'] for row in errors),
    'domain': 'All longitudes; positive native CNRM Gaussian latitudes, near equator to87.86N; no extrapolation to poles.',
    'meaning': 'Raw endpoint samples for per-model temporal variability, not intermodel ensemble significance; UK excluded; no drift correction'}, indent=2), encoding='utf-8')
