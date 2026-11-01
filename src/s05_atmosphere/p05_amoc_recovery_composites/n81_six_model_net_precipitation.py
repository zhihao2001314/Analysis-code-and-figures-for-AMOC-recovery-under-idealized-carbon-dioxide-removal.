# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import false_discovery_control
from threadpoolctl import threadpool_limits


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
base = root / 'output/s05_atmosphere/p05_amoc_recovery_composites'
source = base / 'n72_six_model_annual_circulation_precipitation_samples'
bootstrap_source = base / 'n73_six_model_precipitation_significance'
destination = base / 'n81_six_model_net_precipitation'
assert not destination.exists(), '正式结果存在，须先归档。'
models = ['GFDL-ESM4', 'NorESM2-LM', 'CAS-ESM2-0', 'CESM2', 'CNRM-ESM2-1', 'UKESM1-0-LL']
products = ['recovered', 'not_recovered', 'contrast']
periods = ['early', 'later']
experiments = ['1pctCO2', '1pctCO2-cdr']
starts = {'GFDL-ESM4': [1, 271], 'NorESM2-LM': [1, 271], 'CAS-ESM2-0': [200, 470],
          'CESM2': [1, 131], 'CNRM-ESM2-1': [1850, 2120], 'UKESM1-0-LL': [1850, 2120]}
members = {'GFDL-ESM4': 'r1i1p1f1', 'NorESM2-LM': 'r1i1p1f1', 'CAS-ESM2-0': 'r1i1p1f1',
           'CESM2': 'r1i1p1f1', 'CNRM-ESM2-1': 'r1i1p1f2', 'UKESM1-0-LL': 'r1i1p1f2'}
replicates, block_lengths, alpha = 9999, [2, 3, 4], .05
contrast_coefficients = np.array([-1/3, -1/3, -1/3, 1/3, 1/3, 1/3])
inventory, endpoint_rows, evaporation_parts = [], [], []


# %% 3. 读取现有降水及各模式初末十年月蒸发

with xr.open_dataset(source / 'annual_precipitation_samples.nc') as dataset:
    precipitation = dataset.sel(model=models).load()
assert precipitation.pr.dims == ('model', 'period', 'year_in_window', 'lat', 'lon')
latitude, longitude = precipitation.lat.values, precipitation.lon.values
day_weights = precipitation.day_weight.transpose('model', 'period', 'year_in_window').values
assert list(precipitation.period.values) == periods
for model_index, model in enumerate(models):
    phase_fields = []
    for period_index, period in enumerate(periods):
        folder = Path('F:/sea/cmip6evspsblad' if period_index == 0 else 'F:/sea/cmip6evspsbl')
        first_year = starts[model][period_index]
        paths = []
        for path in sorted(folder.rglob(f'evspsbl_Amon_{model}_{experiments[period_index]}_{members[model]}_*.nc')):
            dates = path.stem.rsplit('_', 1)[1].split('-')
            if int(dates[0][:4]) <= first_year + 9 and int(dates[1][:4]) >= first_year:
                paths.append(path)
        assert paths, (model, period)
        pieces = []
        for path in paths:
            with xr.open_dataset(path, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as dataset:
                assert dataset.attrs['source_id'] == model
                assert dataset.attrs['experiment_id'] == experiments[period_index]
                assert dataset.attrs['variant_label'] == members[model]
                assert dataset.evspsbl.dims == ('time', 'lat', 'lon')
                assert dataset.evspsbl.attrs['units'] == 'kg m-2 s-1'
                expected_name = 'water_evaporation_flux' if model == 'CNRM-ESM2-1' and period == 'early' else 'water_evapotranspiration_flux'
                assert dataset.evspsbl.attrs['standard_name'] == expected_name
                selected = dataset.time.dt.year.isin(range(first_year, first_year + 10))
                # evspsbl 向上为正，包含升华和蒸腾；1 kg/m² 水相当于 1 mm 水深。
                field = dataset.evspsbl.isel(time=selected).load().astype(float) * 86400.
                pieces.append(field.reset_coords(drop=True))
                inventory.append({'model': model, 'period': period, 'path': str(path),
                                  'raw_dimensions': str(dataset.evspsbl.dims), 'raw_shape': str(dataset.evspsbl.shape),
                                  'units': dataset.evspsbl.attrs['units'], 'standard_name': dataset.evspsbl.attrs['standard_name'],
                                  'calendar': dataset.time.dt.calendar, 'months_selected': field.sizes['time'],
                                  'bytes': path.stat().st_size, 'mtime_ns': path.stat().st_mtime_ns})
        monthly = xr.concat(pieces, 'time', join='exact').sortby('time')
        np.testing.assert_array_equal(monthly.time.dt.year.values * 12 + monthly.time.dt.month.values,
                                      np.arange(first_year * 12 + 1, (first_year + 10) * 12 + 1))
        assert monthly.shape[0] == 120 and np.isfinite(monthly).all()
        annual_fields = []
        for year_index, year in enumerate(range(first_year, first_year + 10)):
            one_year = monthly.sel(time=monthly.time.dt.year == year)
            days = one_year.time.dt.days_in_month.astype(float)
            assert float(days.sum()) == day_weights[model_index, period_index, year_index]
            native = one_year.weighted(days).mean('time', skipna=False)
            native = native.assign_coords(lon=(native.lon + 180) % 360 - 180).sortby('lon').sortby('lat')
            periodic = xr.concat([native.isel(lon=[-1]).assign_coords(lon=native.lon.values[[-1]] - 360),
                                  native, native.isel(lon=[0]).assign_coords(lon=native.lon.values[[0]] + 360)], 'lon')
            mapped = periodic.interp(lat=latitude, lon=longitude).transpose('lat', 'lon')
            assert np.isfinite(mapped).all(), (model, year)
            annual_fields.append(mapped.expand_dims(year_in_window=[year_index + 1]))
        phase_fields.append(xr.concat(annual_fields, 'year_in_window').expand_dims(period=[period]))
        endpoint_rows.append({'model': model, 'period': period, 'first_calendar_year': first_year,
                              'last_calendar_year': first_year + 9, 'months': 120,
                              'monthly_values_sha256': hashlib.sha256(np.ascontiguousarray(monthly.values).tobytes()).hexdigest()})
        print(model, period, first_year, '120 months of evaporation complete', flush=True)
    evaporation_parts.append(xr.concat(phase_fields, 'period').expand_dims(model=[model]))


# %% 4. 净降水、三模式等权分组和符号一致性

evaporation = xr.concat(evaporation_parts, 'model').transpose('model', 'period', 'year_in_window', 'lat', 'lon')
samples = xr.Dataset({'pr': precipitation.pr, 'evspsbl': evaporation, 'pe': precipitation.pr - evaporation,
                      'day_weight': precipitation.day_weight})
for variable in ['pr', 'evspsbl', 'pe']:
    samples[variable].attrs['units'] = 'mm day-1'
samples.pe.attrs['definition'] = 'P minus upward surface evapotranspiration E, including sublimation; positive means net surface freshwater gain.'
weights = samples.day_weight / samples.day_weight.sum('year_in_window')
climate = (samples[['pr', 'evspsbl', 'pe']] * weights).sum('year_in_window', skipna=False)
member_changes = climate.sel(period='later', drop=True) - climate.sel(period='early', drop=True)
group_climate = xr.concat([climate.sel(model=models[:3]).mean('model', skipna=False).expand_dims(product=['recovered']),
                           climate.sel(model=models[3:]).mean('model', skipna=False).expand_dims(product=['not_recovered'])], 'product')
group_response = group_climate.sel(period='later', drop=True) - group_climate.sel(period='early', drop=True)
contrast = (group_response.sel(product='not_recovered', drop=True)
            - group_response.sel(product='recovered', drop=True)).expand_dims(product=['contrast'])
changes = xr.concat([group_response, contrast], 'product')
for dataset in [samples, climate, member_changes, group_climate, changes]:
    np.testing.assert_allclose(dataset.pe, dataset.pr - dataset.evspsbl, atol=1e-12, rtol=0)
with xr.open_dataarray(source / 'group_precipitation_changes.nc') as dataset:
    np.testing.assert_allclose(changes.pr, dataset.sel(product=products), atol=1e-10, rtol=0)
values = samples.pe.values.reshape(6, 2, 10, -1)
observed = changes.pe.sel(product='contrast').values.ravel()
probabilities = np.empty((3, values.shape[-1]))


# %% 5. 组间差值重新执行循环年块 bootstrap 与 FDR

# 沿用降水图保存的随机年索引，但对 P-E 样本独立重算统计量和 p/q。
indices_archive = np.load(bootstrap_source / 'bootstrap_year_indices.npz')
for length_index, length in enumerate(block_lengths):
    indices = indices_archive['L' + str(length)]
    assert indices.shape == (6, 2, replicates, 10)
    resampling = np.empty((6, 2, replicates, 10))
    for model_index in range(6):
        for period_index in range(2):
            counts = np.eye(10)[indices[model_index, period_index]].sum(axis=1)
            weighted_counts = counts * day_weights[model_index, period_index][None]
            resampling[model_index, period_index] = weighted_counts / weighted_counts.sum(axis=1, keepdims=True)
    for offset in range(0, values.shape[-1], 256):
        last = min(offset + 256, values.shape[-1])
        with threadpool_limits(limits=4):
            boot = np.zeros((replicates, last - offset))
            for model_index, coefficient in enumerate(contrast_coefficients):
                boot += coefficient * (resampling[model_index, 1] @ values[model_index, 1, :, offset:last]
                                       - resampling[model_index, 0] @ values[model_index, 0, :, offset:last])
            errors = boot - observed[None, offset:last]
            probabilities[length_index, offset:last] = (1 + np.sum(abs(errors) >= abs(observed[None, offset:last]), axis=0)) / (replicates + 1)
    print('P-E contrast bootstrap block', length, 'complete', flush=True)
p_max = probabilities.max(axis=0)
q_values = false_discovery_control(p_max, method='bh')
passed = np.empty(changes.pe.shape, dtype=bool)
for index, selected_models in enumerate([models[:3], models[3:]]):
    individual = member_changes.pe.sel(model=selected_models).values
    passed[index] = np.all(individual > 0, axis=0) | np.all(individual < 0, axis=0)
passed[2] = (q_values < alpha).reshape(changes.pe.shape[1:])
display = xr.Dataset({'delta': changes.pe, 'passed': (changes.pe.dims, passed),
                      'not_passed': (changes.pe.dims, ~passed)})
display['criterion'] = xr.DataArray(['all_three_models_same_sign', 'all_three_models_same_sign',
                                    'two_sided_block_bootstrap_BH_q_lt_0.05'], dims='product', coords={'product': products})
display['p_max'] = xr.DataArray(p_max.reshape(changes.pe.shape[1:]), dims=['lat', 'lon'], coords={'lat': latitude, 'lon': longitude})
display['q'] = xr.DataArray(q_values.reshape(changes.pe.shape[1:]), dims=['lat', 'lon'], coords={'lat': latitude, 'lon': longitude})
display['p_by_block'] = xr.DataArray(probabilities.reshape((3,) + changes.pe.shape[1:]), dims=['block_length', 'lat', 'lon'],
                                      coords={'block_length': block_lengths, 'lat': latitude, 'lon': longitude})


# %% 6. 保存结果、诊断和来源

destination.mkdir(parents=True)
samples.to_netcdf(destination / 'annual_hydrology_samples.nc')
climate.to_netcdf(destination / 'member_climatology.nc')
member_changes.to_netcdf(destination / 'member_changes.nc')
group_climate.to_netcdf(destination / 'group_climatology.nc')
changes.to_netcdf(destination / 'group_changes.nc')
display.to_netcdf(destination / 'pe_display_criteria.nc')
pd.DataFrame(inventory).to_csv(destination / 'evaporation_source_inventory.csv', index=False)
pd.DataFrame(endpoint_rows).to_csv(destination / 'endpoint_qc.csv', index=False)
diagnostics = []
area = np.broadcast_to(np.cos(np.deg2rad(latitude))[:, None], changes.pe.shape[1:]).copy()
area /= area.sum()
for product_index, product in enumerate(products):
    for variable in ['pr', 'evspsbl', 'pe']:
        field = changes[variable].sel(product=product).values
        diagnostics.append({'product': product, 'variable': variable, 'minimum_mm_day': float(field.min()),
                            'maximum_mm_day': float(field.max()), 'global_mean_change_mm_day': float(np.sum(area * field)),
                            'outside_3_mm_day_cells': int(np.sum(abs(field) > 3)),
                            'pe_criterion_passed_cells': int(passed[product_index].sum())})
pd.DataFrame(diagnostics).to_csv(destination / 'group_diagnostics.csv', index=False)
source_paths = [source / 'annual_precipitation_samples.nc', source / 'group_precipitation_changes.nc',
                bootstrap_source / 'bootstrap_year_indices.npz', bootstrap_source / 'metadata.json']
metadata = {'status': 'passed', 'models': models, 'groups': {'recovered': models[:3], 'not_recovered': models[3:]},
            'periods': periods, 'calendar_starts': starts, 'units': 'mm day-1', 'definition': 'pe=pr-evspsbl; E upward positive',
            'time': 'CF month-length-weighted annual means, CF year-length-weighted initial/final decades.',
            'space': 'Native annual evaporation bilinearly mapped periodically to the existing 2-degree precipitation grid, no extrapolation.',
            'drift_correction': False, 'group_weights': 'Equal weight for three fixed models; no model resampling.',
            'raw_dimensions': 'evspsbl(time,lat,lon); units kg m-2 s-1; CNRM ramp-up water_evaporation_flux; other files water_evapotranspiration_flux. All describe upward surface vapor flux.',
            'saved_dimensions': 'samples(model,period,year_in_window,lat,lon); changes/masks(product,lat,lon)',
            'group_criterion': 'All three individual-model P-E changes strictly positive or all strictly negative; zero fails. Not a p-value test.',
            'contrast_criterion': 'Centered two-sided circular block bootstrap,9999 draws, max p over L=2,3,4, then spatial BH q<0.05.',
            'bootstrap_indices': 'Reuse n73 random indices with same six-model order and day weights; recompute P-E p/q; no copied precipitation mask.',
            'bootstrap_limits': 'Temporal sampling conditional on fixed model membership; short ten-year endpoint samples.',
            'initial_final': 'Experiment years1-10 and271-280; individual calendar windows in endpoint_qc.csv.',
            'P_peak_markers': 'Omitted from the new P-E map; precipitation maxima are not net-precipitation maxima.',
            'source_files': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in source_paths]}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print(pd.DataFrame(diagnostics).to_string(index=False), flush=True)
print(destination.resolve(), flush=True)
