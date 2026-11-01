# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import shutil
import netCDF4 as nc
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
destination = root / 'output/s05_atmosphere/p07_toa_energy_balance/n23_figure3_bipolar_ice_edges'
assert not destination.exists(), '正式结果存在，须先归档。'
observation_source = root / 'output/s03_tem/n12_global_temperature_ice_patterns'
land_path = observation_source / 'global_tas_ice_fields.npz'
models = ['CESM2', 'CNRM-ESM2-1', 'UKESM1-0-LL', 'CAS-ESM2-0', 'GFDL-ESM4', 'NorESM2-LM']
members = ['r1i1p1f1', 'r1i1p1f2', 'r1i1p1f2', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1']
period_years = [[[1, 10], [131, 140]], [[1850, 1859], [2120, 2129]],
                [[1850, 1859], [2120, 2129]], [[200, 209], [470, 479]],
                [[6, 15], [271, 280]], [[1, 10], [271, 280]]]
dimensions = [('time', 'nj', 'ni'), ('time', 'y', 'x'), ('time', 'j', 'i'),
              ('time', 'j', 'i'), ('time', 'y', 'x'), ('time', 'j', 'i')]
coordinate_names = [('lat', 'lon'), ('lat', 'lon'), ('latitude', 'longitude'),
                    ('latitude', 'longitude'), ('lat', 'lon'), ('latitude', 'longitude')]
latitude, longitude = np.arange(-89.5, 90., 1.), np.arange(-179.5, 180., 1.)
longitude_grid, latitude_grid = np.meshgrid(longitude, latitude)
target_lat_rad, target_lon_rad = np.deg2rad(latitude_grid.ravel()), np.deg2rad(longitude_grid.ravel())
target_xyz = np.column_stack([np.cos(target_lat_rad) * np.cos(target_lon_rad),
                              np.cos(target_lat_rad) * np.sin(target_lon_rad), np.sin(target_lat_rad)])
model_sic = np.full((6, 2, 180, 360), np.nan)
inventory, diagnostics = [], []


# %% 3. 两个十年全年平均海冰，CAS降低段按时间边界中点识别月份

for model_index, model in enumerate(models):
    for period_index, (experiment, folder) in enumerate([('1pctCO2', 'cmipseaicead'),
                                                        ('1pctCO2-cdr', 'cmipseaice')]):
        first_year, last_year = period_years[model_index][period_index]
        cas_cdr = model == 'CAS-ESM2-0' and period_index == 1
        paths = [Path('F:/sea/cmipseaice/CAS/siconc.nc')] if cas_cdr else sorted(
            Path('F:/sea', folder).rglob(f'siconc_*_{model}_{experiment}_{members[model_index]}_*.nc'))
        month_ids, field_sum, complete = [], None, None
        for path in paths:
            with nc.Dataset(path) as dataset:
                if cas_cdr:
                    variable, lat_name, lon_name = 'aice', 'TLAT', 'TLON'
                    expected_dimensions = ('time', 'nj', 'ni')
                    time = dataset['time']
                    date_values = np.asarray(dataset['time_bnds'][:]).mean(axis=1)
                else:
                    assert dataset.source_id == model and dataset.experiment_id == experiment
                    assert dataset.variant_label == members[model_index]
                    variable = 'siconc'
                    lat_name, lon_name = coordinate_names[model_index]
                    expected_dimensions = dimensions[model_index]
                    time = dataset['time']
                    date_values = time[:]
                field = dataset[variable]
                assert field.dimensions == expected_dimensions and field.units == '%'
                dates = nc.num2date(date_values, time.units, calendar=time.calendar)
                selected = np.flatnonzero([first_year <= date.year <= last_year for date in dates])
                if not len(selected):
                    continue
                native_latitude = np.asarray(dataset[lat_name][:], dtype=float)
                native_longitude = np.asarray(dataset[lon_name][:], dtype=float)
                assert native_latitude.shape == field.shape[1:] == native_longitude.shape
                if field_sum is None:
                    field_sum = np.zeros(field.shape[1:])
                    complete = np.ones(field.shape[1:], dtype=bool)
                for offset in range(0, len(selected), 12):
                    indices = selected[offset:offset + 12]
                    values = np.asarray(np.ma.filled(field[indices], np.nan), dtype=float) / 100.
                    values[(values < 0) | (values > 1)] = np.nan
                    field_sum += np.where(np.isfinite(values), values, 0).sum(axis=0)
                    complete &= np.isfinite(values).all(axis=0)
                month_ids.extend([dates[index].year * 12 + dates[index].month - 1 for index in selected])
                inventory.append({'model': model, 'period': ['initial', 'final'][period_index],
                                  'path': str(path), 'variable': variable, 'dimensions': str(field.dimensions),
                                  'shape': str(field.shape), 'coordinates': str((lat_name, lon_name)),
                                  'units': field.units, 'calendar': time.calendar, 'months': len(selected),
                                  'first_selected': str(dates[selected[0]]), 'last_selected': str(dates[selected[-1]]),
                                  'time_assignment': 'CF bounds midpoint' if cas_cdr else 'CF time coordinate',
                                  'file_bytes': path.stat().st_size, 'mtime_ns': path.stat().st_mtime_ns})
        np.testing.assert_array_equal(np.sort(month_ids), np.arange(first_year * 12, (last_year + 1) * 12))
        assert len(month_ids) == 120
        native_mean = np.where(complete, field_sum / 120., np.nan)
        valid = np.isfinite(native_mean) & np.isfinite(native_latitude) & np.isfinite(native_longitude)
        lat_rad, lon_rad = np.deg2rad(native_latitude[valid]), np.deg2rad(native_longitude[valid])
        native_xyz = np.column_stack([np.cos(lat_rad) * np.cos(lon_rad),
                                      np.cos(lat_rad) * np.sin(lon_rad), np.sin(lat_rad)])
        chord, indices = cKDTree(native_xyz).query(target_xyz, k=8)
        distance = 6371.0088 * 2 * np.arcsin(np.clip(chord / 2, 0, 1))
        weights = 1 / np.maximum(distance, 1e-6) ** 2
        mapped = np.sum(weights * native_mean[valid][indices], axis=1) / weights.sum(axis=1)
        exact = distance[:, 0] <= 1e-6
        mapped[exact] = native_mean[valid][indices[exact, 0]]
        mapped[distance[:, 0] > 300.] = np.nan
        model_sic[model_index, period_index] = mapped.reshape(180, 360)
        diagnostics.append({'model': model, 'period': ['initial', 'final'][period_index],
                            'months': len(month_ids), 'native_complete_cells': int(valid.sum()),
                            'target_valid_cells': int(np.isfinite(mapped).sum())})
        print(model, ['initial', 'final'][period_index], 'SIC complete', flush=True)


# %% 4. 固定三模式浓度组均值，海陆和两时期共同有效范围

with np.load(land_path) as saved:
    np.testing.assert_array_equal(saved['latitude'], latitude)
    np.testing.assert_array_equal(saved['longitude'], longitude)
    land_mask = saved['land_mask'].copy()
assert land_mask.shape == (180, 360) and land_mask[:10].mean() > .9
group_sic = np.stack([model_sic[:3].mean(axis=0), model_sic[3:].mean(axis=0)])
common = np.isfinite(group_sic).all(axis=1) & ~land_mask[None]
ice_for_contour = np.where(common[:, None] & (np.abs(latitude_grid)[None, None] >= 45), group_sic, np.nan)
for group_index in range(2):
    for period_index in range(2):
        for hemisphere in [latitude >= 45, latitude <= -45]:
            values = ice_for_contour[group_index, period_index, hemisphere]
            assert np.nanmin(values) < .15 < np.nanmax(values)
destination.mkdir(parents=True)
np.savez_compressed(destination / 'grouped_bipolar_sic.npz', latitude=latitude, longitude=longitude,
                    models=np.asarray(models), model_sic=model_sic, group_sic=group_sic,
                    ice_for_contour=ice_for_contour, land_mask=land_mask,
                    longitude_cyclic=np.append(longitude, longitude[0] + 360),
                    ice_cyclic=np.concatenate([ice_for_contour, ice_for_contour[:, :, :, :1]], axis=3))
pd.DataFrame(inventory).to_csv(destination / 'ice_source_inventory.csv', index=False)
pd.DataFrame(diagnostics).to_csv(destination / 'ice_window_qc.csv', index=False)
for hemisphere in ['north', 'south']:
    shutil.copy2(observation_source / f'observed_{hemisphere}.npz', destination / f'observed_{hemisphere}.npz')
metadata = {'status': 'passed', 'models': models, 'groups': {'non_recovery': models[:3], 'recovery': models[3:]},
            'period_years': period_years, 'definition': '15% contour of equal-weight three-model annual-mean SIC; initial and final decades; both hemispheres.',
            'missing_data': 'All120months per native cell and all3models required; paired-period common support; no zero-fill or partial mean.',
            'regrid': '8-neighbour spherical IDW power2 on1degree grid, max nearest300km; land excluded for contours, abs(latitude)>=45.',
            'CAS_cdr_time': 'aice labels are interval ends; select by CF time_bnds midpoints, final470-479. File has no CMIP member attributes; use previously supplied project CAS CDR file.',
            'observations': 'NOAA/NSIDC1988-2020 annual-mean SIC, native EPSG3411 north and3412south; reuse validated n12 fields without smoothing.',
            'panel_f': 'No ice edges on difference map; it has no single initial/final SIC state.',
            'styles': {'initial': 'purple dashed', 'final': 'green dashed', 'observed': 'magenta solid'},
            'source_files': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
                             for path in [land_path, observation_source / 'observed_north.npz', observation_source / 'observed_south.npz']]}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print(destination.resolve(), flush=True)
