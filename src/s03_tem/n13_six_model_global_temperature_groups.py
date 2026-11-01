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
source = root / 'output/s03_tem/p05_global_tas_sst_slp/n01'
cartography_source = root / 'output/s03_tem/n12_global_temperature_ice_patterns/coastline_global_50m.geojson'
destination = root / 'output/s03_tem/n13_six_model_global_temperature_groups'
assert not destination.exists() or not any(destination.iterdir()), '正式结果存在，须先归档。'
groups = ['non_recovery', 'recovery']
members = [['CESM2', 'CNRM-ESM2-1', 'UKESM1-0-LL'],
           ['CAS-ESM2-0', 'GFDL-ESM4', 'NorESM2-LM']]
models = members[0] + members[1]
model_fields = np.full((6, 180, 360), np.nan)
source_paths, membership_rows, summary_rows = [], [], []
window_audit = []
aligned_models = ['GFDL-ESM4', 'UKESM1-0-LL']
aligned_years = {'GFDL-ESM4': [[6, 15], [271, 280]],
                 'UKESM1-0-LL': [[1850, 1859], [2120, 2129]]}
aligned_members = {'GFDL-ESM4': 'r1i1p1f1', 'UKESM1-0-LL': 'r1i1p1f2'}


# %% 3. 复用更新后的六模式全球近地面气温残余场

for model_index, model in enumerate(models):
    path = source / f'global_tas_sst_slp_{model}_pattern_1deg.npz'
    with np.load(path) as saved:
        latitude = saved['latitude_deg_north'].copy()
        longitude = saved['longitude_deg_east'].copy()
        np.testing.assert_allclose(latitude, np.arange(-89.5, 90., 1.))
        np.testing.assert_allclose(longitude, np.arange(-179.5, 180., 1.))
        model_fields[model_index] = saved['raw_delta_tas_degC']
    assert np.isfinite(model_fields[model_index]).all()
    source_paths.append(path)
    membership_rows.append({'group': groups[model_index // 3], 'model': model, 'weight': 1 / 3})


# %% 4. 对齐GFDL初期与UKESM末期窗口，只重读这两模式标准月资料

target_lon, target_lat = np.meshgrid(longitude, latitude)
target_lat_rad, target_lon_rad = np.deg2rad(target_lat.ravel()), np.deg2rad(target_lon.ravel())
target_xyz = np.column_stack([np.cos(target_lat_rad) * np.cos(target_lon_rad),
                              np.cos(target_lat_rad) * np.sin(target_lon_rad), np.sin(target_lat_rad)])
for model in aligned_models:
    mapped_periods = []
    for phase_index, (experiment, folder) in enumerate([('1pctCO2', 'cmip6tasad'),
                                                      ('1pctCO2-cdr', 'cmip6tas')]):
        first_year, last_year = aligned_years[model][phase_index]
        field_parts, month_parts = [], []
        for path in sorted(Path('F:/sea', folder).rglob(
                f'tas_Amon_{model}_{experiment}_{aligned_members[model]}_*.nc')):
            with nc.Dataset(path) as dataset:
                assert dataset.source_id == model and dataset.experiment_id == experiment
                assert dataset.variant_label == aligned_members[model]
                field = dataset['tas']
                assert field.dimensions == ('time', 'lat', 'lon') and field.units == 'K'
                time = dataset['time']
                dates = nc.num2date(time[:], time.units, calendar=time.calendar)
                selected = np.flatnonzero([first_year <= date.year <= last_year for date in dates])
                if not len(selected):
                    continue
                values = np.asarray(np.ma.filled(field[selected], np.nan), dtype=float)
                assert np.isfinite(values).all()
                field_parts.append(values)
                month_parts.extend([dates[index].year * 12 + dates[index].month - 1 for index in selected])
                native_latitude = np.asarray(dataset['lat'][:], dtype=float)
                native_longitude = np.asarray(dataset['lon'][:], dtype=float)
                window_audit.append({'model': model, 'experiment': experiment, 'path': str(path),
                                     'dimensions': str(field.dimensions), 'units': field.units,
                                     'first_selected': str(dates[selected[0]]),
                                     'last_selected': str(dates[selected[-1]]), 'months': len(selected),
                                     'file_bytes': path.stat().st_size, 'mtime_ns': path.stat().st_mtime_ns})
        order = np.argsort(month_parts)
        np.testing.assert_array_equal(np.asarray(month_parts)[order],
                                       np.arange(first_year * 12, (last_year + 1) * 12))
        native_mean = np.concatenate(field_parts)[order].mean(axis=0) - 273.15
        native_lon, native_lat = np.meshgrid(native_longitude, native_latitude)
        source_lat_rad, source_lon_rad = np.deg2rad(native_lat.ravel()), np.deg2rad(native_lon.ravel())
        source_xyz = np.column_stack([np.cos(source_lat_rad) * np.cos(source_lon_rad),
                                     np.cos(source_lat_rad) * np.sin(source_lon_rad), np.sin(source_lat_rad)])
        chord, indices = cKDTree(source_xyz).query(target_xyz, k=4)
        distance = 6371.0088 * 2 * np.arcsin(np.clip(chord / 2, 0, 1))
        neighbor_weights = 1 / np.maximum(distance, 1e-6) ** 2
        mapped = np.sum(neighbor_weights * native_mean.ravel()[indices], axis=1) / neighbor_weights.sum(axis=1)
        exact = distance[:, 0] <= 1e-6
        mapped[exact] = native_mean.ravel()[indices[exact, 0]]
        assert np.all(distance[:, 0] <= 450.)
        mapped_periods.append(mapped.reshape(180, 360))
    model_fields[models.index(model)] = mapped_periods[1] - mapped_periods[0]


# %% 5. 固定三成员等权组均值与不恢复减恢复

group_fields = model_fields.reshape(2, 3, 180, 360).mean(axis=1)
panels = np.concatenate([group_fields, (group_fields[0] - group_fields[1])[None]], axis=0)
latitude_bounds = np.deg2rad(np.column_stack([latitude - .5, latitude + .5]))
latitude_weights = np.sin(latitude_bounds[:, 1]) - np.sin(latitude_bounds[:, 0])
weights = np.broadcast_to(latitude_weights[:, None], (180, 360)).copy()
weights /= weights.sum()
for name, field in zip(models + groups + ['non_recovery_minus_recovery'],
                       list(model_fields) + list(panels)):
    summary_rows.append({'series': name, 'global_mean_degC': float(np.sum(field * weights)),
                         'min_degC': float(field.min()), 'max_degC': float(field.max()),
                         'valid_cells': int(np.isfinite(field).sum())})
color_limit = float(np.ceil(np.max(np.abs(panels)) / 2) * 2)


# %% 6. 保存绘图数组与来源说明

destination.mkdir(parents=True, exist_ok=True)
np.savez_compressed(destination / 'global_temperature_groups.npz', latitude=latitude, longitude=longitude,
                    model_fields=model_fields, group_fields=group_fields, temperature_panels=panels)
pd.DataFrame(membership_rows).to_csv(destination / 'membership.csv', index=False)
pd.DataFrame(summary_rows).to_csv(destination / 'summary.csv', index=False)
pd.DataFrame(window_audit).to_csv(destination / 'aligned_temperature_windows.csv', index=False)
shutil.copy2(cartography_source, destination / cartography_source.name)
source_paths.extend([source / 'global_tas_sst_slp_metadata.json', cartography_source])
caption = ('Residual near-surface air temperature for the non-recovery group (CESM2, CNRM-ESM2-1, '
           'UKESM1-0-LL), recovery group (CAS-ESM2-0, GFDL-ESM4, NorESM2-LM), and their difference '
           '(non-recovery minus recovery). Residual changes are experiment years271-280 minus1-10; '
           'each group is an equal-weight mean of three models. All panels share the same color scale.')
metadata = {'status': 'passed', 'groups': dict(zip(groups, members)), 'models': models,
            'n_models_per_group': 3, 'contrast': 'non_recovery minus recovery',
            'variable': 'tas: near-surface (2 m) air temperature', 'units': 'degC',
            'periods': {'initial': [1, 10], 'final': [271, 280]},
            'definition': 'Mean individual-model final-minus-initial maps; fixed membership, no skipna.',
            'drift_correction': 'none', 'temporal_mean': 'Inherited equal-month endpoint means.',
            'CAS_source': 'Updated standard monthly2m tas, not the old lowest-model-level annual proxy.',
            'regrid': 'Reuse four existing1degree fields; GFDL/UK re-evaluated at fixed windows using identical4-neighbor IDW,450km limit; no smoothing.',
            'window_alignment': aligned_years,
            'raw_window_files': window_audit,
            'dimensions': 'model6 xlat180 xlon360; panel3 xlat180 xlon360',
            'color_limit_degC': color_limit, 'caption': caption,
            'significance': 'No new statistical test; supplied endpoint means are not annual replicates.',
            'source_files': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
                             for path in source_paths]}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
(destination / 'caption.txt').write_text(caption + '\n', encoding='utf-8')
print(pd.DataFrame(summary_rows).to_string(index=False))
print(destination.resolve())
