# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay, cKDTree
from shapely import contains_xy
from shapely.geometry import box, mapping, shape


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
comparison_directory = root / 'output/s02_sia/s02_arctic_region/p24_mld_method_comparison/n079_mld_method_comparison'
native_directory = root / 'output/s02sia/s02_arctic_region/p06_mld/n02_na_deep_mixing_three_phase'
geography_directory = root / 'output/s02_sia/s02_arctic_region/p23_stratification_mld_composite/n077_regional_definition_comparison'
destination = root / 'output/s02_sia/s02_arctic_region/p24_mld_method_comparison/n081_mld_method_pattern_inputs'
models = ['ACCESS-ESM1-5', 'CanESM5', 'CESM2', 'CNRM-ESM2-1',
          'GFDL-ESM4', 'NorESM2-LM', 'UKESM1-0-LL']
methods = ['native', 'derived']
periods = ['initial', 'late']
extent = [-65.0, 20.0, 50.0, 85.0]
grid_spacing = 0.25
longitude = np.arange(extent[0], extent[1] + grid_spacing / 2, grid_spacing)
latitude = np.arange(extent[2], extent[3] + grid_spacing / 2, grid_spacing)
longitude_grid, latitude_grid = np.meshgrid(longitude, latitude)
change_limit_m = 1500.0
color_levels = np.arange(-change_limit_m, change_limit_m + 50, 50)
color_ticks = np.arange(-change_limit_m, change_limit_m + 500, 500)
initial_contour_levels_m = np.array([200, 500, 1000, 1500, 2000], dtype=float)
assert not destination.exists(), '正式结果已存在，重跑前先归档。'


# %% 3. 读取配对结果与现有底图

cells_path = comparison_directory / 'paired_native_cells.csv'
metadata_path = comparison_directory / 'metadata.json'
cells = pd.read_csv(cells_path)
comparison_metadata = json.loads(metadata_path.read_text(encoding='utf-8'))
assert comparison_metadata['models'] == models
assert not cells.duplicated(['model', 'latitude_degN', 'longitude_degE']).any()
destination.mkdir(parents=True)
source_paths = [cells_path, metadata_path]
clip_geometry = box(extent[0] - 2, extent[2] - 2, extent[1] + 2, extent[3] + 2)
land_mask = np.zeros(longitude_grid.shape, dtype=bool)
for geography_name in ['land', 'coastline']:
    source_path = geography_directory / ('map_' + geography_name + '_50m.geojson')
    source_paths.append(source_path)
    source_features = json.loads(source_path.read_text(encoding='utf-8'))['features']
    prepared_features = []
    for feature in source_features:
        geometry = shape(feature['geometry']).intersection(clip_geometry)
        if geometry.is_empty:
            continue
        prepared_features.append({'type': 'Feature', 'properties': {}, 'geometry': mapping(geometry)})
        if geography_name == 'land':
            land_mask |= contains_xy(geometry, longitude_grid, latitude_grid)
    (destination / (geography_name + '_50m.geojson')).write_text(
        json.dumps({'type': 'FeatureCollection', 'features': prepared_features}), encoding='utf-8')


# %% 4. 还原原生网格缺测并统一线性显示映射

mapped_fields = np.full((7, 2, 2, latitude.size, longitude.size), np.nan)
diagnostic_rows = []
for model_index, model in enumerate(models):
    paired = cells.loc[cells.model == model].copy()
    native_path = native_directory / (model + '_na_march_mld_three_phase_native.npz')
    source_paths.append(native_path)
    with np.load(native_path) as saved:
        native_latitude = saved['latitude_deg_north'].copy()
        native_longitude = (saved['longitude_deg_east'].copy() + 180) % 360 - 180
        source_initial = saved['up_initial_march_mld_m'].copy()
        source_late = saved['down_final_march_mld_m'].copy()
    assert native_latitude.shape == native_longitude.shape == source_initial.shape
    native_coordinates = np.column_stack([native_longitude.ravel(), native_latitude.ravel()])
    pair_coordinates = paired[['longitude_degE', 'latitude_degN']].to_numpy()
    distance, paired_index = cKDTree(native_coordinates).query(pair_coordinates)
    assert distance.max() < 1e-10
    assert len(np.unique(paired_index)) == len(paired)
    fields = np.full((2, 2, native_latitude.size), np.nan)
    for method_index, method in enumerate(methods):
        for period_index, period in enumerate(periods):
            fields[method_index, period_index, paired_index] = paired[method + '_' + period + '_m'].to_numpy()
    assert np.allclose(fields[1, 0, paired_index], source_initial.ravel()[paired_index], atol=1e-8, rtol=0)
    assert np.allclose(fields[1, 1, paired_index], source_late.ravel()[paired_index], atol=1e-8, rtol=0)
    assert np.array_equal(np.isfinite(fields[0]), np.isfinite(fields[1]))

    # 包含缺测原生节点参与三角剖分，避免仅用有效点跨越内部缺口补值。
    source_domain = np.isfinite(native_coordinates).all(axis=1)
    source_domain &= (native_coordinates[:, 0] >= extent[0] - 1) & (native_coordinates[:, 0] <= extent[1] + 1)
    source_domain &= (native_coordinates[:, 1] >= extent[2] - 1) & (native_coordinates[:, 1] <= extent[3] + 1)
    points = native_coordinates[source_domain]
    values = fields[:, :, source_domain].reshape(4, -1).T
    assert len(np.unique(points, axis=0)) == len(points)
    triangulation = Delaunay(points)
    interpolator = LinearNDInterpolator(triangulation, values, fill_value=np.nan)
    mapped = interpolator(longitude_grid, latitude_grid)
    mapped = mapped.reshape(latitude.size, longitude.size, 2, 2).transpose(2, 3, 0, 1)
    mapped[:, :, land_mask] = np.nan
    assert np.array_equal(np.isfinite(mapped[0]), np.isfinite(mapped[1]))
    mapped_fields[model_index] = mapped
    np.savez_compressed(destination / (model + '_paired_native_grid.npz'),
                        longitude=native_longitude, latitude=native_latitude,
                        fields=fields.reshape((2, 2) + native_latitude.shape),
                        paired_flat_index=paired_index)

    for method_index, method in enumerate(methods):
        change = mapped[method_index, 1] - mapped[method_index, 0]
        native_change = paired[method + '_delta_m'].to_numpy()
        weights = paired.area_m2.to_numpy()
        diagnostic_rows.append({'model': model, 'method': method, 'paired_native_cells': len(paired),
                                'display_valid_cells': int(np.isfinite(change).sum()),
                                'native_change_min_m': float(native_change.min()),
                                'native_change_max_m': float(native_change.max()),
                                'display_change_min_m': float(np.nanmin(change)),
                                'display_change_max_m': float(np.nanmax(change)),
                                'native_area_fraction_outside_color_limits': float(np.average(np.abs(native_change) > change_limit_m, weights=weights)),
                                'display_cells_outside_color_limits': int((np.abs(change) > change_limit_m).sum())})
    print(model, 'display cells', np.isfinite(mapped[0, 0]).sum(), flush=True)


# %% 5. 保存共同显示数组及来源记录

changes = mapped_fields[:, :, 1] - mapped_fields[:, :, 0]
assert np.isfinite(changes).any(axis=(2, 3)).all()
assert np.array_equal(np.isfinite(changes[:, 0]), np.isfinite(changes[:, 1]))
np.savez_compressed(destination / 'paired_mld_pattern_fields.npz',
                    models=np.array(models), methods=np.array(methods),
                    longitude=longitude, latitude=latitude,
                    initial=mapped_fields[:, :, 0], late=mapped_fields[:, :, 1], change=changes,
                    land_mask=land_mask, color_levels=color_levels, color_ticks=color_ticks,
                    contour_levels=initial_contour_levels_m)
pd.DataFrame(diagnostic_rows).to_csv(destination / 'pattern_diagnostics.csv', index=False)
source_records = []
for path in source_paths:
    source_records.append({'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
metadata = {'models': models, 'methods': methods, 'extent': extent, 'grid_spacing_deg': grid_spacing,
            'dimensions': 'initial/late/change(model7,method2,latitude141,longitude341); native fields(method2,period2,row,column)',
            'color': 'RdBu_r; negative=shoaling,positive=deepening; common limits+-1500m with extended tails',
            'shading': 'Late minus initial March MLD in metres',
            'contours': 'Initial March MLD in metres; identical levels for both methods',
            'time_support': comparison_metadata['periods'],
            'native_pairing': 'Identical to n079; no new model selection or raw-data computation',
            'display_mapping': 'Linear interpolation on one shared Delaunay triangulation per model,0.25degree display grid; no smoothing or extrapolation. Missing native nodes retained; interpolated triangles touching missing values remain missing.',
            'geography': 'Existing50m land/coast geometry clipped to regional view; land display mask applied identically to both methods',
            'statistics': 'No scientific statistics recomputed from display grids; native-grid diagnostics remain in n079',
            'pi_control_correction': False, 'sources': source_records}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print('OUTPUT', destination.resolve())
