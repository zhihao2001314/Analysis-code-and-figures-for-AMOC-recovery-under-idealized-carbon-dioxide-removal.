# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json

import netCDF4
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
derived_directory = root / 'output/s02sia/s02_arctic_region/p06_mld/n02_na_deep_mixing_three_phase'
derived_metadata_path = derived_directory / 'na_march_mld_three_phase_metadata.json'
inventory_path = root / 'output/s02sia/s02_arctic_region/p09_mld_mlotst_compare/n00_input_qa/mlotst_file_inventory.csv'
area_registry_path = root / 'output/s02sia/s02_arctic_region/p07_stratification/n01_initial_profiles/na_initial_march_stratification_metadata.json'
destination = root / 'output/s02_sia/s02_arctic_region/p24_mld_method_comparison/n079_mld_method_comparison'
models = ['ACCESS-ESM1-5', 'CanESM5', 'CESM2', 'CNRM-ESM2-1',
          'GFDL-ESM4', 'NorESM2-LM', 'UKESM1-0-LL']
regions = ['LI', 'NS', 'NA_domain']
periods = ['initial', 'late']
experiments = ['1pctCO2', '1pctCO2-cdr']
window_keys = ['up_initial', 'down_final']
derived_keys = ['up_initial_march_mld_m', 'down_final_march_mld_m']
radius_m = 6371000.0
coordinate_tolerance_m = 1.0
material_change_m = 50.0
registry = {
    'ACCESS-ESM1-5': ('j', 'i', 'latitude', 'longitude', 'j', 'i'),
    'CanESM5': ('j', 'i', 'latitude', 'longitude', 'j', 'i'),
    'CESM2': ('nlat', 'nlon', 'lat', 'lon', 'nlat', 'nlon'),
    'CNRM-ESM2-1': ('y', 'x', 'lat', 'lon', 'y', 'x'),
    'GFDL-ESM4': ('y', 'x', 'lat', 'lon', 'yh', 'xh'),
    'NorESM2-LM': ('j', 'i', 'latitude', 'longitude', 'j', 'i'),
    'UKESM1-0-LL': ('j', 'i', 'latitude', 'longitude', 'j', 'i'),
}
assert not destination.exists(), '正式结果已存在，重跑前先归档。'


# %% 3. 读取已保存温盐MLD及准确年月信息

derived_metadata = json.loads(derived_metadata_path.read_text(encoding='utf-8'))
area_registry = json.loads(area_registry_path.read_text(encoding='utf-8'))
inventory = pd.read_csv(inventory_path)
assert set(models).issubset(set(inventory.model))
destination.mkdir(parents=True)
cell_frames, statistics_rows, regional_rows = [], [], []
month_rows, input_rows, grid_rows = [], [], []
source_paths = [derived_metadata_path, area_registry_path, inventory_path]


# %% 4. 原生网格、相同年月及相同成员配对

for model in models:
    print(model, flush=True)
    y_dim, x_dim, latitude_name, longitude_name, area_y_dim, area_x_dim = registry[model]
    derived_path = derived_directory / (model + '_na_march_mld_three_phase_native.npz')
    source_paths.append(derived_path)
    with np.load(derived_path) as saved:
        latitude = saved['latitude_deg_north'].copy()
        longitude = saved['longitude_deg_east'].copy()
        temperature_salinity_mld = np.stack([saved[key].copy() for key in derived_keys])
        derived_valid = saved['common_valid_mask'].copy()
    assert latitude.shape == longitude.shape == derived_valid.shape
    assert temperature_salinity_mld.shape == (2,) + latitude.shape
    longitude = (longitude + 180.0) % 360.0 - 180.0
    target_coordinates = np.column_stack([latitude.ravel(), longitude.ravel()])
    target_lat_rad = np.deg2rad(target_coordinates[:, 0])
    target_lon_rad = np.deg2rad(target_coordinates[:, 1])
    target_xyz = np.column_stack([
        np.cos(target_lat_rad) * np.cos(target_lon_rad),
        np.cos(target_lat_rad) * np.sin(target_lon_rad), np.sin(target_lat_rad),
    ])
    native_mld = np.full((2,) + latitude.shape, np.nan)
    model_info = derived_metadata['model_details'][model]

    for period_index, (period, experiment, window_key) in enumerate(zip(periods, experiments, window_keys)):
        first_year = model_info[window_key + '_target_march_first'] // 100
        last_year = model_info[window_key + '_target_march_last'] // 100
        expected_months = np.arange(first_year, last_year + 1) * 100 + 3
        assert len(expected_months) == 10
        member_key = 'member_up' if period_index == 0 else 'member_cdr'
        member = model_info[member_key]
        files = inventory.loc[(inventory.model == model) & (inventory.experiment == experiment)
                              & (inventory.variant_label == member), 'file'].tolist()
        assert files
        selected_fields, selected_months = [], []

        for file_text in files:
            file_path = Path(file_text)
            with netCDF4.Dataset(file_path) as dataset:
                variable = dataset.variables['mlotst']
                assert variable.dimensions == ('time', y_dim, x_dim)
                assert variable.units == 'm'
                assert dataset.variant_label == member and dataset.experiment_id == experiment
                time = dataset.variables['time']
                dates = netCDF4.num2date(time[:], time.units, calendar=time.calendar)
                month_ids = np.array([date.year * 100 + date.month for date in dates])
                selected = np.flatnonzero(np.isin(month_ids, expected_months))
                if selected.size == 0:
                    continue

                native_latitude = np.asarray(dataset.variables[latitude_name][:], dtype=float)
                native_longitude = (np.asarray(dataset.variables[longitude_name][:], dtype=float) + 180.0) % 360.0 - 180.0
                assert dataset.variables[latitude_name].dimensions == (y_dim, x_dim)
                assert native_latitude.shape == native_longitude.shape
                native_valid_coordinates = np.isfinite(native_latitude) & np.isfinite(native_longitude)
                lat_rad = np.deg2rad(native_latitude[native_valid_coordinates])
                lon_rad = np.deg2rad(native_longitude[native_valid_coordinates])
                native_xyz = np.column_stack([np.cos(lat_rad) * np.cos(lon_rad),
                                              np.cos(lat_rad) * np.sin(lon_rad), np.sin(lat_rad)])
                distance, coordinate_index = cKDTree(native_xyz).query(target_xyz)
                distance_m = distance * radius_m
                assert float(distance_m.max()) < coordinate_tolerance_m, (model, distance_m.max())
                # 只做同一原生网格坐标的索引配对，不进行空间插值。
                flat_index = np.flatnonzero(native_valid_coordinates)[coordinate_index]
                values = np.asarray(np.ma.filled(variable[selected.tolist(), :, :], np.nan), dtype=float)
                assert not np.any(np.isfinite(values) & ((values < 0) | (values >= 1e10)))
                matched_values = values.reshape(selected.size, -1)[:, flat_index]
                selected_fields.append(matched_values)
                selected_months.extend(month_ids[selected].tolist())
                for month_id in month_ids[selected]:
                    month_rows.append({'model': model, 'period': period, 'experiment': experiment,
                                       'source_year_month': int(month_id), 'member': member, 'file': str(file_path)})
                input_rows.append({'model': model, 'period': period, 'file': str(file_path),
                                   'dimensions': ','.join(variable.dimensions), 'member': member,
                                   'calendar': time.calendar, 'cell_methods': getattr(variable, 'cell_methods', ''),
                                   'comment': getattr(variable, 'comment', ''),
                                   'selected_month_count': int(selected.size),
                                   'maximum_coordinate_distance_m': float(distance_m.max())})

        selected_months = np.array(selected_months)
        assert np.array_equal(np.sort(selected_months), expected_months), (model, period, selected_months)
        monthly_values = np.concatenate(selected_fields, axis=0)[np.argsort(selected_months)]
        assert monthly_values.shape == (10, latitude.size)
        # 与温盐产品相同：10个March等权，任何缺年均不作缺测跳过平均。
        native_mld[period_index] = monthly_values.mean(axis=0).reshape(latitude.shape)

    area_path = Path(area_registry['models_metadata'][model]['areacello_file'])
    source_paths.append(area_path)
    with netCDF4.Dataset(area_path) as dataset:
        area_variable = dataset.variables['areacello']
        assert area_variable.dimensions == (area_y_dim, area_x_dim)
        assert area_variable.units.replace(' ', '').replace('^', '') == 'm2'
        area_values = np.asarray(np.ma.filled(area_variable[:], np.nan), dtype=float)
        area_latitude = np.asarray(dataset.variables[latitude_name][:], dtype=float)
        area_longitude = (np.asarray(dataset.variables[longitude_name][:], dtype=float) + 180.0) % 360.0 - 180.0
        assert area_latitude.shape == area_longitude.shape == area_values.shape
        valid_area_coordinates = np.isfinite(area_latitude) & np.isfinite(area_longitude)
        lat_rad = np.deg2rad(area_latitude[valid_area_coordinates])
        lon_rad = np.deg2rad(area_longitude[valid_area_coordinates])
        area_xyz = np.column_stack([np.cos(lat_rad) * np.cos(lon_rad),
                                    np.cos(lat_rad) * np.sin(lon_rad), np.sin(lat_rad)])
        distance, coordinate_index = cKDTree(area_xyz).query(target_xyz)
        assert float(distance.max() * radius_m) < coordinate_tolerance_m
        flat_area_index = np.flatnonzero(valid_area_coordinates)[coordinate_index]
        area = area_values.ravel()[flat_area_index].reshape(latitude.shape)
        grid_rows.append({'model': model, 'mlotst_dimensions': 'time,' + y_dim + ',' + x_dim,
                          'area_dimensions': area_y_dim + ',' + area_x_dim,
                          'derived_shape': str(latitude.shape), 'area_file': str(area_path),
                          'area_source_id': getattr(dataset, 'source_id', ''),
                          'maximum_area_coordinate_distance_m': float(distance.max() * radius_m)})

    common = derived_valid & np.isfinite(native_mld).all(axis=0)
    common &= np.isfinite(temperature_salinity_mld).all(axis=0) & np.isfinite(area) & (area > 0)
    masks = {
        'LI': (latitude >= 50) & (latitude < 65) & (longitude >= -65) & (longitude < -20),
        'NS': (latitude >= 65) & (latitude < 80) & (longitude >= -30) & (longitude < 20),
        'NA_domain': (latitude >= 50) & (latitude < 85) & (longitude >= -65) & (longitude < 20),
    }
    model_cells = pd.DataFrame({
        'model': model, 'latitude_degN': latitude[common], 'longitude_degE': longitude[common],
        'area_m2': area[common], 'native_initial_m': native_mld[0, common], 'native_late_m': native_mld[1, common],
        'derived_initial_m': temperature_salinity_mld[0, common], 'derived_late_m': temperature_salinity_mld[1, common],
    })
    model_cells['native_delta_m'] = model_cells.native_late_m - model_cells.native_initial_m
    model_cells['derived_delta_m'] = model_cells.derived_late_m - model_cells.derived_initial_m
    for region in regions:
        model_cells['in_' + region] = masks[region][common]
    cell_frames.append(model_cells)


    # %% 5. 同支持区域均值与逐格误差

    for region in regions:
        paired = model_cells.loc[model_cells['in_' + region]]
        assert len(paired) >= 3, (model, region)
        weights = paired.area_m2.to_numpy(copy=True)
        weights /= weights.sum()
        initial_reference = float(np.sum(weights * paired.native_initial_m.to_numpy()))
        base_support = masks[region] & derived_valid & np.isfinite(area) & (area > 0)
        area_fraction = float(paired.area_m2.sum() / area[base_support].sum())
        assert 0 < area_fraction <= 1.0 + 1e-10
        regional_record = {'model': model, 'region': region, 'n_cells': len(paired),
                           'paired_area_fraction_of_derived_support': area_fraction}

        for quantity in ['initial', 'late', 'delta']:
            native_values = paired['native_' + quantity + '_m'].to_numpy()
            derived_values = paired['derived_' + quantity + '_m'].to_numpy()
            native_mean = float(np.sum(weights * native_values))
            derived_mean = float(np.sum(weights * derived_values))
            difference = derived_values - native_values
            covariance = np.sum(weights * (native_values - native_mean) * (derived_values - derived_mean))
            native_variance = np.sum(weights * (native_values - native_mean) ** 2)
            derived_variance = np.sum(weights * (derived_values - derived_mean) ** 2)
            correlation = float(covariance / np.sqrt(native_variance * derived_variance)) if native_variance > 0 and derived_variance > 0 else np.nan
            mae = float(np.sum(weights * np.abs(difference)))
            material = (np.abs(native_values) >= material_change_m) | (np.abs(derived_values) >= material_change_m)
            material_area_fraction = float(weights[material].sum())
            sign_agreement = float(np.sum(weights[material] * (np.sign(native_values[material]) == np.sign(derived_values[material]))) / material_area_fraction) if material_area_fraction > 0 else np.nan
            statistics_rows.append({
                'model': model, 'region': region, 'quantity': quantity, 'n_cells': len(paired),
                'native_mean_m': native_mean, 'derived_mean_m': derived_mean,
                'bias_derived_minus_native_m': derived_mean - native_mean,
                'mae_m': mae, 'rmse_m': float(np.sqrt(np.sum(weights * difference ** 2))),
                'spatial_correlation': correlation,
                'mae_pct_native_initial_mean': 100 * mae / initial_reference,
                'area_fraction_abs_difference_le_50m': float(weights[np.abs(difference) <= 50].sum()),
                'area_fraction_abs_difference_le_100m': float(weights[np.abs(difference) <= 100].sum()),
                'material_area_fraction': material_area_fraction,
                'material_sign_agreement_fraction': sign_agreement,
                'paired_area_fraction_of_derived_support': area_fraction,
            })
            regional_record['native_' + quantity + '_m'] = native_mean
            regional_record['derived_' + quantity + '_m'] = derived_mean
        regional_record['delta_bias_m'] = regional_record['derived_delta_m'] - regional_record['native_delta_m']
        regional_record['delta_bias_pct_native_initial'] = 100 * regional_record['delta_bias_m'] / initial_reference
        regional_record['regional_delta_sign_agrees'] = bool(np.sign(regional_record['native_delta_m']) == np.sign(regional_record['derived_delta_m']))
        regional_rows.append(regional_record)


# %% 6. 保存正式结果与绘图参数

statistics = pd.DataFrame(statistics_rows)
regional = pd.DataFrame(regional_rows)
cells = pd.concat(cell_frames, ignore_index=True)
assert len(statistics) == 63 and len(regional) == 21 and len(month_rows) == 140
cells.to_csv(destination / 'paired_native_cells.csv', index=False)
statistics.to_csv(destination / 'paired_statistics.csv', index=False)
regional.to_csv(destination / 'regional_comparison.csv', index=False)
pd.DataFrame(month_rows).to_csv(destination / 'matched_months.csv', index=False)
pd.DataFrame(input_rows).to_csv(destination / 'mlotst_input_inventory.csv', index=False)
pd.DataFrame(grid_rows).to_csv(destination / 'grid_inventory.csv', index=False)

plot_limits = {}
main_regions = regional.loc[regional.region.isin(['LI', 'NS'])]
for quantity in ['initial', 'late', 'delta']:
    values = main_regions[['native_' + quantity + '_m', 'derived_' + quantity + '_m']].to_numpy()
    limit = float(np.ceil(np.max(np.abs(values)) * 1.10 / 100) * 100)
    plot_limits[quantity] = [-limit, limit] if quantity == 'delta' else [0.0, limit]

source_records = []
for source_path in source_paths:
    source_records.append({'path': str(source_path), 'sha256': hashlib.sha256(source_path.read_bytes()).hexdigest()})
metadata = {
    'models': models, 'excluded_models': {'MIROC-ES2L': 'No local forced mlotst files',
                                        'CAS-ESM2-0': 'No local CDR mlotst or existing comparable derived endpoint product'},
    'periods': 'Same ten March calendar dates, experiment and member as saved thetao-so endpoint product; GFDL initial years6-15.',
    'derived_product': str(derived_directory),
    'derived_method': 'TEOS10 sigma0, reference density interpolated to10m, 0.03kg/m3 first crossing; no crossing assigned deepest valid sampled depth.',
    'native_method': 'Native CMIP6 Omon mlotst; model definitions and internal averaging need not be identical to monthly-thetao-so diagnosis.',
    'time_average': 'Ten March samples equally weighted; require all ten values at a cell for both periods.',
    'grid_matching': 'Exact native-coordinate correspondence, maximum permitted separation1m; no spatial interpolation or regridding.',
    'weighting': 'Native areacello, same complete paired support for both methods and both periods.',
    'regions': {'LI': '50<=lat<65,-65<=lon<-20', 'NS': '65<=lat<80,-30<=lon<20', 'NA_domain': '50<=lat<85,-65<=lon<20'},
    'difference_sign': 'Derived thetao-so minus native mlotst',
    'delta': 'late minus initial for each method',
    'material_change_m': material_change_m,
    'magnitude_bins': '50m and100m are descriptive difference thresholds, not significance tests or acceptance criteria.',
    'scope': 'Endpoint climatologies and their residual changes; this does not validate the separate legacy280yr regional time series.',
    'pi_control_correction': False, 'dimensions': 'paired_cells(model,native_cell);statistics(model7,region3,quantity3);months(model7,period2,March10)',
    'plot_limits': plot_limits, 'sources': source_records,
}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print(main_regions[['model', 'region', 'native_initial_m', 'derived_initial_m', 'native_late_m',
                    'derived_late_m', 'native_delta_m', 'derived_delta_m', 'regional_delta_sign_agrees']].round(2).to_string(index=False))
print('OUTPUT', destination.resolve(), flush=True)
