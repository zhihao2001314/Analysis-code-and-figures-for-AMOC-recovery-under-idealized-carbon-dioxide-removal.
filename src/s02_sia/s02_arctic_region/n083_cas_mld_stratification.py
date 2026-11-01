# %% 1. 导入依赖
from pathlib import Path
import hashlib
import json
import re
import shutil
import warnings

import gsw
import netCDF4
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from shapely import contains_xy
from shapely.geometry import shape

warnings.filterwarnings('ignore', message='Mean of empty slice')


# %% 2. 参数与路径设置
root = Path('F:/sea/Antic sea ice constraint')
sea = Path('F:/sea')
destination = root / 'output/s02_sia/s02_arctic_region/p25_cas_hydrography/n083_cas_mld_stratification'
area_path = sea / 'areacello/CAS/areacello_Ofx_CAS-ESM2-0_1pctCO2_r1i1p1f1_gn.nc'
coast_directory = root / 'output/s02sia/s02_arctic_region/p14_mld_fixed_regions_raw'
geometry_path = coast_directory / 'mld_coast_regions_20261012.geojson'
rect_profile_directory = root / 'output/s02sia/s02_arctic_region/p07_stratification/n01_initial_profiles'
rect_observation_directory = root / 'output/s02sia/s02_arctic_region/p07_stratification/n02_iapv4_observation'
rect_history_path = root / 'output/s02_sia/s02_arctic_region/p21_grouped_mld/n061_grouped_mld_timeseries/member_timeseries.csv'
geography_directory = root / 'output/s02_sia/s02_arctic_region/p23_stratification_mld_composite/n077_regional_definition_comparison'
model = 'CAS-ESM2-0'
member = 'r1i1p1f1'
versions = ['rectangles', 'coastal']
regions = ['Labrador-Irminger', 'Nordic Seas']
experiments = ['1pctCO2', '1pctCO2-cdr']
experiment_starts = [200, 1]
temperature_roots = [sea / 'cmipthetaoad/CAS', sea / 'cmipthetao/CAS']
salinity_roots = [sea / 'cmip6_so' / experiment / model for experiment in experiments]
standard_depth = np.arange(0., 1000.1, 10.)
n2_depth = (standard_depth[:-1] + standard_depth[1:]) / 2
reference_depth = 10.0
density_threshold = 0.03
periods = ['up_initial', 'up_final', 'down_initial', 'down_final']
period_slices = [slice(0, 10), slice(130, 140), slice(140, 150), slice(270, 280)]
assert not destination.exists(), '正式CAS结果已存在，重跑前须先归档。'
destination.mkdir(parents=True)
source_records, inventory, month_records, samples = [], [], [], []


# %% 3. 核对原生面积与两段实验的温盐文件
with netCDF4.Dataset(area_path) as dataset:
    assert dataset['areacello'].dimensions == ('lat', 'lon')
    assert dataset['areacello'].units == 'm2'
    latitude_full = np.asarray(dataset['lat'][:], dtype=float)
    longitude_full = np.asarray(dataset['lon'][:], dtype=float)
    area_full = np.ma.filled(dataset['areacello'][:], np.nan).astype(float)
longitude_signed = (longitude_full + 180) % 360 - 180
latitude_indices = np.flatnonzero((latitude_full >= 50) & (latitude_full <= 85))
longitude_indices = np.flatnonzero((longitude_signed >= -65) & (longitude_signed <= 25))
longitude_indices = longitude_indices[np.argsort(longitude_signed[longitude_indices])]
latitude_slice = slice(latitude_indices[0], latitude_indices[-1] + 1)
latitude = latitude_full[latitude_indices]
longitude = longitude_signed[longitude_indices]
longitude_grid, latitude_grid = np.meshgrid(longitude, latitude)
area = area_full[np.ix_(latitude_indices, longitude_indices)]
wet_area = np.isfinite(area) & (area > 0)
ny, nx = area.shape
month_tables = {}
native_depth = None
for experiment, start_year, temperature_root, salinity_root in zip(
        experiments, experiment_starts, temperature_roots, salinity_roots):
    for variable, data_root in [('thetao', temperature_root), ('so', salinity_root)]:
        pattern = rf'{variable}_Omon_CAS-ESM2-0_{experiment}_{member}_gn_\d{{6}}-\d{{6}}\.nc'
        paths = sorted(path for path in data_root.glob('*.nc') if re.fullmatch(pattern, path.name))
        assert paths, (experiment, variable)
        rows = []
        for path in paths:
            with netCDF4.Dataset(path) as dataset:
                field = dataset[variable]
                assert field.dimensions == ('time', 'lev', 'lat', 'lon')
                assert dataset.source_id == model and dataset.experiment_id == experiment
                assert dataset.variant_label == member and dataset.grid_label == 'gn'
                np.testing.assert_array_equal(dataset['lat'][:], latitude_full)
                np.testing.assert_array_equal(dataset['lon'][:], longitude_full)
                z = np.asarray(dataset['lev'][:], dtype=float)
                assert dataset['lev'].units == 'm' and np.all(np.diff(z) > 0)
                if native_depth is None:
                    native_depth = z
                np.testing.assert_array_equal(z, native_depth)
                assert field.units in (['degC'] if variable == 'thetao' else ['0.001', 'psu', '1e-3', 'g kg-1'])
                time = dataset['time']
                dates = netCDF4.num2date(time[:], time.units, time.calendar)
                serial = np.array([date.year * 12 + date.month for date in dates])
                assert np.all(np.diff(serial) == 1)
                for index, date in enumerate(dates):
                    if date.month == 3 and start_year <= date.year < start_year + 140:
                        rows.append({'experiment': experiment, 'variable': variable, 'calendar_year': date.year,
                                     'year_month': date.year * 100 + 3, 'local_index': index, 'path': str(path)})
                inventory.append({'path': str(path), 'variable': variable, 'experiment': experiment,
                                  'dimensions': str(field.dimensions), 'shape': str(field.shape),
                                  'units': field.units, 'first_month': str(dates[0]), 'last_month': str(dates[-1]),
                                  'bytes': path.stat().st_size, 'mtime_ns': path.stat().st_mtime_ns})
        table = pd.DataFrame(rows).sort_values('year_month').reset_index(drop=True)
        np.testing.assert_array_equal(table.calendar_year, np.arange(start_year, start_year + 140))
        assert len(table) == 140 and not table.year_month.duplicated().any()
        month_tables[(experiment, variable)] = table
    np.testing.assert_array_equal(month_tables[(experiment, 'thetao')].year_month,
                                  month_tables[(experiment, 'so')].year_month)
    print(f'{experiment}: matched 140 March temperature/salinity records', flush=True)
assert native_depth[0] == 5 and native_depth[1] == 15


# %% 4. 保存分区并准备旧八模式参照
features = json.loads(geometry_path.read_text(encoding='utf-8'))['features']
coastal_geometry = {item['properties']['region']: shape(item['geometry']) for item in features}
masks = {
    ('rectangles', regions[0]): (latitude_grid >= 50) & (latitude_grid < 65) & (longitude_grid >= -65) & (longitude_grid <= -20),
    ('rectangles', regions[1]): (latitude_grid >= 65) & (latitude_grid <= 80) & (longitude_grid >= -30) & (longitude_grid <= 20),
    ('coastal', regions[0]): contains_xy(coastal_geometry[regions[0]], longitude_grid, latitude_grid),
    ('coastal', regions[1]): contains_xy(coastal_geometry[regions[1]], longitude_grid, latitude_grid),
}
# N²矩形沿用半开边界；旧MLD矩形包含部分外边界，两种掩膜分别保留。
n2_masks = dict(masks)
n2_masks[('rectangles', regions[0])] = (latitude_grid >= 50) & (latitude_grid < 65) & (longitude_grid >= -65) & (longitude_grid < -20)
n2_masks[('rectangles', regions[1])] = (latitude_grid >= 65) & (latitude_grid < 80) & (longitude_grid >= -30) & (longitude_grid < 20)
region_keys = list(masks)
for key in region_keys:
    assert np.any(masks[key] & wet_area) and np.any(n2_masks[key] & wet_area)
reference_profile_frames, reference_history_frames = [], []
for version in versions:
    profile_path = (rect_profile_directory / 'na_initial_march_stratification_profiles.csv' if version == 'rectangles'
                    else coast_directory / 'coast_region_20261012_initial_n2_profiles.csv')
    observation_path = (rect_observation_directory / 'iapv4_march_2005_2024_stratification_profiles.csv'
                        if version == 'rectangles' else coast_directory / 'coast_region_20261012_iapv4_profiles.csv')
    history_path = rect_history_path if version == 'rectangles' else coast_directory / 'coast_region_20261012_plot_data.csv'
    profiles = pd.read_csv(profile_path)
    profiles = profiles.loc[profiles.model != 'MMM', ['model', 'region', 'depth_m', 'initial_march_n2_s2']].copy()
    profiles = profiles.rename(columns={'initial_march_n2_s2': 'n2_s2'})
    observations = pd.read_csv(observation_path)[['region', 'depth_m', 'iapv4_march_mean_n2_s2']].copy()
    observations = observations.rename(columns={'iapv4_march_mean_n2_s2': 'n2_s2'})
    observations['model'] = 'IAP'
    reference = pd.concat([profiles, observations], ignore_index=True)
    reference['version'] = version
    reference['region'] = reference.region.replace({'Labrador-Irminger Seas': regions[0]})
    reference_profile_frames.append(reference)
    history = pd.read_csv(history_path)
    if version == 'coastal':
        history = history.rename(columns={'march_mld_mean_crossing_only_m': 'march_mld_m',
                                          'march_mld_mean_crossing_only_m_11yr': 'march_mld_11yr_m'})
    history = history[['model', 'region', 'simulation_year', 'march_mld_m', 'march_mld_11yr_m']].copy()
    history['version'] = version
    reference_history_frames.append(history)
    for path in [profile_path, observation_path, history_path]:
        source_records.append({'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
reference_profiles = pd.concat(reference_profile_frames, ignore_index=True)
reference_history = pd.concat(reference_history_frames, ignore_index=True)
assert reference_profiles.loc[reference_profiles.model != 'IAP', 'model'].nunique() == 8
for path in [geometry_path, area_path]:
    source_records.append({'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
for name in ['map_land_50m.geojson', 'map_coastline_50m.geojson']:
    shutil.copy2(geography_directory / name, destination / name)
shutil.copy2(geometry_path, destination / geometry_path.name)


# %% 5. 逐年诊断原生密度混合层与共同深度N²
nz = len(native_depth)
depth = native_depth[:, None]
lat_flat = latitude_grid.ravel()[None, :]
lon_flat = longitude_grid.ravel()[None, :]
pressure_native = gsw.p_from_z(-depth, lat_flat)
pressure_standard = gsw.p_from_z(-standard_depth[:, None], lat_flat)
columns = np.arange(ny * nx)
mld_all = np.full((280, ny * nx), np.nan)
mld_crossing = np.full_like(mld_all, np.nan)
mld_nearest_crossing = np.full_like(mld_all, np.nan)
no_cross_flags = np.zeros(mld_all.shape, dtype=bool)
regional_n2 = np.full((280, len(region_keys), len(n2_depth)), np.nan)
regional_n2_coverage = np.full_like(regional_n2, np.nan)
n2_period_sum = np.zeros((4, len(n2_depth), ny * nx))
n2_period_count = np.zeros(n2_period_sum.shape, dtype=np.int16)
history_rows = []
for experiment_index, experiment in enumerate(experiments):
    t_table = month_tables[(experiment, 'thetao')]
    s_table = month_tables[(experiment, 'so')]
    for year_index in range(140):
        simulation_index = experiment_index * 140 + year_index
        t_record, s_record = t_table.iloc[year_index], s_table.iloc[year_index]
        with netCDF4.Dataset(t_record.path) as dataset:
            pt = np.ma.filled(dataset['thetao'][int(t_record.local_index), :, latitude_slice, :], np.nan)
            pt = np.asarray(pt[:, :, longitude_indices], dtype=float).reshape(nz, -1)
        with netCDF4.Dataset(s_record.path) as dataset:
            salinity = np.ma.filled(dataset['so'][int(s_record.local_index), :, latitude_slice, :], np.nan)
            salinity = np.asarray(salinity[:, :, longitude_indices], dtype=float).reshape(nz, -1)
        valid = np.isfinite(pt) & np.isfinite(salinity) & (salinity >= 0) & wet_area.ravel()[None, :]
        # CAS固定z网格的有效层必须从表层连续向下，允许海底以下缺测。
        assert not np.any(np.diff(valid.astype(int), axis=0) > 0), (experiment, year_index, 'internal vertical gap')
        pt = np.where(valid, pt, np.nan)
        salinity = np.where(valid, salinity, np.nan)
        sa = gsw.SA_from_SP(salinity, pressure_native, lon_flat, lat_flat)
        ct = gsw.CT_from_pt(sa, pt)
        sigma0 = gsw.sigma0(sa, ct)
        # 10m正好在CAS的5/15m中心之间；按原空间图方法插值参考密度。
        sigma_reference = (sigma0[0] + sigma0[1]) / 2
        anomaly = sigma0 - sigma_reference[None, :]
        crossing = valid & (depth >= reference_depth) & (anomaly >= density_threshold)
        has_crossing = crossing.any(axis=0)
        k2 = crossing.argmax(axis=0)
        k1 = np.maximum(k2 - 1, 0)
        lower_depth = np.maximum(native_depth[k1], reference_depth)
        lower_anomaly = np.where(native_depth[k1] < reference_depth, 0.0, anomaly[k1, columns])
        upper_anomaly = anomaly[k2, columns]
        good = has_crossing & np.isfinite(sigma_reference) & (upper_anomaly > lower_anomaly)
        mld_crossing[simulation_index, good] = lower_depth[good] + (
            (density_threshold - lower_anomaly[good]) / (upper_anomaly[good] - lower_anomaly[good])
            * (native_depth[k2[good]] - lower_depth[good]))
        deepest = np.max(np.where(valid, depth, -np.inf), axis=0)
        no_cross = np.isfinite(sigma_reference) & ~has_crossing & np.isfinite(deepest)
        mld_all[simulation_index] = mld_crossing[simulation_index]
        mld_all[simulation_index, no_cross] = deepest[no_cross]
        no_cross_flags[simulation_index] = no_cross
        # 旧矩形MLD序列以最近10m原生层为参照；CAS等距时选择第一个5m层。
        nearest_anomaly = sigma0 - sigma0[0][None, :]
        nearest_cross = valid & (nearest_anomaly >= density_threshold)
        nearest_has_cross = nearest_cross.any(axis=0)
        nearest_k2 = nearest_cross.argmax(axis=0)
        nearest_k1 = np.maximum(nearest_k2 - 1, 0)
        denominator = sigma0[nearest_k2, columns] - sigma0[nearest_k1, columns]
        good_nearest = nearest_has_cross & np.isfinite(denominator) & (denominator > 0)
        mld_nearest_crossing[simulation_index, good_nearest] = native_depth[nearest_k1[good_nearest]] + (
            (sigma0[0, good_nearest] + density_threshold - sigma0[nearest_k1[good_nearest], columns[good_nearest]])
            / denominator[good_nearest] * (native_depth[nearest_k2[good_nearest]] - native_depth[nearest_k1[good_nearest]]))
        pt_standard = interp1d(native_depth, pt, axis=0, bounds_error=False, fill_value=np.nan,
                               assume_sorted=True)(standard_depth)
        sp_standard = interp1d(native_depth, salinity, axis=0, bounds_error=False, fill_value=np.nan,
                               assume_sorted=True)(standard_depth)
        sa_standard = gsw.SA_from_SP(sp_standard, pressure_standard, lon_flat, lat_flat)
        ct_standard = gsw.CT_from_pt(sa_standard, pt_standard)
        n2, _ = gsw.Nsquared(sa_standard, ct_standard, pressure_standard, lat=lat_flat, axis=0)
        for period_index, window in enumerate(period_slices):
            if window.start <= simulation_index < window.stop:
                n2_period_sum[period_index] += np.where(np.isfinite(n2), n2, 0)
                n2_period_count[period_index] += np.isfinite(n2)
        for region_index, (version, region) in enumerate(region_keys):
            selected = masks[(version, region)].ravel() & wet_area.ravel()
            weights = area.ravel()[selected]
            record = {'model': model, 'version': version, 'region': region, 'simulation_year': simulation_index + 1,
                      'experiment': experiment, 'experiment_year': year_index + 1, 'year_month': int(t_record.year_month)}
            for name, field in [('all', mld_all[simulation_index]), ('crossing', mld_crossing[simulation_index]),
                                ('nearest_crossing', mld_nearest_crossing[simulation_index])]:
                values = field[selected]
                support = np.isfinite(values)
                record['mld_' + name + '_m'] = np.sum(values[support] * weights[support]) / weights[support].sum()
            record['no_cross_area_fraction'] = np.sum(weights * no_cross[selected]) / weights.sum()
            history_rows.append(record)
            selected_n2 = n2_masks[(version, region)].ravel() & wet_area.ravel()
            weights_n2 = area.ravel()[selected_n2]
            values = n2[:, selected_n2]
            valid_n2 = np.isfinite(values)
            denominator = (valid_n2 * weights_n2).sum(axis=1)
            regional_n2[simulation_index, region_index] = np.divide(
                np.nansum(values * weights_n2, axis=1), denominator,
                out=np.full(len(n2_depth), np.nan), where=denominator > 0)
            regional_n2_coverage[simulation_index, region_index] = denominator / weights_n2.sum()
        month_records.append({'simulation_year': simulation_index + 1, 'experiment': experiment,
                              'year_month': int(t_record.year_month), 'thetao_file': t_record.path,
                              'so_file': s_record.path, 'thetao_index': int(t_record.local_index),
                              'so_index': int(s_record.local_index)})
        if simulation_index in [0, 130, 140, 270]:
            sample_candidates = np.flatnonzero(valid[1] & np.isfinite(sigma_reference))
            selected_columns = sample_candidates[np.linspace(0, len(sample_candidates) - 1, 12, dtype=int)]
            samples.append({'simulation_index': simulation_index, 'columns': selected_columns,
                            'pt': pt[:, selected_columns], 'sp': salinity[:, selected_columns],
                            'n2': n2[:, selected_columns]})
        if year_index % 20 == 0 or year_index == 139:
            print(f'{experiment}: March {year_index + 1}/140; simulation year {simulation_index + 1}', flush=True)


# %% 6. 汇总端点、完整居中平滑与独立样本
history = pd.DataFrame(history_rows)
for (version, region), indices in history.groupby(['version', 'region'], sort=False).groups.items():
    for column in ['mld_all_m', 'mld_crossing_m', 'mld_nearest_crossing_m']:
        history.loc[indices, column + '_11yr'] = history.loc[indices, column].rolling(11, center=True, min_periods=11).mean().to_numpy()
profile_rows, summary_rows, native_mld_periods = [], [], []
for period_index, (period, window) in enumerate(zip(periods, period_slices)):
    valid_count = np.isfinite(mld_all[window]).sum(axis=0)
    native_mld_periods.append(np.where(valid_count == 10, np.mean(mld_all[window], axis=0), np.nan))
    for region_index, (version, region) in enumerate(region_keys):
        matrix = regional_n2[window, region_index]
        complete = np.isfinite(matrix).all(axis=0)
        means = np.where(complete, np.mean(matrix, axis=0), np.nan)
        sd = np.where(complete, np.std(matrix, axis=0, ddof=1), np.nan)
        for depth_index, depth_m in enumerate(n2_depth):
            profile_rows.append({'model': model, 'version': version, 'region': region, 'period': period,
                                 'depth_m': depth_m, 'n2_s2': means[depth_index], 'interannual_sd_s2': sd[depth_index],
                                 'valid_march_count': int(np.isfinite(matrix[:, depth_index]).sum())})
        region_history = history.loc[(history.version == version) & (history.region == region)]
        record = {'version': version, 'region': region, 'period': period,
                  'simulation_year_first': window.start + 1, 'simulation_year_last': window.stop}
        for column in ['mld_all_m', 'mld_crossing_m', 'mld_nearest_crossing_m', 'no_cross_area_fraction']:
            record[column] = region_history.iloc[window][column].mean()
        summary_rows.append(record)
native_mld_periods = np.stack(native_mld_periods).reshape(4, ny, nx)
map_domain = (latitude_grid < 85) & (longitude_grid < 20)
common = np.isfinite(native_mld_periods).all(axis=0) & map_domain
native_mld_periods[:, ~common] = np.nan
n2_native_periods = np.divide(n2_period_sum, n2_period_count, out=np.full_like(n2_period_sum, np.nan), where=n2_period_count == 10)


# %% 7. 保存正式结果及来源
history.to_csv(destination / 'cas_march_mld_timeseries.csv', index=False)
pd.DataFrame(profile_rows).to_csv(destination / 'cas_march_n2_profiles.csv', index=False)
pd.DataFrame(summary_rows).to_csv(destination / 'cas_endpoint_summary.csv', index=False)
pd.DataFrame(inventory).to_csv(destination / 'input_inventory.csv', index=False)
pd.DataFrame(month_records).to_csv(destination / 'matched_marches.csv', index=False)
reference_profiles.to_csv(destination / 'reference_initial_n2_profiles.csv', index=False)
reference_history.to_csv(destination / 'reference_mld_timeseries.csv', index=False)
np.savez_compressed(destination / 'cas_native_mld.npz', latitude=latitude, longitude=longitude, area_m2=area,
                    native_depth_m=native_depth, mld_all_m=mld_all.reshape(280, ny, nx),
                    mld_crossing_m=mld_crossing.reshape(280, ny, nx),
                    mld_nearest_crossing_m=mld_nearest_crossing.reshape(280, ny, nx),
                    no_cross=no_cross_flags.reshape(280, ny, nx),
                    period_mld_m=native_mld_periods, common_map_mask=common,
                    delta_up_m=native_mld_periods[1] - native_mld_periods[0],
                    delta_down_m=native_mld_periods[3] - native_mld_periods[2],
                    delta_late_early_m=native_mld_periods[3] - native_mld_periods[0])
np.savez_compressed(destination / 'cas_regional_n2.npz', depth_m=n2_depth, n2_s2=regional_n2,
                    valid_area_fraction=regional_n2_coverage,
                    native_period_n2_s2=n2_native_periods.reshape(4, len(n2_depth), ny, nx),
                    region_masks=np.stack([masks[key] for key in region_keys]),
                    n2_region_masks=np.stack([n2_masks[key] for key in region_keys]))
sample_arrays = {}
for sample in samples:
    for name in ['columns', 'pt', 'sp', 'n2']:
        sample_arrays[f'year{sample["simulation_index"] + 1}_{name}'] = sample[name]
np.savez_compressed(destination / 'independent_raw_samples.npz', **sample_arrays)
metadata = {'model': model, 'member': member, 'analysis_month': 3, 'years': 280,
            'source_calendar_windows': {'1pctCO2': [200, 339], '1pctCO2-cdr': [1, 140]},
            'dimensions': {'thetao_so': '(time,lev,lat,lon); CAS native 30x196x360',
                           'area': '(lat,lon)', 'saved_mld': [280, ny, nx],
                           'saved_regional_n2': [280, 4, 100], 'saved_native_period_n2': [4, 100, ny, nx]},
            'area_coordinate_match': 'exact at every native coordinate in both experiments and both variables',
            'mld_map_and_coastal_reference': 'sigma0 threshold0.03 kg/m3; sigma0 interpolated to10m; first downward crossing interpolated; no-cross deepest valid depth retained only in all/mapped MLD and flagged',
            'mld_rectangular_legacy_reference': 'nearest valid native level to10m; CAS5m; no-cross excluded; used only for comparison with the old rectangle time series',
            'n2_method': 'T/S linearly interpolated to0..1000m every10m, no extrapolation; SA/CT; gsw.Nsquared at5..995m; local computation precedes native-area averaging',
            'n2_missing': '5m output lacks vertical support and remains missing; negative N2 retained',
            'periods': dict(zip(periods, [[1, 10], [131, 140], [141, 150], [271, 280]])),
            'region_order': [list(key) for key in region_keys],
            'rectangular_mask_difference': 'N2 uses half-open65N/80N and longitude boundaries; MLD matches legacy inclusive outer nodes; saved separately',
            'coastal_regions': 'unchanged user-defined mld_coast_regions_20261012.geojson polygons',
            'smoothing': '11-year centered MLD display only; full11 samples; no endpoint fill or N2 smoothing',
            'drift_correction': 'none', 'other_models': 'saved references only, not recalculated or regrouped',
            'reference_time_caution': 'GFDL initial N2 years6-15 retained; legacy8-model MLD alignment limitations retained; CAS actual March labels fully verified',
            'source_files': source_records, 'raw_inventory': inventory}
(destination / 'metadata.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
print(destination, flush=True)
