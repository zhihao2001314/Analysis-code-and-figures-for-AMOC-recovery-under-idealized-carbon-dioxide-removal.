# %% 1. 导入依赖
from pathlib import Path
import hashlib
import json
import shutil

import numpy as np
import pandas as pd
from scipy.interpolate import griddata


# %% 2. 参数与路径设置
root = Path('F:/sea/Antic sea ice constraint')
cas_directory = root / 'output/s02_sia/s02_arctic_region/p25_cas_hydrography/n083_cas_mld_stratification'
destination = root / 'output/s02_sia/s02_arctic_region/p25_cas_hydrography/n085_nine_model_hydrography_panels'
legacy_map_directory = root / 'output/s02sia/s02_arctic_region/p06_mld/n058_observed_ice_overlay'
legacy_profile_directory = root / 'output/s02sia/s02_arctic_region/p07_stratification/n01_initial_profiles'
observation_directory = root / 'output/s02sia/s02_arctic_region/p07_stratification/n02_iapv4_observation'
coast_directory = root / 'output/s02sia/s02_arctic_region/p14_mld_fixed_regions_raw'
ice_directory = root / 'output/s02sia/s02_arctic_region/p04_pattern_unadjusted/n01_arctic_sia_pattern'
old_models = ['ACCESS-ESM1-5', 'CanESM5', 'CESM2', 'CNRM-ESM2-1', 'GFDL-ESM4', 'MIROC-ES2L', 'UKESM1-0-LL', 'NorESM2-LM']
models = old_models + ['CAS-ESM2-0']
profile_path = legacy_profile_directory / 'na_initial_march_stratification_profiles.csv'
observation_path = observation_directory / 'iapv4_march_2005_2024_stratification_profiles.csv'
history_path = coast_directory / 'coast_region_20261012_plot_data.csv'
co2_path = coast_directory / 'na_march_mld_coast_regions_20261012_regions_co2.csv'
map_path = legacy_map_directory / 'mld_ice_display.npz'
ice_path = ice_directory / 'arctic_sia_pattern_cells.csv'
assert not destination.exists(), '已有正式九模式绘图输入须先归档。'
destination.mkdir(parents=True)
sources = [profile_path, observation_path, history_path, co2_path, map_path, ice_path,
           ice_directory / 'arctic_sia_pattern_metadata.json',
           legacy_map_directory / 'metadata.json',
           cas_directory / 'cas_march_n2_profiles.csv', cas_directory / 'cas_march_mld_timeseries.csv',
           cas_directory / 'cas_native_mld.npz', cas_directory / 'metadata.json']


# %% 3. 九个单模式初始剖面与同一IAP参照
profiles = pd.read_csv(profile_path)
profiles = profiles.loc[(profiles.region == 'Labrador-Irminger Seas') & profiles.model.isin(old_models),
                         ['model', 'depth_m', 'initial_march_n2_s2']].copy()
profiles = profiles.rename(columns={'initial_march_n2_s2': 'n2_s2'})
cas_profiles = pd.read_csv(cas_directory / 'cas_march_n2_profiles.csv')
cas_profiles = cas_profiles.loc[(cas_profiles.version == 'rectangles') &
                               (cas_profiles.region == 'Labrador-Irminger') &
                               (cas_profiles.period == 'up_initial'), ['model', 'depth_m', 'n2_s2']]
profiles = pd.concat([profiles, cas_profiles], ignore_index=True)
profiles['n2_1e5_s2'] = profiles.n2_s2 * 1e5
observations = pd.read_csv(observation_path)
observations = observations.loc[observations.region == 'Labrador-Irminger Seas'].copy()
assert len(profiles) == 900 and len(observations) == 100
assert not profiles.duplicated(['model', 'depth_m']).any()
profiles.to_csv(destination / 'initial_n2_profiles.csv', index=False)
observations.to_csv(destination / 'iap_initial_n2_reference.csv', index=False)


# %% 4. 九模式海岸分区绝对MLD序列
history = pd.read_csv(history_path)
history = history.loc[history.model.isin(old_models),
                      ['model', 'region', 'simulation_year', 'march_mld_mean_crossing_only_m',
                       'march_mld_mean_crossing_only_m_11yr']].copy()
cas_history = pd.read_csv(cas_directory / 'cas_march_mld_timeseries.csv')
cas_history = cas_history.loc[cas_history.version == 'coastal',
                              ['model', 'region', 'simulation_year', 'mld_crossing_m', 'mld_crossing_m_11yr']].copy()
cas_history = cas_history.rename(columns={'mld_crossing_m': 'march_mld_mean_crossing_only_m',
                                          'mld_crossing_m_11yr': 'march_mld_mean_crossing_only_m_11yr'})
history = pd.concat([history, cas_history], ignore_index=True)
assert len(history) == 9 * 2 * 280
assert not history.duplicated(['model', 'region', 'simulation_year']).any()
history.to_csv(destination / 'coastal_mld_timeseries.csv', index=False)
shutil.copy2(co2_path, destination / 'atmospheric_co2.csv')
plot_settings = json.loads((coast_directory / 'coast_region_20261012_plot_metadata.json').read_text(encoding='utf-8'))
profile_settings = json.loads((observation_directory / 'iapv4_march_2005_2024_stratification_metadata.json').read_text(encoding='utf-8'))


# %% 5. 原样保留八模式地图显示数组并加入CAS
with np.load(map_path) as saved:
    payload = {name: saved[name].copy() for name in saved.files}
with np.load(cas_directory / 'cas_native_mld.npz') as saved:
    cas_longitude, cas_latitude = np.meshgrid(saved['longitude'], saved['latitude'])
    valid = saved['common_map_mask'] & np.isfinite(saved['period_mld_m'][0]) & np.isfinite(saved['delta_late_early_m'])
    points = np.column_stack([cas_longitude[valid], cas_latitude[valid]])
    targets = np.column_stack([payload['lon2_plot'].ravel(), payload['lat2_plot'].ravel()])
    for field_name, values in [('up_initial_march_mld_m', saved['period_mld_m'][0]),
                                ('delta_later_early_march_mld_m', saved['delta_late_early_m'])]:
        payload['CAS-ESM2-0_' + field_name] = griddata(points, values[valid], targets, method='linear').reshape(payload['lon2_plot'].shape)
    np.savez_compressed(destination / 'cas_native_map_source.npz', latitude=cas_latitude, longitude=cas_longitude,
                        valid=valid, initial_mld_m=saved['period_mld_m'][0], delta_mld_m=saved['delta_late_early_m'])


# %% 6. 加入CAS既有初末期年平均海冰场，保留原有缺测
ice = pd.read_csv(ice_path)
cas_ice = ice.loc[ice.model == 'CAS-ESM2-0'].copy()
ice_latitude = np.arange(50.5, 85, 1.)
ice_longitude = np.arange(-179.5, 180, 1.)
ice_longitude_grid, ice_latitude_grid = np.meshgrid(ice_longitude, ice_latitude)
payload['CAS-ESM2-0_ice_lon'] = ice_longitude_grid
payload['CAS-ESM2-0_ice_lat'] = ice_latitude_grid
for column, period in [('edge_raw_co2_sic_fraction', 'initial'), ('edge_raw_cdr_sic_fraction', 'later')]:
    field = cas_ice.pivot(index='latitude_deg_north', columns='longitude_deg_east', values=column)
    field = field.reindex(index=ice_latitude, columns=ice_longitude).to_numpy()
    assert field.shape == (35, 360) and np.nanmin(field) >= 0 and np.nanmax(field) <= 1
    assert np.nanmin(field) < .15 < np.nanmax(field)
    payload['CAS-ESM2-0_ice_' + period] = field
np.savez_compressed(destination / 'nine_model_mld_ice_display.npz', **payload)
for name in ['map_land_50m.geojson', 'map_coastline_50m.geojson']:
    source = cas_directory / name
    shutil.copy2(source, destination / name)
    sources.append(source)


# %% 7. 保存口径与来源
sources.extend([coast_directory / 'coast_region_20261012_plot_metadata.json',
                observation_directory / 'iapv4_march_2005_2024_stratification_metadata.json'])
metadata = {
    'models': models, 'n2_region': 'rectangular Labrador-Irminger:50<=lat<65,-65<=lon<-20',
    'n2_initial': 'model initial10 March samples; original GFDL6-15 exception retained; CAS200-209',
    'observation': 'IAP March2005-2024, repeated unchanged in every individual model panel; not contemporaneous validation',
    'n2_x_limits': profile_settings['plot_x_limits_n2_1e5_s2'],
    'mld_region': 'original coastal polygons;10m interpolated reference sigma0;0.03kg/m3 crossing-only;11-year centered means; no endpoint filling',
    'mld_y_limits': [plot_settings['plot_y_min_m'], plot_settings['plot_y_max_m']],
    'co2_y_limits': [plot_settings['co2_plot_y_min_ppm'], plot_settings['co2_plot_y_max_ppm']],
    'co2_ticks': plot_settings['co2_plot_ticks_ppm'],
    'map': 'late10 minus initial10 March MLD; contours initial MLD; old8 display arrays and color scale unchanged; CAS linear griddata mapped to identical0.25degree display grid',
    'map_ice': 'initial/final annual mean15percent SIC; old8 exact saved fields retained, including legacy sparse zero restoration; newCAS reads current complete1degree plot-ready product without replacing NaN with zero',
    'cas_ice_source_note': 'current saved SIC product has native0:120 up and1560:1680 down monthly windows,8-neighbour spherical IDW display; old8 inherited display fields retained to match source figure',
    'observed_ice': 'unchanged original EPSG3411 25km1988-2020 annual mean display field; inherited sigma0.8cell display smoothing, no new smoothing',
    'dimensions': {'model_n2': [9, 100], 'observation_n2': [100], 'mld_timeseries': [9, 2, 280],
                   'map_display': list(payload['lon2_plot'].shape), 'cas_ice': [35, 360]},
    'preservation': 'other8 scientific and display arrays unchanged; no regrouping or raw data recomputation; legacy8 GFDL MLD time-alignment limitation remains',
    'source_files': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in sources]
}
(destination / 'metadata.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
print(destination, flush=True)
