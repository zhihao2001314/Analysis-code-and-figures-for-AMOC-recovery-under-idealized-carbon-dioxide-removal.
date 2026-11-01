# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import shutil

import numpy as np
import pandas as pd


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
stage = Path('s02_sia/s02_arctic_region/p23_stratification_mld_composite')
destination = root / 'output' / stage / 'n093_compact_regional_hydrography'
li_source = root / 'output/s02_sia/s02_arctic_region/p25_cas_hydrography/n085_nine_model_hydrography_panels'
ns_source = root / 'output/s02_sia/s02_arctic_region/p25_cas_hydrography/n087_nordic_initial_n2_nine_models'
history_source = root / 'output/s02_sia/s02_arctic_region/p25_cas_hydrography/n089_rectangular_mld_nine_models'
prior_source = root / 'output' / stage / 'n091_grouped_hydrography_with_cas'
color_source = root / 'output/s01_amoc/n79_seven_panel_hysteresis_40n/metadata.json'
assert not destination.exists() or not any(destination.iterdir()), '正式输出存在，须先归档。'
regions = ['Labrador-Irminger', 'Nordic Seas']
groups = ['non_recovery', 'recovery', 'larger_bias']
group_names = ['Non-recovery', 'Recovery', 'Larger bias']
members = [['CESM2', 'CNRM-ESM2-1', 'UKESM1-0-LL'],
           ['NorESM2-LM', 'GFDL-ESM4', 'CAS-ESM2-0'],
           ['ACCESS-ESM1-5', 'CanESM5', 'MIROC-ES2L']]
observation_depth_limits = [10., 1000.]
n2_x_ticks = [0, 1, 2, 4, 6, 8]
n2_scale_values = [-1.6, 0., 1., 2., 4., 6., 8., 9.2]
n2_scale_positions = [-.22, 0., 1., 2., 3., 4., 5., 5.6]
contour_levels_m = np.arange(200, 1001, 200)
color_levels_m = np.arange(-1000., 1001., 50.)
color_ticks_m = np.arange(-1000., 1001., 500.)
profile_parts, observation_parts, membership_records, diagnostics = [], [], [], []
source_paths = []


# %% 3. 合并两区九模式剖面及截取观测深度

color_metadata = json.loads(color_source.read_text(encoding='utf-8'))
models, model_colors = color_metadata['models'], color_metadata['model_colors']
for region, directory in zip(regions, [li_source, ns_source]):
    profiles = pd.read_csv(directory / 'initial_n2_profiles.csv', float_precision='round_trip')
    profiles['region'] = region
    profile_parts.append(profiles)
    observations = pd.read_csv(directory / 'iap_initial_n2_reference.csv', float_precision='round_trip')
    observations = observations.loc[observations.depth_m.between(*observation_depth_limits)].copy()
    observations['region'] = region
    observations['lower_1sd_n2_1e5_s2'] = observations.iapv4_march_mean_n2_1e5_s2 - observations.iapv4_march_interannual_sd_n2_1e5_s2
    observations['upper_1sd_n2_1e5_s2'] = observations.iapv4_march_mean_n2_1e5_s2 + observations.iapv4_march_interannual_sd_n2_1e5_s2
    observation_parts.append(observations)
    source_paths.extend([directory / 'initial_n2_profiles.csv', directory / 'iap_initial_n2_reference.csv'])
profiles = pd.concat(profile_parts, ignore_index=True)
observations = pd.concat(observation_parts, ignore_index=True)
assert len(profiles) == 2 * 9 * 100 and len(observations) == 2 * 99
assert not profiles.duplicated(['region', 'model', 'depth_m']).any()
assert observations.depth_m.min() >= 10 and observations.depth_m.max() <= 1000
# 不在10m/1000m人为补点；原观测中间层深度筛选后为15–995m。
assert np.isfinite(observations[['iapv4_march_mean_n2_1e5_s2', 'lower_1sd_n2_1e5_s2', 'upper_1sd_n2_1e5_s2']]).all().all()


# %% 4. 保留两区矩形MLD序列并固定三组成员

history = pd.read_csv(history_source / 'rectangular_mld_timeseries.csv', float_precision='round_trip')
assert len(history) == 2 * 9 * 280
assert set(history.region) == set(regions) and set(history.model) == set(models)
for group, label, group_members in zip(groups, group_names, members):
    for order, model in enumerate(group_members):
        membership_records.append({'group': group, 'label': label, 'model': model,
                                    'member_order': order, 'group_size': 3})
membership = pd.DataFrame(membership_records)
assert membership.model.nunique() == 9


# %% 5. 组内等权平均MLD空间场和SIC场

with np.load(li_source / 'nine_model_mld_ice_display.npz') as saved:
    maps = {name: saved[name].copy() for name in saved.files}
payload = {'longitude': maps['lon2_plot'], 'latitude': maps['lat2_plot'],
           'color_levels_m': color_levels_m, 'color_ticks_m': color_ticks_m,
           'contour_levels_m': contour_levels_m,
           'ice_longitude': maps['CESM2_ice_lon'], 'ice_latitude': maps['CESM2_ice_lat'],
           'observed_x_m': maps['observed_x_m'], 'observed_y_m': maps['observed_y_m'],
           'observed_sic_fraction': maps['observed_sic_fraction']}
for group, group_members in zip(groups, members):
    initial = np.stack([maps[model + '_up_initial_march_mld_m'] for model in group_members])
    delta = np.stack([maps[model + '_delta_later_early_march_mld_m'] for model in group_members])
    valid = np.isfinite(initial) & np.isfinite(delta)
    count = valid.sum(axis=0)
    common = count == 3
    payload[group + '_initial_mld_m'] = np.where(common, initial.mean(axis=0), np.nan)
    payload[group + '_delta_mld_m'] = np.where(common, delta.mean(axis=0), np.nan)
    payload[group + '_later_mld_m'] = payload[group + '_initial_mld_m'] + payload[group + '_delta_mld_m']
    payload[group + '_mld_valid_count'] = count
    for model in group_members:
        np.testing.assert_array_equal(maps[model + '_ice_lon'], payload['ice_longitude'])
        np.testing.assert_array_equal(maps[model + '_ice_lat'], payload['ice_latitude'])
    for period in ['initial', 'later']:
        ice = np.stack([maps[model + '_ice_' + period] for model in group_members])
        ice_count = np.isfinite(ice).sum(axis=0)
        # 先等权平均浓度，再由绘图取15%等值线，不平均冰缘坐标。
        # CAS缺测保留，禁止用0或skipna补成看似完整的三模式场。
        ice_mean = np.where(ice_count == 3, ice.mean(axis=0), np.nan)
        payload[group + '_ice_' + period] = ice_mean
        payload[group + '_ice_' + period + '_valid_count'] = ice_count
        assert np.nanmin(ice_mean) < .15 < np.nanmax(ice_mean)
        assert np.nanmin(ice_mean) >= 0 and np.nanmax(ice_mean) <= 1
    diagnostics.append({'group': group, 'model_count': 3, 'models': '; '.join(group_members),
                        'common_mld_cells': int(common.sum()),
                        'delta_min_m': float(np.nanmin(payload[group + '_delta_mld_m'])),
                        'delta_max_m': float(np.nanmax(payload[group + '_delta_mld_m'])),
                        'initial_ice_common_cells': int(np.isfinite(payload[group + '_ice_initial']).sum()),
                        'later_ice_common_cells': int(np.isfinite(payload[group + '_ice_later']).sum())})


# %% 6. 保存正式输入和底图

destination.mkdir(parents=True, exist_ok=True)
profiles.to_csv(destination / 'member_n2_profiles.csv', index=False)
observations.to_csv(destination / 'iap_profiles_10_1000m.csv', index=False)
history.to_csv(destination / 'member_mld_timeseries.csv', index=False)
membership.to_csv(destination / 'map_membership.csv', index=False)
pd.DataFrame(diagnostics).to_csv(destination / 'map_diagnostics.csv', index=False)
np.savez_compressed(destination / 'grouped_mld_sic_fields.npz', **payload)
for name in ['map_land_50m.geojson', 'map_coastline_50m.geojson', 'mld_rectangular_vertices.csv']:
    shutil.copy2(prior_source / name, destination / name)
    source_paths.append(prior_source / name)


# %% 7. 保存口径与来源

caption = ('Initial March stratification and March mixed-layer depth in the Labrador-Irminger Seas '
           '(a,c) and Nordic Seas (b,d), with nine individual models shown in matching colors and no '
           'model-group mean curves. The piecewise-linear horizontal axes in a,b give equal visual spacing '
           'to0,1,2,4,6,8 without modifying the plotted physical N2 values. Black lines and gray bands show IAP March2005-2024 mean and '
           'interannual standard deviation using only observations within10-1000m; actual retained '
           'mid-depth levels are15-995m, without boundary extrapolation. The MLD curves retain the '
           'saved rectangular-region11yr centred means and missing end windows. The vertical dotted '
           'lines mark CO2 forcing reversal. Panels e-g show residual March MLD (final minus initial '
           'decades) for three fixed groups of3models: non-recovery (CESM2,CNRM-ESM2-1,UKESM1-0-LL), '
           'recovery (NorESM2-LM,GFDL-ESM4,CAS-ESM2-0), and larger bias '
           '(ACCESS-ESM1-5,CanESM5,MIROC-ES2L). Black contours show initial group-mean MLD at '
           '200,400,600,800 and1000m. The filled scale spans-1000 to1000m, with extensions for '
           'values beyond that range and labelled ticks every500m. '
           'Initial and final ice edges are15percent contours of equal-weight group-mean annual SIC '
           'fields; the observed edge is the same NOAA/NSIDC1988-2020 annual climatology in all maps. '
           'All fixed members must have valid data for each averaged cell. No new smoothing, '
           'piControl correction or raw hydrographic computation is introduced.')
source_paths.extend([history_source / 'rectangular_mld_timeseries.csv', li_source / 'nine_model_mld_ice_display.npz',
                     li_source / 'metadata.json', color_source])
metadata = {'models': models, 'model_colors': model_colors, 'regions': regions,
            'groups': groups, 'group_names': group_names, 'members': members, 'group_sizes': [3, 3, 3],
            'panel_count': 7, 'layout': 'Row1 LI/NS N2; row2 LI/NS MLD; row3 three grouped MLD maps with ice edges.',
            'no_group_mean_curves': True, 'n2_x_ticks': n2_x_ticks, 'n2_x_limits': [-1.6, 9.2],
            'n2_axis_transform': {'type': 'piecewise_linear', 'raw_knots': n2_scale_values,
                                  'display_knots': n2_scale_positions, 'display_ticks_equal_spacing': True,
                                  'negative_domain_compressed': True, 'on_panel_axis_note': False},
            'mld_y_limits_m': [0, 700], 'mld_y_ticks_m': list(range(0, 701, 100)),
            'observation_depth_selection_m': observation_depth_limits, 'observation_actual_retained_depth_m': [15., 995.],
            'model_depth_unchanged': True, 'observation_period': 'March2005-2024',
            'regions_definition': 'Rectangles: LI50<=lat<65,-65<=lon<-20; NS65<=lat<80,-30<=lon<20. MLD saved legacy rectangle edge-cell convention retained.',
            'map_extent': [-65., 20., 50., 85.], 'ice_threshold': .15,
            'ice_average': 'Equal-weight3member SIC field first, then15percent contour; all3must be finite, no skipna, no new zero imputation.',
            'observed_ice': 'Reuse saved NOAA/NSIDC annual1988-2020 EPSG3411 SIC, including inherited display smoothing; same field in each panel.',
            'source_limitations': 'Legacy8model ice fields include inherited sparse-table zero restoration; CAS source keeps missing cells. No new repair or regridding. Map versus regional MLD reference treatments remain those of their saved source products.',
            'dimensions': {'n2': 'region2 xmodel9 xdepth100', 'iap': 'region2 xdepth99 after selection',
                           'MLD': 'region2 xmodel9 xyear280', 'map': 'group3 xlat141 xlon341',
                           'SIC': 'group3 xperiod2 xlat35 xlon360', 'observed_SIC': 'projected_y448 xprojected_x304'},
            'caption': caption, 'drift_correction': 'none',
            'source_files': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in source_paths]}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
(destination / 'caption.txt').write_text(caption + '\n', encoding='utf-8')
print(pd.DataFrame(diagnostics).to_string(index=False), flush=True)
print(destination.resolve(), flush=True)
