# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import pandas as pd


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
stage = Path('s02_sia/s02_arctic_region')
n2_directory = root / 'output' / stage / 'p20_grouped_initial_stratification/n059_grouped_initial_n2_profiles'
mld_directory = root / 'output' / stage / 'p21_grouped_mld/n061_grouped_mld_timeseries'
map_directory = root / 'output' / stage / 'p06_mld/n071_model_observation_panels'
registry_path = root / 'output/s05_atmosphere/p05_amoc_recovery_composites/n49_group_naming/group_labels.json'
destination = root / 'output' / stage / 'p23_stratification_mld_composite/n074_grouped_mld_patterns'
assert not destination.exists(), '正式计算结果存在，重跑前先归档。'
groups = ['better_nonrecovery', 'better_recovery', 'larger_bias']
map_extent = [-65., 20., 50., 85.]
region_n2, region_mld = 'Labrador-Irminger Seas', 'Labrador-Irminger'
map_shape = (141, 341)
source_files = [registry_path]


# %% 3. 读取既有剖面、序列和共同显示网格

membership = pd.read_csv(n2_directory / 'membership.csv')
pd.testing.assert_frame_equal(membership, pd.read_csv(mld_directory / 'membership.csv'))
registry = json.loads(registry_path.read_text(encoding='utf-8'))
with np.load(map_directory / 'mld_model_observation_display.npz') as saved:
    longitude, latitude = saved['lon2_plot'].copy(), saved['lat2_plot'].copy()
    assert longitude.shape == latitude.shape == map_shape
    initial_fields, later_fields, delta_fields = {}, {}, {}
    for model in membership.model:
        initial_fields[model] = saved[model + '_up_initial_march_mld_m'].copy()
        later_fields[model] = saved[model + '_down_final_march_mld_m'].copy()
        delta_fields[model] = saved[model + '_delta_later_early_march_mld_m'].copy()
        assert initial_fields[model].shape == later_fields[model].shape == delta_fields[model].shape == map_shape
        np.testing.assert_allclose(later_fields[model] - initial_fields[model], delta_fields[model],
                                   rtol=1e-10, atol=1e-9, equal_nan=True)
    color_levels, color_ticks = saved['color_levels'].copy(), saved['color_ticks'].copy()
assert json.loads((n2_directory / 'qc_report.json').read_text(encoding='utf-8'))['status'] == 'pass'
source_files.extend([map_directory / 'mld_model_observation_display.npz', map_directory / 'metadata.json'])


# %% 4. 固定组成员等权平均，不随缺测改变成员数

payload = {'longitude': longitude, 'latitude': latitude, 'color_levels_m': color_levels,
           'color_ticks_m': color_ticks, 'contour_levels_m': np.array([200, 500, 1000, 1500, 2000])}
map_diagnostics, group_records = [], []
for group in groups:
    members = membership.loc[membership.group == group].sort_values('member_order').model.tolist()
    assert set(members) == set(registry['members'][group])
    initial_stack = np.stack([initial_fields[model] for model in members])
    later_stack = np.stack([later_fields[model] for model in members])
    delta_stack = np.stack([delta_fields[model] for model in members])
    valid_members = np.isfinite(initial_stack) & np.isfinite(later_stack) & np.isfinite(delta_stack)
    valid_count = valid_members.sum(axis=0)
    common_mask = valid_count == len(members)
    # 普通mean传播NaN；仅保留全部固定成员共同有效格点。
    initial_mean = np.where(common_mask, initial_stack.mean(axis=0), np.nan)
    later_mean = np.where(common_mask, later_stack.mean(axis=0), np.nan)
    delta_mean = np.where(common_mask, delta_stack.mean(axis=0), np.nan)
    np.testing.assert_allclose(later_mean - initial_mean, delta_mean, rtol=1e-10, atol=1e-9, equal_nan=True)
    payload[group + '_initial_mld_m'] = initial_mean
    payload[group + '_later_mld_m'] = later_mean
    payload[group + '_delta_mld_m'] = delta_mean
    payload[group + '_valid_member_count'] = valid_count
    payload[group + '_common_mask'] = common_mask
    group_records.append({'group': group, 'members': members, 'n_models': len(members)})
    map_diagnostics.append({'group': group, 'n_models': len(members), 'valid_common_cells': int(common_mask.sum()),
                            'incomplete_membership_cells': int(((valid_count > 0) & ~common_mask).sum()),
                            'minimum_delta_m': float(np.nanmin(delta_mean)), 'maximum_delta_m': float(np.nanmax(delta_mean)),
                            'outside_color_limits_cells': int(np.sum(np.abs(delta_mean) > color_levels[-1])),
                            'max_budget_identity_error_m': float(np.nanmax(np.abs(later_mean - initial_mean - delta_mean)))})


# %% 5. 原样保存曲线输入，计算观测均值加减年际标准差

destination.mkdir(parents=True)
curve_tables = [(n2_directory, 'membership.csv', None),
                (n2_directory, 'member_profiles.csv', region_n2),
                (n2_directory, 'group_mean_profiles.csv', region_n2),
                (n2_directory, 'observation_profiles.csv', region_n2),
                (mld_directory, 'member_timeseries.csv', region_mld),
                (mld_directory, 'group_mean_timeseries.csv', region_mld)]
for directory, name, region in curve_tables:
    source_path = directory / name
    table = pd.read_csv(source_path)
    if region is not None:
        table = table.loc[table.region == region].copy()
    if name == 'observation_profiles.csv':
        observation_mean = table.iapv4_march_mean_n2_1e5_s2
        observation_sd = table.iapv4_march_interannual_sd_n2_1e5_s2
        assert np.all(np.isfinite(observation_sd)) and np.all(observation_sd >= 0)
        table['lower_1sd_n2_1e5_s2'] = observation_mean - observation_sd
        table['upper_1sd_n2_1e5_s2'] = observation_mean + observation_sd
    table.to_csv(destination / name, index=False, encoding='utf-8-sig')
    source_files.append(source_path)
for feature in ['land', 'coastline']:
    source_path = map_directory / f'map_{feature}_50m.geojson'
    shutil.copy2(source_path, destination / source_path.name)
    source_files.append(source_path)
np.savez_compressed(destination / 'grouped_mld_patterns.npz', **payload)
pd.DataFrame(map_diagnostics).to_csv(destination / 'map_diagnostics.csv', index=False)


# %% 6. 保存处理口径与来源

metadata = {
    'status': 'passed', 'groups': group_records, 'region_profiles': region_n2, 'region_timeseries': region_mld,
    'scope': 'Labrador-Irminger nine-panel replacement only; Nordic figure unchanged.',
    'dimensions': {'model_mld_maps': 'model8, latitude141, longitude341',
                   'group_mld_maps': 'group3, latitude141, longitude341; saved as separate named2D arrays',
                   'profiles': 'model8/group3/observation1, depth100', 'timeseries': 'model8/group3, year280'},
    'map_extent': map_extent,
    'map_method': 'Arithmetic equal-model mean of existing n058/n0710.25degree linearly mapped display fields; no new interpolation or smoothing. All fixed group members must be finite in both periods. No skipna or time-varying membership.',
    'map_fill': 'Mean of each member March MLD later-minus-early; positive=deepening, negative=shoaling.',
    'map_contours': 'Equal-model mean initial March absolute MLD,200/500/1000/1500/2000m; not contours of change.',
    'map_color_limits_m': [float(color_levels[0]), float(color_levels[-1])],
    'map_ice_edges': 'none; no sea-ice data averaged or plotted',
    'curve_policy': 'Original individual-model curves and group means copied unchanged. Eleven-year centered MLD smoothing and endpoint gaps retained. Middle row has no Argo/ELG reference; xlabel Year.',
    'observation': 'Top row only: IAPv4 March2005-2024 mean plus/minus1sample interannual SD from existing profile table; not measurement error, SE or CI. Same reference repeated across groups, not independent observations.',
    'n2_viewports': 'Panels a-b use -0.5..3.4 instead of old zoom to include full IAP mean+/-SD. Panel c retains -1..7.1. Depth0..1000m. Different x scales must be read from ticks.',
    'limitations': ['GFDL initial N2 years6-15 retained; MLD inherited index-based timing caveat retained; no new exact recovery timing inference.',
                    'N2 coastal polygons, MLD legacy rectangles and spatial map differ in geographical support; map area is broader North Atlantic, not a direct area reconstruction of the middle-row curves.',
                    'Earlier MLD source regridding may interpolate across coastal gaps; no higher effective native resolution or significance claimed.',
                    'Existing March endpoint windows and absence of piControl drift correction retained.'],
    'source_files': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in source_files],
    'software_skill': 'scientific-visualization: fixed members, common map scale, explicit observation SD, actual exports to be inspected; software reference retained in upstream n059 readme.'}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print(pd.DataFrame(map_diagnostics).to_string(index=False), flush=True)
print(destination, flush=True)
