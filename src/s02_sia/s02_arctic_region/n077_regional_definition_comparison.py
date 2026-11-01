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
source_base = root / 'output' / stage / 'n074_grouped_mld_patterns'
locator_directory = root / 'output' / stage / 'n076_region_locator_inputs'
destination = root / 'output' / stage / 'n077_regional_definition_comparison'
assert not destination.exists(), '正式对照结果存在，重跑前先归档。'
rect_n2 = root / 'output/s02sia/s02_arctic_region/p07_stratification/n01_initial_profiles'
rect_obs = root / 'output/s02sia/s02_arctic_region/p07_stratification/n02_iapv4_observation'
coastal_directory = root / 'output/s02sia/s02_arctic_region/p14_mld_fixed_regions_raw'
region_n2, region_mld = 'Labrador-Irminger Seas', 'Labrador-Irminger'
groups = ['better_nonrecovery', 'better_recovery', 'larger_bias']
versions = ['rectangles', 'coastal']
source_records, diagnostics = [], []


# %% 3. 准备共同成员、原样保留第三行地图与海区边界

membership = pd.read_csv(source_base / 'membership.csv')
assert membership.groupby('group').size().to_dict() == {'better_nonrecovery': 2, 'better_recovery': 2, 'larger_bias': 4}
destination.mkdir(parents=True)
common_sources = [(source_base / name, name) for name in ['membership.csv', 'grouped_mld_patterns.npz',
                   'map_land_50m.geojson', 'map_coastline_50m.geojson']]
common_sources.extend([(locator_directory / name, name) for name in
                       ['n2_coastal_regions.geojson', 'mld_rectangular_vertices.csv', 'mld_rectangular_regions.csv']])
for source, name in common_sources:
    shutil.copy2(source, destination / name)
    source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    assert hashlib.sha256((destination / name).read_bytes()).hexdigest() == source_hash
    source_records.append({'path': str(source), 'sha256': source_hash})


# %% 4. 整理两套既有N²与各自对应观测，不混用区域

for version in versions:
    version_directory = destination / version
    version_directory.mkdir()
    if version == 'rectangles':
        profile_path = rect_n2 / 'na_initial_march_stratification_profiles.csv'
        observation_path = rect_obs / 'iapv4_march_2005_2024_stratification_profiles.csv'
        metadata_paths = [rect_n2 / 'na_initial_march_stratification_metadata.json',
                          rect_obs / 'iapv4_march_2005_2024_stratification_metadata.json']
        history_path = source_base / 'member_timeseries.csv'
    else:
        profile_path = coastal_directory / 'coast_region_20261012_initial_n2_profiles.csv'
        observation_path = coastal_directory / 'coast_region_20261012_iapv4_profiles.csv'
        metadata_paths = [coastal_directory / 'coast_region_20261012_initial_n2_metadata.json',
                          coastal_directory / 'coast_region_20261012_iapv4_metadata.json',
                          coastal_directory / 'coast_region_20261012_verified_metadata.json',
                          coastal_directory / 'coast_region_20261012_plot_metadata.json']
        history_path = coastal_directory / 'coast_region_20261012_plot_data.csv'
    profiles = pd.read_csv(profile_path)
    profiles = profiles.loc[(profiles.region == region_n2) & profiles.model.isin(membership.model)].copy()
    profiles = profiles.drop(columns=['group', 'panel_index'], errors='ignore')
    profiles = profiles.merge(membership[['model', 'group']], on='model', validate='many_to_one')
    assert len(profiles) == 800 and profiles.model.nunique() == 8
    assert profiles.loc[profiles.depth_m == 5, 'initial_march_n2_1e5_s2'].isna().all()
    observations = pd.read_csv(observation_path)
    observations = observations.loc[observations.region == region_n2].sort_values('depth_m').copy()
    assert len(observations) == 100
    observations['lower_1sd_n2_1e5_s2'] = observations.iapv4_march_mean_n2_1e5_s2 - observations.iapv4_march_interannual_sd_n2_1e5_s2
    observations['upper_1sd_n2_1e5_s2'] = observations.iapv4_march_mean_n2_1e5_s2 + observations.iapv4_march_interannual_sd_n2_1e5_s2
    assert observations.lower_1sd_n2_1e5_s2.min() >= -.5 and observations.upper_1sd_n2_1e5_s2.max() <= 3.4

    profile_mean_rows = []
    for group in groups:
        members = membership.loc[membership.group == group].sort_values('member_order').model.tolist()
        matrix = profiles.loc[profiles.group == group].pivot(index='depth_m', columns='model', values='initial_march_n2_1e5_s2')[members]
        means = np.mean(matrix.to_numpy(), axis=1)
        profile_mean_rows.append(pd.DataFrame({'group': group, 'region': region_n2, 'depth_m': matrix.index,
                                              'group_mean_n2_1e5_s2': means, 'n_models_expected': len(members),
                                              'n_models_valid': matrix.notna().sum(axis=1)}))
    profile_means = pd.concat(profile_mean_rows, ignore_index=True)

    # 整理该分区既有MLD序列，不伪称已统一算法与物理年份。

    history = pd.read_csv(history_path)
    history = history.loc[(history.region == region_mld) & history.model.isin(membership.model)].copy()
    if version == 'coastal':
        history = history.rename(columns={'march_mld_mean_crossing_only_m': 'march_mld_m',
                                          'march_mld_mean_crossing_only_m_11yr': 'march_mld_11yr_m'})
        history = history.merge(membership[['model', 'group']], on='model', validate='many_to_one')
    history = history[['model', 'group', 'region', 'simulation_year', 'march_mld_m', 'march_mld_11yr_m']].copy()
    assert len(history) == 2240 and history.march_mld_m.notna().all()
    assert history.march_mld_11yr_m.max() < 900
    history_mean_rows = []
    for group in groups:
        members = membership.loc[membership.group == group].sort_values('member_order').model.tolist()
        matrix = history.loc[history.group == group].pivot(index='simulation_year', columns='model', values='march_mld_m')[members]
        assert matrix.index.tolist() == list(range(1, 281))
        means = np.mean(matrix.to_numpy(), axis=1)
        smooth = pd.Series(means).rolling(11, center=True, min_periods=11).mean().to_numpy()
        history_mean_rows.append(pd.DataFrame({'group': group, 'region': region_mld, 'simulation_year': matrix.index,
                                              'group_mean_m': means, 'group_mean_11yr_m': smooth}))
    history_means = pd.concat(history_mean_rows, ignore_index=True)
    for name, table in [('member_profiles.csv', profiles), ('group_mean_profiles.csv', profile_means),
                        ('observation_profiles.csv', observations), ('member_timeseries.csv', history),
                        ('group_mean_timeseries.csv', history_means)]:
        table.to_csv(version_directory / name, index=False, encoding='utf-8-sig')
    for source in [profile_path, observation_path, history_path] + metadata_paths:
        source_records.append({'path': str(source), 'sha256': hashlib.sha256(source.read_bytes()).hexdigest()})
    for member in membership.itertuples(index=False):
        sub = history.loc[history.model == member.model].sort_values('simulation_year')
        diagnostics.append({'version': version, 'model': member.model, 'group': member.group,
                            'initial_10_index_mean_mld_m': float(sub.iloc[:10].march_mld_m.mean()),
                            'final_10_index_mean_mld_m': float(sub.iloc[-10:].march_mld_m.mean()),
                            'maximum_initial_n2_1e5_s2': float(profiles.loc[profiles.model == member.model, 'initial_march_n2_1e5_s2'].max())})


# %% 5. 保存明确的对照范围与局限

pd.DataFrame(diagnostics).to_csv(destination / 'regional_diagnostics.csv', index=False)
metadata = {'status': 'passed', 'versions': versions, 'groups': groups, 'map_extent': [-65., 20., 50., 85.],
            'scope': 'Two existing regional diagnostic versions for visual comparison. Not a controlled region-only sensitivity test.',
            'first_row': 'Each definition uses its own8model N2 profiles and matching IAP March2005-2024 mean+/-sample interannual SD. Model method common10m geometric-depth TEOS-10 N2; first10March years, GFDL6-15 in both. No reclassification of fixed groups.',
            'second_row': 'Existing resolved-threshold-crossing regional mean MLD and saved11year centered member curves; fixed2/2/4equal-model mean recomputed from saved raw annual values. No Argo/ELG.',
            'method_confounders': ['Legacy rectangle MLD uses nearest valid native level to10m for reference density, coastal MLD interpolates reference density to10m; both use sigma0+0.03kg/m3 resolved crossings.',
                                   'GFDL legacy rectangle and coastal MLD use inherited index/calendar alignments that are not a matched physical-year experiment. Coastal source labels up calendar6-145 as index1-140, CDR calendar141-280 as index141-280. No timing correction or padding added.',
                                   'Rectangular N2 uses upper-exclusive bounds; saved legacy MLD includes some east/north edges. Locator shows geographical boxes, not exact edge-cell equality.',
                                   'Coastal regions are user-drawn/offshore-boundary coastline-clipped sectors, not standard basin boundaries; Nordic sector extends into Iceland Basin.',
                                   'Differences in MLD cannot be attributed solely to region boundaries without harmonizing the above methods and time axes. No claim that stronger recovery-group separation establishes the scientifically preferred region.'],
            'third_row': 'Bit-for-bit copy of n074group maps, identical in both figures: original later-minus-early fill and initial group-mean contours, no region masking, no sea ice, same color scale. No new MLD scientific calculation.',
            'dimensions': {'profiles': 'version2,model8,depth100', 'observations': 'version2,depth100',
                           'histories': 'version2,model8,year280', 'group_maps': 'group3,latitude141,longitude341'},
            'display': {'n2_x_limits': [[-.5, 3.4], [-.5, 3.4], [-1., 8.7]], 'mld_y_limits': [0, 900],
                        'inset': 'Inside a blank region; each version shows its own region outlines; does not count as a new main panel.'},
            'drift_correction': 'none', 'source_files': source_records,
            'software_skill': 'scientific-visualization: disclose non-region confounders, match viewports, freeze maps and verify inset boundaries and actual exports; existing software citation retained in upstream readme.'}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print(pd.DataFrame(diagnostics).round(3).to_string(index=False), flush=True)
print(destination, flush=True)
