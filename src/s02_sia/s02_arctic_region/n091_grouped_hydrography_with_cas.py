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
previous = root / 'output' / stage / 'n077_regional_definition_comparison'
cas_source = root / 'output/s02_sia/s02_arctic_region/p25_cas_hydrography/n083_cas_mld_stratification'
history_source = root / 'output/s02_sia/s02_arctic_region/p25_cas_hydrography/n089_rectangular_mld_nine_models'
map_source = root / 'output/s02_sia/s02_arctic_region/p25_cas_hydrography/n085_nine_model_hydrography_panels'
old_map_source = root / 'output/s02_sia/s02_arctic_region/p06_mld/n071_model_observation_panels'
assessment_source = root / 'output/s06_manuscript/n08_recovery_group_extension_assessment'
recovery26_source = root / 'output/s01_amoc/n75_hysteresis_morphology_26p5'
recovery40_source = root / 'output/s01_amoc/n22_amoc_recovery_definition_audit'
destination = root / 'output' / stage / 'n091_grouped_hydrography_with_cas'
assert not destination.exists(), '正式输出存在，覆盖前先归档。'
groups = ['better_nonrecovery', 'better_recovery', 'larger_bias']
region_n2, region_mld = 'Labrador-Irminger Seas', 'Labrador-Irminger'
groups_expected = {'better_nonrecovery': 2, 'better_recovery': 3, 'larger_bias': 4}
profile_mean_parts, history_mean_parts, map_diagnostics = [], [], []
source_paths = []


# %% 3. 保持既有成员并将CAS加入基本恢复组

membership = pd.read_csv(previous / 'membership.csv', keep_default_na=False)
cas_member = {name: '' for name in membership.columns}
cas_member.update({'model': 'CAS-ESM2-0', 'group': 'better_recovery', 'panel_index': 1,
                   'panel_title': 'Better-simulated stratification\nSubstantial recovery',
                   'group_size': 3, 'member_order': 2, 'member': 'r1i1p1f1',
                   'color': '#B54F2B', 'line_style': 'solid', 'calendar_first': 200,
                   'calendar_last': 209, 'experiment_year_first': 1, 'experiment_year_last': 10,
                   'n_marches': 10, 'grouping_basis': 'CAS initial LI stratification close to IAP; late/initial AMOC exceeds0.8 at both26.5N and40N.'})
membership = pd.concat([membership, pd.DataFrame([cas_member])], ignore_index=True)
membership.loc[membership.group.eq('better_recovery'), 'group_size'] = 3
membership.loc[membership.model.eq('GFDL-ESM4'), ['experiment_year_first', 'experiment_year_last']] = [1, 10]
membership.loc[membership.model.eq('GFDL-ESM4'), 'baseline_exception'] = 'Author convention: source year6 is analysis year1; initial source6-15 retained.'
assert membership.groupby('group').size().to_dict() == groups_expected
assert membership.model.nunique() == 9


# %% 4. 合并矩形区域剖面并重算等权组平均

profiles = pd.read_csv(previous / 'rectangles/member_profiles.csv', float_precision='round_trip')
cas = pd.read_csv(cas_source / 'cas_march_n2_profiles.csv', float_precision='round_trip')
cas = cas.loc[cas.version.eq('rectangles') & cas.region.eq(region_mld) & cas.period.eq('up_initial')].copy()
cas = cas[['model', 'depth_m', 'n2_s2', 'interannual_sd_s2', 'valid_march_count']]
cas = cas.rename(columns={'n2_s2': 'initial_march_n2_s2', 'interannual_sd_s2': 'interannual_sd_n2_s2'})
cas['region'] = region_n2
cas['group'] = 'better_recovery'
cas['initial_march_n2_1e5_s2'] = cas.initial_march_n2_s2 * 1e5
profiles = pd.concat([profiles, cas], ignore_index=True)
assert len(profiles) == 900 and not profiles.duplicated(['model', 'depth_m']).any()
assert profiles.loc[profiles.depth_m.eq(5), 'initial_march_n2_s2'].isna().all()
for group in groups:
    models = membership.loc[membership.group.eq(group)].sort_values('member_order').model.tolist()
    matrix = profiles.loc[profiles.group.eq(group)].pivot(index='depth_m', columns='model', values='initial_march_n2_1e5_s2')[models]
    profile_mean_parts.append(pd.DataFrame({'group': group, 'region': region_n2, 'depth_m': matrix.index,
                                            'group_mean_n2_1e5_s2': np.mean(matrix.to_numpy(), axis=1),
                                            'n_models_expected': len(models), 'n_models_valid': matrix.notna().sum(axis=1)}))
profile_means = pd.concat(profile_mean_parts, ignore_index=True)


# %% 5. 读取统一矩形MLD序列并更新组平均

history = pd.read_csv(history_source / 'rectangular_mld_timeseries.csv', float_precision='round_trip')
history = history.loc[history.region.eq(region_mld)].copy()
history = history.merge(membership[['model', 'group']], on='model', validate='many_to_one')
history = history[['model', 'group', 'region', 'simulation_year', 'march_mld_m', 'march_mld_11yr_m']]
assert len(history) == 9 * 280 and history.march_mld_m.notna().all()
original_history = pd.read_csv(previous / 'rectangles/member_timeseries.csv', float_precision='round_trip')
for model in original_history.model.unique():
    old = original_history.loc[original_history.model.eq(model)].sort_values('simulation_year')
    new = history.loc[history.model.eq(model)].sort_values('simulation_year')
    np.testing.assert_allclose(new[['march_mld_m', 'march_mld_11yr_m']], old[['march_mld_m', 'march_mld_11yr_m']], atol=1e-9, rtol=0, equal_nan=True)
for group in groups:
    models = membership.loc[membership.group.eq(group)].sort_values('member_order').model.tolist()
    matrix = history.loc[history.group.eq(group)].pivot(index='simulation_year', columns='model', values='march_mld_m')[models]
    mean = np.mean(matrix.to_numpy(), axis=1)
    smooth = pd.Series(mean).rolling(11, center=True, min_periods=11).mean().to_numpy()
    history_mean_parts.append(pd.DataFrame({'group': group, 'region': region_mld, 'simulation_year': matrix.index,
                                            'group_mean_m': mean, 'group_mean_11yr_m': smooth}))
history_means = pd.concat(history_mean_parts, ignore_index=True)


# %% 6. 仅重算加入CAS的中列空间平均

with np.load(previous / 'grouped_mld_patterns.npz') as saved:
    payload = {name: saved[name].copy() for name in saved.files}
with np.load(old_map_source / 'mld_model_observation_display.npz') as saved:
    initial_fields = [saved[model + '_up_initial_march_mld_m'].copy() for model in ['NorESM2-LM', 'GFDL-ESM4']]
    later_fields = [saved[model + '_down_final_march_mld_m'].copy() for model in ['NorESM2-LM', 'GFDL-ESM4']]
    delta_fields = [saved[model + '_delta_later_early_march_mld_m'].copy() for model in ['NorESM2-LM', 'GFDL-ESM4']]
with np.load(map_source / 'nine_model_mld_ice_display.npz') as saved:
    np.testing.assert_array_equal(saved['lon2_plot'], payload['longitude'])
    np.testing.assert_array_equal(saved['lat2_plot'], payload['latitude'])
    cas_initial = saved['CAS-ESM2-0_up_initial_march_mld_m'].copy()
    cas_delta = saved['CAS-ESM2-0_delta_later_early_march_mld_m'].copy()
initial_fields.append(cas_initial)
delta_fields.append(cas_delta)
later_fields.append(cas_initial + cas_delta)
initial_stack, later_stack, delta_stack = np.stack(initial_fields), np.stack(later_fields), np.stack(delta_fields)
valid = np.isfinite(initial_stack) & np.isfinite(later_stack) & np.isfinite(delta_stack)
valid_count = valid.sum(axis=0)
common_mask = valid_count == 3
# 固定三个模式共同有效，禁止skipna造成空间上成员数变化。
payload['better_recovery_initial_mld_m'] = np.where(common_mask, initial_stack.mean(axis=0), np.nan)
payload['better_recovery_later_mld_m'] = np.where(common_mask, later_stack.mean(axis=0), np.nan)
payload['better_recovery_delta_mld_m'] = np.where(common_mask, delta_stack.mean(axis=0), np.nan)
payload['better_recovery_valid_member_count'] = valid_count
payload['better_recovery_common_mask'] = common_mask
for group in groups:
    np.testing.assert_allclose(payload[group + '_later_mld_m'] - payload[group + '_initial_mld_m'],
                               payload[group + '_delta_mld_m'], atol=1e-9, rtol=0, equal_nan=True)
    map_diagnostics.append({'group': group, 'n_models': groups_expected[group],
                            'valid_common_cells': int(payload[group + '_common_mask'].sum()),
                            'minimum_delta_m': float(np.nanmin(payload[group + '_delta_mld_m'])),
                            'maximum_delta_m': float(np.nanmax(payload[group + '_delta_mld_m']))})


# %% 7. 保存诊断依据与正式绘图输入

rmsd = pd.read_csv(assessment_source / 'initial_n2_iap_comparison.csv')
recovery26 = pd.read_csv(recovery26_source / 'endpoint_diagnostics.csv').set_index('model')
recovery40 = pd.read_csv(recovery40_source / 'nine_model_endpoint_diagnostics.csv').set_index('model')
cas_rmsd = rmsd.loc[rmsd.model.eq('CAS-ESM2-0') & rmsd.region.eq('LI'), 'rmsd_1e5_s2'].iloc[0]
assert recovery26.loc['CAS-ESM2-0', 'endpoint_ratio'] >= .8
assert recovery40.loc['CAS-ESM2-0', 'endpoint_ratio'] >= .8
destination.mkdir(parents=True)
version_directory = destination / 'rectangles'
version_directory.mkdir()
membership.to_csv(destination / 'membership.csv', index=False)
for name, table in [('member_profiles.csv', profiles), ('group_mean_profiles.csv', profile_means),
                    ('member_timeseries.csv', history), ('group_mean_timeseries.csv', history_means)]:
    table.to_csv(version_directory / name, index=False)
shutil.copy2(previous / 'rectangles/observation_profiles.csv', version_directory / 'observation_profiles.csv')
for name in ['map_land_50m.geojson', 'map_coastline_50m.geojson', 'n2_coastal_regions.geojson',
             'mld_rectangular_vertices.csv', 'mld_rectangular_regions.csv']:
    shutil.copy2(previous / name, destination / name)
np.savez_compressed(destination / 'grouped_mld_patterns.npz', **payload)
pd.DataFrame(map_diagnostics).to_csv(destination / 'map_diagnostics.csv', index=False)
decision = {'model': 'CAS-ESM2-0', 'assigned_group': 'better_recovery',
            'initial_LI_n2_rmsd_1e5_s2': float(cas_rmsd), 'rmsd_depth_range_m': [15, 995],
            'AMOC_recovery_ratio_26p5N': float(recovery26.loc['CAS-ESM2-0', 'endpoint_ratio']),
            'AMOC_recovery_ratio_40N': float(recovery40.loc['CAS-ESM2-0', 'endpoint_ratio']),
            'recovery_threshold': .8,
            'scope': 'Initial LI stratification quality plus late AMOC recovery; not the hysteresis-shape grouping from Figure1. Other8models retain their previous groups.'}
(destination / 'cas_grouping_basis.json').write_text(json.dumps(decision, indent=2), encoding='utf-8')


# %% 8. 保存口径与来源

source_paths.extend([previous / name for name in ['membership.csv', 'grouped_mld_patterns.npz', 'metadata.json',
                     'rectangles/member_profiles.csv', 'rectangles/member_timeseries.csv', 'rectangles/observation_profiles.csv']])
source_paths.extend([cas_source / 'cas_march_n2_profiles.csv', history_source / 'rectangular_mld_timeseries.csv',
                     old_map_source / 'mld_model_observation_display.npz', map_source / 'nine_model_mld_ice_display.npz',
                     assessment_source / 'initial_n2_iap_comparison.csv', recovery26_source / 'endpoint_diagnostics.csv',
                     recovery40_source / 'nine_model_endpoint_diagnostics.csv'])
caption = ('Initial March stratification and mixed-layer depth responses. Columns show CESM2/CNRM-ESM2-1 '
           '(Non-recovery, n=2), NorESM2-LM/GFDL-ESM4/CAS-ESM2-0 (Substantial recovery, n=3), and '
           'ACCESS-ESM1-5/CanESM5/MIROC-ES2L/UKESM1-0-LL (larger initial stratification bias, n=4). '
           'Profiles and March MLD time series use the rectangular Labrador-Irminger region. '
           'The same IAP March2005-2024 mean and interannual SD are shown in all top panels. '
           'Member MLD lines retain full-window11yr centred running means and missing end windows. '
           'Maps show fixed-member equal-weight late10-minus-initial10 March MLD, with initial '
           'group-mean contours; all group members must be valid at each mapped cell. '
           'CAS is added only to the middle column; original8member data and two outer map groups are unchanged. '
           'The right column is defined by stratification bias and contains differing AMOC recovery outcomes.')
metadata = {'status': 'passed', 'groups': groups, 'group_sizes': [2, 3, 4], 'map_extent': [-65., 20., 50., 85.],
            'scope': 'Add CAS only to the original better-recovery column; preserve other8model group membership. Entire composite remains rectangular LI diagnostics.',
            'dimensions': {'profiles': 'model9 xdepth100; observation depth100', 'history': 'model9 xyear280',
                           'group_maps': 'group3 xlatitude141 xlongitude341'},
            'regional_support': 'N2:50<=lat<65,-65<=lon<-20; MLD matches saved legacy rectangle inclusive outer-edge nodes. NS rectangle appears only in locator.',
            'mld_method': 'Reuse n089 unified rectangles, nearest valid native10m reference; sigma0 threshold0.03kg/m3 resolved crossing only; native-area means. No raw reprocessing.',
            'GFDL_time': 'Author-confirmed source6 maps analysis1; initial6-15. Existing source values unchanged.',
            'map_support': 'Full North Atlantic saved gridded MLD, not regional means; same interpolation/color scale/initial contours as original. CAS already-saved map added with fixed3member common validity.',
            'observation': 'IAP March2005-2024 mean+/-sample interannual SD; file copied unchanged.',
            'drift_correction': 'none', 'cas_grouping_basis': decision, 'caption': caption,
            'source_files': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in source_paths]}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
(destination / 'caption.txt').write_text(caption + '\n', encoding='utf-8')
print(json.dumps(decision, ensure_ascii=True, indent=2), flush=True)
print(destination.resolve(), flush=True)
