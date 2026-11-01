# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import shutil

import numpy as np
import pandas as pd


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
source = root / 'output/s01_amoc/n73_hysteresis_morphology_groups'
annual_source = root / 'output/s01_amoc/n22_amoc_recovery_definition_audit'
latitude_source = root / 'output/amoc_s01_n01/l_40'
destination = root / 'output/s01_amoc/n79_seven_panel_hysteresis_40n'
assert not destination.exists() or not any(destination.iterdir()), '正式结果存在，须先归档。'
initial_source = root / 'output/amoc_s01_n01/l_40/CMIP6_AMOC_9models_annual_series.csv'
group_names = ['Recovery with small A', 'Recovery with large A', 'Non-recovery']
schematic_mapping = {'c': 'b', 'f': 'd', 'j': 'f'}
model_panels = ['c', 'e', 'g']
colors = ['#AA4499', '#117733', '#CC6677', '#332288', '#44AA99',
          '#88CCEE', '#882255', '#999933', '#DDCC77']
curve_parts, schematic_parts, model_arrows, schematic_arrows = [], [], [], []
panel_records = [{'panel': 'a', 'kind': 'overview', 'row': 0, 'models': 'all9', 'label': ''}]
annotation_records = []


# %% 3. 读取已核对的40.0度结果

metadata = json.loads((source / 'metadata.json').read_text(encoding='utf-8'))
assert metadata['overview_target_latitude_degN'] == metadata['pathway_target_latitude_degN'] == 40.0
old_panels = pd.read_csv(source / 'panels.csv', keep_default_na=False)
paths = pd.read_csv(source / 'path_coordinates.csv', float_precision='round_trip')
arrows = pd.read_csv(source / 'path_arrows.csv', float_precision='round_trip')
colors_by_model = dict(zip(metadata['models'], colors))
# 图例采用未平滑的强迫实验初始1–10年平均，与异常和归一化参考一致。
initial_annual = pd.read_csv(initial_source, float_precision='round_trip')
initial_decade = initial_annual.loc[initial_annual['Simulation year'].between(1, 10)]
assert initial_decade.groupby('Model').size().eq(10).all()
initial_means = initial_decade.groupby('Model')['Forced AMOC strength (Sv)'].mean()
initial_records = []
for model in metadata['models']:
    value = float(initial_means.loc[model])
    initial_records.append({'model': model, 'initial_mean_sv': value,
                            'first_simulation_year': 1, 'last_simulation_year': 10,
                            'label': f'{model} ({value:.1f})'})
mean_initial = float(initial_means.loc[metadata['models']].mean())
initial_records.append({'model': 'mmm', 'initial_mean_sv': mean_initial,
                        'first_simulation_year': 1, 'last_simulation_year': 10,
                        'label': f'9-model mean ({mean_initial:.1f})'})
legend_values = pd.DataFrame(initial_records)
normalized = pd.read_csv(annual_source / 'nine_model_normalized_series.csv', float_precision='round_trip')
normalized['relative_amoc'] = normalized.smooth_11yr_percent_initial / 100.
latitude_mapping = pd.read_csv(latitude_source / 'CMIP6_AMOC_9models_processing_metadata.csv')


# %% 4. 每组模式合并到同一面板并串联完整往返轨迹

for row, (members, label, model_panel) in enumerate(zip(metadata['members'], group_names, model_panels), 1):
    schematic_panel = ['b', 'd', 'f'][row - 1]
    panel_records.append({'panel': schematic_panel, 'kind': 'schematic', 'row': row,
                          'models': '; '.join(members), 'label': label})
    panel_records.append({'panel': model_panel, 'kind': 'models', 'row': row,
                          'models': '; '.join(members), 'label': label})
    for model in members:
        old_panel = old_panels.loc[old_panels.model.eq(model), 'panel'].iloc[0]
        selected = paths.loc[paths.panel.eq(old_panel)].sort_values(['simulation_year', 'node'])
        # 峰值140年在两分支各出现一次；仅去掉同一个转折点的重复副本。
        selected = selected.drop_duplicates('simulation_year').sort_values('simulation_year')
        assert selected.simulation_year.astype(int).tolist() == list(range(1, 281))
        reference = normalized.loc[normalized.model.eq(model)].sort_values('simulation_year')
        np.testing.assert_allclose(selected.relative_amoc, reference.relative_amoc, atol=1e-12, rtol=0)
        curve_parts.append(pd.DataFrame({'panel': model_panel, 'model': model,
                                         'simulation_year': selected.simulation_year.astype(int).to_numpy(),
                                         'co2_ppm': selected.co2_ppm.to_numpy(),
                                         'relative_amoc': selected.relative_amoc.to_numpy(),
                                         'color': colors_by_model[model]}))
        selected_arrows = arrows.loc[arrows.panel.eq(old_panel)].copy()
        selected_arrows['panel'] = model_panel
        selected_arrows['model'] = model
        selected_arrows['color'] = colors_by_model[model]
        model_arrows.extend(selected_arrows.to_dict('records'))


# %% 5. 移除原b概念图并把面积标注移至两幅示意

for old_panel, new_panel in schematic_mapping.items():
    branches = {}
    for phase in ['increase', 'decrease']:
        selected = paths.loc[paths.panel.eq(old_panel) & paths.phase.eq(phase)].sort_values('node').copy()
        if new_panel == 'b':
            # 概念小面积闭合滞回；仅修改示意，不拟合或更改实际模式数据。
            endpoint = float(paths.loc[paths.panel.eq(old_panel) & paths.phase.eq('increase')]
                             .sort_values('node').relative_amoc.iloc[-1])
            fraction = (selected.co2_ppm.to_numpy() - selected.co2_ppm.iloc[0]) / (selected.co2_ppm.iloc[-1] - selected.co2_ppm.iloc[0])
            centre = 1. + (endpoint - 1.) * fraction + .012 * np.sin(2. * np.pi * fraction)
            separation = .044 * np.sin(np.pi * fraction) ** 2
            branch = centre + (.5 * separation if phase == 'increase' else -.5 * separation)
            branch[[0, -1]] = [1., endpoint]
            selected['relative_amoc'] = branch
        selected['panel'] = new_panel
        schematic_parts.append(selected)
        branches[phase] = selected
        arrow = arrows.loc[arrows.panel.eq(old_panel) & arrows.phase.eq(phase)].iloc[0]
        schematic_arrows.append({'panel': new_panel, 'phase': phase,
                                  'start_x': arrow.start_x, 'end_x': arrow.end_x,
                                  'start_y': float(np.interp(arrow.start_x, selected.co2_ppm, selected.relative_amoc)),
                                  'end_y': float(np.interp(arrow.end_x, selected.co2_ppm, selected.relative_amoc))})
    if new_panel == 'd':
        target_x = metadata['reference_co2_ppm'] + .45 * (metadata['maximum_co2_ppm'] - metadata['reference_co2_ppm'])
        target_y = .5 * (np.interp(target_x, branches['increase'].co2_ppm, branches['increase'].relative_amoc)
                         + np.interp(target_x, branches['decrease'].co2_ppm, branches['decrease'].relative_amoc))
        annotation_records.append({'panel': new_panel, 'text_axes_x': .65, 'text_axes_y': .22,
                                    'target_co2_ppm': target_x, 'target_relative_amoc': target_y})


# %% 6. 保存七面板输入和来源

curves = pd.concat(curve_parts, ignore_index=True)
schematics = pd.concat(schematic_parts, ignore_index=True)
assert len(curves) == 9 * 280 and curves.model.nunique() == 9
assert len(panel_records) == 7 and len(annotation_records) == 1
destination.mkdir(parents=True, exist_ok=True)
legend_values.to_csv(destination / 'initial_legend_values.csv', index=False)
curves.to_csv(destination / 'model_paths.csv', index=False)
schematics.to_csv(destination / 'schematic_paths.csv', index=False)
pd.DataFrame(model_arrows).to_csv(destination / 'model_arrows.csv', index=False)
pd.DataFrame(schematic_arrows).to_csv(destination / 'schematic_arrows.csv', index=False)
pd.DataFrame(panel_records).to_csv(destination / 'panels.csv', index=False)
pd.DataFrame(annotation_records).to_csv(destination / 'area_annotations.csv', index=False)
for name in ['co2_forcing.csv', 'model_hysteresis_metrics.csv']:
    shutil.copy2(source / name, destination / name)
shutil.copy2(annual_source / 'nine_model_endpoint_diagnostics.csv', destination / 'endpoint_diagnostics.csv')
latitude_mapping.to_csv(destination / 'latitude_mapping.csv', index=False)
normalized.to_csv(destination / 'nine_model_normalized_series.csv', index=False)
# 图a与路径统一为分析1–280年内的同一40N强度及11年平均，不使用后续稳定阶段。
strength = normalized.pivot(index='simulation_year', columns='model', values='raw_sv').loc[1:280, metadata['models']]
smoothed = normalized.pivot(index='simulation_year', columns='model', values='smooth_11yr_sv').loc[1:280, metadata['models']]
initial = strength.loc[1:10].mean()
np.testing.assert_allclose(smoothed, strength.rolling(11, center=True, min_periods=1).mean(), atol=1e-12, rtol=0)
overview = smoothed - initial
overview.index.name = 'simulation_year'
forcing = pd.read_csv(source / 'co2_forcing.csv').set_index('simulation_year')
overview['co2_ppm'] = forcing.co2_ppm
overview['mmm'] = overview[metadata['models']].mean(axis=1)
overview['std'] = overview[metadata['models']].std(axis=1, ddof=0)
overview['lower'] = overview.mmm - overview['std']
overview['upper'] = overview.mmm + overview['std']
previous_overview = pd.read_csv(source / 'overview_anomaly.csv').set_index('simulation_year')
np.testing.assert_allclose(overview.loc[1:275, metadata['models']], previous_overview.loc[1:275, metadata['models']], atol=1e-12, rtol=0)
overview.to_csv(destination / 'overview_anomaly.csv')
caption = ('AMOC evolution and hysteresis pathways at the native latitude nearest40N. '
           'Panel a shows nine-model initial-decade anomalies, the equal-weight mean and intermodel population SD, '
           'with prescribed CO2 on a separate right axis. Panels b,d,f are idealized schematics of recovery '
           'with small A, recovery with large A, and non-recovery. Only panel d labels the shaded hysteresis '
           'area A, with an arrow from its blank lower portion; panel f retains open-path shading without an area label. The small-A schematic shows gently curved, separated branches and '
           'encloses a narrow closed area. Panels c,e,g overlay2,3,4actual model trajectories, respectively. '
           'Each model keeps its panel-a color for its complete CO2 rise-and-fall trajectory; all paths are '
           'solid and arrows show time direction. Actual trajectories remain the saved40N11yr centred '
           'means normalized by each model\'s raw years1-10 strength. No actual endpoint is forced to close '
           'and no actual small nonzero area is replaced by zero; the group names are descriptive author labels. '
           'Schematics are conceptual rather than fitted model curves or equilibrium branches. '
           'For panel a, the final five smoothed years now exclude post280 stabilization data, consistently with the model paths. '
           'The author designates this40N version as the primary manuscript figure;26.5N is a comparison version.')
caption += ' Parentheses in the panel-a legend give mean AMOC strength over experiment years1-10 in Sv at the same native latitude; the9-model mean is the equal-weight average of those means.'
metadata.update({'panel_count': 7, 'group_names': group_names, 'model_colors': colors_by_model,
                 'initial_legend': 'Parentheses: raw forced AMOC averaged over experiment years1-10, Sv; same reference as anomalies and normalization.',
                 'primary_manuscript_version': True,
                 'caption': caption, 'layout': 'One centred moderate-width overview above three rows of schematic plus overlaid member trajectories.',
                 'revision': 'Remove oldb concept; collapse nine member panels into3model overlays; renumbera-g. All actual40N values retained.',
                 'small_area_schematic': 'Panel b shows a narrow closed hysteresis loop with gentle curvature; other2schematic shapes retained.',
                 'overview_data_range_sv': [float(overview[metadata['models']].min().min()), float(overview[metadata['models']].max().max())],
                 'smoothing': 'Same40N11yr centred min_periods1 within years1-280 for panel a and all paths; initial reference raw1-10.',
                 'primary_manuscript_version': True,
                 'dimensions': 'model_paths:model9 xyear280; schematic_paths:panel3 xphase2 xnode601; no raw NetCDF/MAT read.'})
source_names = ['metadata.json', 'panels.csv', 'path_coordinates.csv', 'path_arrows.csv',
                'overview_anomaly.csv', 'co2_forcing.csv', 'model_hysteresis_metrics.csv']
additional_sources = [annual_source / 'nine_model_normalized_series.csv',
                      annual_source / 'nine_model_endpoint_diagnostics.csv',
                      latitude_source / 'CMIP6_AMOC_9models_processing_metadata.csv', initial_source]
for key in ['layout_revision', 'concept_source_scope', 'picontrol_legend']:
    metadata.pop(key, None)
metadata['sources'] = [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
                       for path in [source / name for name in source_names] + additional_sources]
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
(destination / 'caption.txt').write_text(caption + '\n', encoding='utf-8')
print('Seven panels; nine complete40N trajectories preserved.', flush=True)
print(destination.resolve(), flush=True)
