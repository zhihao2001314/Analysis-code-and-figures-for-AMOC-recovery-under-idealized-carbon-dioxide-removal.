# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import shutil

import numpy as np
import pandas as pd


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
fig1_source = root / 'output/s01_amoc/n66_amoc_overview_26p5'
annual_source = root / 'output/s01_amoc/n22_amoc_recovery_definition_audit'
destination = root / 'output/s01_amoc/n68_four_row_pathway_preview'
assert not destination.exists(), '正式结果已存在，须先归档。'
members = [['GFDL-ESM4', 'NorESM2-LM', 'CAS-ESM2-0'],
           ['CESM2', 'CNRM-ESM2-1', 'UKESM1-0-LL'],
           ['ACCESS-ESM1-5', 'CanESM5', 'MIROC-ES2L']]
group_names = ['Substantial recovery', 'Non-recovery', 'Larger initial biases']
row_panels = [['c', 'd', 'e', 'f'], ['g', 'h', 'i', 'j'], ['k', 'l', 'm', 'n']]
anchor_x = np.array([0., .08, .16, .30, .50, .70, .84, 1.])
drawing_x = np.linspace(0, 1, 601)
coordinate_parts, arrow_rows, panel_rows, summary_rows, anchor_rows, actual_mean_parts = [], [], [], [], [], []
recovery_shapes = {}


# %% 3. 读取既有首行、九模式路径和初末强度

annual = pd.read_csv(annual_source / 'nine_model_normalized_series.csv')
endpoints = pd.read_csv(annual_source / 'nine_model_endpoint_diagnostics.csv').set_index('model')
old_paths = pd.read_csv(fig1_source / 'path_coordinates.csv')
old_arrows = pd.read_csv(fig1_source / 'path_arrows.csv')
old_metadata = json.loads((fig1_source / 'metadata.json').read_text(encoding='utf-8'))
reference_co2 = old_metadata['reference_co2_ppm']
maximum_co2 = old_metadata['maximum_co2_ppm']
coordinate_parts.append(old_paths.loc[old_paths.panel.eq('b')].copy())
arrow_rows.extend(old_arrows.loc[old_arrows.panel.eq('b')].to_dict('records'))
panel_rows.append({'panel': 'b', 'row': 0, 'column': 2, 'kind': 'concept', 'model': '', 'label': 'Conceptual'})


# %% 4. 三行示意和九条实际模式路径

for row_index, (models, group_name, letters) in enumerate(zip(members, group_names, row_panels), 1):
    records = annual.loc[annual.model.isin(models)].copy()
    assert len(records) == 840
    group_mean = records.groupby('simulation_year', sort=True)[['co2_ppm', 'smooth_11yr_percent_initial']].mean()
    group_mean['relative_amoc'] = group_mean.smooth_11yr_percent_initial / 100
    group_mean['normalized_co2'] = (group_mean.co2_ppm - reference_co2) / (maximum_co2 - reference_co2)
    saved_mean = group_mean.reset_index()
    saved_mean['group'] = group_name
    actual_mean_parts.append(saved_mean)
    final_ratio = float(endpoints.loc[models, 'endpoint_ratio'].mean())
    peak = float(group_mean.loc[140, 'relative_amoc'])
    summary_rows.append({'row': row_index, 'group': group_name, 'models': '; '.join(models),
                         'mean_member_endpoint_ratio': final_ratio,
                         'minimum_member_endpoint_ratio': float(endpoints.loc[models, 'endpoint_ratio'].min()),
                         'maximum_member_endpoint_ratio': float(endpoints.loc[models, 'endpoint_ratio'].max()),
                         'n_recovered': int((endpoints.loc[models, 'endpoint_ratio'] >= .8).sum()),
                         'peak_relative_group_mean': peak})
    schematic_panel = letters[0]
    panel_rows.append({'panel': schematic_panel, 'row': row_index, 'column': 0,
                       'kind': 'schematic', 'model': '', 'label': group_name})
    for phase in ['increase', 'decrease']:
        if row_index <= 2:
            old_panel = 'c' if row_index == 1 else 'f'
            original = old_paths.loc[old_paths.panel.eq(old_panel) & old_paths.phase.eq(phase)].sort_values('node')
            old_y = original.relative_amoc.to_numpy()
            old_peak = float(old_y[-1])
            if row_index == 1 or phase == 'increase':
                values = peak + (old_y - old_peak) * (1 - peak) / (1 - old_peak)
            else:
                old_final = float(old_y[0])
                values = old_y + (final_ratio - old_final) * (1 - drawing_x) + (peak - old_peak) * drawing_x
            for x, value in zip(anchor_x, np.interp(anchor_x, drawing_x, values)):
                anchor_rows.append({'panel': schematic_panel, 'phase': phase, 'normalized_co2': x, 'relative_amoc': value})
        else:
            # k复用c的平滑轮廓；只映射概念端点，保留低CO₂端开口。
            reference_shape = recovery_shapes[phase]
            reference_peak = float(reference_shape[-1])
            start_value = 1. if phase == 'increase' else final_ratio
            values = peak + (reference_shape - reference_peak) * (start_value - peak) / (1 - reference_peak)
            assert np.all(np.diff(values) <= 1e-12)
            for x, value in zip(anchor_x, np.interp(anchor_x, drawing_x, values)):
                anchor_rows.append({'panel': schematic_panel, 'phase': phase, 'normalized_co2': x, 'relative_amoc': value})
        if row_index == 1:
            recovery_shapes[phase] = values.copy()
        part = pd.DataFrame({'panel': schematic_panel, 'phase': phase, 'node': np.arange(601),
                             'simulation_year': np.nan,
                             'co2_ppm': reference_co2 + drawing_x * (maximum_co2 - reference_co2),
                             'relative_amoc': values})
        coordinate_parts.append(part)
        start, end = [270, 320] if phase == 'increase' else [330, 280]
        arrow_rows.append({'panel': schematic_panel, 'phase': phase,
                           'start_x': float(part.co2_ppm.iloc[start]), 'end_x': float(part.co2_ppm.iloc[end]),
                           'start_y': float(values[start]), 'end_y': float(values[end])})
    for column, model in enumerate(models, 1):
        letter = letters[column]
        panel_rows.append({'panel': letter, 'row': row_index, 'column': column,
                           'kind': 'model', 'model': model, 'label': ''})
        history = annual.loc[annual.model.eq(model)].sort_values('simulation_year').set_index('simulation_year')
        for phase in ['increase', 'decrease']:
            branch = history.loc[history.index <= 140 if phase == 'increase' else history.index >= 140]
            part = pd.DataFrame({'panel': letter, 'phase': phase, 'node': np.arange(len(branch)),
                                 'simulation_year': branch.index, 'co2_ppm': branch.co2_ppm.to_numpy(),
                                 'relative_amoc': branch.smooth_11yr_percent_initial.to_numpy() / 100})
            coordinate_parts.append(part)
            points = history.loc[[65, 75] if phase == 'increase' else [210, 220]]
            arrow_rows.append({'panel': letter, 'phase': phase,
                               'start_x': float(points.co2_ppm.iloc[0]), 'end_x': float(points.co2_ppm.iloc[1]),
                               'start_y': float(points.smooth_11yr_percent_initial.iloc[0] / 100),
                               'end_y': float(points.smooth_11yr_percent_initial.iloc[1] / 100)})


# %% 5. 逐点核对并保存试排输入

paths = pd.concat(coordinate_parts, ignore_index=True)
panels = pd.DataFrame(panel_rows)
assert len(panels) + 1 == 14 and panels.kind.eq('model').sum() == 9
assert np.isfinite(paths[['co2_ppm', 'relative_amoc']]).all().all()
assert paths.relative_amoc.between(0, 1.2).all()
for row in panels.loc[panels.kind.eq('model')].itertuples(index=False):
    history = annual.loc[annual.model.eq(row.model)].set_index('simulation_year')
    part = paths.loc[paths.panel.eq(row.panel)]
    np.testing.assert_allclose(part.relative_amoc, history.loc[part.simulation_year.astype(int), 'smooth_11yr_percent_initial'].to_numpy() / 100, atol=1e-12)
np.testing.assert_array_equal(paths.loc[paths.panel.eq('b'), ['co2_ppm', 'relative_amoc']],
                              old_paths.loc[old_paths.panel.eq('b'), ['co2_ppm', 'relative_amoc']])
destination.mkdir(parents=True)
for name in ['overview_anomaly.csv', 'co2_forcing.csv']:
    shutil.copy2(fig1_source / name, destination / name)
paths.to_csv(destination / 'path_coordinates.csv', index=False)
pd.DataFrame(arrow_rows).to_csv(destination / 'path_arrows.csv', index=False)
panels.to_csv(destination / 'panels.csv', index=False)
pd.DataFrame(summary_rows).to_csv(destination / 'group_summary.csv', index=False)
pd.DataFrame(anchor_rows).to_csv(destination / 'schematic_anchors.csv', index=False)
pd.concat(actual_mean_parts, ignore_index=True).to_csv(destination / 'actual_normalized_group_means.csv', index=False)
caption = ('Four-row layout preview. The first row retains the existing nine-model26.5°N AMOC-anomaly data and conceptual hysteresis-area diagram. '
           'The remaining rows contain one schematic and three actual40°N AMOC–CO₂ model pathways each. '
           'Substantial recovery: GFDL-ESM4, NorESM2-LM, CAS-ESM2-0. Non-recovery: CESM2, CNRM-ESM2-1, UKESM1-0-LL. '
           'Larger initial biases: ACCESS-ESM1-5, CanESM5, MIROC-ES2L; this is an initial-state grouping containing two recovery outcomes, not a third homogeneous recovery class. '
           'The third schematic reuses the smooth outline of the substantial-recovery guide, with its peak and open final endpoint set by the three-member reference values; its final conceptual value is the mean member recovery ratio. '
           'The first recovery schematic remains an idealized closed-loop guide; the non-recovery guide retains its previous outline with updated three-member endpoint references. '
           'Schematics are not simulations, equilibrium branches or numerical confidence intervals. All nine real model paths retain the original initial-decade normalization and11-year centred running means without extra smoothing, drift correction or forced closure. '
           'Blue and orange solid lines represent increasing and decreasing CO₂, with arrows indicating time direction. '
           'Panels a and b each span two column positions; there are14 lettered panels. This preview does not replace the original Fig1 or downstream statistical composites.')
(destination / 'caption.txt').write_text(caption + '\n', encoding='utf-8')
sources = [fig1_source / name for name in ['overview_anomaly.csv', 'co2_forcing.csv', 'path_coordinates.csv', 'path_arrows.csv', 'metadata.json']]
sources.extend([annual_source / 'nine_model_normalized_series.csv', annual_source / 'nine_model_endpoint_diagnostics.csv'])
metadata = {'preview_only': True, 'panel_count': 14, 'layout': '4rows x4positions; a spans columns1–2 and b columns3–4 with equal widths; three lower rows schematic+3models.',
            'models': old_metadata['models'], 'members': members, 'group_names': group_names,
            'overview_y_limits_sv': old_metadata['overview_y_limits_sv'], 'overview_y_ticks_sv': old_metadata['overview_y_ticks_sv'],
            'reference_co2_ppm': reference_co2, 'maximum_co2_ppm': maximum_co2,
            'normalization': old_metadata['path_normalization'], 'caption': caption,
            'dimensions': 'overview simulation_year280 x9models; path long table panel13 xphase2 xnode; actual path node140/141, schematics601.',
            'third_schematic': 'Reuse panel c smooth outlines with affine vertical endpoint mapping. Initial1, decreasing-branch endpoint meanR0.8607, peak0.49946; monotone branches, retained low-CO2 opening, not a fitted model path or homogeneous recovery claim.',
            'sources': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in sources]}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print(pd.DataFrame(summary_rows).to_string(index=False), flush=True)
