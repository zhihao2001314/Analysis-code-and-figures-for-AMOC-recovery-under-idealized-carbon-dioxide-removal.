# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
destination = root / 'output/s01_amoc/n62_amoc_pathway_composite'
assert not destination.exists() or not any(destination.iterdir()), '正式结果存在，须先归档。'
initial_directory = root / 'output/amoc_s01_n01/l_40'
path_directory = root / 'output/s01_amoc/n57_amoc_schematic_model_panels'
annual_directory = root / 'output/s01_amoc/n22_amoc_recovery_definition_audit'
models = ['ACCESS-ESM1-5', 'CanESM5', 'CESM2', 'CNRM-ESM2-1', 'GFDL-ESM4',
          'MIROC-ES2L', 'NorESM2-LM', 'UKESM1-0-LL', 'CAS-ESM2-0']
reference_co2 = 284.7
maximum_co2 = reference_co2 * 1.01 ** 140
panel_map = {1: 'd', 2: 'e', 4: 'g', 5: 'h'}
selected_members = [['GFDL-ESM4', 'NorESM2-LM'], ['CESM2', 'CNRM-ESM2-1']]
coordinate_parts, arrow_parts, anchor_rows = [], [], []


# %% 3. 读取并核验已有异常与实际模式路径

overview = pd.read_csv(initial_directory / 'amoc_initial_reference_plot_data.csv')
initial_metadata = json.loads((initial_directory / 'amoc_initial_reference_metadata.json').read_text(encoding='utf-8'))
coordinates = pd.read_csv(path_directory / 'coordinates.csv')
arrows = pd.read_csv(path_directory / 'arrows.csv')
panels = pd.read_csv(path_directory / 'panels.csv', keep_default_na=False)
source_metadata = json.loads((path_directory / 'metadata.json').read_text(encoding='utf-8'))
annual = pd.read_csv(annual_directory / 'nine_model_normalized_series.csv')
overview = overview.loc[overview.simulation_year.between(1, 280)].copy()
assert overview.simulation_year.to_list() == list(range(1, 281))
assert np.isfinite(overview[models + ['mmm', 'std', 'lower', 'upper']]).all().all()
assert initial_metadata['reference'] == 'each model initial years1-10 mean'
np.testing.assert_allclose(overview[models].mean(axis=1), overview.mmm, atol=1e-12, rtol=0)
np.testing.assert_allclose(overview[models].std(axis=1, ddof=0), overview['std'], atol=1e-12, rtol=0)

for old_panel, label in panel_map.items():
    selected = coordinates.loc[coordinates.panel.eq(old_panel)].copy()
    selected['panel'] = label
    coordinate_parts.append(selected)
    selected_arrows = arrows.loc[arrows.panel.eq(old_panel)].copy()
    selected_arrows['panel'] = label
    arrow_parts.append(selected_arrows)
    if old_panel not in [0, 3]:
        model = panels.loc[panels.panel.eq(old_panel), 'model'].iloc[0]
        history = annual.loc[annual.model.eq(model)].set_index('simulation_year')
        for phase in ['increase', 'decrease']:
            branch = selected.loc[selected.phase.eq(phase)]
            expected = history.loc[branch.simulation_year.astype(int)]
            np.testing.assert_allclose(branch.relative_amoc, expected.smooth_11yr_percent_initial / 100, atol=1e-14, rtol=0)
            np.testing.assert_allclose(branch.co2_ppm, expected.co2_ppm, atol=1e-10, rtol=0)


# %% 4. 按手绘重建菱形、窄斜环和平台开口环，仅为概念示意

concept_x = np.linspace(0, 1, 601)
recovery_reference = coordinates.loc[coordinates.panel.eq(0)].sort_values('node')
incomplete_reference = coordinates.loc[coordinates.panel.eq(3)].sort_values('node')
recovery_peak = recovery_reference.loc[recovery_reference.phase.eq('increase'), 'relative_amoc'].iloc[-1]
incomplete_final = incomplete_reference.loc[incomplete_reference.phase.eq('decrease'), 'relative_amoc'].iloc[0]
incomplete_peak = incomplete_reference.loc[incomplete_reference.phase.eq('increase'), 'relative_amoc'].iloc[-1]
schematic_anchors = {
    'b': {'x': [0, .08, .20, .35, .50, .65, .80, .92, 1],
          'increase': [1, 1, 1, .96, .85, .64, .45, .34, .30],
          'decrease': [1, .96, .85, .66, .45, .34, .30, .30, .30]},
    'c': {'x': [0, .08, .16, .30, .50, .70, .84, 1],
          'increase': [1, .985, .94, .84, .69, .54, .44, recovery_peak],
          'decrease': [1, .91, .80, .59, .43, .40, .39, recovery_peak]},
    'f': {'x': [0, .08, .16, .30, .50, .70, .84, 1],
          'increase': [1, 1, .94, .78, .57, .40, .32, incomplete_peak],
          'decrease': [incomplete_final, .34, .325, .31, .30, .29, .285, incomplete_peak]}}
for label, anchors in schematic_anchors.items():
    for phase in ['increase', 'decrease']:
        values = PchipInterpolator(anchors['x'], anchors[phase])(concept_x)
        branch = pd.DataFrame({'panel': label, 'phase': phase, 'node': np.arange(601),
                               'simulation_year': np.nan,
                               'co2_ppm': reference_co2 + concept_x * (maximum_co2 - reference_co2),
                               'relative_amoc': values})
        coordinate_parts.append(branch)
        for anchor_index, (anchor_x, anchor_y) in enumerate(zip(anchors['x'], anchors[phase])):
            anchor_rows.append({'panel': label, 'phase': phase, 'anchor': anchor_index,
                                'normalized_co2': anchor_x, 'relative_amoc': anchor_y})
        indices = [270, 320] if phase == 'increase' else [330, 280]
        arrow_parts.append(pd.DataFrame([{'panel': label, 'phase': phase,
                                          'start_x': branch.co2_ppm.iloc[indices[0]],
                                          'end_x': branch.co2_ppm.iloc[indices[1]],
                                          'start_y': values[indices[0]], 'end_y': values[indices[1]]}]))

paths = pd.concat(coordinate_parts, ignore_index=True)
path_arrows = pd.concat(arrow_parts, ignore_index=True)
path_rows = [
    {'panel': 'b', 'kind': 'concept', 'row': 0, 'column': 2, 'label': 'Conceptual', 'model': ''},
    {'panel': 'c', 'kind': 'schematic', 'row': 1, 'column': 0, 'label': 'Substantial recovery', 'model': ''},
    {'panel': 'd', 'kind': 'model', 'row': 1, 'column': 1, 'label': 'GFDL-ESM4', 'model': 'GFDL-ESM4'},
    {'panel': 'e', 'kind': 'model', 'row': 1, 'column': 2, 'label': 'NorESM2-LM', 'model': 'NorESM2-LM'},
    {'panel': 'f', 'kind': 'schematic', 'row': 2, 'column': 0, 'label': 'Incomplete recovery', 'model': ''},
    {'panel': 'g', 'kind': 'model', 'row': 2, 'column': 1, 'label': 'CESM2', 'model': 'CESM2'},
    {'panel': 'h', 'kind': 'model', 'row': 2, 'column': 2, 'label': 'CNRM-ESM2-1', 'model': 'CNRM-ESM2-1'}]


# %% 5. 保存绘图输入、面积定义及来源

destination.mkdir(parents=True, exist_ok=True)
overview.to_csv(destination / 'overview_anomaly.csv', index=False)
co2_forcing = annual.loc[annual.model.eq(models[0]), ['simulation_year', 'co2_ppm']].sort_values('simulation_year')
assert co2_forcing.simulation_year.to_list() == list(range(1, 281))
for model in models[1:]:
    other_forcing = annual.loc[annual.model.eq(model)].sort_values('simulation_year')
    np.testing.assert_allclose(other_forcing.co2_ppm, co2_forcing.co2_ppm, atol=1e-10, rtol=0)
co2_forcing.to_csv(destination / 'co2_forcing.csv', index=False)
paths.to_csv(destination / 'path_coordinates.csv', index=False)
path_arrows.to_csv(destination / 'path_arrows.csv', index=False)
pd.DataFrame(path_rows).to_csv(destination / 'panels.csv', index=False)
pd.DataFrame(anchor_rows).to_csv(destination / 'schematic_anchors.csv', index=False)
sources = [initial_directory / 'amoc_initial_reference_plot_data.csv',
           initial_directory / 'amoc_initial_reference_metadata.json',
           annual_directory / 'nine_model_normalized_series.csv']
sources.extend(path_directory / name for name in ['coordinates.csv', 'arrows.csv', 'panels.csv', 'metadata.json'])
caption = ('(a) Nine-model AMOC anomalies relative to each model\'s initial years 1–10 at the grid latitude nearest 40°N, '
           'using the Atlantic–Arctic overturning maximum below 500 m. Lines retain the existing 11-year centred running means; '
           'the black line and grey shading are the equal-model mean and ±1 inter-model population standard deviation, not a confidence interval. '
           'The grey dashed line uses the right axis for the prescribed atmospheric CO2 concentration. '
           '(b) Conceptual transient hysteresis in an AMOC–CO2 diagram; Initial CO2 and Peak CO2 identify the forcing endpoints. '
           'The shaded hysteresis area A is the absolute separation of the increasing- '
           'and decreasing-CO2 branches integrated over the common linear CO2 range, not a time lag or proof of equilibrium bistability. '
           '(c, f) Idealized schematic guides for substantial and incomplete recovery, with the corresponding model trajectories '
           'in (d, e) and (g, h). Recovery groups retain the original late-decade/initial-decade criterion R ≥ 0.8 versus R < 0.8. '
           'All AMOC–CO2 curves are solid; arrows indicate temporal direction. Schematic starting strength is set to 1; '
           'real model curves retain their instantaneous smoothed endpoints. '
           'The idealized path in (c) closes at both forcing endpoints to illustrate a recoverable trajectory; its shape is not fitted to the model curves. '
           'The schematic in (f) retains its group-mean final-to-initial strength ratio. '
           'Schematic and actual panels share physical CO2 units; schematic endpoint C0 is conceptual. '
           'Only the four previously designated better-simulated models appear in the pathway rows; (a) retains all nine models.')
metadata = {'models': models, 'selected_members': selected_members,
            'layout': 'Eight lettered panels: a spans top-left two columns; b top-right; c/d/e and f/g/h below.',
            'reference_co2_ppm': reference_co2, 'maximum_co2_ppm': maximum_co2,
            'overview_reference': initial_metadata, 'overview_years': [1, 280],
            'dimensions': 'overview: simulation_year280 x model9 (wide CSV); paths: panel7 x phase2 x node (long CSV), schematic601 each; model140/141 each.',
            'no_added_correction_or_smoothing': True,
            'path_normalization': source_metadata['normalization'],
            'schematics': 'User-sketch-inspired PCHIP guides: b has exact central symmetry M_up(x)+M_down(1-x)=1.30 and offset short plateaus; c is a closed recoverable loop with a steep early decreasing-CO2 branch following the user sketch; f is an open declining/low-return loop with its high plateau shortened to normalized CO2 x=0.08. c closes at initial and peak CO2; f preserves its original group-mean endpoint ratio. Drawing knots are not fitted or simulated data, equilibrium branches, or tipping thresholds.',
            'area_name': 'Hysteresis area, A; explicitly transient branch separation, not a recovery timescale or evidence of equilibrium bistability.',
            'area_definition': 'Integral |M_increase(C)-M_decrease(C)| dC over the common linear CO2 range; relative-AMOC times ppm, unless both axes are separately normalized. No numerical empirical area statistic is newly estimated in this layout task.',
            'recovery_names': ['Substantial recovery', 'Incomplete recovery'],
            'scope': 'Two recovery outcomes illustrated by the original four better-simulated members; no grouping registry changes or claim that these represent all nine models.',
            'caption': caption,
            'sources': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in sources]}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
(destination / 'caption.txt').write_text(caption + '\n', encoding='utf-8')
print(destination, flush=True)
