# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import shutil

import numpy as np
import pandas as pd


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
source = root / 'output/s01_amoc/n70_four_row_pathway_40n'
annual_source = root / 'output/s01_amoc/n22_amoc_recovery_definition_audit'
destination = root / 'output/s01_amoc/n73_hysteresis_morphology_groups'
assert not destination.exists() or not any(destination.iterdir()), '正式输出存在，须先归档。'
members = [['ACCESS-ESM1-5', 'CAS-ESM2-0'],
           ['GFDL-ESM4', 'NorESM2-LM', 'CanESM5'],
           ['CESM2', 'CNRM-ESM2-1', 'UKESM1-0-LL', 'MIROC-ES2L']]
group_names = ['Closed, small area', 'Near-closed, larger area', 'Open, larger area']
letters = [['c', 'd', 'e'], ['f', 'g', 'h', 'i'], ['j', 'k', 'l', 'm', 'n']]
schematic_sources = ['k', 'c', 'g']
schematic_x = np.linspace(0., 1., 601)
anchor_x = np.array([0., .08, .16, .3, .5, .7, .84, 1.])
small_gap_amplitude = .016
coordinate_parts, arrow_records, panel_records = [], [], []
anchor_records, metric_records, summary_records, mean_parts = [], [], [], []


# %% 3. 读取已存序列和路径

metadata = json.loads((source / 'metadata.json').read_text(encoding='utf-8'))
old_paths = pd.read_csv(source / 'path_coordinates.csv')
old_arrows = pd.read_csv(source / 'path_arrows.csv')
old_panels = pd.read_csv(source / 'panels.csv', keep_default_na=False)
annual = pd.read_csv(annual_source / 'nine_model_normalized_series.csv')
endpoints = pd.read_csv(annual_source / 'nine_model_endpoint_diagnostics.csv').set_index('model')
reference_co2, maximum_co2 = metadata['reference_co2_ppm'], metadata['maximum_co2_ppm']
coordinate_parts.append(old_paths.loc[old_paths.panel.eq('b')].copy())
arrow_records.extend(old_arrows.loc[old_arrows.panel.eq('b')].to_dict('records'))
panel_records.append({'panel': 'b', 'row': 0, 'column': 2, 'kind': 'concept',
                      'model': '', 'label': 'Conceptual', 'old_panel': 'b'})


# %% 4. 调组并逐点保留九模式轨迹

for row, (models, group_name, panel_letters) in enumerate(zip(members, group_names, letters), 1):
    selected = annual.loc[annual.model.isin(models)].copy()
    assert len(selected) == len(models) * 280
    mean = selected.groupby('simulation_year', sort=True)[['co2_ppm', 'smooth_11yr_percent_initial']].mean()
    mean['relative_amoc'] = mean.smooth_11yr_percent_initial / 100.
    mean['group'] = group_name
    mean_parts.append(mean.reset_index())
    peak = float(mean.loc[140, 'relative_amoc'])
    final_ratio = float(endpoints.loc[models, 'endpoint_ratio'].mean())
    summary_records.append({'row': row, 'group': group_name, 'model_count': len(models),
                            'models': '; '.join(models), 'mean_endpoint_ratio': final_ratio,
                            'peak_relative_group_mean': peak})
    panel_records.append({'panel': panel_letters[0], 'row': row, 'column': 0,
                          'kind': 'schematic', 'model': '', 'label': group_name,
                          'old_panel': schematic_sources[row - 1]})
    for column, model in enumerate(models, 1):
        original_panel = old_panels.loc[old_panels.model.eq(model), 'panel'].iloc[0]
        letter = panel_letters[column]
        panel_records.append({'panel': letter, 'row': row, 'column': column, 'kind': 'model',
                              'model': model, 'label': '', 'old_panel': original_panel})
        paths = old_paths.loc[old_paths.panel.eq(original_panel)].copy()
        paths['panel'] = letter
        coordinate_parts.append(paths)
        arrows = old_arrows.loc[old_arrows.panel.eq(original_panel)].copy()
        arrows['panel'] = letter
        arrow_records.extend(arrows.to_dict('records'))
        up = paths.loc[paths.phase.eq('increase')].sort_values('co2_ppm')
        down = paths.loc[paths.phase.eq('decrease')].sort_values('co2_ppm')
        # 仅在两阶段共同浓度范围内线性配对；积分绝对差，避免交叉抵消。
        matched_x = np.unique(np.r_[up.co2_ppm.to_numpy(), down.co2_ppm.to_numpy()])
        matched_x = matched_x[(matched_x >= max(up.co2_ppm.min(), down.co2_ppm.min()))
                              & (matched_x <= min(up.co2_ppm.max(), down.co2_ppm.max()))]
        differences = (np.interp(matched_x, up.co2_ppm, up.relative_amoc)
                       - np.interp(matched_x, down.co2_ppm, down.relative_amoc))
        # 增加分支交叉点，使绝对差在线性段上的积分也严格分段。
        crossings = np.flatnonzero(differences[:-1] * differences[1:] < 0.)
        crossing_x = (matched_x[crossings] - differences[crossings]
                      * np.diff(matched_x)[crossings] / np.diff(differences)[crossings])
        matched_x = np.unique(np.r_[matched_x, crossing_x])
        separation = np.abs(np.interp(matched_x, up.co2_ppm, up.relative_amoc)
                            - np.interp(matched_x, down.co2_ppm, down.relative_amoc))
        area = float(np.trapezoid(separation, matched_x))
        metric_records.append({'model': model, 'group': group_name, 'panel': letter,
                               'area_relative_amoc_ppm': area,
                               'area_divided_by_co2_span': area / (matched_x[-1] - matched_x[0]),
                               'co2_min_ppm': matched_x[0], 'co2_max_ppm': matched_x[-1],
                               'endpoint_ratio_raw_decades': float(endpoints.loc[model, 'endpoint_ratio']),
                               'increase_low_co2_smooth': float(up.relative_amoc.iloc[0]),
                               'decrease_low_co2_smooth': float(down.relative_amoc.iloc[0])})
    # 示意图单独建模，实际模式轨迹不做强制闭合或额外平滑。
    for phase in ['increase', 'decrease']:
        if row == 1:
            # 近零面积组：近直线返回、两端闭合；只改示意图，不改实际模式。
            values = 1. - (1. - peak) * schematic_x
            if phase == 'decrease':
                values = values - small_gap_amplitude * np.sin(np.pi * schematic_x)
        else:
            original = old_paths.loc[old_paths.panel.eq(schematic_sources[row - 1])
                                     & old_paths.phase.eq(phase)].sort_values('node')
            old_values = original.relative_amoc.to_numpy()
            old_peak, old_start = float(old_values[-1]), float(old_values[0])
            start = 1. if phase == 'increase' else final_ratio
            values = peak + (old_values - old_peak) * (start - peak) / (old_start - old_peak)
        assert len(values) == 601
        if row == 1:
            assert np.all(np.diff(values) <= 1e-10)
        letter = panel_letters[0]
        paths = pd.DataFrame({'panel': letter, 'phase': phase, 'node': np.arange(601),
                              'simulation_year': np.nan,
                              'co2_ppm': reference_co2 + schematic_x * (maximum_co2 - reference_co2),
                              'relative_amoc': values})
        coordinate_parts.append(paths)
        for x, y in zip(anchor_x, np.interp(anchor_x, schematic_x, values)):
            anchor_records.append({'panel': letter, 'phase': phase, 'normalized_co2': x,
                                   'relative_amoc': y})
        start_node, end_node = ([200, 250] if phase == 'increase' else [410, 360]) if row == 1 else ([270, 320] if phase == 'increase' else [330, 280])
        arrow_records.append({'panel': letter, 'phase': phase,
                              'start_x': float(paths.co2_ppm.iloc[start_node]),
                              'end_x': float(paths.co2_ppm.iloc[end_node]),
                              'start_y': float(values[start_node]), 'end_y': float(values[end_node])})


# %% 5. 核对并保存分组绘图输入

panels = pd.DataFrame(panel_records)
coordinates = pd.concat(coordinate_parts, ignore_index=True)
metrics = pd.DataFrame(metric_records)
assert len(panels) + 1 == 14 and panels.kind.eq('model').sum() == 9
assert panels.loc[panels.kind.eq('model')].groupby('row').size().tolist() == [2, 3, 4]
assert panels.loc[panels.kind.eq('model'), 'model'].nunique() == 9
for record in panels.loc[panels.kind.eq('model')].itertuples(index=False):
    actual = coordinates.loc[coordinates.panel.eq(record.panel)].drop(columns='panel').reset_index(drop=True)
    original = old_paths.loc[old_paths.panel.eq(record.old_panel)].drop(columns='panel').reset_index(drop=True)
    pd.testing.assert_frame_equal(actual, original)
pd.testing.assert_frame_equal(coordinates.loc[coordinates.panel.eq('b')].reset_index(drop=True),
                              old_paths.loc[old_paths.panel.eq('b')].reset_index(drop=True))
assert np.isfinite(coordinates[['co2_ppm', 'relative_amoc']]).all().all()
assert coordinates.relative_amoc.between(0, 1.2).all()
destination.mkdir(parents=True, exist_ok=True)
for name in ['overview_anomaly.csv', 'co2_forcing.csv']:
    shutil.copy2(source / name, destination / name)
coordinates.to_csv(destination / 'path_coordinates.csv', index=False)
pd.DataFrame(arrow_records).to_csv(destination / 'path_arrows.csv', index=False)
panels.to_csv(destination / 'panels.csv', index=False)
pd.DataFrame(anchor_records).to_csv(destination / 'schematic_anchors.csv', index=False)
pd.DataFrame(summary_records).to_csv(destination / 'group_summary.csv', index=False)
metrics.to_csv(destination / 'model_hysteresis_metrics.csv', index=False)
pd.concat(mean_parts, ignore_index=True).to_csv(destination / 'actual_normalized_group_means.csv', index=False)


# %% 6. 保存定义和来源

caption = ('AMOC evolution and transient hysteresis pathways at40N. Panel a retains the nine-model anomaly '
           'series, equal-weight mean and intermodel population SD, with prescribed CO2 on the right axis; '
           'b explains the interbranch hysteresis area. The three lower rows are a descriptive morphology '
           'preview: closed with small area (ACCESS-ESM1-5, CAS-ESM2-0; c-e), near-closed with larger '
           'area (GFDL-ESM4, NorESM2-LM, CanESM5; f-i), and open with larger area '
           '(CESM2, CNRM-ESM2-1, UKESM1-0-LL, MIROC-ES2L; j-n). '
           'Panels c,f,j are idealized schematics, not model averages or equilibrium branches. '
           'The first schematic uses almost coincident straight return paths; actual ACCESS and CAS areas '
           'are small relative to the other models, not zero. Near-closed denotes an approximate return, '
           'not exact endpoint identity. Larger-area categories overlap in area and are also distinguished '
           'by endpoint recovery; no area cutoff or new inferential classification is claimed. '
           'All nine actual model trajectories and their arrows are unchanged from the previous40N preview, '
           'using the original initial-decade normalization and11-year running means. '
           'Blue and orange solid paths indicate increasing and decreasing CO2, respectively. '
           'This visual regrouping does not change downstream recovery composites or statistical analyses.')
metadata.update({'members': members, 'group_names': group_names, 'caption': caption,
                 'layout': 'Top a/b retain equal two-column spans; lower rows have schematic+2/3/4models. Common schematic width; model widths adapt to row count.',
                 'revision': 'User morphology regrouping: old f/m exchanged, old n moved into the open group, old third group moved to first lower row.',
                 'grouping_basis': 'Descriptive trajectory separation and endpoint return; replaces initial-bias grouping in this preview only. No universal area cutoff.',
                 'small_area_schematic': 'Straight mean slope with0.016*sin(pi*x) branch separation, closed endpoints; idealized, not fitted.',
                 'metric_definition': 'Integral abs(increase-decrease) over common CO2 interval, piecewise-linear pairing on union of original CO2 nodes and all branch crossings. Also divided by CO2 span for dimensionless comparison; uses saved normalized11yr paths.',
                 'dimensions': 'overview simulation_year280 x9model; path long panel13 xphase2 xnode; actual paths140/141nodes, schematics601nodes. No raw NetCDF/MAT read.'})
metadata.pop('third_schematic', None)
source_paths = [source / name for name in ['overview_anomaly.csv', 'co2_forcing.csv', 'path_coordinates.csv', 'path_arrows.csv', 'panels.csv', 'metadata.json']]
source_paths.extend([annual_source / 'nine_model_normalized_series.csv', annual_source / 'nine_model_endpoint_diagnostics.csv'])
metadata['sources'] = [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in source_paths]
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
(destination / 'caption.txt').write_text(caption + '\n', encoding='utf-8')
print(metrics[['model', 'group', 'area_divided_by_co2_span']].to_string(index=False), flush=True)
print(destination.resolve(), flush=True)
