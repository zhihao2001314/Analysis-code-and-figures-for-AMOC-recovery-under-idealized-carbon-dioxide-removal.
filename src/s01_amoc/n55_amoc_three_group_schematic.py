# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
source = root / 'output/s01_amoc/n22_amoc_recovery_definition_audit'
registry_path = root / 'output/s05_atmosphere/p05_amoc_recovery_composites/n49_group_naming/group_labels.json'
destination = root / 'output/s01_amoc/n55_amoc_three_group_schematic'
assert not destination.exists(), '正式输出已存在，需先归档。'
group_order = ['better_recovery', 'better_nonrecovery', 'larger_bias']
reference_co2 = 284.7
maximum_co2 = reference_co2 * 1.01 ** 140
anchor_x = np.array([0., .12, .28, .48, .72, 1.])
drawing_x = np.linspace(0, 1, 601)
arrow_indices = {'increase': [280, 330], 'decrease': [310, 260]}
coordinate_parts, anchor_rows, arrow_rows, summary_rows, member_parts = [], [], [], [], []


# %% 3. 读取既有八模式分组和实际路径

registry = json.loads(registry_path.read_text(encoding='utf-8'))
annual = pd.read_csv(source / 'nine_model_normalized_series.csv')
endpoints = pd.read_csv(source / 'nine_model_endpoint_diagnostics.csv')
assert annual.simulation_year.min() == 1 and annual.simulation_year.max() == 280
assert list(registry['members']['better_recovery']) == ['NorESM2-LM', 'GFDL-ESM4']
assert list(registry['members']['better_nonrecovery']) == ['CESM2', 'CNRM-ESM2-1']
assert len(registry['members']['larger_bias']) == 4


# %% 4. 从实际组路径提取稀疏锚点，生成简化示意线

for panel, group in enumerate(group_order):
    members = registry['members'][group]
    records = annual.loc[annual.model.isin(members)].copy()
    assert len(records) == len(members) * 280
    # 先按各模式初始十年归一化，再作模式等权平均；不重算上游11年平滑。
    group_series = records.groupby('simulation_year', sort=True)[['co2_ppm', 'smooth_11yr_percent_initial']].mean()
    group_series['relative_amoc'] = group_series.smooth_11yr_percent_initial / 100
    group_series['co2_coordinate'] = (group_series.co2_ppm - reference_co2) / (maximum_co2 - reference_co2)
    up = group_series.loc[1:140]
    down = group_series.loc[140:280].sort_values('co2_coordinate')
    member_endpoints = endpoints.loc[endpoints.model.isin(members)].copy()
    member_endpoints['group'] = group
    member_parts.append(member_endpoints[['model', 'member', 'group', 'initial_mean_sv', 'late_mean_sv', 'endpoint_ratio']])
    final_ratio = float(member_endpoints.endpoint_ratio.mean())

    # 只有概念端点：x=0上升线取归一化基准1，降低线取实际十年恢复率均值。
    # x=0并非观测到的模式第0年；其他锚点由实际浓度域内的平滑组路径线性取样。
    increase_y = np.interp(anchor_x, up.co2_coordinate, up.relative_amoc)
    increase_y[0] = 1.
    minimum_year = int(group_series.loc[140:280, 'relative_amoc'].idxmin())
    minimum_x = float(group_series.loc[minimum_year, 'co2_coordinate'])
    decrease_x = np.sort(np.unique(np.r_[anchor_x, minimum_x]))
    decrease_y = np.interp(decrease_x, down.co2_coordinate, down.relative_amoc)
    decrease_y[0] = final_ratio
    decrease_y[-1] = increase_y[-1]
    increase = PchipInterpolator(anchor_x, increase_y)(drawing_x)
    decrease = PchipInterpolator(decrease_x, decrease_y)(drawing_x)
    assert np.all(np.diff(increase) < 0)
    assert np.all(increase <= 1 + 1e-12)
    assert np.all(decrease > 0)
    assert np.all(decrease <= increase + 1e-12)

    coordinate_parts.append(pd.DataFrame({'panel': panel, 'group': group, 'co2_coordinate': drawing_x,
                                           'increase_amoc_coordinate': increase, 'decrease_amoc_coordinate': decrease}))
    for phase, x_values, y_values in [('increase', anchor_x, increase_y), ('decrease', decrease_x, decrease_y)]:
        for x_value, y_value in zip(x_values, y_values):
            anchor_rows.append({'panel': panel, 'group': group, 'phase': phase,
                                'co2_coordinate': x_value, 'relative_amoc': y_value})
        start, end = arrow_indices[phase]
        y_curve = increase if phase == 'increase' else decrease
        arrow_rows.append({'panel': panel, 'phase': phase, 'start_x': drawing_x[start],
                           'end_x': drawing_x[end], 'start_y': y_curve[start], 'end_y': y_curve[end]})
    straight_chord = 1 + (increase_y[-1] - 1) * drawing_x
    summary_rows.append({'panel': panel, 'group': group, 'n_models': len(members),
                         'endpoint_ratio_mean': final_ratio,
                         'endpoint_ratio_min': float(member_endpoints.endpoint_ratio.min()),
                         'endpoint_ratio_max': float(member_endpoints.endpoint_ratio.max()),
                         'n_above_recovery_threshold': int(np.sum(member_endpoints.endpoint_ratio >= .8)),
                         'increase_end_coordinate': float(increase[-1]),
                         'minimum_actual_group_year': minimum_year, 'minimum_co2_coordinate': minimum_x,
                         'max_increase_above_straight_chord': float(np.max(increase - straight_chord))})


# %% 5. 保存坐标、真实成员范围和口径

destination.mkdir(parents=True)
pd.concat(coordinate_parts, ignore_index=True).to_csv(destination / 'schematic_coordinates.csv', index=False)
pd.DataFrame(anchor_rows).to_csv(destination / 'schematic_anchors.csv', index=False)
pd.DataFrame(arrow_rows).to_csv(destination / 'schematic_arrows.csv', index=False)
pd.concat(member_parts, ignore_index=True).to_csv(destination / 'membership.csv', index=False)
pd.DataFrame(summary_rows).to_csv(destination / 'group_summary.csv', index=False)
metadata = {
    'type': 'Data-informed three-group schematic; simplified guide curves, not the complete actual trajectories or equilibrium branches',
    'groups': group_order, 'titles': [registry['three_group_titles'][group] for group in group_order],
    'members': {group: registry['members'][group] for group in group_order},
    'source_outputs': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
                       for path in [registry_path, source / 'nine_model_normalized_series.csv', source / 'nine_model_endpoint_diagnostics.csv']],
    'dimensions': 'Input model x simulation_year; schematic panel3 x drawing_node601; anchors panel x phase x sparse_node; membership8 rows. No new NetCDF or MAT read.',
    'normalization': 'AMOC for each model divided by its own raw initial1-10 mean, then equal-model mean; no added drift correction.',
    'co2_axis': 'Linear normalized concentration x=(C-284.7)/(284.7*1.01**140-284.7); C0 and Cmax conceptual axis labels, not elapsed time.',
    'shape': 'Replace ellipse quadrant with shape-preserving cubic interpolation through actual-path-informed sparse knots. Blue paths lie predominantly below the initial-to-peak chord: earlier weakening followed by flattening, not an imposed global mathematical second-derivative constraint.',
    'anchor_x': anchor_x.tolist(),
    'anchor_method': 'Use existing11-year-smoothed relative group means; linear sampling at fixed CO2 anchor positions; add actual group minimum position to decrease anchors. PCHIP creates the simplified drawing only, not new simulated data.',
    'endpoint_convention': 'Conceptual increase endpoint at x0 is exactly1. Decrease endpoint at x0 is mean member R from raw years271-280/1-10. These window means are not exact instantaneous AMOC at C0. Peak shared from actual smoothed year140.',
    'threshold': .8, 'third_group': 'Larger-bias selection is not a third recovery class. Its members include both recovery outcomes; show actual member endpoint range and individual endpoint ticks, not a confidence interval.',
    'range_is_not_CI': True,
    'cautions': ['Gray interbranch fill is qualitative path separation, not statistical uncertainty.',
                 'No new class boundary, significance, fitted timescale, loop-area statistic or permanent irreversibility claim.',
                 'Old two-class schematic output remains historical; only its superseded figure is archived after new figure review.',
                 'Other actual model trajectory products and grouping registries are unchanged.'],
    'skill': 'scientific-visualization; actual anchor provenance, explicit schematic label, no invented homogeneous third recovery outcome, common axes and line-style redundancy. Software citation retained in stage documentation.'}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print(pd.DataFrame(summary_rows).to_string(index=False))
