# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
schematic_directory = root / 'output/s01_amoc/n55_amoc_three_group_schematic'
annual_directory = root / 'output/s01_amoc/n22_amoc_recovery_definition_audit'
destination = root / 'output/s01_amoc/n57_amoc_schematic_model_panels'
assert not destination.exists(), '正式结果存在，需先归档。'
reference_co2 = 284.7
maximum_co2 = reference_co2 * 1.01 ** 140
groups = ['better_recovery', 'better_nonrecovery', 'larger_bias']
row_titles = ['Better-simulated · Substantial recovery · n = 2',
              'Better-simulated · Incomplete recovery · n = 2',
              'Larger-bias models · Substantial recovery · n = 4']
fixed_members = [['GFDL-ESM4', 'NorESM2-LM'], ['CESM2', 'CNRM-ESM2-1']]
phase_names = ['increase', 'decrease']
panel_rows, coordinate_parts, arrow_rows, similarity_rows = [], [], [], []


# %% 3. 读取既有示意路径与模式序列

schematic = pd.read_csv(schematic_directory / 'schematic_coordinates.csv')
schematic_arrows = pd.read_csv(schematic_directory / 'schematic_arrows.csv')
summary = pd.read_csv(schematic_directory / 'group_summary.csv')
membership = pd.read_csv(schematic_directory / 'membership.csv')
annual = pd.read_csv(annual_directory / 'nine_model_normalized_series.csv')
assert annual.groupby('model').simulation_year.nunique().eq(280).all()


# %% 4. 按完整增碳和降碳路径筛选第三组两个示例

third_curve = schematic.loc[schematic.group == 'larger_bias'].sort_values('co2_coordinate')
for model in membership.loc[membership.group == 'larger_bias', 'model']:
    history = annual.loc[annual.model == model].sort_values('simulation_year')
    phase_errors = []
    for phase in phase_names:
        # 转折年140在两段共用；两个阶段等权，不根据末期是否恢复挑选。
        branch = history.loc[history.simulation_year.le(140) if phase == 'increase' else history.simulation_year.ge(140)]
        normalized_co2 = (branch.co2_ppm.to_numpy() - reference_co2) / (maximum_co2 - reference_co2)
        schematic_y = np.interp(normalized_co2, third_curve.co2_coordinate,
                                third_curve[f'{phase}_amoc_coordinate'])
        differences = branch.smooth_11yr_percent_initial.to_numpy() / 100 - schematic_y
        phase_errors.append(float(np.mean(differences ** 2)))
    similarity_rows.append({'model': model, 'increase_rmse': np.sqrt(phase_errors[0]),
                            'decrease_rmse': np.sqrt(phase_errors[1]),
                            'equal_phase_rmse': np.sqrt(np.mean(phase_errors))})
similarity = pd.DataFrame(similarity_rows).sort_values(['equal_phase_rmse', 'model']).reset_index(drop=True)
similarity['rank'] = np.arange(1, 5)
similarity['selected'] = similarity['rank'].le(2)
chosen_members = fixed_members + [similarity.loc[similarity.selected, 'model'].tolist()]


# %% 5. 保存九面板绘图坐标及沿时间的箭头

for row, group in enumerate(groups):
    curve = schematic.loc[schematic.group == group].sort_values('co2_coordinate')
    info = summary.loc[summary.group == group].iloc[0]
    for column in range(3):
        panel = row * 3 + column
        model = '' if column == 0 else chosen_members[row][column - 1]
        panel_rows.append({'panel': panel, 'row': row, 'column': column, 'group': group,
                           'row_title': row_titles[row], 'kind': 'schematic' if column == 0 else 'model',
                           'model': model, 'title': 'Schematic' if column == 0 else model,
                           'schematic_final_mean': info.endpoint_ratio_mean,
                           'schematic_final_min': info.endpoint_ratio_min,
                           'schematic_final_max': info.endpoint_ratio_max})
        for phase in phase_names:
            if column == 0:
                part = pd.DataFrame({'co2_ppm': reference_co2 + curve.co2_coordinate.to_numpy() * (maximum_co2 - reference_co2),
                                     'relative_amoc': curve[f'{phase}_amoc_coordinate'].to_numpy(),
                                     'simulation_year': np.nan})
                arrow = schematic_arrows.loc[(schematic_arrows.panel == row) & (schematic_arrows.phase == phase)].iloc[0]
                start_x = reference_co2 + arrow.start_x * (maximum_co2 - reference_co2)
                end_x = reference_co2 + arrow.end_x * (maximum_co2 - reference_co2)
                start_y, end_y = arrow.start_y, arrow.end_y
            else:
                history = annual.loc[annual.model == model].sort_values('simulation_year')
                branch = history.loc[history.simulation_year.le(140) if phase == 'increase' else history.simulation_year.ge(140)]
                part = branch[['co2_ppm', 'simulation_year']].copy()
                part['relative_amoc'] = branch.smooth_11yr_percent_initial.to_numpy() / 100
                # 箭头两端直接取实际序列，降碳段从高浓度指向低浓度。
                arrow_years = [65, 75] if phase == 'increase' else [210, 220]
                arrow_points = history.set_index('simulation_year').loc[arrow_years]
                start_x, end_x = arrow_points.co2_ppm.to_numpy()
                start_y, end_y = arrow_points.smooth_11yr_percent_initial.to_numpy() / 100
            part['panel'], part['phase'] = panel, phase
            part['node'] = np.arange(len(part))
            coordinate_parts.append(part[['panel', 'phase', 'node', 'simulation_year', 'co2_ppm', 'relative_amoc']])
            arrow_rows.append({'panel': panel, 'phase': phase, 'start_x': start_x, 'end_x': end_x,
                               'start_y': start_y, 'end_y': end_y})
coordinates = pd.concat(coordinate_parts, ignore_index=True)
assert np.isfinite(coordinates[['co2_ppm', 'relative_amoc']]).all().all()
input_files = [schematic_directory / name for name in ['schematic_coordinates.csv', 'schematic_arrows.csv',
                'group_summary.csv', 'membership.csv', 'metadata.json']]
input_files.append(annual_directory / 'nine_model_normalized_series.csv')
metadata = {'layout': '3 rows x 3 columns: group schematic, actual model1, actual model2.',
            'groups': groups, 'selected_members': chosen_members,
            'reference_co2_ppm': reference_co2, 'maximum_co2_ppm': maximum_co2,
            'dimensions': 'Saved plotting coordinates: panel x phase x node; model node counts140/141, schematic601/601.',
            'axis': 'Physical atmospheric CO2 concentration in ppm, not simulation year or log concentration.',
            'normalization': 'Reuse n22 smooth_11yr_percent_initial/100; each model normalized by its own raw years1-10 mean; no new smoothing or drift correction.',
            'selection': 'Rank all4 larger-bias members by sqrt((MSE_increase+MSE_decrease)/2) against the existing third schematic at each actual annual CO2. Lowest2 selected. Equal phase and equal annual weights within phase. Descriptive example selection, not validation.',
            'display_naming': 'User-requested third-row name: Larger-bias models / Substantial recovery. Display label only;members,paths,endpoint ratios and original per-model R>=0.8 criterion unchanged.',
            'cautions': ['Left column remains sparse-anchor schematic, not annual model data or equilibrium branches.',
                         'Model panels retain every point in the existing smoothed sequence, including inherited shortened-window endpoints; no forced endpoint or curve reshaping.',
                         'Third schematic represents all4 members, while right panels show only2 examples; neither pair nor schematic spans the entire group diversity.',
                         'Third-row Substantial recovery is the user-requested group label;it does not assert that all4 individual endpoint ratios exceed0.8 or that their trajectories are monotonic.',
                         'CO2 reverses at year140; this point is shared by both phases. Actual paths may cross or have wiggles.',
                         'Schematic final square uses mean member late-decade recovery ratios at conceptual C0; model panels omit endpoint squares to avoid conflating instantaneous smoothed endpoints and late-decade means.'],
            'sources': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in input_files]}
destination.mkdir(parents=True)
pd.DataFrame(panel_rows).to_csv(destination / 'panels.csv', index=False)
coordinates.to_csv(destination / 'coordinates.csv', index=False)
pd.DataFrame(arrow_rows).to_csv(destination / 'arrows.csv', index=False)
similarity.to_csv(destination / 'third_group_similarity.csv', index=False)
metadata['member_endpoint_ratios'] = membership.set_index('model').endpoint_ratio.to_dict()
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print(similarity.to_string(index=False))
print('Selected:', chosen_members)
print('AMOC range:', coordinates.relative_amoc.min(), coordinates.relative_amoc.max())
