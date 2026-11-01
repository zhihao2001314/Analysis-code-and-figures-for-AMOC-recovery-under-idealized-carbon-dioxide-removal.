# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
source = root / 'output/s01_amoc/n73_hysteresis_morphology_groups'
latitude_source = root / 'output/amoc_s01_n01/l_26.5'
destination = root / 'output/s01_amoc/n75_hysteresis_morphology_26p5'
assert not destination.exists() or not any(destination.iterdir()), '正式输出存在，须先归档。'
annual_path = latitude_source / 'CMIP6_AMOC_9models_annual_series.csv'
mapping_path = latitude_source / 'CMIP6_AMOC_9models_processing_metadata.csv'
years = np.arange(1, 281)
initial_years, final_years = np.arange(1, 11), np.arange(271, 281)
peak_year, smoothing_window = 140, 11
schematic_x = np.linspace(0., 1., 601)
anchor_x = np.array([0., .08, .16, .3, .5, .7, .84, 1.])
small_gap_amplitude = .016
coordinate_parts, arrow_records, normalized_parts = [], [], []
metric_records, summary_records, mean_parts, anchor_records = [], [], [], []


# %% 3. 读取26.5度强迫原值和概念图形状

metadata = json.loads((source / 'metadata.json').read_text(encoding='utf-8'))
models, members = metadata['models'], metadata['members']
annual = pd.read_csv(annual_path)
mapping = pd.read_csv(mapping_path).set_index('Model').loc[models]
strength = annual.pivot(index='Simulation year', columns='Model', values='Forced AMOC strength (Sv)').loc[years, models]
forcing = annual.pivot(index='Simulation year', columns='Model', values='Atmospheric CO2 (ppm)').loc[years, models]
assert strength.notna().all().all()
np.testing.assert_allclose(forcing.to_numpy(), np.repeat(forcing.iloc[:, :1].to_numpy(), len(models), axis=1), atol=1e-7, rtol=0)
co2 = forcing.iloc[:, 0]
panels = pd.read_csv(source / 'panels.csv', keep_default_na=False)
concept_paths = pd.read_csv(source / 'path_coordinates.csv')
concept_arrows = pd.read_csv(source / 'path_arrows.csv')
coordinate_parts.append(concept_paths.loc[concept_paths.panel.eq('b')].copy())
arrow_records.extend(concept_arrows.loc[concept_arrows.panel.eq('b')].to_dict('records'))


# %% 4. 同纬度统一计算异常与初始强度归一化

initial = strength.loc[initial_years].mean()
late = strength.loc[final_years].mean()
# 仅平滑被分析的1–280年，不借用后续稳定阶段；端点采用缩短窗口。
smoothed = strength.rolling(smoothing_window, center=True, min_periods=1).mean()
relative = smoothed / initial
overview = smoothed - initial
overview.index.name = 'simulation_year'
overview['co2_ppm'] = co2
overview['mmm'] = overview[models].mean(axis=1)
overview['std'] = overview[models].std(axis=1, ddof=0)
overview['lower'] = overview.mmm - overview['std']
overview['upper'] = overview.mmm + overview['std']
endpoint_rows = []
for model in models:
    normalized_parts.append(pd.DataFrame({'model': model, 'simulation_year': years,
                                          'co2_ppm': co2.to_numpy(), 'raw_sv': strength[model].to_numpy(),
                                          'smooth_11yr_sv': smoothed[model].to_numpy(),
                                          'relative_amoc': relative[model].to_numpy(),
                                          'initial_mean_sv': initial[model]}))
    endpoint_rows.append({'model': model, 'initial_mean_sv': initial[model], 'late_mean_sv': late[model],
                          'endpoint_ratio': late[model] / initial[model],
                          'late_minus_initial_sv': late[model] - initial[model],
                          'minimum_11yr_sv': smoothed[model].min(),
                          'minimum_11yr_year': int(smoothed[model].idxmin()),
                          'actual_latitude_degN': mapping.loc[model, 'Actual latitude (degN)'],
                          'member': mapping.loc[model, 'Selected member']})
endpoints = pd.DataFrame(endpoint_rows).set_index('model')


# %% 5. 重建九条模式路径及滞回面积

for record in panels.loc[panels.kind.eq('model')].itertuples(index=False):
    for phase in ['increase', 'decrease']:
        phase_years = years[years <= peak_year] if phase == 'increase' else years[years >= peak_year]
        coordinate_parts.append(pd.DataFrame({'panel': record.panel, 'phase': phase,
                                              'node': np.arange(len(phase_years)),
                                              'simulation_year': phase_years,
                                              'co2_ppm': co2.loc[phase_years].to_numpy(),
                                              'relative_amoc': relative.loc[phase_years, record.model].to_numpy()}))
        arrow_years = [65, 75] if phase == 'increase' else [210, 220]
        arrow_records.append({'panel': record.panel, 'phase': phase,
                              'start_x': float(co2.loc[arrow_years[0]]), 'end_x': float(co2.loc[arrow_years[1]]),
                              'start_y': float(relative.loc[arrow_years[0], record.model]),
                              'end_y': float(relative.loc[arrow_years[1], record.model])})
    up = pd.DataFrame({'co2': co2.loc[1:140], 'relative': relative.loc[1:140, record.model]}).sort_values('co2')
    down = pd.DataFrame({'co2': co2.loc[140:280], 'relative': relative.loc[140:280, record.model]}).sort_values('co2')
    matched_x = np.unique(np.r_[up.co2.to_numpy(), down.co2.to_numpy()])
    matched_x = matched_x[(matched_x >= max(up.co2.min(), down.co2.min()))
                          & (matched_x <= min(up.co2.max(), down.co2.max()))]
    differences = np.interp(matched_x, up.co2, up.relative) - np.interp(matched_x, down.co2, down.relative)
    crossings = np.flatnonzero(differences[:-1] * differences[1:] < 0.)
    crossing_x = matched_x[crossings] - differences[crossings] * np.diff(matched_x)[crossings] / np.diff(differences)[crossings]
    matched_x = np.unique(np.r_[matched_x, crossing_x])
    separation = np.abs(np.interp(matched_x, up.co2, up.relative) - np.interp(matched_x, down.co2, down.relative))
    area = float(np.trapezoid(separation, matched_x))
    metric_records.append({'model': record.model, 'group': metadata['group_names'][record.row - 1],
                           'panel': record.panel, 'target_latitude_degN': 26.5,
                           'area_relative_amoc_ppm': area, 'area_sv_ppm': area * initial[record.model],
                           'area_divided_by_co2_span': area / (matched_x[-1] - matched_x[0]),
                           'co2_min_ppm': matched_x[0], 'co2_max_ppm': matched_x[-1],
                           'endpoint_ratio_raw_decades': late[record.model] / initial[record.model],
                           'increase_low_co2_smooth': float(up.relative.iloc[0]),
                           'decrease_low_co2_smooth': float(down.relative.iloc[0])})


# %% 6. 更新26.5度参考的三幅示意图

for record in panels.loc[panels.kind.eq('schematic')].itertuples(index=False):
    group_models = members[record.row - 1]
    mean = relative[group_models].mean(axis=1)
    peak = float(mean.loc[peak_year])
    endpoint_ratio = float(endpoints.loc[group_models, 'endpoint_ratio'].mean())
    mean_parts.append(pd.DataFrame({'simulation_year': years, 'co2_ppm': co2.to_numpy(),
                                    'relative_amoc': mean.to_numpy(), 'group': record.label}))
    summary_records.append({'row': record.row, 'group': record.label, 'model_count': len(group_models),
                            'models': '; '.join(group_models), 'mean_endpoint_ratio': endpoint_ratio,
                            'peak_relative_group_mean': peak})
    for phase in ['increase', 'decrease']:
        if record.row == 1:
            values = 1. - (1. - peak) * schematic_x
            if phase == 'decrease':
                values = values - small_gap_amplitude * np.sin(np.pi * schematic_x)
        else:
            shape = concept_paths.loc[concept_paths.panel.eq(record.panel) & concept_paths.phase.eq(phase)].sort_values('node')
            old_values = shape.relative_amoc.to_numpy()
            old_peak, old_start = float(old_values[-1]), float(old_values[0])
            start = 1. if phase == 'increase' else endpoint_ratio
            values = peak + (old_values - old_peak) * (start - peak) / (old_start - old_peak)
        coordinate_parts.append(pd.DataFrame({'panel': record.panel, 'phase': phase, 'node': np.arange(601),
                                              'simulation_year': np.nan,
                                              'co2_ppm': metadata['reference_co2_ppm'] + schematic_x * (metadata['maximum_co2_ppm'] - metadata['reference_co2_ppm']),
                                              'relative_amoc': values}))
        for x, y in zip(anchor_x, np.interp(anchor_x, schematic_x, values)):
            anchor_records.append({'panel': record.panel, 'phase': phase, 'normalized_co2': x, 'relative_amoc': y})
        start_node, end_node = ([200, 250] if phase == 'increase' else [410, 360]) if record.row == 1 else ([270, 320] if phase == 'increase' else [330, 280])
        arrow_records.append({'panel': record.panel, 'phase': phase,
                              'start_x': float(metadata['reference_co2_ppm'] + schematic_x[start_node] * (metadata['maximum_co2_ppm'] - metadata['reference_co2_ppm'])),
                              'end_x': float(metadata['reference_co2_ppm'] + schematic_x[end_node] * (metadata['maximum_co2_ppm'] - metadata['reference_co2_ppm'])),
                              'start_y': float(values[start_node]), 'end_y': float(values[end_node])})


# %% 7. 核对同纬度一致性并保存

coordinates = pd.concat(coordinate_parts, ignore_index=True)
metrics = pd.DataFrame(metric_records)
for record in panels.loc[panels.kind.eq('model')].itertuples(index=False):
    part = coordinates.loc[coordinates.panel.eq(record.panel)]
    lookup = part.simulation_year.astype(int)
    np.testing.assert_allclose(part.relative_amoc, relative.loc[lookup, record.model], atol=1e-12, rtol=0)
    np.testing.assert_allclose(part.relative_amoc * initial[record.model] - initial[record.model],
                              overview.loc[lookup, record.model], atol=1e-12, rtol=0)
assert len(panels) + 1 == 14 and panels.kind.eq('model').sum() == 9
assert coordinates.relative_amoc.between(0, 1.2).all(), '路径超出纵轴范围，需要统一调整范围。'
assert np.isfinite(coordinates[['co2_ppm', 'relative_amoc']]).all().all()
destination.mkdir(parents=True, exist_ok=True)
overview.to_csv(destination / 'overview_anomaly.csv')
pd.DataFrame({'simulation_year': years, 'co2_ppm': co2.to_numpy()}).to_csv(destination / 'co2_forcing.csv', index=False)
coordinates.to_csv(destination / 'path_coordinates.csv', index=False)
pd.DataFrame(arrow_records).to_csv(destination / 'path_arrows.csv', index=False)
panels.to_csv(destination / 'panels.csv', index=False)
pd.DataFrame(anchor_records).to_csv(destination / 'schematic_anchors.csv', index=False)
metrics.to_csv(destination / 'model_hysteresis_metrics.csv', index=False)
endpoints.to_csv(destination / 'endpoint_diagnostics.csv')
pd.DataFrame(summary_records).to_csv(destination / 'group_summary.csv', index=False)
pd.concat(normalized_parts, ignore_index=True).to_csv(destination / 'nine_model_normalized_series.csv', index=False)
pd.concat(mean_parts, ignore_index=True).to_csv(destination / 'actual_normalized_group_means.csv', index=False)
mapping.reset_index().to_csv(destination / 'latitude_mapping.csv', index=False)


# %% 8. 保存修正口径和来源

caption = ('AMOC evolution and transient hysteresis pathways at the native grid latitude nearest26.5N, '
           'using the maximum Atlantic-Arctic overturning below500m. Panel a shows initial-decade anomalies '
           'and the equal-weight nine-model mean with intermodel population SD; b is a conceptual area diagram. '
           'All nine actual trajectories in d,e,g-i,k-n are newly calculated from26.5N forced annual strength, '
           'normalized by each model\'s raw years1-10 mean. Both the overview and paths use the same11-year '
           'centred means within years1-280, shortening the end windows and excluding the stabilization period. '
           'First-row differences from the previous version are limited to the final five smoothed years. '
           'Rows retain the requested descriptive2/3/4model arrangement and labels; schematics c,f,j are idealized '
           'illustrations with26.5N peak references, not fits or group-mean trajectories. The first schematic '
           'is closed and has small area; actual ACCESS/CAS endpoints are not forced to close. '
           'Hysteresis area integrates absolute interbranch separation over common CO2 concentrations, '
           'including linear crossing points. Blue/orange solid paths indicate increasing/decreasing CO2 '
           'and arrows indicate time direction. This corrected version replaces the earlier mixed-latitude preview.')
metadata.update({'caption': caption, 'overview_target_latitude_degN': 26.5,
                 'pathway_target_latitude_degN': 26.5,
                 'overview_actual_latitudes': mapping['Actual latitude (degN)'].to_dict(),
                 'pathway_actual_latitudes': mapping['Actual latitude (degN)'].to_dict(),
                 'normalization': 'Each model26.5N forced strength divided by its own raw years1-10 mean; no piControl subtraction or drift correction.',
                 'smoothing': '11yr centred min_periods1, strictly analysis years1-280 for both a and all nine paths.',
                 'revision': 'Corrected all nine actual paths, arrows, endpoints, areas and schematic reference values to26.5N. No40N model trajectory is copied.',
                 'concept_source_scope': 'n73 contributes only layout, user group labels and abstract concept geometry, never actual model path values.',
                 'dimensions': 'Input CSV Model9 xSimulation year301; selected1-280. Plot overview time280 x9model; normalized long table2520rows; model branches140/141nodes. No raw NetCDF/MAT read.'})
sources = [annual_path, mapping_path]
sources.extend([source / name for name in ['path_coordinates.csv', 'path_arrows.csv', 'panels.csv', 'metadata.json']])
metadata['sources'] = [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in sources]
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
(destination / 'caption.txt').write_text(caption + '\n', encoding='utf-8')
print(metrics[['model', 'endpoint_ratio_raw_decades', 'area_divided_by_co2_span']].to_string(index=False), flush=True)
print(destination.resolve(), flush=True)
