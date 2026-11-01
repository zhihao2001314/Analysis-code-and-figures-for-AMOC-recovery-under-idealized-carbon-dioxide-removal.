# %% 1. 导入依赖
from pathlib import Path
import hashlib
import json
import shutil

import numpy as np
import pandas as pd


# %% 2. 参数与路径设置
root = Path('F:/sea/Antic sea ice constraint')
prior_directory = root / 'output/s01_amoc/n62_amoc_pathway_composite'
latitude_directory = root / 'output/amoc_s01_n01/l_26.5'
destination = root / 'output/s01_amoc/n66_amoc_overview_26p5'
annual_path = latitude_directory / 'CMIP6_AMOC_9models_annual_series.csv'
saved_overview_path = latitude_directory / 'amoc_initial_reference_plot_data.csv'
latitude_mapping_path = latitude_directory / 'CMIP6_AMOC_9models_processing_metadata.csv'
reference_metadata_path = latitude_directory / 'amoc_initial_reference_metadata.json'
unchanged_files = ['co2_forcing.csv', 'path_coordinates.csv', 'path_arrows.csv',
                   'panels.csv', 'schematic_anchors.csv']
assert not destination.exists(), '正式结果已存在，覆盖前须先归档。'
destination.mkdir(parents=True)


# %% 3. 从既有26.5N原始年度强度重建图a
metadata = json.loads((prior_directory / 'metadata.json').read_text(encoding='utf-8'))
models = metadata['models']
annual = pd.read_csv(annual_path)
strength = annual.pivot(index='Simulation year', columns='Model', values='Forced AMOC strength (Sv)')
strength = strength.reindex(columns=models).sort_index()
assert strength.notna().all().all()
initial = strength.loc[1:10].mean()
anomaly = strength - initial
# 与原图相同：在完整年度来源上平滑，再显示1–280年；边缘缩短窗口。
smoothed = anomaly.rolling(11, center=True, min_periods=1).mean()
overview = smoothed.loc[1:280].copy()
overview.index.name = 'simulation_year'
forcing = pd.read_csv(prior_directory / 'co2_forcing.csv').set_index('simulation_year')
overview['co2_ppm'] = forcing.co2_ppm
overview['mmm'] = overview[models].mean(axis=1)
overview['std'] = overview[models].std(axis=1, ddof=0)
overview['lower'] = overview.mmm - overview['std']
overview['upper'] = overview.mmm + overview['std']
assert list(overview.index) == list(range(1, 281))
saved_overview = pd.read_csv(saved_overview_path).set_index('simulation_year').loc[1:280]
np.testing.assert_allclose(overview[models + ['mmm', 'std', 'lower', 'upper']],
                           saved_overview[models + ['mmm', 'std', 'lower', 'upper']], atol=1e-12, rtol=0)


# %% 4. 保留原图b至h的所有数据与CO2强迫
for name in unchanged_files:
    shutil.copy2(prior_directory / name, destination / name)
    assert (prior_directory / name).read_bytes() == (destination / name).read_bytes()
overview.to_csv(destination / 'overview_anomaly.csv')
latitude_mapping = pd.read_csv(latitude_mapping_path)
latitude_mapping = latitude_mapping.set_index('Model').loc[models].reset_index()
latitude_mapping.to_csv(destination / 'panel_a_latitude_mapping.csv', index=False)
endpoints = pd.DataFrame({'initial_mean_sv': initial,
                          'late_mean_sv': strength.loc[271:280].mean(),
                          'late_minus_initial_sv': anomaly.loc[271:280].mean()})
endpoints.index.name = 'model'
endpoints.to_csv(destination / 'panel_a_endpoint_summary.csv')


# %% 5. 保存纬度、坐标范围与图注
data_min = float(overview[models + ['lower']].min().min())
data_max = float(overview[models + ['upper']].max().max())
y_limits = [min(-21., float(np.floor(data_min / 5) * 5)),
            max(3., float(np.ceil(data_max / 5) * 5))]
# 数据未超过旧坐标范围时不改变纵轴。
if data_max <= 3:
    y_limits[1] = 3.
if data_min >= -21:
    y_limits[0] = -21.
caption = metadata['caption'].replace('at the grid latitude nearest 40°N,',
                                       'at the grid latitude nearest 26.5°N,', 1)
caption = caption.replace('in (d, e) and (g, h).',
                           'in (d, e) and (g, h), which retain the original 40°N AMOC definition.')
reference_metadata = json.loads(reference_metadata_path.read_text(encoding='utf-8'))
metadata['overview_reference'] = reference_metadata
metadata['overview_target_latitude_degrees_north'] = 26.5
metadata['model_pathway_target_latitude_degrees_north'] = 40.
metadata['overview_y_limits_sv'] = y_limits
metadata['overview_y_ticks_sv'] = list(np.arange(-20, y_limits[1] + 1, 5))
metadata['overview_data_range_sv'] = [data_min, data_max]
metadata['overview_actual_latitudes'] = dict(zip(latitude_mapping.Model,
                                                latitude_mapping['Actual latitude (degN)']))
metadata['revision'] = 'Only panel a changes from40N to26.5N; b-h and prescribedCO2 retained byte-for-byte. Recovery groups still use the original40N definition.'
metadata['caption'] = caption
metadata['dimensions'] = 'Saved CSV indexed by simulation_year280, nine model columns plusCO2/MMM/populationSD/bounds; original annual input long table indexed by Model and Simulation year. No raw NetCDF/MAT read.'
sources = [prior_directory / name for name in unchanged_files + ['overview_anomaly.csv', 'metadata.json']]
sources.extend([annual_path, saved_overview_path, latitude_mapping_path, reference_metadata_path])
metadata['sources'] = [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in sources]
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
(destination / 'caption.txt').write_text(caption + '\n', encoding='utf-8')
print(json.dumps({'output': str(destination), 'panel_a': '26.5N', 'other_panels': 'unchanged',
                 'y_limits': y_limits, 'data_range': [data_min, data_max]}, indent=2), flush=True)
