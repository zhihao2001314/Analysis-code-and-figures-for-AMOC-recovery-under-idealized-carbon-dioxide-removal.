# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import shutil

import numpy as np
import pandas as pd


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
previous = root / 'output/s01_amoc/n68_four_row_pathway_preview'
latitude_directory = root / 'output/amoc_s01_n01/l_40'
destination = root / 'output/s01_amoc/n70_four_row_pathway_40n'
assert not destination.exists(), '正式结果已存在，须先归档。'
annual_path = latitude_directory / 'CMIP6_AMOC_9models_annual_series.csv'
reference_path = latitude_directory / 'amoc_initial_reference_plot_data.csv'
mapping_path = latitude_directory / 'CMIP6_AMOC_9models_processing_metadata.csv'
unchanged_files = ['co2_forcing.csv', 'path_coordinates.csv', 'path_arrows.csv', 'panels.csv',
                   'group_summary.csv', 'schematic_anchors.csv', 'actual_normalized_group_means.csv']


# %% 3. 以40°N原始强度重建图a

metadata = json.loads((previous / 'metadata.json').read_text(encoding='utf-8'))
models = metadata['models']
annual = pd.read_csv(annual_path)
strength = annual.pivot(index='Simulation year', columns='Model', values='Forced AMOC strength (Sv)')
strength = strength.reindex(columns=models).sort_index()
assert strength.notna().all().all()
initial = strength.loc[1:10].mean()
# 沿用原Fig1：完整已存年度序列先平滑，再显示1–280年。
smoothed = (strength - initial).rolling(11, center=True, min_periods=1).mean()
overview = smoothed.loc[1:280].copy()
overview.index.name = 'simulation_year'
forcing = pd.read_csv(previous / 'co2_forcing.csv').set_index('simulation_year')
overview['co2_ppm'] = forcing.co2_ppm
overview['mmm'] = overview[models].mean(axis=1)
overview['std'] = overview[models].std(axis=1, ddof=0)
overview['lower'] = overview.mmm - overview['std']
overview['upper'] = overview.mmm + overview['std']
saved_reference = pd.read_csv(reference_path).set_index('simulation_year').loc[1:280]
np.testing.assert_allclose(overview[models + ['mmm', 'std', 'lower', 'upper']],
                           saved_reference[models + ['mmm', 'std', 'lower', 'upper']], atol=1e-12, rtol=0)
assert overview[models + ['lower', 'upper']].min().min() >= -21
assert overview[models + ['lower', 'upper']].max().max() <= 3


# %% 4. 保留其余十三个面板的全部坐标

destination.mkdir(parents=True)
for name in unchanged_files:
    shutil.copy2(previous / name, destination / name)
    assert (previous / name).read_bytes() == (destination / name).read_bytes()
overview.to_csv(destination / 'overview_anomaly.csv')
mapping = pd.read_csv(mapping_path).set_index('Model').loc[models].reset_index()
mapping.to_csv(destination / 'panel_a_latitude_mapping.csv', index=False)
summary = pd.DataFrame({'initial_mean_sv': initial, 'late_mean_sv': strength.loc[271:280].mean()})
summary['late_minus_initial_sv'] = summary.late_mean_sv - summary.initial_mean_sv
summary.index.name = 'model'
summary.to_csv(destination / 'panel_a_endpoint_summary.csv')


# %% 5. 保存40°N版本说明和来源

caption = metadata['caption'].replace('existing nine-model26.5°N AMOC-anomaly data', 'nine-model40°N AMOC-anomaly data')
caption += ' This version changes only panel a from26.5°N to40°N; all14-panel layout settings and panels b–n remain unchanged.'
metadata['caption'] = caption
metadata['overview_target_latitude_degN'] = 40.0
metadata['pathway_target_latitude_degN'] = 40.0
metadata['overview_actual_latitudes'] = dict(zip(mapping.Model, mapping['Actual latitude (degN)']))
metadata['revision'] = 'Additional40N overview version for comparison; only panel a science changes, b–n coordinates preserved byte-for-byte.'
metadata['overview_data_range_sv'] = [float(overview[models].min().min()), float(overview[models].max().max())]
metadata['dimensions'] = 'Annual input Model xSimulation year; overview simulation_year280 x9model wide CSV; path dimensions unchanged fromn68.'
sources = [previous / name for name in unchanged_files + ['metadata.json']]
sources.extend([annual_path, reference_path, mapping_path])
metadata['sources'] = [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in sources]
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
(destination / 'caption.txt').write_text(caption + '\n', encoding='utf-8')
print('40N overview source independently verified; other13panels byte-identical.', flush=True)
