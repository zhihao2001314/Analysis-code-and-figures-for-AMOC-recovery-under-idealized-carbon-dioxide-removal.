# %% 1. 导入依赖
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd


# %% 2. 参数与路径设置
root = Path('F:/sea/Antic sea ice constraint')
destination = root / 'output/s02_sia/s02_arctic_region/p25_cas_hydrography/n087_nordic_initial_n2_nine_models'
profile_path = root / 'output/s02sia/s02_arctic_region/p07_stratification/n01_initial_profiles/na_initial_march_stratification_profiles.csv'
iap_directory = root / 'output/s02sia/s02_arctic_region/p07_stratification/n02_iapv4_observation'
observation_path = iap_directory / 'iapv4_march_2005_2024_stratification_profiles.csv'
cas_path = root / 'output/s02_sia/s02_arctic_region/p25_cas_hydrography/n083_cas_mld_stratification/cas_march_n2_profiles.csv'
reference_metadata_path = root / 'output/s02_sia/s02_arctic_region/p25_cas_hydrography/n085_nine_model_hydrography_panels/metadata.json'
sources = [profile_path, observation_path, cas_path, reference_metadata_path,
           iap_directory / 'iapv4_march_2005_2024_stratification_metadata.json']
assert not destination.exists(), '正式结果已存在，覆盖前须先归档。'
destination.mkdir(parents=True)


# %% 3. 提取既有八模式和CAS矩形分区初始剖面
reference_metadata = json.loads(reference_metadata_path.read_text(encoding='utf-8'))
models = reference_metadata['models']
original = pd.read_csv(profile_path)
profiles = original.loc[(original.region == 'Nordic Seas') & original.model.isin(models[:-1]),
                         ['model', 'depth_m', 'initial_march_n2_s2']].copy()
profiles = profiles.rename(columns={'initial_march_n2_s2': 'n2_s2'})
cas_original = pd.read_csv(cas_path)
cas = cas_original.loc[(cas_original.version == 'rectangles') &
                       (cas_original.region == 'Nordic Seas') &
                       (cas_original.period == 'up_initial'), ['model', 'depth_m', 'n2_s2']].copy()
profiles = pd.concat([profiles, cas], ignore_index=True)
profiles['n2_1e5_s2'] = profiles.n2_s2 * 1e5
original_observation = pd.read_csv(observation_path)
observation = original_observation.loc[original_observation.region == 'Nordic Seas'].copy()


# %% 4. 保存并核对绘图输入
assert len(profiles) == 900 and len(observation) == 100
assert set(profiles.model) == set(models)
assert not profiles.duplicated(['model', 'depth_m']).any()
for model in models:
    selected = profiles.loc[profiles.model == model].sort_values('depth_m')
    assert len(selected) == 100
    np.testing.assert_array_equal(selected.depth_m, np.arange(5, 1000, 10))
    assert np.nanmax(selected.n2_1e5_s2) < reference_metadata['n2_x_limits'][1]
    assert np.nanmin(selected.n2_1e5_s2) > reference_metadata['n2_x_limits'][0]
profiles.to_csv(destination / 'initial_n2_profiles.csv', index=False)
observation.to_csv(destination / 'iap_initial_n2_reference.csv', index=False)
saved = pd.read_csv(destination / 'initial_n2_profiles.csv')
np.testing.assert_allclose(saved.n2_s2, profiles.n2_s2, rtol=1e-14, atol=1e-18)
metadata = {
    'models': models, 'region': 'Nordic Seas',
    'region_mask': '65<=latitude<80 degrees N; -30<=longitude<20 degrees E; ocean cells only',
    'region_shape': 'rectangle, not the coastal polygon used for the MLD time series',
    'model_period': 'initial 10 March samples, inherited GFDL model years6-15 exception; CAS years1-10, calendar200-209',
    'observation_period': 'IAP March2005-2024 climatology, identical profile in every panel; not contemporaneous evaluation',
    'dimensions': {'model_depth_long_table': [9, 100], 'iap_depth_long_table': [100]},
    'n2_x_limits': reference_metadata['n2_x_limits'], 'depth_limits_m': [0, 1000],
    'processing': 'saved initial profiles only; no recomputation, vertical smoothing or new piControl correction',
    'source_files': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in sources],
    'qc': {'nine_models_and_iap': True, 'csv_roundtrip_passed': True, 'common_x_limits_do_not_clip_data': True}
}
(destination / 'metadata.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
print(destination, flush=True)
