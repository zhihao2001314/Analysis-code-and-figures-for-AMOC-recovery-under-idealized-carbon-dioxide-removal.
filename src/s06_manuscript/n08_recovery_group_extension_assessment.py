# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import sys

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding='utf-8')


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
destination = root / 'output/s06_manuscript/n08_recovery_group_extension_assessment'
assert not destination.exists(), '正式结果已存在，须先归档。'
amoc_paths = {40.0: root / 'output/amoc_s01_n01/l_40/CMIP6_AMOC_9models_annual_series.csv',
              26.5: root / 'output/amoc_s01_n01/l_26.5/CMIP6_AMOC_9models_annual_series.csv'}
profile_base = root / 'output/s02_sia/s02_arctic_region/p25_cas_hydrography'
regions = {'LI': 'n085_nine_model_hydrography_panels', 'NS': 'n087_nordic_initial_n2_nine_models'}
threshold = .8
models = ['CESM2', 'CNRM-ESM2-1', 'NorESM2-LM', 'GFDL-ESM4',
          'ACCESS-ESM1-5', 'CanESM5', 'MIROC-ES2L', 'UKESM1-0-LL', 'CAS-ESM2-0']
source_paths = list(amoc_paths.values())


# %% 3. 原始年值的恢复比例与纬度及末期窗口敏感性

endpoint_rows = []
for latitude, path in amoc_paths.items():
    annual = pd.read_csv(path)
    for model in models:
        selected = annual.loc[annual.Model.eq(model) & annual['Simulation year'].between(1, 280)]
        series = selected.set_index('Simulation year')['Forced AMOC strength (Sv)'].sort_index()
        np.testing.assert_array_equal(series.index, np.arange(1, 281))
        assert np.isfinite(series.to_numpy()).all()
        initial = float(series.loc[1:10].mean())
        for window in [10, 20, 30]:
            final = float(series.loc[281 - window:280].mean())
            ratio = final / initial
            endpoint_rows.append({'model': model, 'latitude_degN': latitude,
                                  'initial_years': '1–10', 'final_window_years': window,
                                  'initial_amoc_sv': initial, 'final_amoc_sv': final,
                                  'residual_amoc_sv': final - initial, 'recovery_ratio': ratio,
                                  'recovery_percent': 100 * ratio,
                                  'category': 'Substantial recovery' if ratio >= threshold else 'Non-recovery'})
endpoints = pd.DataFrame(endpoint_rows)
for model, category in [('CAS-ESM2-0', 'Substantial recovery'), ('UKESM1-0-LL', 'Non-recovery')]:
    assert endpoints.loc[endpoints.model.eq(model), 'category'].eq(category).all()


# %% 4. 相同深度层的初始层结与IAP描述性对比

profile_rows = []
for region, folder in regions.items():
    profile_path = profile_base / folder / 'initial_n2_profiles.csv'
    observation_path = profile_base / folder / 'iap_initial_n2_reference.csv'
    source_paths.extend([profile_path, observation_path])
    profiles = pd.read_csv(profile_path)
    observation = pd.read_csv(observation_path)
    paired = profiles.merge(observation[['depth_m', 'iapv4_march_mean_n2_1e5_s2']], on='depth_m', validate='many_to_one')
    # 所有模式统一15–995m、每10m一层；不以CAS缺失的最浅5m填值。
    paired = paired.loc[paired.depth_m.between(15, 995)]
    for model in models:
        values = paired.loc[paired.model.eq(model)].sort_values('depth_m')
        assert len(values) == 99 and values[['n2_1e5_s2', 'iapv4_march_mean_n2_1e5_s2']].notna().all().all()
        error = values.n2_1e5_s2.to_numpy() - values.iapv4_march_mean_n2_1e5_s2.to_numpy()
        profile_rows.append({'region': region, 'model': model, 'depth_min_m': 15, 'depth_max_m': 995,
                             'sampled_depth_count': 99, 'rmsd_1e5_s2': float(np.sqrt(np.mean(error ** 2))),
                             'mean_bias_1e5_s2': float(np.mean(error))})
profile_diagnostics = pd.DataFrame(profile_rows)
profile_diagnostics['rmsd_rank_of9'] = profile_diagnostics.groupby('region').rmsd_1e5_s2.rank(method='min').astype(int)


# %% 5. 保存评估结果，不修改现有分组或复合图

candidate_groups = pd.DataFrame([
    ['Substantial recovery', 'GFDL-ESM4', 'existing'],
    ['Substantial recovery', 'NorESM2-LM', 'existing'],
    ['Substantial recovery', 'CAS-ESM2-0', 'proposed addition'],
    ['Non-recovery', 'CESM2', 'existing'],
    ['Non-recovery', 'CNRM-ESM2-1', 'existing'],
    ['Non-recovery', 'UKESM1-0-LL', 'proposed addition'],
], columns=['group', 'model', 'status'])
destination.mkdir(parents=True)
endpoints.to_csv(destination / 'amoc_recovery_and_window_sensitivity.csv', index=False)
profile_diagnostics.to_csv(destination / 'initial_n2_iap_comparison.csv', index=False)
candidate_groups.to_csv(destination / 'proposed_six_model_groups.csv', index=False)
metadata = {
    'purpose': 'Evaluate CAS and UK addition, not bulk reassign existing figures or membership.',
    'primary_definition': 'Existing40N, depth>=500m Atlantic–Arctic index; final271–280 mean / initial1–10 mean, threshold0.8.',
    'data': 'Forced AMOC strength only, raw saved annual values, no smoothing or new piControl drift correction.',
    'sensitivity': 'Also26.5N and final10/20/30year windows with fixed initial10year reference.',
    'profiles': 'Initial March regional N2 vs IAP2005–2024 March mean.15–995m on99 uniform10m samples, equal-depth RMSD.',
    'scope_limit': 'Modern observations vs idealized initial climate; RMSD is descriptive, not a formal validation or new calibrated skill threshold. Recovery outcome and initial regional bias are separate criteria.',
    'proposed_names': 'Substantial recovery / Non-recovery; expanded group must not be universally labelled better-simulated.',
    'dimensions': {'annual_csv': 'Model x Simulation year, named columns', 'n2_csv': 'model x depth_m, named columns'},
    'sources': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in source_paths],
}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print(endpoints.loc[endpoints.model.isin(['CAS-ESM2-0', 'UKESM1-0-LL'])].to_string(index=False), flush=True)
print(profile_diagnostics.loc[profile_diagnostics.model.isin(['CAS-ESM2-0', 'UKESM1-0-LL'])].to_string(index=False), flush=True)
