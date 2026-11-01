# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import shutil

import numpy as np
import pandas as pd


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
previous = root / 'output/s06_manuscript/n04_tables_peak_and_variables'
destination = root / 'output/s06_manuscript/n06_gfdl_mld_peak_time_audit'
qa_path = root / 'output/s02sia/s02_arctic_region/p06_mld/n04_na_region_timeseries/na_march_mld_region_qa.csv'
history_path = root / 'output/s02_sia/s02_arctic_region/p21_grouped_mld/n061_grouped_mld_timeseries/member_timeseries.csv'
assert not destination.exists(), '正式结果已存在，须先归档。'


# %% 3. 核对源年份与旧索引

qa = pd.read_csv(qa_path)
up = qa.loc[qa.model.eq('GFDL-ESM4') & qa.experiment.eq('1pctCO2')].iloc[0]
down = qa.loc[qa.model.eq('GFDL-ESM4') & qa.experiment.eq('1pctCO2-cdr')].iloc[0]
assert up.thetao_first_march_year == up.so_first_march_year == 6
assert up.thetao_last_march_year == up.so_last_march_year == 145
assert down.thetao_first_march_year == down.so_first_march_year == 141
assert down.thetao_last_march_year == down.so_last_march_year == 280
peak_legacy_index = 140 - int(up.thetao_first_march_year) + 1
assert peak_legacy_index == 135
history = pd.read_csv(history_path)
table = pd.read_csv(previous / 'table1_peak_state_residual.csv')
original = table.copy(deep=True)
mapping_rows = []
for region, key in [('Labrador-Irminger', 'li'), ('Nordic Seas', 'ns')]:
    series = history.loc[history.model.eq('GFDL-ESM4') & history.region.eq(region)].set_index('simulation_year').march_mld_m
    assert series.index.tolist() == list(range(1, 281))
    value = float(series.loc[peak_legacy_index])
    assert np.isfinite(value) and value > 0
    table.loc[table.model.eq('GFDL-ESM4'), 'peak_mld_' + key + '_m'] = value
    mapping_rows.append({'region': region, 'physical_year': 140, 'legacy_index': peak_legacy_index,
                         'march_mld_m': value, 'old_index140_physical_year': 145,
                         'old_index140_mld_m': float(series.loc[140])})
pd.testing.assert_frame_equal(table.loc[table.model.ne('GFDL-ESM4')], original.loc[original.model.ne('GFDL-ESM4')])
assert table.loc[table.model.eq('GFDL-ESM4'), ['delta_mld_li_m', 'delta_mld_ns_m']].isna().all().all()


# %% 4. 保存修正版独立结果及来源

destination.mkdir(parents=True)
table.to_csv(destination / 'table1_peak_state_residual.csv', index=False)
for name in ['table2_variables_periods.csv', 'table2_source_calendar_windows.csv', 'cas_standard_tas_selected_years.csv']:
    shutil.copy2(previous / name, destination / name)
pd.DataFrame(mapping_rows).to_csv(destination / 'gfdl_peak_mapping.csv', index=False)
metadata = json.loads((previous / 'metadata.json').read_text(encoding='utf-8'))
metadata['GFDL_MLD'] = 'Verified legacy up source6–145 labelled1–140; physical year140=legacy index135. Peak actual MLD now reported. Initial source1–5 absent, so residual against1–10 still not reported. Original series and open Word files unchanged.'
metadata['audit_sources'] = [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
                             for path in [qa_path, history_path, previous / 'table1_peak_state_residual.csv']]
metadata['word_paths'] = [str(root / 'docs/table1_cmip6_peak_and_residual_gfdl140.docx'),
                          str(root / 'docs/table2_variables_and_periods_gfdl140.docx')]
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print(pd.DataFrame(mapping_rows).to_string(index=False), flush=True)
