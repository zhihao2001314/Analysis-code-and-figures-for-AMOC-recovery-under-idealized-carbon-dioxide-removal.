# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import shutil

import numpy as np
import pandas as pd


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
previous = root / 'output/s06_manuscript/n06_gfdl_mld_peak_time_audit'
destination = root / 'output/s06_manuscript/n07_gfdl_start_year6'
history_path = root / 'output/s02_sia/s02_arctic_region/p21_grouped_mld/n061_grouped_mld_timeseries/member_timeseries.csv'
assert not destination.exists(), '正式结果已存在，须先归档。'


# %% 3. 按作者指定的第6年起点计算

history = pd.read_csv(history_path)
table = pd.read_csv(previous / 'table1_peak_state_residual.csv')
original = table.copy(deep=True)
summary = []
for region, key in [('Labrador-Irminger', 'li'), ('Nordic Seas', 'ns')]:
    values = history.loc[history.model.eq('GFDL-ESM4') & history.region.eq(region)].set_index('simulation_year').march_mld_m
    assert values.index.tolist() == list(range(1, 281))
    # 作者确认：升碳源6–145年对应分析1–140年，不再减去5年。
    initial = float(values.loc[1:10].mean())
    peak = float(values.loc[140])
    late = float(values.loc[271:280].mean())
    delta = late - initial
    for name, value in [('initial', initial), ('peak', peak), ('late', late), ('delta', delta)]:
        table.loc[table.model.eq('GFDL-ESM4'), name + '_mld_' + key + '_m'] = value
    summary.append({'region': region, 'initial_source_years': '6–15', 'peak_source_year': 145,
                    'peak_analysis_year': 140, 'late_source_years': '271–280',
                    'initial_mld_m': initial, 'peak_mld_m': peak, 'late_mld_m': late, 'delta_mld_m': delta})
pd.testing.assert_frame_equal(table.loc[table.model.ne('GFDL-ESM4')], original.loc[original.model.ne('GFDL-ESM4')])
assert np.isfinite(table[['peak_mld_li_m', 'peak_mld_ns_m', 'delta_mld_li_m', 'delta_mld_ns_m']]).all().all()
calendar = pd.read_csv(previous / 'table2_source_calendar_windows.csv')
calendar = pd.concat([calendar, pd.DataFrame([{
    'model': 'GFDL-ESM4', 'variable_scope': 'thetao, so; source year 6 = analysis year 1',
    'initial_years': '0006–0015', 'peak_year': '0145', 'final_years': '0271–0280',
}])], ignore_index=True)


# %% 4. 保存最终数据及作者时间口径

destination.mkdir(parents=True)
table.to_csv(destination / 'table1_peak_state_residual.csv', index=False)
calendar.to_csv(destination / 'table2_source_calendar_windows.csv', index=False)
for name in ['table2_variables_periods.csv', 'cas_standard_tas_selected_years.csv']:
    shutil.copy2(previous / name, destination / name)
pd.DataFrame(summary).to_csv(destination / 'gfdl_mld_summary.csv', index=False)
metadata = json.loads((previous / 'metadata.json').read_text(encoding='utf-8'))
metadata['GFDL_MLD'] = 'Author explicitly confirms source year6 as analysis year1. Initial source6–15, peak source145=analysis140; final source271–280. Residual reported using this existing phase-relative convention. Prior missing-first5 explanation and physical140 extraction are superseded. No raw series or figures changed.'
metadata['audit_sources'] = [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
                             for path in [history_path, previous / 'table1_peak_state_residual.csv']]
metadata['word_paths'] = [str(root / 'docs/table1_cmip6_peak_and_residual_gfdl_start6.docx'),
                          str(root / 'docs/table2_variables_and_periods_gfdl_start6.docx')]
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print(pd.DataFrame(summary).to_string(index=False), flush=True)
