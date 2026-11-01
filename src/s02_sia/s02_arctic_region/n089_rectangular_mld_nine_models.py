# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import shutil

import numpy as np
import pandas as pd


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
destination = root / 'output/s02_sia/s02_arctic_region/p25_cas_hydrography/n089_rectangular_mld_nine_models'
rectangle_path = root / 'output/s02_sia/s02_arctic_region/p21_grouped_mld/n061_grouped_mld_timeseries/member_timeseries.csv'
group_figure_path = root / 'output/s02_sia/s02_arctic_region/p23_stratification_mld_composite/n077_regional_definition_comparison/rectangles/member_timeseries.csv'
cas_path = root / 'output/s02_sia/s02_arctic_region/p25_cas_hydrography/n083_cas_mld_stratification/cas_march_mld_timeseries.csv'
coast_directory = root / 'output/s02_sia/s02_arctic_region/p25_cas_hydrography/n085_nine_model_hydrography_panels'
table_path = root / 'output/s06_manuscript/n07_gfdl_start_year6/table1_peak_state_residual.csv'
models = ['ACCESS-ESM1-5', 'CanESM5', 'CESM2', 'CNRM-ESM2-1', 'GFDL-ESM4',
          'MIROC-ES2L', 'UKESM1-0-LL', 'NorESM2-LM', 'CAS-ESM2-0']
regions = ['Labrador-Irminger', 'Nordic Seas']
assert not destination.exists(), '正式结果已存在，须先归档。'


# %% 3. 复用首张分组图的矩形序列并加入同口径CAS

rectangles = pd.read_csv(rectangle_path)
group_figure = pd.read_csv(group_figure_path)
columns = ['model', 'region', 'simulation_year', 'march_mld_m', 'march_mld_11yr_m']
rectangle_table = rectangles[columns].sort_values(['model', 'region', 'simulation_year']).reset_index(drop=True)
group_table = group_figure[columns].sort_values(['model', 'region', 'simulation_year']).reset_index(drop=True)
matched_rectangle = rectangle_table.loc[rectangle_table.region.isin(group_table.region.unique())].reset_index(drop=True)
pd.testing.assert_frame_equal(matched_rectangle, group_table)
cas = pd.read_csv(cas_path)
cas = cas.loc[cas.version.eq('rectangles'), ['model', 'region', 'simulation_year',
                                           'mld_nearest_crossing_m', 'mld_nearest_crossing_m_11yr']].copy()
cas = cas.rename(columns={'mld_nearest_crossing_m': 'march_mld_m',
                          'mld_nearest_crossing_m_11yr': 'march_mld_11yr_m'})
history = pd.concat([rectangles[columns], cas[columns]], ignore_index=True)
assert len(history) == 9 * 2 * 280 and not history.duplicated(['model', 'region', 'simulation_year']).any()
coast = pd.read_csv(coast_directory / 'coastal_mld_timeseries.csv')
table1 = pd.read_csv(table_path).set_index('model')
summary_rows, comparison_rows = [], []
for model in models:
    for region, key in zip(regions, ['li', 'ns']):
        selected = history.loc[history.model.eq(model) & history.region.eq(region)].sort_values('simulation_year')
        np.testing.assert_array_equal(selected.simulation_year, np.arange(1, 281))
        np.testing.assert_allclose(selected.march_mld_11yr_m,
                                   selected.march_mld_m.rolling(11, center=True, min_periods=11).mean(),
                                   atol=1e-7, equal_nan=True)
        initial = float(selected.loc[selected.simulation_year.between(1, 10), 'march_mld_m'].mean())
        late = float(selected.loc[selected.simulation_year.between(271, 280), 'march_mld_m'].mean())
        delta = late - initial
        np.testing.assert_allclose([initial, late, delta],
                                   [table1.loc[model, 'initial_mld_' + key + '_m'],
                                    table1.loc[model, 'late_mld_' + key + '_m'],
                                    table1.loc[model, 'delta_mld_' + key + '_m']], atol=1e-7)
        summary_rows.append({'model': model, 'region': region, 'initial_march_mld_m': initial,
                             'late_march_mld_m': late, 'delta_march_mld_m': delta})
        old = coast.loc[coast.model.eq(model) & coast.region.eq(region)]
        old_initial = float(old.loc[old.simulation_year.between(1, 10), 'march_mld_mean_crossing_only_m'].mean())
        old_late = float(old.loc[old.simulation_year.between(271, 280), 'march_mld_mean_crossing_only_m'].mean())
        comparison_rows.append({'model': model, 'region': region,
                                'rectangle_initial_m': initial, 'rectangle_late_m': late,
                                'coastal_initial_m': old_initial, 'coastal_late_m': old_late,
                                'late_coastal_minus_rectangle_m': old_late - late})


# %% 4. 保存统一后的绘图输入及口径审计

destination.mkdir(parents=True)
history.to_csv(destination / 'rectangular_mld_timeseries.csv', index=False)
pd.DataFrame(summary_rows).to_csv(destination / 'endpoint_summary.csv', index=False)
pd.DataFrame(comparison_rows).to_csv(destination / 'rectangular_vs_coastal.csv', index=False)
shutil.copy2(coast_directory / 'atmospheric_co2.csv', destination / 'atmospheric_co2.csv')
prior_metadata = json.loads((coast_directory / 'metadata.json').read_text(encoding='utf-8'))
metadata = {
    'models': models, 'mld_y_limits': [0, 900],
    'co2_y_limits': prior_metadata['co2_y_limits'], 'co2_ticks': prior_metadata['co2_ticks'],
    'regions': {'Labrador-Irminger': '50<=lat<65N,65W<=lon<=20W',
                'Nordic Seas': '65<=lat<=80N,30W<=lon<=20E'},
    'method': 'Same rectangle definition and saved8-model regional values as selected n078 group figure; CAS nearest10m reference crossing diagnostic fromn083. TEOS10 sigma0+0.03kg/m3 first resolved crossing; no-cross excluded; native areacello weighted.',
    'smoothing': 'Reuse saved11-year centred means with full11years; first/last5years missing. No new smoothing or drift correction.',
    'GFDL': 'Author convention: source6 is analysis1; up6–145 maps1–140; CDR141–280 maps141–280. No index shift introduced.',
    'comparison_scope': 'Old individual panels used coastal polygons (Nordic–Iceland sector) and exactly interpolated10m density reference; rectangles use nearest valid native layer. Differences combine region and reference treatment, not new forced data or drift correction.',
    'dimensions': 'Long CSV model9 xregion2 xsimulation_year280; co2 simulation_year280. No raw NetCDF/MAT read.',
    'checks': 'All2240 original8-model LI rectangle values/curves match selected first figure; Nordic uses the same legacy rectangle product. All18 initial/late/changes match current Table1.',
    'sources': [{'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
                for p in [rectangle_path, group_figure_path, cas_path, table_path,
                          coast_directory / 'coastal_mld_timeseries.csv', coast_directory / 'atmospheric_co2.csv',
                          coast_directory / 'metadata.json']]}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
caption = ('March mixed-layer depth in nine models for the fixed rectangular Labrador–Irminger Seas (orange;50–65°N,65–20°W) '
           'and Nordic Seas (blue;65–80°N,30°W–20°E). Curves are native-area regional means shown as11-year centred running means, '
           'using the same rectangle diagnostics as the selected group figure and Table1. The grey dashed line denotes atmospheric CO₂ on the right axis. '
           'GFDL-ESM4 source year6 is analysis year1. No piControl drift correction is introduced.')
(destination / 'caption.txt').write_text(caption + '\n', encoding='utf-8')
print(pd.DataFrame(comparison_rows).loc[pd.DataFrame(comparison_rows).model.eq('GFDL-ESM4')].to_string(index=False), flush=True)
