# %% 1. 导入依赖

from pathlib import Path
import json
import hashlib
import numpy as np
import pandas as pd
from lxml import html


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
destination = root / 'output/s06_manuscript/n01_table1_quantitative'
assert not destination.exists(), '正式结果存在，须先归档。'
configuration_path = root / 'docs/table1_cmip6_word.html'
temperature_path = root / 'output/s01_amoc/n52_nine_model_amoc_temperature/endpoint_summary.csv'
ice_path = root / 'output/s02sia/s02_arctic_region/p04_pattern_unadjusted/n01_arctic_sia_pattern/arctic_sia_pattern_summary.csv'
mld_path = root / 'output/s02_sia/s02_arctic_region/p21_grouped_mld/n061_grouped_mld_timeseries/member_timeseries.csv'
models = ['CESM2', 'CNRM-ESM2-1', 'NorESM2-LM', 'GFDL-ESM4',
          'ACCESS-ESM1-5', 'CanESM5', 'MIROC-ES2L', 'UKESM1-0-LL']


# %% 3. 读取原表配置和未经漂移订正的定量端点

source_html = configuration_path.read_text(encoding='utf-8').replace('<br>', '\n')
tree = html.fromstring(source_html)
original_rows = tree.xpath('//tbody/tr')
configuration = [[cell.text_content().strip() for cell in row.xpath('./td')[:5]] for row in original_rows]
temperature = pd.read_csv(temperature_path).set_index('model')
ice = pd.read_csv(ice_path)
ice = ice[ice.region.eq('Arctic total (45-90N)')].set_index('model')
mld = pd.read_csv(mld_path)
print('MLD columns:', mld.columns.to_list(), flush=True)
rows = []
for model, settings in zip(models, configuration):
    assert settings[0].replace('\n', '') == model
    assert temperature.loc[model, 'temperature_kind'] == 'standard_2m_tas'
    assert temperature.loc[model, 'standard_tas_valid_years'] == 280
    row = dict(zip(['model_label', 'institution', 'member', 'atmosphere', 'ocean'], settings))
    row['model'] = model
    row['delta_global_tas_degC'] = float(temperature.loc[model, 'late_warming_degC'])
    row['delta_arctic_sia_million_km2'] = float(ice.loc[model, 'raw_delta_sia_million_km2'])
    np.testing.assert_allclose(row['delta_arctic_sia_million_km2'],
                               ice.loc[model, 'cdr_year131_140_sia_million_km2'] - ice.loc[model, 'co2_initial10_sia_million_km2'], atol=1e-8)
    for region, key in [('Labrador-Irminger', 'li'), ('Nordic Seas', 'ns')]:
        selected = mld[mld.model.eq(model) & mld.region.eq(region)]
        early = selected[selected.simulation_year.between(1, 10)]
        late = selected[selected.simulation_year.between(271, 280)]
        assert len(early) == len(late) == 10
        row[f'initial_mld_{key}_m'] = float(early.march_mld_m.mean())
        row[f'late_mld_{key}_m'] = float(late.march_mld_m.mean())
        row[f'delta_mld_{key}_m'] = float(late.march_mld_m.mean() - early.march_mld_m.mean()) if model != 'GFDL-ESM4' else np.nan
    rows.append(row)


# %% 4. 保存正式表格数据及来源口径

destination.mkdir(parents=True)
table = pd.DataFrame(rows)
table.to_csv(destination / 'table1_quantitative.csv', index=False)
metadata = {'period': 'Later years271-280 minus initial years1-10; no added piControl correction.',
            'temperature': 'Global area-weighted near-surface air-temperature change from n52; standard tas only.',
            'ice': 'Annual Arctic sea-ice area45-90N, native-area raw endpoint product; not sea-ice extent or ice-edge displacement.',
            'mld': 'March regional mean diagnostic MLD from the current rectangular-region figure data; unsmoothed annual values, ten-year means.',
            'mld_regions': {'LI': '50<=latitude<65N,65W-20W', 'NS': '65-80N,30W-20E'},
            'GFDL_mld': 'Not reported: existing rectangular MLD series retains unresolved calendar alignment and earliest up thetao begins at year6. No proxy value silently substituted.',
            'configuration': 'First five columns retained from user-approved original Table1 HTML.',
            'rounding': 'Temperature2 decimals; SIA2 decimals in10^6km2; MLDnearest1m with sign.',
            'sources': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
                        for path in [configuration_path, temperature_path, ice_path, mld_path]]}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print(table[['model', 'delta_global_tas_degC', 'delta_arctic_sia_million_km2', 'delta_mld_li_m', 'delta_mld_ns_m']].to_string(index=False), flush=True)
