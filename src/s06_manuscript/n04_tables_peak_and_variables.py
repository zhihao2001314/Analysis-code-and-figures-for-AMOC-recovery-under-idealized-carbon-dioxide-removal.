# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import sys

import numpy as np
import pandas as pd
import xarray as xr

sys.stdout.reconfigure(encoding='utf-8')


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
destination = root / 'output/s06_manuscript/n04_tables_peak_and_variables'
assert not destination.exists(), '正式结果已存在，须先归档。'
models = ['CESM2', 'CNRM-ESM2-1', 'NorESM2-LM', 'GFDL-ESM4',
          'ACCESS-ESM1-5', 'CanESM5', 'MIROC-ES2L', 'UKESM1-0-LL', 'CAS-ESM2-0']
configuration_path = root / 'output/s06_manuscript/n01_table1_quantitative/table1_quantitative.csv'
temperature_path = root / 'output/s01_amoc/n52_nine_model_amoc_temperature/annual_paths.csv'
temperature_months_path = root / 'output/s01_amoc/n52_nine_model_amoc_temperature/standard_monthly_tas.csv'
ice_path = root / 'output/s02sia/s02_arctic_region/p01_regional_sia/n01/arctic_regional_sia_annual.csv'
old_ice_path = root / 'output/s02sia/s02_arctic_region/p04_pattern_unadjusted/n01_arctic_sia_pattern/arctic_sia_pattern_summary.csv'
mld_path = root / 'output/s02_sia/s02_arctic_region/p21_grouped_mld/n061_grouped_mld_timeseries/member_timeseries.csv'
cas_mld_path = root / 'output/s02_sia/s02_arctic_region/p25_cas_hydrography/n083_cas_mld_stratification/cas_march_mld_timeseries.csv'
legacy_time_qa_path = root / 'output/s02sia/s02_arctic_region/p06_mld/n04_na_region_timeseries/na_march_mld_region_qa.csv'
cas_tas_paths = [
    Path('F:/sea/cmip6tasad/CAS/tas_Amon_CAS-ESM2-0_1pctCO2_r1i1p1f1_gn_020001-035012.nc'),
    Path('F:/sea/cmip6tas/CAS/tas_Amon_CAS-ESM2-0_1pctCO2-cdr_r1i1p1f1_gn_034001-056912.nc'),
]
source_paths = [configuration_path, temperature_path, temperature_months_path, ice_path,
                old_ice_path, mld_path, cas_mld_path, legacy_time_qa_path]


# %% 3. 读取已保存结果并核对CAS标准气温

configurations = pd.read_csv(configuration_path).set_index('model')
temperature = pd.read_csv(temperature_path)
temperature_months = pd.read_csv(temperature_months_path)
ice = pd.read_csv(ice_path)
ice = ice.loc[ice.region.eq('Arctic total 45-90N')].copy()
old_ice = pd.read_csv(old_ice_path)
old_ice = old_ice.loc[old_ice.region.eq('Arctic total (45-90N)')].set_index('model')
mld = pd.read_csv(mld_path)
cas_mld = pd.read_csv(cas_mld_path)
legacy_time_qa = pd.read_csv(legacy_time_qa_path)
gfdl_up_qa = legacy_time_qa.loc[legacy_time_qa.model.eq('GFDL-ESM4') & legacy_time_qa.experiment.eq('1pctCO2')].iloc[0]
assert gfdl_up_qa.thetao_first_march_year == gfdl_up_qa.so_first_march_year == 6
assert gfdl_up_qa.thetao_last_march_year == gfdl_up_qa.so_last_march_year == 145
# 作者确认源第6年为分析第1年；分析峰值140年对应源145年。
gfdl_peak_legacy_index = 140
cas_rows, inventory = [], []
for path, experiment, start_year in zip(cas_tas_paths, ['1pctCO2', '1pctCO2-cdr'], [200, 340]):
    with xr.open_dataset(path, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as dataset:
        assert dataset.tas.dims == ('time', 'lat', 'lon')
        assert dataset.attrs['source_id'] == 'CAS-ESM2-0'
        assert dataset.attrs['experiment_id'] == experiment
        assert dataset.attrs['variant_label'] == 'r1i1p1f1' and dataset.tas.attrs['units'] == 'K'
        bounds = dataset[dataset.lat.attrs['bounds']].values
        latitude_weights = np.abs(np.diff(np.sin(np.deg2rad(bounds)), axis=1)[:, 0])
        np.testing.assert_allclose(latitude_weights.sum(), 2, atol=1e-6)
        np.testing.assert_allclose(np.diff(dataset.lon.values), 360 / dataset.sizes['lon'], atol=1e-4)
        weights = latitude_weights[:, None] * np.ones((1, dataset.sizes['lon']))
        weights /= weights.sum()
        years = list(range(200, 210)) + [339] if experiment == '1pctCO2' else list(range(470, 480))
        for calendar_year in years:
            take = np.flatnonzero(dataset.time.dt.year.values == calendar_year)
            assert len(take) == 12
            field = dataset.tas.isel(time=take).load()
            assert np.isfinite(field.values).all()
            monthly = np.einsum('tyx,yx->t', field.values.astype(float), weights)
            days = field.time.dt.days_in_month.values.astype(float)
            annual_k = float(np.average(monthly, weights=days))
            phase_index = int(experiment == '1pctCO2-cdr')
            experiment_year = calendar_year - start_year + 1
            cas_rows.append({'model': 'CAS-ESM2-0', 'experiment': experiment,
                             'calendar_year': calendar_year, 'simulation_year': experiment_year + phase_index * 140,
                             'tas_K': annual_k, 'tas_degC': annual_k - 273.15})
        inventory.append({'path': str(path), 'variable': 'tas', 'experiment': experiment,
                          'dimensions': list(dataset.tas.dims), 'shape': list(dataset.tas.shape),
                          'selected_calendar_years': years, 'units': 'K',
                          'source': dataset.attrs['source'], 'institution_id': dataset.attrs['institution_id'],
                          'bytes': path.stat().st_size, 'mtime_ns': path.stat().st_mtime_ns})
cas_temperature = pd.DataFrame(cas_rows)
print('CAS standard tas annual means read: 21 years', flush=True)


# %% 4. 汇总第140年实际值及末期残余变化

rows, qc_rows = [], []
for model in models:
    if model == 'CAS-ESM2-0':
        record = {'model': model, 'model_label': model, 'institution': 'CAS\n(China)',
                  'member': 'r1i1p1f1', 'atmosphere': 'IAP AGCM 5.0\n128 × 256 grid',
                  'ocean': 'LICOM2.0 / ~1°'}
        series_t = cas_temperature.set_index('simulation_year').tas_K
    else:
        record = configurations.loc[model, ['model_label', 'institution', 'member', 'atmosphere', 'ocean']].to_dict()
        record['model'] = model
        selected_t = temperature.loc[temperature.model.eq(model)]
        assert selected_t.temperature_kind.eq('standard_2m_tas').all()
        series_t = selected_t.set_index('simulation_year').gmst_standard_K
    record['initial_tas_degC'] = float(series_t.loc[1:10].mean() - 273.15)
    record['peak_tas_degC'] = float(series_t.loc[140] - 273.15)
    record['late_tas_degC'] = float(series_t.loc[271:280].mean() - 273.15)
    record['delta_tas_degC'] = record['late_tas_degC'] - record['initial_tas_degC']
    selected_ice = ice.loc[ice.model.eq(model)].set_index('simulation_year')
    assert len(selected_ice) == 280 and selected_ice.index.is_unique
    series_ice = selected_ice.raw_sia_million_km2
    record['initial_sia_million_km2'] = float(series_ice.loc[1:10].mean())
    record['peak_sia_million_km2'] = float(series_ice.loc[140])
    record['late_sia_million_km2'] = float(series_ice.loc[271:280].mean())
    record['delta_sia_million_km2'] = record['late_sia_million_km2'] - record['initial_sia_million_km2']
    for region, key in [('Labrador-Irminger', 'li'), ('Nordic Seas', 'ns')]:
        if model == 'CAS-ESM2-0':
            selected_mld = cas_mld.loc[cas_mld.version.eq('rectangles') & cas_mld.region.eq(region)]
            # 与原八模式矩形序列保持同一最近10m参考层口径。
            series_mld = selected_mld.set_index('simulation_year').mld_nearest_crossing_m
        else:
            selected_mld = mld.loc[mld.model.eq(model) & mld.region.eq(region)]
            series_mld = selected_mld.set_index('simulation_year').march_mld_m
        assert len(series_mld) == 280 and series_mld.index.is_unique
        record['initial_mld_' + key + '_m'] = float(series_mld.loc[1:10].mean())
        record['late_mld_' + key + '_m'] = float(series_mld.loc[271:280].mean())
        peak_index = gfdl_peak_legacy_index if model == 'GFDL-ESM4' else 140
        record['peak_mld_' + key + '_m'] = float(series_mld.loc[peak_index])
        record['delta_mld_' + key + '_m'] = float(series_mld.loc[271:280].mean() - series_mld.loc[1:10].mean())
    rows.append(record)
    qc_rows.append({'model': model, 'old_sia_delta_million_km2': float(old_ice.loc[model, 'raw_delta_sia_million_km2']),
                    'new_sia_delta_million_km2': record['delta_sia_million_km2'],
                    'difference_from_old_spatial_product_million_km2': record['delta_sia_million_km2'] - float(old_ice.loc[model, 'raw_delta_sia_million_km2']),
                    'mld_peak_time_verified': True,
                    'mld_initial_analysis_1_10_available': True})
table = pd.DataFrame(rows)
assert np.isfinite(table[['peak_tas_degC', 'peak_sia_million_km2', 'delta_tas_degC', 'delta_sia_million_km2']]).all().all()
assert table.loc[table.model.ne('GFDL-ESM4'), ['peak_mld_li_m', 'peak_mld_ns_m', 'delta_mld_li_m', 'delta_mld_ns_m']].notna().all().all()
for model in models[:-1]:
    np.testing.assert_allclose(table.set_index('model').loc[model, 'delta_tas_degC'],
                               configurations.loc[model, 'delta_global_tas_degC'], atol=1e-10)


# %% 5. 整理变量、观测时段和实际文件年份

variable_records = [
    ['CMIP6', 'msftmz / msftyz', '9 models', '1–280', 'Annual', 'AMOC strength and regional overturning'],
    ['CMIP6', 'thetao, so', '9 models', '1–280; initial profiles: 1–10*', 'March', 'Mixed-layer depth and N²; ocean temperature'],
    ['CMIP6', 'tas', '9 models', '1–10, 140, 271–280; full series where used', 'Annual', 'Near-surface air temperature'],
    ['CMIP6', 'siconc', '9 models', '1–280', 'Annual', 'Sea-ice area and 15% ice edge'],
    ['CMIP6', 'pr', '4 models†', '1–10 and 271–280', 'Annual', 'Precipitation response'],
    ['CMIP6', 'va', '4 models†', '1–10 and 271–280', 'Annual', 'Atmospheric meridional circulation'],
    ['CMIP6', 'rsdt, rsut, rlut', '4 models†', '1–280; piControl reference', 'Annual', 'TOA energy imbalance and SW/LW components'],
    ['CMIP6', 'mlotst', '7 models‡', '1–10* and 271–280', 'March', 'Comparison with temperature–salinity-derived MLD'],
    ['CMIP6', 'areacello', 'As required', 'Time invariant', 'Fixed', 'Native ocean/ice grid-cell areas'],
    ['IAP derived product', 'N²', 'Observations', '2005–2024', 'March', 'Upper-ocean stratification reference'],
    ['NOAA/NSIDC SIC CDR v6', 'Sea-ice concentration', 'Observations', '1988–2020', 'Monthly → annual', 'Observed 15% ice-edge reference'],
]
variables = pd.DataFrame(variable_records, columns=['source', 'variables', 'coverage', 'period_used', 'sampling', 'application'])
calendar_records = []
for model in models:
    selected_t = temperature_months.loc[temperature_months.model.eq(model)]
    if model == 'CAS-ESM2-0':
        calendar_records.extend([
            [model, 'tas', '0200–0209', '0339', '0470–0479'],
            [model, 'thetao, so', '0200–0209', '0339', '0131–0140'],
            [model, 'siconc', '0200–0209', '0339', 'CDR months 1561–1680§'],
        ])
    else:
        early_years = selected_t.loc[selected_t.simulation_year.between(1, 10), 'calendar_year'].astype(int)
        peak_years = selected_t.loc[selected_t.simulation_year.eq(140), 'calendar_year'].astype(int)
        late_years = selected_t.loc[selected_t.simulation_year.between(271, 280), 'calendar_year'].astype(int)
        assert early_years.nunique() == late_years.nunique() == 10 and peak_years.nunique() == 1
        calendar_records.append([model, 'tas; standard forcing labels',
                                 f'{early_years.min():04d}–{early_years.max():04d}',
                                 f'{peak_years.iloc[0]:04d}', f'{late_years.min():04d}–{late_years.max():04d}'])
calendar = pd.DataFrame(calendar_records, columns=['model', 'variable_scope', 'initial_years', 'peak_year', 'final_years'])


# %% 6. 保存表格数据和核查记录

destination.mkdir(parents=True, exist_ok=True)
table.to_csv(destination / 'table1_peak_state_residual.csv', index=False)
variables.to_csv(destination / 'table2_variables_periods.csv', index=False)
calendar.to_csv(destination / 'table2_source_calendar_windows.csv', index=False)
cas_temperature.to_csv(destination / 'cas_standard_tas_selected_years.csv', index=False)
pd.DataFrame(qc_rows).to_csv(destination / 'comparison_with_previous_table.csv', index=False)
metadata = {
    'peak': 'Actual state in experiment year140; NOT a change or a131–140 mean. TAS and SIA annual; MLD March.',
    'residual': 'Mean of years271–280 minus mean of years1–10. Unsmoothed, no newly introduced piControl drift correction.',
    'tas': 'Native-area global TAS in degreesC; annual actual-month-day weights. Original8 values preserved. CAS from standard2m monthly TAS, no proxy.',
    'sia': 'Annual mean raw SIA45–90N from internally consistent nine-model saved annual series; sum SIC times native area. Some last decimals differ from old endpoint-map product; comparisons saved.',
    'mld': 'March native-area rectangular means, nearest valid layer to10m sigma0 reference, threshold0.03kg/m3, first crossing; no-cross excluded. CAS uses matched nearest-layer diagnostic.',
    'mld_regions': {'LI': '50<=lat<65;65W<=lon<=20W', 'NS': '65<=lat<=80;30W<=lon<=20E'},
    'GFDL_MLD': 'Author confirms source6–145 corresponds to analysis1–140. Initial source6–15, peak source145, final source271–280. Residual reported under this phase-relative convention; do not treat source1–5 as missing analysis years.',
    'CAS_configuration': {'official_source': 'https://wcrp-cmip.github.io/CMIP6_CVs/docs/CMIP6_source_id.html',
                          'atmosphere': 'IAP AGCM5.0,128latitude x256longitude; do not confuse grid count with degrees',
                          'ocean': 'LICOM2.0 nominal1degree,30levels', 'member': 'r1i1p1f1'},
    'dimensions': {'CAS_tas': 'time,lat,lon;128x256', 'saved_tables': 'long-form named columns; no positional scientific axes'},
    'input_inventory': inventory,
    'sources': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in source_paths],
    'word_paths': [str(root / 'docs/table1_cmip6_peak_and_residual.docx'), str(root / 'docs/table2_variables_and_periods.docx')],
}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print(table[['model', 'peak_tas_degC', 'peak_sia_million_km2', 'peak_mld_li_m', 'peak_mld_ns_m', 'delta_tas_degC', 'delta_sia_million_km2']].to_string(index=False), flush=True)
print(destination.resolve(), flush=True)
