# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import xarray as xr


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
destination = root / 'output/s01_amoc/n52_nine_model_amoc_temperature'
assert not destination.exists(), '正式结果已存在，须先归档。'
models = ['ACCESS-ESM1-5', 'CanESM5', 'CAS-ESM2-0', 'CESM2', 'CNRM-ESM2-1',
          'GFDL-ESM4', 'MIROC-ES2L', 'NorESM2-LM', 'UKESM1-0-LL']
starts = {'ACCESS-ESM1-5': [101, 241], 'CanESM5': [1850, 1991], 'CAS-ESM2-0': [200, 340],
          'CESM2': [1, 1], 'CNRM-ESM2-1': [1850, 1990], 'GFDL-ESM4': [1, 141],
          'MIROC-ES2L': [1850, 1989], 'NorESM2-LM': [1, 141], 'UKESM1-0-LL': [1850, 1990]}
members = {'CanESM5': 'r1i1p2f1', 'CNRM-ESM2-1': 'r1i1p1f2',
           'MIROC-ES2L': 'r1i1p1f2', 'UKESM1-0-LL': 'r1i1p1f2'}
reused_models = ['CESM2', 'CNRM-ESM2-1', 'GFDL-ESM4', 'NorESM2-LM']
previous = root / 'output/s01_amoc/n24_amoc_temperature_paths'
amoc_path = root / 'output/amoc_s01_n01/l_40/CMIP6_AMOC_9models_annual_series.csv'
definition_path = root / 'output/amoc_s01_n01/l_40/CMIP6_AMOC_9models_processing_metadata.csv'
proxy_path = Path('F:/sea/cmip6tas/CAS/tas200-499.surface.annual.nc')
years = np.arange(1, 281)
monthly_parts, audit_records, coverage_records, annual_parts, summaries, arrow_records = [], [], [], [], [], []


# %% 3. 复用已审计四模式月序列，读取其余模式标准气温

amoc_source = pd.read_csv(amoc_path)
definitions = pd.read_csv(definition_path).set_index('Model')
old_metadata = json.loads((previous / 'metadata.json').read_text(encoding='utf-8'))
assert hashlib.sha256(amoc_path.read_bytes()).hexdigest() == old_metadata['source_amoc_sha256']
old_monthly = pd.read_csv(previous / 'global_monthly_tas.csv')
old_audit = pd.read_csv(previous / 'input_audit.csv')
for model in models:
    for phase_index, (phase, experiment, folder) in enumerate([
            ('increase', '1pctCO2', 'cmip6tasad'), ('removal', '1pctCO2-cdr', 'cmip6tas')]):
        first_year = starts[model][phase_index]
        member = members.get(model, 'r1i1p1f1')
        assert definitions.loc[model, 'Selected member'] == member
        if model in reused_models:
            monthly = old_monthly.loc[old_monthly.model.eq(model) & old_monthly.phase.eq(phase)].copy()
            monthly['origin'] = 'audited_n24_monthly_series'
            audit_records.extend(old_audit.loc[old_audit.model.eq(model) & old_audit.phase.eq(phase)].to_dict('records'))
        else:
            parts = []
            paths = sorted(Path('F:/sea', folder).rglob(f'tas_Amon_{model}_{experiment}_{member}_*.nc'))
            for path in paths:
                with xr.open_dataset(path, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as dataset:
                    assert dataset.tas.dims == ('time', 'lat', 'lon') and dataset.tas.attrs['units'] == 'K'
                    assert dataset.attrs['source_id'] == model and dataset.attrs['experiment_id'] == experiment
                    assert dataset.attrs['variant_label'] == member
                    take = np.flatnonzero((dataset.time.dt.year.values >= first_year)
                                          & (dataset.time.dt.year.values < first_year + 140))
                    if not take.size:
                        continue
                    bounds = dataset[dataset.lat.attrs['bounds']]
                    assert bounds.dims == ('lat', 'bnds') and bounds.shape == (dataset.sizes['lat'], 2)
                    latitude_weights = np.abs(np.diff(np.sin(np.deg2rad(bounds.values)), axis=1)[:, 0])
                    np.testing.assert_allclose(latitude_weights.sum(), 2, rtol=0, atol=1e-6)
                    np.testing.assert_allclose(np.diff(dataset.lon.values), 360 / dataset.sizes['lon'], rtol=0, atol=1e-4)
                    weights = latitude_weights[:, None] * np.ones((1, dataset.sizes['lon']))
                    weights /= weights.sum()
                    for offset in range(0, len(take), 120):
                        indices = take[offset:offset + 120]
                        field = dataset.tas.isel(time=indices).load()
                        values = field.values.astype(float)
                        assert np.isfinite(values).all(), '不以零替代全球场缺测。'
                        mean_tas = np.einsum('tyx,yx->t', values, weights, optimize=True)
                        part = pd.DataFrame({'calendar_year': field.time.dt.year.values,
                                             'month': field.time.dt.month.values,
                                             'days': field.time.dt.days_in_month.values,
                                             'gmst_K': mean_tas, 'source_file': str(path), 'source_index': indices})
                        parts.append(part)
                    audit_records.append({'model': model, 'phase': phase, 'file': str(path), 'variable': 'tas',
                                          'dimensions': str(dataset.tas.dims), 'shape': str(dataset.tas.shape),
                                          'selected_months': len(take), 'units': 'K', 'member': member,
                                          'area_weighting': 'native latitude-bound spherical cell areas; uniform longitude',
                                          'file_bytes': path.stat().st_size, 'mtime_ns': path.stat().st_mtime_ns})
            if parts:
                monthly = pd.concat(parts, ignore_index=True).sort_values(['calendar_year', 'month'])
                monthly['model'], monthly['phase'] = model, phase
                monthly['simulation_year'] = monthly.calendar_year - first_year + 1 + phase_index * 140
                monthly['origin'] = 'raw_standard_tas'
            else:
                assert model == 'CAS-ESM2-0' and phase == 'removal'
                monthly = pd.DataFrame(columns=['calendar_year', 'month', 'days', 'gmst_K', 'model', 'phase', 'simulation_year'])
        if not monthly.empty:
            ids = monthly.calendar_year.to_numpy() * 12 + monthly.month.to_numpy()
            assert len(np.unique(ids)) == len(ids), '月份重复必须先审计，不能任意覆盖。'
            assert np.all(np.diff(ids) > 0)
            monthly_parts.append(monthly)
        available = np.asarray(monthly.calendar_year * 12 + monthly.month, dtype=int)
        expected = np.arange(first_year * 12 + 1, (first_year + 140) * 12 + 1)
        missing = np.setdiff1d(expected, available)
        coverage_records.append({'model': model, 'phase': phase, 'standard_months': len(available),
                                 'expected_months': 1680, 'missing_months': len(missing),
                                 'missing_calendar_month_ids': missing.tolist()})
        print(model, phase, 'standard months', len(available), flush=True)
monthly_all = pd.concat(monthly_parts, ignore_index=True)


# %% 4. CAS最低模式层代理单独处理，禁止与标准2米气温拼接

cas_native_path = next(Path('F:/sea/cmip6tasad').rglob('tas_Amon_CAS-ESM2-0_1pctCO2_r1i1p1f1_*.nc'))
with xr.open_dataset(cas_native_path, decode_times=False) as dataset:
    cas_lat = dataset.lat.values
    cas_lon = dataset.lon.values
    bounds = dataset.lat_bnds.values
    cas_weights = np.abs(np.diff(np.sin(np.deg2rad(bounds)), axis=1)[:, 0])[:, None] * np.ones((1, len(cas_lon)))
    cas_weights /= cas_weights.sum()
with xr.open_dataset(proxy_path, decode_times=False) as dataset:
    assert dataset.tas.dims == ('time', 'lev', 'lat', 'lon') and dataset.tas.shape == (300, 1, 128, 256)
    assert dataset.tas.attrs['units'] == 'K' and 'time' not in dataset.coords
    # 派生文件将规则坐标重写为双精度；回到原网格单精度后须逐值完全一致。
    np.testing.assert_array_equal(dataset.lat.values.astype(np.float32), cas_lat.astype(np.float32))
    np.testing.assert_array_equal(dataset.lon.values.astype(np.float32), cas_lon.astype(np.float32))
    proxy_values = dataset.tas.isel(time=slice(0, 280), lev=0).load().values.astype(float)
    assert np.isfinite(proxy_values).all()
    cas_proxy = np.einsum('tyx,yx->t', proxy_values, cas_weights, optimize=True)
    proxy_level = float(dataset.lev.values[0])
    audit_records.append({'model': 'CAS-ESM2-0', 'phase': 'both_proxy', 'file': str(proxy_path), 'variable': 'tas',
                          'dimensions': str(dataset.tas.dims), 'shape': str(dataset.tas.shape), 'selected_records': 280,
                          'units': 'K', 'member': 'not recorded in source', 'area_weighting': 'same CAS native latitude bounds',
                          'file_bytes': proxy_path.stat().st_size, 'mtime_ns': proxy_path.stat().st_mtime_ns})


# %% 5. 对齐逐年AMOC、计算异常和分阶段滑动均值

for model in models:
    data = pd.DataFrame({'model': model, 'simulation_year': years,
                         'phase': np.where(years <= 140, 'increase', 'removal')}).set_index('simulation_year')
    annual_standard = np.full(280, np.nan)
    model_monthly = monthly_all.loc[monthly_all.model.eq(model)]
    for year, group in model_monthly.groupby('simulation_year'):
        if len(group) == 12:
            np.testing.assert_array_equal(np.sort(group.month), np.arange(1, 13))
            annual_standard[int(year) - 1] = np.average(group.gmst_K, weights=group.days)
    data['gmst_standard_K'] = annual_standard
    # CAS始终用同一种低层代理，不将两种温度变量接成一条路径。
    path_temperature = cas_proxy if model == 'CAS-ESM2-0' else annual_standard
    data['temperature_path_K'] = path_temperature
    data['temperature_kind'] = 'model_level_proxy_provisional' if model == 'CAS-ESM2-0' else 'standard_2m_tas'
    baseline = float(np.mean(path_temperature[:10]))
    data['temperature_initial10_K'] = baseline
    data['warming_raw_degC'] = path_temperature - baseline
    data['gmst_standard_warming_degC'] = annual_standard - np.mean(annual_standard[:10])
    original = amoc_source.loc[amoc_source.Model.eq(model)].set_index('Simulation year').loc[years]
    data['amoc_raw_sv'] = original['Forced AMOC strength (Sv)'].to_numpy()
    data['amoc_calendar_year'] = original['Calendar year'].to_numpy()
    assert np.isfinite(data.amoc_raw_sv).all()
    for first, last in [(1, 140), (141, 280)]:
        for raw, smooth in [('warming_raw_degC', 'warming_smooth11_degC'), ('amoc_raw_sv', 'amoc_smooth11_sv')]:
            part = data.loc[first:last, raw]
            # 分阶段平滑、不跨强迫转折；缺测中心保持缺测，不能把缺口向两侧填平。
            smoothed = part.rolling(11, center=True, min_periods=6).mean().where(part.notna())
            data.loc[first:last, smooth] = smoothed
    initial_amoc = float(data.loc[1:10, 'amoc_raw_sv'].mean())
    late_amoc = float(data.loc[271:280, 'amoc_raw_sv'].mean())
    summaries.append({'model': model, 'temperature_kind': data.temperature_kind.iloc[0],
                      'standard_tas_valid_years': int(np.isfinite(annual_standard).sum()),
                      'path_valid_years': int(np.isfinite(path_temperature).sum()),
                      'amoc_initial10_sv': initial_amoc, 'amoc_later10_sv': late_amoc,
                      'amoc_later_over_initial': late_amoc / initial_amoc,
                      'late_warming_degC': float(data.loc[271:280, 'warming_raw_degC'].mean()),
                      'max_smooth_warming_degC': float(data.warming_smooth11_degC.max()),
                      'temperature_peak_simulation_year': int(data.warming_smooth11_degC.idxmax()),
                      'amoc_min_smooth_sv': float(data.amoc_smooth11_sv.min()),
                      'amoc_min_simulation_year': int(data.amoc_smooth11_sv.idxmin()),
                      'amoc_actual_latitude_degN': float(definitions.loc[model, 'Actual latitude (degN)'])})
    for first, last in [(55, 65), (105, 115), (175, 185), (240, 250)]:
        xy = data.loc[[first, last], ['warming_smooth11_degC', 'amoc_smooth11_sv']].to_numpy()
        if np.isfinite(xy).all():
            arrow_records.append({'model': model, 'phase': 'increase' if first < 140 else 'removal',
                                  'first_year': first, 'last_year': last, 'x0': xy[0, 0], 'y0': xy[0, 1],
                                  'x1': xy[1, 0], 'y1': xy[1, 1]})
    annual_parts.append(data.reset_index())
annual_all = pd.concat(annual_parts, ignore_index=True)
cas_check = annual_all.loc[annual_all.model.eq('CAS-ESM2-0')].iloc[:140]
cas_offset = float(np.mean(cas_check.temperature_path_K - cas_check.gmst_standard_K))
cas_anomaly_rmse = float(np.sqrt(np.mean((cas_check.warming_raw_degC - cas_check.gmst_standard_warming_degC) ** 2)))


# %% 6. 保存数值、来源和适用范围

destination.mkdir(parents=True)
annual_all.to_csv(destination / 'annual_paths.csv', index=False)
monthly_all.to_csv(destination / 'standard_monthly_tas.csv', index=False)
pd.DataFrame(summaries).to_csv(destination / 'endpoint_summary.csv', index=False)
pd.DataFrame(arrow_records).to_csv(destination / 'direction_arrows.csv', index=False)
pd.DataFrame(audit_records).to_csv(destination / 'input_audit.csv', index=False)
definitions.loc[models].to_csv(destination / 'amoc_definitions.csv')
(destination / 'coverage.json').write_text(json.dumps(coverage_records, ensure_ascii=False, indent=2), encoding='utf-8')
source_paths = [amoc_path, definition_path, previous / 'global_monthly_tas.csv', previous / 'input_audit.csv', proxy_path]
metadata = {'status': 'complete_with_explicit_qualifications', 'models': models, 'years': [1, 280],
            'starts': starts, 'members': definitions.loc[models, 'Selected member'].to_dict(),
            'amoc_definition': 'Existing nearest40N Atlantic-Arctic maximum below500m; raw Forced AMOC strength (Sv) unchanged.',
            'tas_dimensions': ['time', 'lat', 'lon'], 'cas_proxy_dimensions': ['time', 'lev', 'lat', 'lon'],
            'csv_dimensions': ['model', 'simulation_year'], 'monthly_dimensions': ['model', 'phase', 'calendar_year', 'month'],
            'temperature_reference': 'Per-model unsmoothed first10year mean of the SAME temperature variable; not piControl or1850-1900.',
            'annual_weighting': 'Standard TAS: actual month-day weights; CAS proxy: externally provided annual means, within-year averaging not independently recoverable.',
            'smoothing': '11year centered mean separately within each CO2 phase, min6; missing centers remain NaN; no gap interpolation.',
            'phase_labels': 'CO2 increase/removal, not assumed actual warming/cooling; paths remain in time order.',
            'UKESM': 'Newly supplied standard TAS CDR1990-2049 verified; complete1680months in both phases, no calendar shift.',
            'CAS': {'variable': 'lowest model-level Temperature proxy, not standard2m TAS', 'level': proxy_level,
                    'calendar': 'No decodable time/member/experiment attributes; exploratory records1-280 follow existing project assumption of combined ramp-up/CDR. Not independently confirmed calendar.',
                    'path': 'Same proxy used for both phases and initial10 baseline; standard up TAS retained separately for comparison.',
                    'mean_proxy_minus_standard_up_K': cas_offset, 'warming_anomaly_up_rmse_K': cas_anomaly_rmse,
                    'user_authorization': 'User explicitly requested using current CAS proxy for this exploratory nine-model figure.',
                    'quantitative_multimodel_inclusion': False},
            'scope': 'Nine-panel exploratory overview; do not modify existing eight-model groups; one member/model, no ensemble CI.',
            'drift_correction': False, 'piControl_used': False,
            'sources': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in source_paths]}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print(pd.DataFrame(summaries).to_string(index=False), flush=True)
print('CAS proxy comparison', cas_offset, cas_anomaly_rmse, flush=True)
