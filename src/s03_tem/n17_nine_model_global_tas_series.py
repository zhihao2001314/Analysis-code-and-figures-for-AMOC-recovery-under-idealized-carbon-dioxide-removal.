# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import netCDF4 as nc
import numpy as np
import pandas as pd


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
destination = root / 'output/s03_tem/n17_nine_model_global_tas_series'
assert not destination.exists(), '正式结果存在，须先归档。'
color_source = root / 'output/s01_amoc/n79_seven_panel_hysteresis_40n/metadata.json'
models = ['ACCESS-ESM1-5', 'CanESM5', 'CESM2', 'CNRM-ESM2-1', 'GFDL-ESM4',
          'MIROC-ES2L', 'NorESM2-LM', 'UKESM1-0-LL', 'CAS-ESM2-0']
starts = {'ACCESS-ESM1-5': [101, 241], 'CanESM5': [1850, 1991], 'CESM2': [1, 1],
          'CNRM-ESM2-1': [1850, 1990], 'GFDL-ESM4': [1, 141], 'MIROC-ES2L': [1850, 1989],
          'NorESM2-LM': [1, 141], 'UKESM1-0-LL': [1850, 1990], 'CAS-ESM2-0': [200, 340]}
members = {'ACCESS-ESM1-5': 'r1i1p1f1', 'CanESM5': 'r1i1p2f1', 'CESM2': 'r1i1p1f1',
           'CNRM-ESM2-1': 'r1i1p1f2', 'GFDL-ESM4': 'r1i1p1f1', 'MIROC-ES2L': 'r1i1p1f2',
           'NorESM2-LM': 'r1i1p1f1', 'UKESM1-0-LL': 'r1i1p1f2', 'CAS-ESM2-0': 'r1i1p1f1'}
phases = [('increase', '1pctCO2', Path('F:/sea/cmip6tasad')),
          ('removal', '1pctCO2-cdr', Path('F:/sea/cmip6tas'))]
chunk_months, smooth_window = 24, 11
inventory, grid_records, monthly_parts, annual_parts, summary_rows = [], [], [], [], []
colors = json.loads(color_source.read_text(encoding='utf-8'))['model_colors']


# %% 3. 原生大气网格全球面积权重与标准月平均气温

for model in models:
    reference = sorted(phases[0][2].rglob(f'tas_Amon_{model}_1pctCO2_{members[model]}_*.nc'))[0]
    with nc.Dataset(reference) as dataset:
        assert dataset['tas'].dimensions == ('time', 'lat', 'lon')
        latitude = np.asarray(dataset['lat'][:], dtype=float)
        longitude = np.asarray(dataset['lon'][:], dtype=float)
        if model == 'CNRM-ESM2-1':
            nodes, latitude_weights = np.polynomial.legendre.leggauss(len(latitude))
            np.testing.assert_allclose(latitude, np.rad2deg(np.arcsin(nodes)), atol=1e-4, rtol=0)
            area_method = 'Gaussian quadrature weights; uniform periodic longitude.'
        else:
            latitude_bounds = np.asarray(dataset[dataset['lat'].bounds][:], dtype=float)
            assert latitude_bounds.shape == (len(latitude), 2)
            latitude_weights = np.abs(np.diff(np.sin(np.deg2rad(latitude_bounds)), axis=1)[:, 0])
            area_method = 'Native spherical latitude-bound areas; uniform periodic longitude.'
        np.testing.assert_allclose(latitude_weights.sum(), 2., atol=1e-6, rtol=0)
        np.testing.assert_allclose(np.diff(longitude), 360. / len(longitude), atol=1e-4, rtol=0)
        weights = latitude_weights[:, None] * np.ones((1, len(longitude)))
        weights /= weights.sum()
        grid_records.append({'model': model, 'latitude_count': len(latitude),
                             'longitude_count': len(longitude), 'method': area_method,
                             'weight_sum': float(weights.sum()), 'minimum_weight': float(weights.min())})

    for phase_index, (phase, experiment, folder) in enumerate(phases):
        first_year = starts[model][phase_index]
        paths = sorted(folder.rglob(f'tas_Amon_{model}_{experiment}_{members[model]}_*.nc'))
        assert paths, (model, experiment)
        for path in paths:
            with nc.Dataset(path) as dataset:
                assert dataset.source_id == model and dataset.experiment_id == experiment
                assert dataset.variant_label == members[model]
                field, time = dataset['tas'], dataset['time']
                assert field.dimensions == ('time', 'lat', 'lon') and field.units == 'K'
                np.testing.assert_array_equal(dataset['lat'][:], latitude)
                np.testing.assert_array_equal(dataset['lon'][:], longitude)
                dates = nc.num2date(time[:], time.units, calendar=time.calendar)
                source_years = np.array([date.year for date in dates])
                source_months = np.array([date.month for date in dates])
                selected = np.flatnonzero((source_years >= first_year) & (source_years < first_year + 140))
                if not len(selected):
                    continue
                bounds = nc.num2date(dataset[time.bounds][selected].reshape(-1), time.units,
                                    calendar=time.calendar).reshape(-1, 2)
                days = np.array([(pair[1] - pair[0]).total_seconds() / 86400. for pair in bounds])
                assert np.all(days > 0)
                global_tas = np.empty(len(selected))
                for offset in range(0, len(selected), chunk_months):
                    indices = selected[offset:offset + chunk_months]
                    values = np.ma.asarray(field[indices], dtype=np.float64)
                    assert not np.ma.getmaskarray(values).any()
                    values = np.asarray(values)
                    assert np.isfinite(values).all()
                    global_tas[offset:offset + len(indices)] = np.einsum('tyx,yx->t', values, weights, optimize=True)
                monthly_parts.append(pd.DataFrame({'model': model, 'phase': phase,
                    'simulation_year': source_years[selected] - first_year + 1 + phase_index * 140,
                    'source_year': source_years[selected], 'month': source_months[selected],
                    'days': days, 'global_tas_K': global_tas, 'source_file': str(path)}))
                inventory.append({'model': model, 'phase': phase, 'path': str(path),
                    'member': members[model], 'experiment': experiment,
                    'dimensions': str(field.dimensions), 'shape': str(field.shape),
                    'units': field.units, 'calendar': time.calendar, 'time_units': time.units,
                    'first_selected_month': str(dates[selected[0]]),
                    'last_selected_month': str(dates[selected[-1]]), 'selected_months': len(selected),
                    'file_bytes': path.stat().st_size, 'mtime_ns': path.stat().st_mtime_ns,
                    'branch_time_in_parent': dataset.__dict__.get('branch_time_in_parent'),
                    'branch_time_in_child': dataset.__dict__.get('branch_time_in_child'),
                    'parent_time_units': dataset.__dict__.get('parent_time_units')})
        print(model, phase, 'global monthly tas complete', flush=True)


# %% 4. 核对月份并计算年平均、初期异常和连续滑动均值

monthly = pd.concat(monthly_parts, ignore_index=True).sort_values(['model', 'simulation_year', 'month'])
assert len(monthly) == 9 * 280 * 12
assert not monthly.duplicated(['model', 'simulation_year', 'month']).any()
for model in models:
    selected = monthly.loc[monthly.model.eq(model)].copy()
    np.testing.assert_array_equal(selected.simulation_year * 12 + selected.month - 1, np.arange(12, 281 * 12))
    selected['weighted_tas_K'] = selected.global_tas_K * selected.days
    annual = selected.groupby('simulation_year').agg(source_year=('source_year', 'first'),
        phase=('phase', 'first'), months=('month', 'count'), days=('days', 'sum'),
        weighted_tas_K=('weighted_tas_K', 'sum')).reset_index()
    assert annual.months.eq(12).all()
    annual['model'] = model
    annual['global_tas_degC'] = annual.weighted_tas_K / annual.days - 273.15
    initial = float(annual.loc[annual.simulation_year.between(1, 10), 'global_tas_degC'].mean())
    annual['initial10_tas_degC'] = initial
    annual['tas_anomaly_degC'] = annual.global_tas_degC - initial
    # 全部280年使用同一连续窗口，不对转折点平移、插值或强行拼接。
    annual['tas_anomaly_smooth11_degC'] = annual.tas_anomaly_degC.rolling(
        smooth_window, center=True, min_periods=6).mean()
    assert np.isfinite(annual[['global_tas_degC', 'tas_anomaly_smooth11_degC']]).all().all()
    np.testing.assert_allclose(annual.loc[:9, 'tas_anomaly_degC'].mean(), 0., atol=1e-12)
    peak = annual.loc[annual.tas_anomaly_smooth11_degC.idxmax()]
    summary_rows.append({'model': model, 'initial10_tas_degC': initial,
        'year140_anomaly_degC': float(annual.loc[139, 'tas_anomaly_degC']),
        'year141_minus140_degC': float(annual.loc[140, 'global_tas_degC'] - annual.loc[139, 'global_tas_degC']),
        'final10_anomaly_degC': float(annual.loc[270:279, 'tas_anomaly_degC'].mean()),
        'maximum_smooth11_anomaly_degC': float(peak.tas_anomaly_smooth11_degC),
        'maximum_smooth11_year': int(peak.simulation_year)})
    annual_parts.append(annual[['model', 'simulation_year', 'source_year', 'phase', 'months', 'days',
        'global_tas_degC', 'initial10_tas_degC', 'tas_anomaly_degC', 'tas_anomaly_smooth11_degC']])
annual = pd.concat(annual_parts, ignore_index=True)


# %% 5. 保存数据、来源与处理说明

destination.mkdir(parents=True)
monthly.to_csv(destination / 'global_monthly_tas.csv', index=False)
annual.to_csv(destination / 'global_annual_tas.csv', index=False)
pd.DataFrame(summary_rows).to_csv(destination / 'temperature_summary.csv', index=False)
pd.DataFrame(grid_records).to_csv(destination / 'global_area_weights_qc.csv', index=False)
pd.DataFrame(inventory).to_csv(destination / 'source_inventory.csv', index=False)
caption = ('Global-mean near-surface air-temperature anomalies for nine models over experiment years1–280, '
    'relative to each model mean over years1–10. Solid coloured lines show11-year centred running means, '
    'computed continuously across the forcing reversal, using at least6 years at the two experiment endpoints. '
    'Annual global means use native atmospheric-grid area weights and CF month-duration weights. '
    'The vertical dotted line marks the transition from increasing to decreasing CO2 after year140. '
    'One ensemble member is used for each model. No multimodel mean or significance inference is shown.')
metadata = {'models': models, 'model_colors': colors, 'starts': starts, 'members': members,
    'quantity': 'Global mean standard near-surface air temperature (tas), over land and ocean.',
    'source_dimensions': 'All raw tas(time,lat,lon), unitsK; one-dimensional latitude/longitude.',
    'output_dimensions': 'CSV:model,simulation_year; monthly:model,simulation_year,month.',
    'reference': 'Each model unsmoothed experiment years1–10 mean; anomalies in degreesC.',
    'spatial_weighting': 'Native spherical latitude-bound cell areas, with Gaussian quadrature for CNRM; no regridding.',
    'annual_weighting': 'Actual CF monthly time-bound duration; equal weights among annual means in baseline.',
    'smoothing': 'Continuous11year centred means across280years, min6 at experiment endpoints; no filling or offset adjustment.',
    'piControl_used': False, 'drift_correction': 'none', 'global_month_count': len(monthly),
    'GFDL_alignment': 'Source1–140 plus141–280, continuous; no five-year overlap from old regional diagnostics.',
    'CAS_source': 'Standard monthly Amon tas in both phases; no older model-layer annual proxy.',
    'CAS_branch_caveat': 'Dates339→340 and complete months verified. Existing branch-time metadata inconsistency remains unresolved; no artificial offset correction.',
    'shared_y_limits_degC': [-.4, float(np.ceil(annual.tas_anomaly_smooth11_degC.max() * 2) / 2 + .3)],
    'caption': caption, 'color_source': str(color_source),
    'color_source_sha256': hashlib.sha256(color_source.read_bytes()).hexdigest()}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print(pd.DataFrame(summary_rows).to_string(index=False), flush=True)
print(destination.resolve(), flush=True)
