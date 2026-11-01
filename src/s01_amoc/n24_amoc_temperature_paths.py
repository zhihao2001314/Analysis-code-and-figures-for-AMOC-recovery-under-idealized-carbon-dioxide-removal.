# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import re
import numpy as np
import pandas as pd
import xarray as xr


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
destination = root / 'output/s01_amoc/n24_amoc_temperature_paths'
assert not destination.exists(), '正式结果已存在，须先归档。'
models = ['CESM2', 'CNRM-ESM2-1', 'GFDL-ESM4', 'NorESM2-LM']
groups = {'CESM2': 'not_recovered', 'CNRM-ESM2-1': 'not_recovered',
          'GFDL-ESM4': 'recovered', 'NorESM2-LM': 'recovered'}
starts = {'CESM2': [1, 1], 'CNRM-ESM2-1': [1850, 1990],
          'GFDL-ESM4': [1, 141], 'NorESM2-LM': [1, 141]}
years = np.arange(1, 281)
source_path = root / 'output/amoc_s01_n01/l_40/CMIP6_AMOC_9models_annual_series.csv'
source = pd.read_csv(source_path)
amoc = source.pivot(index='Simulation year', columns='Model', values='Forced AMOC strength (Sv)').loc[years, models]
annual_rows, monthly_rows, audit_rows, event_rows, sensitivity_rows, arrow_rows = [], [], [], [], [], []


# %% 3. 原生全球面积加权近地面气温，逐块读取完整280年

for model in models:
    model_annual = []
    member = 'r1i1p1f2' if model == 'CNRM-ESM2-1' else 'r1i1p1f1'
    for phase_index, phase in enumerate(['increase', 'removal']):
        experiment = '1pctCO2' if phase_index == 0 else '1pctCO2-cdr'
        first_year = starts[model][phase_index]
        folder = Path('F:/sea/cmip6tas' + ('ad' if phase_index == 0 else ''))
        paths = sorted(folder.rglob(f'tas_Amon_{model}_{experiment}_{member}_*.nc'))
        phase_months = []
        for path in paths:
            dates = re.search(r'_(\d{6})-(\d{6})', path.name)
            assert dates is not None
            if int(dates[1][:4]) > first_year + 139 or int(dates[2][:4]) < first_year:
                continue
            with xr.open_dataset(path, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as dataset:
                assert dataset.tas.dims == ('time', 'lat', 'lon') and dataset.tas.attrs['units'] == 'K'
                assert dataset.attrs['source_id'] == model and dataset.attrs['experiment_id'] == experiment
                assert dataset.attrs['variant_label'] == member
                take = np.flatnonzero((dataset.time.dt.year.values >= first_year) & (dataset.time.dt.year.values <= first_year + 139))
                if not take.size:
                    continue
                if model == 'CNRM-ESM2-1':
                    gaussian_nodes, latitude_weight = np.polynomial.legendre.leggauss(dataset.sizes['lat'])
                    np.testing.assert_allclose(dataset.lat.values, np.rad2deg(np.arcsin(gaussian_nodes)), rtol=0, atol=1e-4)
                    area_method = 'Gaussian quadrature latitude weights and uniform longitude spacing'
                else:
                    latitude_bounds = dataset[dataset.lat.attrs['bounds']]
                    expected_bound_dims = ('lat', 'nbnd') if model == 'CESM2' else ('lat', 'bnds')
                    assert latitude_bounds.dims == expected_bound_dims
                    latitude_weight = abs(np.diff(np.sin(np.deg2rad(latitude_bounds.values)), axis=1)[:, 0])
                    area_method = 'Native latitude-bound spherical cell areas and uniform longitude spacing'
                np.testing.assert_allclose(np.diff(dataset.lon.values), 360 / dataset.sizes['lon'], rtol=0, atol=1e-4)
                np.testing.assert_allclose(latitude_weight.sum(), 2, rtol=0, atol=1e-6)
                weights = latitude_weight[:, None] * np.ones((1, dataset.sizes['lon']))
                weights /= weights.sum()
                for offset in range(0, len(take), 120):
                    indices = take[offset:offset + 120]
                    field = dataset.tas.isel(time=indices).load()
                    values = field.values.astype(float)
                    assert np.isfinite(values).all(), '全球场缺测，禁止改变分母或填零。'
                    temperature = np.einsum('tyx,yx->t', values, weights, optimize=True)
                    part = pd.DataFrame({'calendar_year': field.time.dt.year.values, 'month': field.time.dt.month.values,
                                         'days': field.time.dt.days_in_month.values, 'gmst_K': temperature})
                    phase_months.append(part)
                audit_rows.append({'model': model, 'phase': phase, 'file': str(path), 'variable': 'tas',
                                   'dimensions': str(dataset.tas.dims), 'shape': str(dataset.tas.shape),
                                   'selected_months': len(take), 'units': 'K', 'member': member, 'area_weighting': area_method})
        monthly = pd.concat(phase_months).sort_values(['calendar_year', 'month']).reset_index(drop=True)
        ids = monthly.calendar_year.to_numpy() * 12 + monthly.month.to_numpy()
        np.testing.assert_array_equal(ids, np.arange(first_year * 12 + 1, (first_year + 140) * 12 + 1))
        annual_temperature = np.sum(monthly.gmst_K.to_numpy().reshape(140, 12) * monthly.days.to_numpy().reshape(140, 12), axis=1) / monthly.days.to_numpy().reshape(140, 12).sum(axis=1)
        model_annual.extend(annual_temperature)
        monthly['model'], monthly['phase'] = model, phase
        monthly['simulation_year'] = monthly.calendar_year - first_year + 1 + 140 * phase_index
        monthly_rows.append(monthly)
        print(model, phase, '140 complete global annual means', flush=True)
    gmst = np.asarray(model_annual)
    baseline = gmst[:10].mean()
    warming = gmst - baseline
    for smoothing in [5, 11, 21]:
        smoothed_amoc = amoc[model].rolling(smoothing, center=True, min_periods=smoothing // 2 + 1).mean().to_numpy()
        smoothed_warming = pd.Series(warming).rolling(smoothing, center=True, min_periods=smoothing // 2 + 1).mean().to_numpy()
        for change_window in [10, 20, 30]:
            half = change_window // 2
            for phase, bounds, direction in [('increase', (1, 140), 'weakening'), ('removal', (141, 280), 'strengthening')]:
                centers = np.arange(bounds[0] + half + smoothing // 2, bounds[1] - half - smoothing // 2 + 1)
                rates = (smoothed_amoc[centers - 1 + half] - smoothed_amoc[centers - 1 - half]) / change_window
                position = np.argmin(rates) if direction == 'weakening' else np.argmax(rates)
                year = int(centers[position])
                record = {'model': model, 'group': groups[model], 'phase': phase, 'event': 'fastest_' + direction,
                          'smoothing_years': smoothing, 'change_window_years': change_window,
                          'simulation_year': year, 'interval_first_year': year - half, 'interval_last_year': year + half,
                          'rate_sv_per_year': float(rates[position]), 'amoc_sv': float(smoothed_amoc[year - 1]),
                          'warming_degC': float(smoothed_warming[year - 1]),
                          'warming_interval_min_degC': float(np.min(smoothed_warming[year - half - 1:year + half])),
                          'warming_interval_max_degC': float(np.max(smoothed_warming[year - half - 1:year + half]))}
                sensitivity_rows.append(record)
                if smoothing == 11 and change_window == 20:
                    event_rows.append(record)
        if smoothing == 11:
            display_amoc, display_warming = smoothed_amoc.copy(), smoothed_warming.copy()
    for year in years:
        annual_rows.append({'model': model, 'group': groups[model], 'simulation_year': year,
                            'phase': 'increase' if year <= 140 else 'removal',
                            'gmst_K': gmst[year - 1], 'gmst_initial10_K': baseline,
                            'warming_raw_degC': warming[year - 1], 'warming_smooth11_degC': display_warming[year - 1],
                            'amoc_raw_sv': amoc.loc[year, model], 'amoc_smooth11_sv': display_amoc[year - 1]})
    for first, last in [(35, 45), (90, 100), (165, 175), (230, 240)]:
        arrow_rows.append({'model': model, 'first_year': first, 'last_year': last,
                           'x0': display_warming[first - 1], 'x1': display_warming[last - 1],
                           'y0': display_amoc[first - 1], 'y1': display_amoc[last - 1],
                           'phase': 'increase' if last <= 140 else 'removal'})


# %% 4. 保存温升路径和变化最快区间，禁止将其直接命名为临界温度

destination.mkdir(parents=True)
pd.DataFrame(annual_rows).to_csv(destination / 'annual_paths.csv', index=False)
pd.concat(monthly_rows, ignore_index=True).to_csv(destination / 'global_monthly_tas.csv', index=False)
pd.DataFrame(audit_rows).to_csv(destination / 'input_audit.csv', index=False)
pd.DataFrame(event_rows).to_csv(destination / 'fast_change_intervals.csv', index=False)
pd.DataFrame(sensitivity_rows).to_csv(destination / 'change_window_sensitivity.csv', index=False)
pd.DataFrame(arrow_rows).to_csv(destination / 'direction_arrows.csv', index=False)
metadata = {'models': models, 'groups': groups, 'source_amoc': str(source_path),
            'source_amoc_sha256': hashlib.sha256(source_path.read_bytes()).hexdigest(), 'amoc_column': 'Forced AMOC strength (Sv)',
            'amoc_definition': 'Existing nearest40N Atlantic-Arctic maximum below500m; unchanged from n20.',
            'tas_dimensions': ['time', 'lat', 'lon'], 'csv_dimensions': ['model', 'simulation_year'],
            'warming_reference': 'Per-model global annual TAS minus unsmoothed initial years1-10 global annual mean. Not exactly piControl or1850-1900.',
            'annual_weighting': 'Monthly GMST weighted by actual days in each calendar month; native global cell-area weights.',
            'time': 'Full140year1pctCO2 +140yearCDR; identical combined index to existing AMOC. Raw calendar starts recorded.',
            'starts': starts, 'smoothing': '11yearcentered mean, min6 values; 5/11/21 sensitivity retained.',
            'events': 'Largest negative20year AMOC change rate in CO2-increase phase; largest positive20year rate in removal phase. Event center is descriptive, not detected abrupt tipping or recovery onset.',
            'event_sensitivity': 'All5/11/21smoothing x10/20/30yearwindows retained. Centers and smoothing support stay within each CO2 phase.',
            'path': 'Always connected in time order, never sorted by temperature. CO2 phase is not necessarily a warming/cooling phase.',
            'limitation': 'Fastest change is not proof of an abrupt change, universal temperature threshold, equilibrium bistability or causality.', 'drift_correction': False}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print(pd.DataFrame(event_rows)[['model', 'event', 'simulation_year', 'warming_degC', 'rate_sv_per_year']].to_string(index=False), flush=True)
