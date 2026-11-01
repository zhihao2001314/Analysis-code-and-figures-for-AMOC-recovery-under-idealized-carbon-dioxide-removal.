# %% 1. 导入依赖

from pathlib import Path
import argparse
import hashlib
import json
import netCDF4 as nc
import numpy as np
import pandas as pd


# %% 2. 参数与路径设置

parser = argparse.ArgumentParser()
parser.add_argument('--model', required=True, choices=['CESM2', 'CNRM-ESM2-1', 'NorESM2-LM', 'GFDL-ESM4'])
arguments = parser.parse_args()
model = arguments.model
root = Path('F:/sea/Antic sea ice constraint')
stage = Path('s05_atmosphere/p07_toa_energy_balance')
destination = root / 'output' / stage / 'n10_picontrol_toa' / model
assert not destination.exists(), '正式结果存在，须先归档。'
member = 'r1i1p1f2' if model == 'CNRM-ESM2-1' else 'r1i1p1f1'
fit_windows = {'CESM2': [501, 780], 'CNRM-ESM2-1': [1850, 2129],
               'NorESM2-LM': [1600, 1879], 'GFDL-ESM4': [101, 380]}
variables = ['rsdt', 'rsut', 'rlut']
chunk_months = 24
inventory, monthly_parts = [], []
weight_path = root / 'output' / stage / 'n01_toa_energy_balance' / f'weights_{model}.npz'
with np.load(weight_path) as saved:
    latitude_reference, longitude_reference, weights = saved['lat'], saved['lon'], saved['weight']
assert weights.shape == (len(latitude_reference), len(longitude_reference))
np.testing.assert_allclose(weights.sum(), 1., atol=1e-12)
destination.mkdir(parents=True)


# %% 3. 逐文件检查真实维度并读取控制辐射

for variable in variables:
    paths = sorted(Path('F:/sea/cmip6' + variable + '-picontrol').rglob(
        f'{variable}_Amon_{model}_piControl_{member}_*.nc'))
    assert paths, (model, variable)
    for path in paths:
        with nc.Dataset(path) as dataset:
            assert dataset.source_id == model and dataset.experiment_id == 'piControl'
            assert dataset.variant_label == member and dataset.variable_id == variable
            field, time = dataset.variables[variable], dataset.variables['time']
            assert field.dimensions == ('time', 'lat', 'lon') and field.units == 'W m-2'
            np.testing.assert_allclose(dataset.variables['lat'][:], latitude_reference, atol=1e-12, rtol=0)
            np.testing.assert_allclose(dataset.variables['lon'][:], longitude_reference, atol=1e-12, rtol=0)
            dates = nc.num2date(time[:], time.units, calendar=time.calendar)
            years = np.array([date.year for date in dates])
            months = np.array([date.month for date in dates])
            np.testing.assert_array_equal(np.diff(years * 12 + months), 1)
            time_bounds = np.asarray(dataset.variables[time.bounds][:], dtype=float)
            assert time_bounds.shape == (len(time), 2)
            bounds_dates = nc.num2date(time_bounds.reshape(-1), time.units,
                                      calendar=time.calendar).reshape(-1, 2)
            days = np.array([(pair[1] - pair[0]).total_seconds() / 86400 for pair in bounds_dates])
            for index, pair in enumerate(bounds_dates):
                assert (pair[0].year, pair[0].month, pair[0].day) == (years[index], months[index], 1)
                assert pair[1].year * 12 + pair[1].month == years[index] * 12 + months[index] + 1
            assert np.all(days > 0)
            global_mean = np.empty(len(time))
            for offset in range(0, len(time), chunk_months):
                values = np.ma.asarray(field[offset:offset + chunk_months, :, :], dtype=np.float64)
                assert not np.ma.getmaskarray(values).any(), '完整全球通量不能跳过缺测或缩小分母。'
                values = np.asarray(values)
                assert np.isfinite(values).all()
                global_mean[offset:offset + len(values)] = np.einsum('tyx,yx->t', values, weights, optimize=True)
            monthly_parts.append(pd.DataFrame({'model': model, 'variable': variable,
                                               'source_year': years, 'month': months, 'days': days,
                                               'global_Wm2': global_mean, 'path': str(path),
                                               'file_index': np.arange(len(time))}))
            digest = hashlib.sha256()
            with path.open('rb') as handle:
                while True:
                    block = handle.read(4 * 1024 * 1024)
                    if not block:
                        break
                    digest.update(block)
            file_hash = digest.hexdigest()
            inventory.append({'model': model, 'variable': variable, 'path': str(path),
                              'sha256': file_hash, 'bytes': path.stat().st_size,
                              'member': dataset.variant_label, 'grid': dataset.grid_label,
                              'dimensions': list(field.dimensions), 'shape': list(field.shape),
                              'units': field.units, 'calendar': time.calendar, 'time_units': time.units,
                              'time_bounds_variable': time.bounds, 'first_year': int(years[0]),
                              'last_year': int(years[-1]), 'months': len(time)})
        print(model, variable, path.name, len(time_bounds), 'months complete', flush=True)


# %% 4. 按真实月界时长计算完整控制年

monthly = pd.concat(monthly_parts, ignore_index=True).sort_values(['variable', 'source_year', 'month'])
assert not monthly.duplicated(['variable', 'source_year', 'month']).any()
annual_parts = []
coverage = []
for variable in variables:
    part = monthly[monthly.variable == variable].copy()
    serial_month = part.source_year.to_numpy() * 12 + part.month.to_numpy()
    np.testing.assert_array_equal(np.diff(serial_month), 1)
    part['weighted_flux'] = part.global_Wm2 * part.days
    annual_variable = part.groupby('source_year').agg(
        months=('month', 'count'), days=('days', 'sum'), weighted=('weighted_flux', 'sum')).reset_index()
    assert (annual_variable.months == 12).all()
    annual_variable[variable] = annual_variable.weighted / annual_variable.days
    annual_parts.append(annual_variable[['source_year', variable]])
    coverage.append({'model': model, 'variable': variable,
                     'first_year': int(annual_variable.source_year.min()),
                     'last_year': int(annual_variable.source_year.max()),
                     'years': len(annual_variable), 'unique_continuous_months': len(part)})
annual = annual_parts[0].merge(annual_parts[1], on='source_year', how='outer').merge(
    annual_parts[2], on='source_year', how='outer').sort_values('source_year').reset_index(drop=True)
annual.insert(0, 'model', model)
annual['asr'] = annual.rsdt - annual.rsut
annual['negative_olr'] = -annual.rlut
annual['eei'] = annual.asr - annual.rlut
annual['in_matched_fit_window'] = annual.source_year.between(*fit_windows[model])
selected = annual[annual.in_matched_fit_window]
assert len(selected) == 280 and selected[variables + ['eei']].notna().all().all()
for name in variables + ['asr', 'negative_olr', 'eei']:
    annual[name + '_smooth11'] = annual[name].rolling(11, center=True, min_periods=11).mean()


# %% 5. 保存结果与处理口径

monthly.to_csv(destination / 'monthly_global_long.csv', index=False)
annual.to_csv(destination / 'annual_control.csv', index=False)
pd.DataFrame(inventory).to_csv(destination / 'source_inventory.csv', index=False)
pd.DataFrame(coverage).to_csv(destination / 'coverage.csv', index=False)
metadata = {'status': 'passed', 'model': model, 'member': member,
            'fit_window': fit_windows[model], 'complete_fit_years': len(selected),
            'formula': 'N=global(rsdt-rsut-rlut), positive downward, W m-2',
            'dimensions': 'raw(time,lat,lon); monthly(variable,time); annual(source_year)',
            'spatial_weights': 'Exactly same native global weights as n01, after checking actual coordinates.',
            'annual_weights': 'CF month-bound durations, all12months and all global cells required.',
            'smoothing': '11year centered complete window for display only, not used in fitting.',
            'source_files': inventory, 'coverage': coverage,
            'weight_file': str(weight_path), 'weight_sha256': hashlib.sha256(weight_path.read_bytes()).hexdigest(),
            'code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print('PASSED', model, len(annual), 'control years;', len(selected), 'fit years', flush=True)
