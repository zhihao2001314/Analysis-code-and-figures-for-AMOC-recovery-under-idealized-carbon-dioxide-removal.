# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import xarray as xr
from scipy.integrate import trapezoid


# %% 2. 参数与强度参照

root = Path('F:/sea/Antic sea ice constraint')
base = root / 'output/s05_atmosphere/p05_amoc_recovery_composites'
destination = base / 'check_circulation_precipitation_strength'
assert not destination.exists(), '正式核验结果已存在，须先归档。'
models = ['GFDL-ESM4', 'NorESM2-LM', 'CESM2', 'CNRM-ESM2-1']
starts = [[1, 271], [1, 271], [1, 131], [1850, 2120]]
periods = ['early', 'later']
pressure_levels = np.array([1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 25000, 30000], dtype=float)
targets = [-10., 10.]
radius, gravity = 6371000., 9.80665
audit = pd.read_csv(base / 'n35_global_meridional_samples/input_audit.csv')
with xr.open_dataset(base / 'n35_global_meridional_samples/member_climatology.nc') as dataset:
    climate = dataset.load()
with xr.open_dataset(base / 'n35_global_meridional_samples/group_changes.nc') as dataset:
    changes = dataset.load()
raw_rows, strength_rows, source_rows = [], [], []


# %% 3. 四模式初末首年原始v独立重建300hPa质量流函数

for model_index, model in enumerate(models):
    for period_index, period in enumerate(periods):
        year = starts[model_index][period_index]
        paths = audit.loc[audit.model.eq(model) & audit.period.eq(period) & audit.year.eq(year) & audit.variable.eq('va'), 'file'].unique()
        assert len(paths) == 1
        path = Path(paths[0])
        with xr.open_dataset(path, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as dataset:
            assert dataset.attrs['source_id'] == model
            expected_experiment = '1pctCO2' if period == 'early' else '1pctCO2-cdr'
            assert dataset.attrs['experiment_id'] == expected_experiment
            assert dataset.va.attrs['units'] == 'm s-1'
            source = dataset.sel(lat=slice(-14, 14)).isel(time=dataset.time.dt.year == year)
            assert source.sizes['time'] == 12
            np.testing.assert_array_equal(source.time.dt.month, np.arange(1, 13))
            days = source.time.dt.days_in_month.values.astype(float)
            if model == 'NorESM2-LM':
                assert source.va.dims == ('time', 'lev', 'lat', 'lon')
                native = source.va.load().values.astype(float)
                native_ps = source.ps.load().values.astype(float)
                pressure = source.a.values[None, :, None, None] * float(source.p0) + source.b.values[None, :, None, None] * native_ps[:, None]
                if np.mean(pressure[:, 0]) > np.mean(pressure[:, -1]):
                    pressure, native = pressure[:, ::-1], native[:, ::-1]
                assert np.all(np.diff(pressure, axis=1) > 0)
                assert np.max(pressure[:, 0]) <= pressure_levels[0] and np.min(pressure[:, -1]) >= pressure_levels[-1]
                # 使用逐柱np.interp独立于原计算的矢量化层索引实现。
                values = np.empty((12, len(pressure_levels), source.sizes['lat'], source.sizes['lon']))
                for month in range(12):
                    for latitude_index in range(source.sizes['lat']):
                        for longitude_index in range(source.sizes['lon']):
                            values[month, :, latitude_index, longitude_index] = np.interp(
                                np.log(pressure_levels), np.log(pressure[month, :, latitude_index, longitude_index]),
                                native[month, :, latitude_index, longitude_index])
            else:
                assert source.va.dims == ('time', 'plev', 'lat', 'lon') and source.plev.attrs['units'] == 'Pa'
                values = source.va.sel(plev=pressure_levels).load().values.astype(float)
            assert np.isfinite(values).all()
            latitude = source.lat.values
            transport = trapezoid(values.mean(axis=3), x=pressure_levels, axis=1)
            psi_monthly = transport * (2 * np.pi * radius / gravity) * np.cos(np.deg2rad(latitude))[None] / 1e9
            psi_annual = np.average(psi_monthly, weights=days, axis=0)
            with xr.open_dataset(base / 'n35_global_meridional_samples' / (model + '_samples.nc')) as saved:
                stored = saved.psi.sel(model=model, period=period, season='Annual', year_in_window=1, plev=30000)
                for target in targets:
                    expected = np.interp(target, latitude, psi_annual)
                    actual = float(stored.sel(lat=target))
                    error = actual - expected
                    assert abs(error) < 1e-8
                    raw_rows.append({'model': model, 'period': period, 'source_year': year,
                                     'latitude': target, 'pressure_Pa': 30000,
                                     'independent_psi_1e9kgs': expected, 'saved_psi_1e9kgs': actual,
                                     'error_1e9kgs': error})
            source_rows.append({'model': model, 'period': period, 'path': str(path), 'year': year,
                                'dimensions': str(source.va.dims), 'units': source.va.attrs['units'],
                                'member': dataset.attrs['variant_label'], 'calendar': source.time.dt.calendar})
        print(model + ' ' + period + ': raw streamfunction verified', flush=True)


# %% 4. 文献量级、正负号与底部残余审计

for model in models:
    for period in periods:
        field = climate.psi.sel(model=model, period=period, season='Annual')
        nh = field.sel(plev=50000, lat=slice(0, 30))
        sh = field.sel(plev=50000, lat=slice(-30, 0))
        tropical = field.sel(plev=slice(20000, 70000), lat=slice(-30, 30))
        record = {'model': model, 'period': period,
                  'NH500_max_1e9kgs': float(nh.max()), 'SH500_min_1e9kgs': float(sh.min()),
                  'tropical_max_1e9kgs': float(tropical.max()), 'tropical_min_1e9kgs': float(tropical.min()),
                  'surface_residual_absmax_1e9kgs': float(abs(field.sel(plev=100000)).max())}
        assert record['NH500_max_1e9kgs'] > 0 and record['SH500_min_1e9kgs'] < 0
        if period == 'early':
            # 仅量级筛查，不是模式技能/偏差验收阈值；历史再分析与理想化初期不能逐值比。
            assert 30 < record['tropical_max_1e9kgs'] < 300
            assert 30 < -record['tropical_min_1e9kgs'] < 300
        strength_rows.append(record)
assert float(abs(climate.psi.sel(plev=1000)).max()) == 0.
np.testing.assert_allclose(changes.psi.sel(product='contrast'),
                           changes.psi.sel(product='not_recovered') - changes.psi.sel(product='recovered'), atol=1e-12, rtol=0)


# %% 5. 保存验收门槛及不可外推的结论

destination.mkdir(parents=True)
pd.DataFrame(raw_rows).to_csv(destination / 'independent_raw_streamfunction.csv', index=False)
pd.DataFrame(source_rows).to_csv(destination / 'raw_source_inventory.csv', index=False)
pd.DataFrame(strength_rows).to_csv(destination / 'annual_strength_audit.csv', index=False)
report = {'status': 'passed', 'scope': 'Independent computation/units/signs/order-of-magnitude gate, not validation of model circulation biases.',
          'raw_months_checked': 96, 'raw_comparisons': len(raw_rows),
          'raw_max_error_1e9kgs': float(pd.DataFrame(raw_rows).error_1e9kgs.abs().max()),
          'formula': '2*pi*a*cos(lat)/g integral_10hPa^p zonal_mean(v) dp, then divide1e9.',
          'literature_context': 'Annual HC extrema typically order1e11kg/s; seasonal winter extrema may exceed2e11. Compare units and broad scale, not equality across epochs/models.',
          'references': ['https://doi.org/10.1175/JCLI-D-12-00224.1', 'https://doi.org/10.1007/s13143-021-00256-z'],
          'initial_tropical_order_of_magnitude_screen_1e9kgs': [30, 300],
          'interpretation': 'Initial contour values about75-115 in1e9kg/s mean0.75-1.15e11kg/s, not75kg/s and not ocean Sv.',
          'limits': '10hPa finite top; no column mass adjustment; nonzero near-surface residual documented. This is not an observational match, model-skill certification, or causal/significance test.',
          'code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(destination / 'metadata.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
print(pd.DataFrame(strength_rows).to_string(index=False), flush=True)
print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
