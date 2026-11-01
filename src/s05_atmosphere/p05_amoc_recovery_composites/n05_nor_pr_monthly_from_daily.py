# %% 1. 导入依赖

from pathlib import Path
import json
import numpy as np
import pandas as pd
import xarray as xr
import cftime


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
parent = root / 'output/s05_atmosphere/p05_amoc_recovery_composites'
download = parent / 'n04_nor_precipitation_input'
output = parent / 'n05_nor_pr_monthly_from_daily'
assert not output.exists(), '正式衍生产品已存在，先归档。'
verified = json.loads((download / 'download_verification.json').read_text(encoding='utf-8'))
records = [record for record in verified['files'] if record['experiment'] == '1pctCO2-cdr']
early_path = Path('F:/sea/cmip6prad/NORESM2-LM2/pr_Amon_NorESM2-LM_1pctCO2_r1i1p1f1_gn_000101-001012.nc')
monthly_fields = []
monthly_bounds = []
monthly_checks = []
annual_daily_integrals = {}
target_years = range(271, 281)
time_units = 'days since 0001-01-01 00:00:00'
calendar = 'noleap'


# %% 3. 核对原始月值和网格

with xr.open_dataset(early_path, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as source:
    assert source.pr.dims == ('time', 'lat', 'lon')
    assert source.pr.attrs['units'] == 'kg m-2 s-1'
    early_ids = source.time.dt.year.values * 12 + source.time.dt.month.values
    assert np.array_equal(early_ids, np.arange(13, 133))
    target_lat = source.lat.load()
    target_lon = source.lon.load()
    latitude_bounds = source.lat_bnds.load()
    longitude_bounds = source.lon_bnds.load()


# %% 4. 逐月积分原始日均降水，不填补缺日或无效格点

for record in records:
    assert record['sha256_matches_esgf']
    with xr.open_dataset(record['file'], decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as source:
        assert source.pr.dims == ('time', 'lat', 'lon')
        assert source.time_bnds.dims == ('time', 'bnds')
        assert source.time.encoding['calendar'] == calendar
        assert np.array_equal(source.lat, target_lat) and np.array_equal(source.lon, target_lon)
        source_ids = source.time.dt.year.values * 10000 + source.time.dt.month.values * 100 + source.time.dt.day.values
        assert len(source_ids) == len(np.unique(source_ids))
        for year in target_years:
            if not bool((source.time.dt.year == year).any()):
                continue
            # 独立年度日值积分，与逐月积分相加核对守恒。
            annual = source.pr.sel(time=source.time.dt.year == year).load().astype('float64')
            annual_bounds = source.time_bnds.sel(time=source.time.dt.year == year).values
            annual_duration = np.array([(end - start).total_seconds() / 86400. for start, end in annual_bounds])
            assert len(annual_duration) == 365 and np.all(annual_duration == 1.)
            annual_daily_integrals[year] = (annual * xr.DataArray(annual_duration, dims='time', coords={'time': annual.time})).sum('time', skipna=False).values
            for month in range(1, 13):
                selected = (source.time.dt.year == year) & (source.time.dt.month == month)
                field = source.pr.isel(time=selected).load().astype('float64')
                bounds = source.time_bnds.isel(time=selected).values
                durations = np.array([(end - start).total_seconds() / 86400. for start, end in bounds])
                days = cftime.DatetimeNoLeap(year, month, 1).daysinmonth
                assert field.sizes['time'] == days and np.all(durations == 1.)
                assert np.array_equal(field.time.dt.day, np.arange(1, days + 1))
                expected_start = cftime.DatetimeNoLeap(year, month, 1)
                expected_end = cftime.DatetimeNoLeap(year + int(month == 12), month % 12 + 1, 1)
                assert bounds[0, 0] == expected_start and bounds[-1, 1] == expected_end
                assert np.all(bounds[:-1, 1] == bounds[1:, 0])
                weights = xr.DataArray(durations, dims='time', coords={'time': field.time})
                value = (field * weights).sum('time', skipna=False) / weights.sum()
                assert np.isfinite(value).all() and float(value.min()) >= 0.
                midpoint = bounds[0, 0] + (bounds[-1, 1] - bounds[0, 0]) / 2
                monthly_fields.append(value.expand_dims(time=[midpoint]))
                monthly_bounds.append([bounds[0, 0], bounds[-1, 1]])
                monthly_checks.append({'year': year, 'month': month, 'days_present': days, 'expected_days': days, 'duration_days': float(durations.sum()), 'minimum_kg_m2_s': float(value.min()), 'maximum_kg_m2_s': float(value.max())})
            print('Daily-to-monthly verified', year, flush=True)


# %% 5. 独立守恒核对及保存，明确标为日值聚合衍生产品

monthly = xr.concat(monthly_fields, 'time').sortby('time')
assert monthly.dims == ('time', 'lat', 'lon') and monthly.sizes['time'] == 120
ids = monthly.time.dt.year.values * 12 + monthly.time.dt.month.values
assert np.array_equal(ids, np.arange(271 * 12 + 1, 281 * 12 + 1))
conservation_errors = []
for year in target_years:
    selected = monthly.sel(time=monthly.time.dt.year == year)
    reconstructed = (selected * selected.time.dt.days_in_month).sum('time', skipna=False).values
    error = float(np.max(abs(reconstructed - annual_daily_integrals[year])))
    assert error < 1.e-10
    conservation_errors.append({'year': year, 'max_annual_integral_error_kg_m2_s_times_day': error})
monthly.attrs = {'standard_name': 'precipitation_flux', 'long_name': 'Monthly mean precipitation derived from published daily mean pr', 'units': 'kg m-2 s-1', 'cell_methods': 'area: mean time: mean', 'source_frequency': 'day'}
result = monthly.to_dataset(name='pr')
result['lat_bnds'] = latitude_bounds
result['lon_bnds'] = longitude_bounds
# 拼接顺序已按文件时间递增，独立核对月界后以统一单位保存，避免混用原文件time原点。
result['time_bnds'] = xr.DataArray(cftime.date2num(np.array(monthly_bounds).ravel().tolist(), units=time_units, calendar=calendar).reshape(120, 2), dims=('time', 'bnds'))
result.time_bnds.attrs = {'units': time_units, 'calendar': calendar}
result.time.attrs['bounds'] = 'time_bnds'
result.attrs = {'Conventions': 'CF-1.7', 'source_id': 'NorESM2-LM', 'experiment_id': '1pctCO2-cdr', 'variant_label': 'r1i1p1f1', 'grid_label': 'gn',
                'table_id': 'derived_monthly_from_day', 'frequency': 'mon', 'variable_id': 'pr', 'is_derived': 1,
                'source_table_id': 'day', 'source_tracking_ids': ';'.join(record['tracking_id'] for record in records),
                'history': 'Duration-weighted monthly means of complete original daily records, 0271-0280;no interpolation,missing-day fill,drift correction,or regridding.',
                'not_original_Amon': 'This is a local derived monthly product,not an ESGF-published Amon file.'}
output.mkdir(parents=True)
result.to_netcdf(output / 'nor_pr_cdr_027101-028012_monthly_derived.nc', encoding={'pr': {'zlib': True, 'complevel': 4, 'dtype': 'float64'}, 'time': {'units': time_units, 'calendar': calendar}})
pd.DataFrame(monthly_checks).to_csv(output / 'monthly_day_coverage.csv', index=False)
(output / 'aggregation_qc.json').write_text(json.dumps({'sources': records, 'dimensions': list(monthly.dims), 'months': 120, 'days': 3650, 'calendar': calendar, 'raw_early_monthly_file': str(early_path), 'early_months': 120, 'grid_matches_early': True, 'yearly_conservation': conservation_errors, 'interpretation': 'Available for the same later271-280 monthly window as other models;derived from day,not original Amon;rawdaily retained;no group figures recalculated.'}, indent=2), encoding='utf-8')
print(output, flush=True)
