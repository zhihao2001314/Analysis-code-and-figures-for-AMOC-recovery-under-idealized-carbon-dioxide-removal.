# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import xarray as xr


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
base = root / 'output/s05_atmosphere/p05_amoc_recovery_composites'
destination = base / 'n65_global_precipitation_rainband'
assert not destination.exists(), '正式计算结果存在，需先归档。'
models = ['GFDL-ESM4', 'NorESM2-LM', 'CESM2', 'CNRM-ESM2-1']
members = ['r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f2']
starts = [[1, 271], [1, 271], [1, 131], [1850, 2120]]
periods, seasons = ['early', 'later'], ['Annual', 'DJF']
latitude, longitude = np.arange(-88., 90., 2.), np.arange(-180., 180., 2.)
rainband_bounds = [-20., 20.]
input_registry = base / 'n01_endpoint_groups/qc.json'
records = json.loads(input_registry.read_text(encoding='utf-8'))['inputs']
with xr.open_dataset(base / 'n35_global_meridional_samples/member_climatology.nc') as dataset:
    old_members = dataset.pr.load()
with xr.open_dataset(base / 'n35_global_meridional_samples/group_climatology.nc') as dataset:
    old_groups = dataset.pr.load()
audit_rows, reproduction_rows, model_fields, rainband_rows = [], [], [], []


# %% 3. 原生全球降水按实际日数平均，独立复现旧全球纬向曲线

for model_index, model in enumerate(models):
    period_fields = []
    for period_index, period in enumerate(periods):
        first_year = starts[model_index][period_index]
        paths = sorted(set(row['file'] for row in records if row['model'] == model
                           and row['period'] == period and row['variable'] == 'pr'))
        weighted_sums, total_days = {}, {}
        for year in range(first_year, first_year + 10):
            pieces = []
            for filename in paths:
                with xr.open_dataset(filename, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as dataset:
                    assert dataset.pr.dims == ('time', 'lat', 'lon')
                    assert dataset.pr.attrs['units'] == 'kg m-2 s-1'
                    assert dataset.attrs['source_id'] == model and dataset.attrs['variant_label'] == members[model_index]
                    take = dataset.time.dt.year == year
                    if not bool(take.any()):
                        continue
                    raw = dataset.pr.isel(time=take).load()
                    pieces.append(raw.astype(float) * 86400.)
                    source_path = Path(filename)
                    audit_rows.append({'model': model, 'period': period, 'year': year, 'file': filename,
                                       'dimensions': 'time,lat,lon', 'months': int(take.sum()),
                                       'file_bytes': source_path.stat().st_size,
                                       'selected_values_sha256': hashlib.sha256(raw.values.tobytes()).hexdigest()})
            monthly = xr.concat(pieces, 'time').sortby('time')
            assert np.array_equal(monthly.time.dt.month.values, np.arange(1, 13))
            assert monthly.time.dt.year.values.tolist() == [year] * 12
            assert np.isfinite(monthly.values).all()
            assert np.allclose(np.diff(monthly.lon), 360 / monthly.sizes['lon'], atol=1e-4)
            for season in seasons:
                take = np.ones(12, dtype=bool) if season == 'Annual' else np.isin(monthly.time.dt.month, [1, 2, 12])
                branch = monthly.isel(time=take)
                weights = branch.time.dt.days_in_month.astype(float)
                year_sum = (branch * weights).sum('time', skipna=False)
                day_sum = float(weights.sum())
                if year == first_year:
                    weighted_sums[season], total_days[season] = year_sum, day_sum
                else:
                    weighted_sums[season] = weighted_sums[season] + year_sum
                    total_days[season] += day_sum
        season_fields = []
        for season in seasons:
            native = weighted_sums[season] / total_days[season]
            # 曲线继续使用旧n35值；本次从原始全球降水独立复现，核对窗口及单位。
            reproduction = native.mean('lon').interp(lat=latitude)
            previous = old_members.sel(model=model, period=period, season=season)
            error = float(abs(reproduction - previous).max())
            assert error < 1e-10, (model, period, season, error)
            reproduction_rows.append({'model': model, 'period': period, 'season': season,
                                      'maximum_old_zonal_error_mm_day': error, 'days': total_days[season]})
            # 统一到2°公共格点，周期经度插值；不外推到极点，不以该地图算全球水量守恒。
            native = native.assign_coords(lon=((native.lon + 180) % 360) - 180).sortby('lon').sortby('lat')
            periodic = xr.concat([native.isel(lon=[-1]).assign_coords(lon=native.lon.values[[-1]] - 360),
                                  native, native.isel(lon=[0]).assign_coords(lon=native.lon.values[[0]] + 360)], 'lon')
            common = periodic.interp(lat=latitude, lon=longitude).transpose('lat', 'lon')
            assert np.isfinite(common.values).all()
            season_fields.append(common.expand_dims(season=[season]))
        period_fields.append(xr.concat(season_fields, 'season').expand_dims(period=[period]))
        print(model, period, 'complete', flush=True)
    model_fields.append(xr.concat(period_fields, 'period').expand_dims(model=[model]))


# %% 4. 模式等权组平均、差值与逐经度热带最大降水格点

member_maps = xr.concat(model_fields, 'model').rename('pr').transpose('model', 'period', 'season', 'lat', 'lon')
group_maps = xr.concat([member_maps.sel(model=models[:2]).mean('model', skipna=False).expand_dims(product=['recovered']),
                       member_maps.sel(model=models[2:]).mean('model', skipna=False).expand_dims(product=['not_recovered'])], 'product')
group_maps = group_maps.transpose('product', 'period', 'season', 'lat', 'lon')
changes = group_maps.sel(period='later', drop=True) - group_maps.sel(period='early', drop=True)
contrast = (changes.sel(product='not_recovered', drop=True) - changes.sel(product='recovered', drop=True)).expand_dims(product=['contrast'])
changes = xr.concat([changes, contrast], 'product').transpose('product', 'season', 'lat', 'lon')
for product in ['recovered', 'not_recovered']:
    for season in seasons:
        for period in periods:
            tropical = group_maps.sel(product=product, season=season, period=period, lat=slice(*rainband_bounds))
            values = tropical.values
            indices = np.argmax(values, axis=0)
            maxima = np.max(values, axis=0)
            for longitude_index, lon in enumerate(longitude):
                peak_lat = float(tropical.lat.values[indices[longitude_index]])
                rainband_rows.append({'product': product, 'season': season, 'period': period, 'lon': lon,
                                      'peak_lat': peak_lat, 'peak_pr_mm_day': maxima[longitude_index],
                                      'boundary_peak': abs(peak_lat) == 20.,
                                      'n_exact_ties': int(np.sum(values[:, longitude_index] == maxima[longitude_index]))})
rainband = pd.DataFrame(rainband_rows)
paired = rainband.pivot(index=['product', 'season', 'lon'], columns='period', values='peak_lat').reset_index()
paired['later_minus_early_degN'] = paired.later - paired.early


# %% 5. 保存地图、标记与原曲线，记录实际维序和限制

destination.mkdir(parents=True)
for field in [member_maps, group_maps, changes]:
    field.attrs['units'] = 'mm day-1'
member_maps.to_netcdf(destination / 'member_precipitation_maps.nc')
group_maps.to_netcdf(destination / 'group_precipitation_maps.nc')
changes.to_netcdf(destination / 'group_precipitation_changes.nc')
# 保留旧曲线全纬度数据和既有ITCZ中位纬度，不按新地图重新诊断它们。
old_groups.to_netcdf(destination / 'global_zonal_precipitation.nc')
positions = pd.read_csv(base / 'n36_global_circulation_positions/endpoint_positions.csv')
positions.loc[positions['product'].isin(['recovered', 'not_recovered']) &
              positions.metric.eq('ITCZ_precipitation_median')].to_csv(destination / 'itcz_median_positions.csv', index=False)
rainband.to_csv(destination / 'rainband_peak_gridpoints.csv', index=False)
paired.to_csv(destination / 'rainband_peak_gridpoint_changes.csv', index=False)
pd.DataFrame(audit_rows).to_csv(destination / 'input_audit.csv', index=False)
pd.DataFrame(reproduction_rows).to_csv(destination / 'zonal_reproduction.csv', index=False)
source_paths = [input_registry, base / 'n35_global_meridional_samples/member_climatology.nc',
                base / 'n35_global_meridional_samples/group_climatology.nc',
                base / 'n35_global_meridional_samples/group_changes.nc',
                base / 'n36_global_circulation_positions/endpoint_positions.csv']
metadata = {'groups': {'recovered': models[:2], 'not_recovered': models[2:]},
            'periods': 'Initial1-10;Later271-280 in combined experiment;actual source starts recorded below.',
            'starts': starts, 'models': models, 'seasons': seasons,
            'season_definition': 'Annual12months;DJF Jan/Feb/Dec of same calendar year, preserving n35;actual day weighting over10years.',
            'source_dimensions': 'All8 pr sources confirmed(time,lat,lon);kg m-2 s-1 multiplied86400 to mm/day. Native arrays checked before selection, all monthly values finite.',
            'output_dimensions': {'member_maps': ['model', 'period', 'season', 'lat', 'lon'],
                                  'group_maps': ['product', 'period', 'season', 'lat', 'lon'],
                                  'changes': ['product', 'season', 'lat', 'lon'],
                                  'zonal_profiles': ['product', 'period', 'season', 'lat']},
            'map_grid': '2degree:lat-88..88,lon-180..178;bilinear linear lat/lon interpolation after dayweighted native climatology;periodic longitude;not conservative regridding;polar caps not extrapolated.',
            'peak_definition': 'User-approved single largest precipitation gridpoint at each longitude within20S..20N, diagnosed from each2-model group map. All180 longitudes shown per period. In exact tie choose southernmost gridpoint. No peak smoothing, no significance dots.',
            'median_definition': 'Old n36 precipitation area-weighted median latitude within20S..20N, all longitudes;unchanged. This differs from the map single-peak tracer.',
            'contrast': 'Delta incomplete minus Delta substantial. No single initial/later rainband exists for this contrast, so no invented peak markers on contrast map.',
            'cautions': ['Peak tracers show wettest tropical grid cells, not a dynamical ITCZ diagnosis or rigid displacement of the entire rainbelt.',
                         'A single maximum can jump between northern/southern branches or land monsoon systems;ties, boundary maxima and raw gridpoint shifts are saved.',
                         'Maps have2degree latitude resolution;do not infer sub-grid peak shifts or significance. Unchanged peaks can coexist with changed precipitation median.',
                         'Curves are old global zonal averages with actual support88S..88N;axes90S..90N do not extrapolate poles.',
                         'No statistical or AMOC-only causal inference; groups each n=2.'],
            'source_outputs': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in source_paths]}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print('PR change absolute maxima by product:', abs(changes).max(['season', 'lat', 'lon']).values)
print('Peak boundary counts:', rainband.groupby(['product', 'season', 'period']).boundary_peak.sum().to_string())
print(destination)
