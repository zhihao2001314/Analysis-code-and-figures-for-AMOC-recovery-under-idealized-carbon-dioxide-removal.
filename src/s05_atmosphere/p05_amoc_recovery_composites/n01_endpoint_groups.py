# %% 1. 导入依赖

from pathlib import Path
import json
import numpy as np
import pandas as pd
import xarray as xr
from scipy.spatial import cKDTree
import shapely
import cartopy.io.shapereader as shp


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
output = root / 'output/s05_atmosphere/p05_amoc_recovery_composites/n01_endpoint_groups'
assert not output.exists(), '正式产品已存在；先按项目规则归档。'
ocean = root / 'output/s02sia/s02_arctic_region/p14_mld_fixed_regions_raw'
jet = root / 'output/s05_atmosphere/p02_ua'
models = ['ACCESS-ESM1-5', 'CanESM5', 'CESM2', 'CNRM-ESM2-1', 'GFDL-ESM4', 'MIROC-ES2L', 'NorESM2-LM', 'UKESM1-0-LL']
members = dict(zip(models, ['r1i1p1f1', 'r1i1p2f1', 'r1i1p1f1', 'r1i1p1f2', 'r1i1p1f1', 'r1i1p1f2', 'r1i1p1f1', 'r1i1p1f2']))
starts = dict(zip(models, [(101, 371), (1850, 2121), (1, 131), (1850, 2120), (1, 271), (1850, 2119), (1, 271), (1850, 2120)]))
ice_dims = dict(zip(models, [('j', 'i'), ('j', 'i'), ('nj', 'ni'), ('y', 'x'), ('y', 'x'), ('y', 'x'), ('j', 'i'), ('j', 'i')]))
ice_coords = dict(zip(models, [('latitude', 'longitude'), ('latitude', 'longitude'), ('lat', 'lon'), ('lat', 'lon'), ('lat', 'lon'), ('latitude', 'longitude'), ('latitude', 'longitude'), ('latitude', 'longitude')]))
periods = ['early', 'later']
initial_bias_ratio = 1.5
amoc_recovery_threshold = 0.8
excluded_models = ['UKESM1-0-LL']
radius_m = 6371000.
mapping_cutoff_m = 150000.
nor_monthly_pr = root / 'output/s05_atmosphere/p05_amoc_recovery_composites/n05_nor_pr_monthly_from_daily/nor_pr_cdr_027101-028012_monthly_derived.nc'
assert nor_monthly_pr.exists()
provenance = []
availability = []


# %% 3. 保留海洋指标分类，单独记录用户指定的成员排除

profiles = pd.read_csv(ocean / 'coast_region_20261012_initial_n2_profiles.csv')
observations = pd.read_csv(ocean / 'coast_region_20261012_n2_1000m_amoc_recovery_observations.csv')
obs = observations.set_index('region').loc['Labrador-Irminger Seas']
n2_metadata = json.loads((ocean / 'coast_region_20261012_initial_n2_metadata.json').read_text(encoding='utf-8'))
amoc = pd.read_csv(root / 'output/amoc_s01_n01/l_40/CMIP6_AMOC_9models_annual_series.csv')
mld = pd.read_csv(ocean / 'na_march_mld_coast_regions_20261012_regions_annual.csv')
# GFDL旧MLD表含物理年错位，仅读已纠正的连续年份序列；缺1–5年不补。
gfdl = pd.read_csv(ocean / 'gfdl_coast_stabilization_20261013_annual.csv')
gfdl = gfdl.rename(columns={'mld_mean_m': 'march_mld_mean_crossing_only_m', 'deep_area_1000m_percent': 'resolved_mld_ge_1000m_fraction'})
gfdl['resolved_mld_ge_1000m_fraction'] /= 100.
gfdl['model'] = 'GFDL-ESM4'
mld = pd.concat([mld[mld.model != 'GFDL-ESM4'], gfdl], ignore_index=True)
records = []
normalized_amoc = []
for model in models:
    profile = profiles[(profiles.model == model) & (profiles.region == 'Labrador-Irminger Seas') & profiles.depth_m.between(15, 995)].sort_values('depth_m')
    assert np.array_equal(profile.depth_m.values, np.arange(15., 1000., 10.))
    assert profile.initial_march_n2_1e5_s2.notna().all()
    mean_n2 = float(profile.initial_march_n2_1e5_s2.mean())
    series = amoc[(amoc.Model == model) & amoc['Simulation year'].between(1, 280)].sort_values('Simulation year')
    assert np.array_equal(series['Simulation year'], np.arange(1, 281))
    strength = series['Forced AMOC strength (Sv)'].to_numpy()
    early_amoc = float(strength[:10].mean())
    later_amoc = float(strength[-10:].mean())
    smoothed = pd.Series(strength).rolling(11, center=True, min_periods=11).mean().to_numpy()
    normalized_amoc.append(strength / early_amoc)
    record = {'model': model, 'member': members[model], 'initial_n2_1e5_s2': mean_n2,
              'n2_over_obs_mean': mean_n2 / float(obs.mean_n2_1e5_s2),
              'amoc_early_sv': early_amoc, 'amoc_later_sv': later_amoc,
              'amoc_later_over_early': later_amoc / early_amoc,
              'amoc_min_11yr_over_early': float(np.nanmin(smoothed) / early_amoc),
              'amoc_min_11yr_year': int(np.nanargmin(smoothed) + 1),
              'n2_calendar_first': n2_metadata['models_metadata'][model]['initial_first_year'],
              'n2_calendar_last': n2_metadata['models_metadata'][model]['initial_last_year'],
              'n2_baseline_caveat': 'March years6-15, not1-10' if model == 'GFDL-ESM4' else ''}
    for region, prefix in [('Labrador-Irminger', 'li'), ('Nordic Seas', 'nordic')]:
        regional = mld[(mld.model == model) & (mld.region == region)]
        early = regional[regional.simulation_year.between(1, 10)].dropna(subset=['march_mld_mean_crossing_only_m'])
        later = regional[regional.simulation_year.between(271, 280)].dropna(subset=['march_mld_mean_crossing_only_m'])
        assert len(later) == 10 and len(early) == (5 if model == 'GFDL-ESM4' else 10)
        for name, frame in [('early', early), ('later', later)]:
            record[prefix + '_' + name + '_mld_m'] = float(frame.march_mld_mean_crossing_only_m.mean())
            record[prefix + '_' + name + '_deep_fraction'] = float(frame.resolved_mld_ge_1000m_fraction.mean())
            record[prefix + '_' + name + '_n_years'] = len(frame)
        record[prefix + '_mld_ratio'] = record[prefix + '_later_mld_m'] / record[prefix + '_early_mld_m']
        record[prefix + '_deep_area_ratio'] = record[prefix + '_later_deep_fraction'] / record[prefix + '_early_deep_fraction']
        record[prefix + '_both_80pct'] = record[prefix + '_mld_ratio'] >= 0.8 and record[prefix + '_deep_area_ratio'] >= 0.8
    record['group'] = 'high_initial_stratification' if record['n2_over_obs_mean'] > initial_bias_ratio else ('amoc_recovered' if record['amoc_later_over_early'] >= amoc_recovery_threshold else 'amoc_not_recovered')
    record['group_before_user_exclusion'] = record['group']
    record['exclusion_reason'] = 'User-requested exclusion from composite; not an initial-stratification rejection' if model in excluded_models else ''
    if model in excluded_models:
        record['group'] = 'user_excluded'
    records.append(record)
membership = pd.DataFrame(records)
recovered = membership.loc[membership.group == 'amoc_recovered', 'model'].tolist()
not_recovered = membership.loc[membership.group == 'amoc_not_recovered', 'model'].tolist()
eligible = recovered + not_recovered
assert len(recovered) >= 2 and len(not_recovered) >= 2
print(membership[['model', 'n2_over_obs_mean', 'amoc_later_over_early', 'li_mld_ratio', 'li_deep_area_ratio', 'group']].to_string(index=False), flush=True)
sensitivity = []
for n2_limit in [1.25, 1.5, 2.0, 2.25]:
    for recovery_limit in [0.7, 0.8, 0.9]:
        pool = membership[(membership.n2_over_obs_mean <= n2_limit) & ~membership.model.isin(excluded_models)]
        group_r = pool.loc[pool.amoc_later_over_early >= recovery_limit, 'model'].tolist()
        group_n = pool.loc[pool.amoc_later_over_early < recovery_limit, 'model'].tolist()
        sensitivity.append({'n2_ratio_limit': n2_limit, 'amoc_ratio_limit': recovery_limit, 'recovered': ';'.join(group_r), 'not_recovered': ';'.join(group_n), 'n_recovered': len(group_r), 'n_not_recovered': len(group_n), 'at_least_two_each': len(group_r) >= 2 and len(group_n) >= 2})


# %% 4. 复用同一10年窗口的温度、风及同层梯度正式产品

with xr.open_dataset(jet / 'n12_endpoint10_temperature_gradient_jet/temperature_gradient_jet.nc') as source:
    climate = source[['tas', 'u850']].sel(model=models).load()
with xr.open_dataset(jet / 'n23_eight_model_pressure_thermal/pressure_matched_patterns.nc') as source:
    climate['ta850'] = source.ta.sel(model=models, plev=85000, drop=True).load()
    climate['gradient850'] = source.gradient.sel(model=models, plev=85000, drop=True).load()
for variable in climate.data_vars:
    assert climate[variable].dims == ('model', 'period', 'lat', 'lon')
latitude = climate.lat.values
longitude = climate.lon.values
target_lon, target_lat = np.meshgrid(longitude, latitude)
target_xyz = np.column_stack([(np.cos(np.deg2rad(target_lat)) * np.cos(np.deg2rad(target_lon))).ravel(), (np.cos(np.deg2rad(target_lat)) * np.sin(np.deg2rad(target_lon))).ravel(), np.sin(np.deg2rad(target_lat)).ravel()])
land_geometry = shapely.union_all(list(shp.Reader(shp.natural_earth(resolution='110m', category='physical', name='land')).geometries()))
target_ocean = ~shapely.contains_xy(land_geometry, (target_lon + 180) % 360 - 180, target_lat)
ocean_mask = xr.DataArray(target_ocean, dims=('lat', 'lon'), coords={'lat': latitude, 'lon': longitude})
region_geojson = json.loads((ocean / 'mld_coast_regions_20261012.geojson').read_text(encoding='utf-8'))
li_geometry = shapely.from_geojson(json.dumps(next(feature for feature in region_geojson['features'] if feature['properties']['region'] == 'Labrador-Irminger')['geometry']))
li_mask = xr.DataArray(shapely.contains_xy(li_geometry, (target_lon + 180) % 360 - 180, target_lat) & target_ocean, dims=('lat', 'lon'), coords={'lat': latitude, 'lon': longitude})


# %% 5. 重新读取海冰和水文，严格统一为每窗30个DJF日历月

for variable in ['siconc', 'pr', 'evspsbl']:
    model_fields = []
    for model in models:
        phase_fields = []
        native_grids = []
        for index, period in enumerate(periods):
            folder = Path('F:/sea/' + ('cmipseaice' if variable == 'siconc' else 'cmip6' + variable) + ('ad' if index == 0 else ''))
            experiment = '1pctCO2' if index == 0 else '1pctCO2-cdr'
            table = 'SImon' if variable == 'siconc' else 'Amon'
            paths = sorted(folder.rglob(f'{variable}_{table}_{model}_{experiment}_{members[model]}_*.nc'))
            paths = [path for path in paths if 'yearmean' not in path.name.lower() and '(' not in path.name]
            start = starts[model][index]
            paths = [path for path in paths if int(path.stem.split('_')[-1].split('-')[-1][:4]) >= start and int(path.stem.split('_')[-1].split('-')[0][:4]) <= start + 9]
            # 已明确Nor CDR无发布Amon/pr；仅读经过完整日覆盖及守恒核验的衍生月值。
            derived_monthly = variable == 'pr' and model == 'NorESM2-LM' and period == 'later'
            if derived_monthly:
                paths = [nor_monthly_pr]
            assert paths, (variable, model, period, start)
            parts = []
            for path in paths:
                with xr.open_dataset(path, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as source:
                    selected = source.time.dt.year.isin(range(start, start + 10)) & source.time.dt.month.isin([1, 2, 12])
                    if not bool(selected.any()):
                        continue
                    assert source.attrs['source_id'] == model and source.attrs['variant_label'] == members[model]
                    assert source.attrs['experiment_id'] == experiment
                    if derived_monthly:
                        assert source.attrs['is_derived'] == 1 and source.attrs['source_table_id'] == 'day'
                        assert source.attrs['table_id'] == 'derived_monthly_from_day'
                    field = source[variable].isel(time=selected).load()
                    expected_dims = ('time',) + ice_dims[model] if variable == 'siconc' else ('time', 'lat', 'lon')
                    assert field.dims == expected_dims
                    assert field.attrs['units'] == ('%' if variable == 'siconc' else 'kg m-2 s-1')
                    parts.append(field)
                    if variable == 'siconc':
                        native_lat = source[ice_coords[model][0]].values
                        native_lon = source[ice_coords[model][1]].values
                    provenance.append({'model': model, 'variable': variable, 'period': period, 'file': str(path), 'dimensions': list(field.dims), 'months': field.sizes['time'], 'member': members[model], 'source_kind': 'duration_weighted_monthly_from_complete_daily' if derived_monthly else 'published_monthly'})
            monthly = xr.concat(parts, 'time').sortby('time')
            ids = monthly.time.dt.year.values * 12 + monthly.time.dt.month.values
            expected = np.array([year * 12 + month for year in range(start, start + 10) for month in [1, 2, 12]])
            assert np.array_equal(ids, expected) and len(np.unique(ids)) == 30, (variable, model, period, ids)
            field_mean = monthly.weighted(monthly.time.dt.days_in_month).mean('time').where(monthly.notnull().all('time'))
            if variable == 'siconc':
                assert float(field_mean.min()) >= -0.001 and float(field_mean.max()) <= 100.01
                phase_fields.append(field_mean.values)
                native_grids.append((native_lat.copy(), native_lon.copy()))
            else:
                # kg m-2 s-1 = mm s-1；E使用向上为正的原evspsbl，绝不反号。
                field_mean = field_mean * 86400.
                field_mean = field_mean.assign_coords(lon=field_mean.lon % 360).sortby('lon')
                left = field_mean.isel(lon=[-1]).assign_coords(lon=[float(field_mean.lon[-1]) - 360])
                right = field_mean.isel(lon=[0]).assign_coords(lon=[float(field_mean.lon[0]) + 360])
                field_mean = xr.concat([left, field_mean, right], 'lon').interp(lat=latitude, lon=longitude)
                phase_fields.append(field_mean)
            availability.append({'model': model, 'variable': variable, 'period': period, 'complete': True, 'reason': '', 'months': 30, 'calendar_first': start, 'calendar_last': start + 9})
        if variable == 'siconc':
            # Nor部分文件陆地坐标为1e30，另一些为NaN；仅在两期有效海洋坐标核对。
            valid = np.isfinite(phase_fields[0]) & np.isfinite(phase_fields[1]) & np.isfinite(native_lat) & np.isfinite(native_lon) & (abs(native_lat) <= 90) & (abs(native_lon) <= 360)
            assert np.allclose(native_grids[0][0][valid], native_grids[1][0][valid])
            assert np.allclose((native_grids[0][1][valid] - native_grids[1][1][valid] + 180) % 360 - 180, 0.)
            xyz = np.column_stack([np.cos(np.deg2rad(native_lat[valid])) * np.cos(np.deg2rad(native_lon[valid])), np.cos(np.deg2rad(native_lat[valid])) * np.sin(np.deg2rad(native_lon[valid])), np.sin(np.deg2rad(native_lat[valid]))])
            distance, nearest = cKDTree(xyz).query(target_xyz)
            mapped_fields = []
            for values in phase_fields:
                mapped = values[valid][nearest].reshape(target_lat.shape)
                mapped[distance.reshape(target_lat.shape) * radius_m > mapping_cutoff_m] = np.nan
                mapped_fields.append(xr.DataArray(mapped, dims=('lat', 'lon'), coords={'lat': latitude, 'lon': longitude}).where(ocean_mask))
            phase_fields = mapped_fields
        model_field = xr.concat(phase_fields, pd.Index(periods, name='period'))
        model_field = model_field.where(model_field.notnull().all('period'))
        model_fields.append(model_field.expand_dims(model=[model]))
        print(variable, model, 'completed', flush=True)
    climate[variable] = xr.concat(model_fields, 'model').transpose('model', 'period', 'lat', 'lon')
climate['pe'] = climate.pr - climate.evspsbl
units = {'tas': 'K', 'u850': 'm s-1', 'ta850': 'K', 'gradient850': 'K per 1000 km', 'siconc': 'percentage points', 'pr': 'mm day-1', 'evspsbl': 'mm day-1', 'pe': 'mm day-1'}
for variable in units:
    climate[variable].attrs['units'] = units[variable]
delta = climate.sel(period='later', drop=True) - climate.sel(period='early', drop=True)
assert float(abs(delta.pe - (delta.pr - delta.evspsbl)).max()) < 1.e-10


# %% 6. 固定成员等权合成、配对符号及逐一剔除检验

groups = []
for name, names in [('recovered', recovered), ('not_recovered', not_recovered)]:
    groups.append(delta.sel(model=names).mean('model', skipna=False).expand_dims(group=[name]))
group_fields = xr.concat(groups, 'group')
contrast = group_fields.sel(group='not_recovered', drop=True) - group_fields.sel(group='recovered', drop=True)
# 每个变量在四个入组成员共同有效支持上比较；不让有效成员数随格点变化。
common = delta.sel(model=eligible).notnull().all('model')
group_fields = group_fields.where(common)
contrast = contrast.where(common)
loo_fields = []
pair_fields = []
for omitted in eligible:
    group_r = [model for model in recovered if model != omitted]
    group_n = [model for model in not_recovered if model != omitted]
    value = delta.sel(model=group_n).mean('model', skipna=False) - delta.sel(model=group_r).mean('model', skipna=False)
    loo_fields.append(value.where(common).expand_dims(omitted_model=[omitted]))
for model_n in not_recovered:
    for model_r in recovered:
        pair_fields.append((delta.sel(model=model_n, drop=True) - delta.sel(model=model_r, drop=True)).where(common).expand_dims(pair=[model_n + ' minus ' + model_r]))
loo = xr.concat(loo_fields, 'omitted_model')
pairs = xr.concat(pair_fields, 'pair')
pair_agreement = ((np.sign(pairs) == np.sign(contrast)) & (abs(pairs) > 1.e-12)).mean('pair').where(common)
loo_sign_stable = ((np.sign(loo) == np.sign(contrast)) & (abs(loo) > 1.e-12)).all('omitted_model').where(common)
group_availability = []
for variable in units:
    for name, names in [('recovered', recovered), ('not_recovered', not_recovered)]:
        complete = [model for model in names if bool(np.isfinite(delta[variable].sel(model=model)).any())]
        group_availability.append({'variable': variable, 'group': name, 'expected_models': ';'.join(names), 'available_models': ';'.join(complete), 'full_group_available': len(names) == len(complete)})


# %% 7. 固定区域统计及独立恒等式复核

longitude_signed = (climate.lon + 180) % 360 - 180
region_masks = {'LI_coast': li_mask,
                'SPNA_ocean_45_65N_60W_0': ocean_mask & (climate.lat >= 45) & (climate.lat <= 65) & (longitude_signed >= -60) & (longitude_signed <= 0),
                'NA_ocean_30_60N_60W_0': ocean_mask & (climate.lat >= 30) & (climate.lat <= 60) & (longitude_signed >= -60) & (longitude_signed <= 0),
                'Europe_land_35_65N_10W_30E': (~ocean_mask) & (climate.lat >= 35) & (climate.lat <= 65) & (longitude_signed >= -10) & (longitude_signed <= 30)}
# 固定地理分带，在查看新降水结果之前确定；不按异常中心移动边界。
for first, last in [(30, 45), (45, 60)]:
    region_masks[f'NA_ocean_{first}_{last}N_60W_0'] = ocean_mask & (climate.lat >= first) & (climate.lat < last) & (longitude_signed >= -60) & (longitude_signed <= 0)
for first, last in [(35, 50), (50, 65)]:
    region_masks[f'Europe_land_{first}_{last}N_10W_30E'] = (~ocean_mask) & (climate.lat >= first) & (climate.lat < last) & (longitude_signed >= -10) & (longitude_signed <= 30)
regional_records = []
for region, mask in region_masks.items():
    weights = np.cos(np.deg2rad(climate.lat)) * xr.ones_like(climate.lon)
    for variable in units:
        support = common[variable] & mask
        for model in models:
            field = delta[variable].sel(model=model).where(support)
            value = float(field.weighted(weights).mean()) if bool(field.notnull().any()) else np.nan
            initial_field = climate[variable].sel(model=model, period='early').where(mask)
            initial_value = float(initial_field.weighted(weights).mean()) if bool(initial_field.notnull().any()) else np.nan
            regional_records.append({'region': region, 'variable': variable, 'model': model, 'group': membership.set_index('model').loc[model, 'group'], 'early_mean': initial_value, 'delta_on_common_model_support': value, 'unit': units[variable]})
regional = pd.DataFrame(regional_records)
checks = {}
for variable in units:
    array = delta[variable].sel(model=not_recovered).values
    other = delta[variable].sel(model=recovered).values
    manual = np.mean(array, axis=0) - np.mean(other, axis=0)
    expected = contrast[variable].values
    finite = np.isfinite(expected)
    checks[variable] = float(np.max(abs(manual[finite] - expected[finite]))) if finite.any() else None
    assert checks[variable] is None or checks[variable] < 1.e-10
assert bool(contrast.pr.notnull().any()) and bool(contrast.pe.notnull().any())
assert all(record['complete'] for record in availability)
assert all(record['full_group_available'] for record in group_availability)


# %% 8. 保存计算结果、成员与完整限制

output.mkdir(parents=True)
climate.attrs = {'early': '1pctCO2 years1-10', 'later': 'CDR131-140;continuous271-280', 'season': '30calendar DJF months Jan/Feb/Dec per10year window;dayweighted', 'grid': 'Can64Gaussian latitude and128longitude;bilinear atmosphere;nearest native common-valid SIC <=150km and landmask', 'group_rule': 'Exploratory:high initial LI N2 >1.5*modern IAP mean;remaining split by AMOC late/early >=0.8', 'interpretation': 'AMOC endpoint groups, not universal deep-convection recovery groups;no causal identification or proven tipping', 'drift_correction': 'none'}
climate.attrs['user_excluded_models'] = ';'.join(excluded_models)
climate.attrs['selection_caveat'] = 'UKESM omitted at user request after examining five-model results; do not interpret stronger separation as independent confirmation.'
climate.attrs['common_support'] = 'Common valid cells of the four included models; excluded model does not determine support.'
climate.to_netcdf(output / 'endpoint_climatology.nc')
xr.concat([mask.astype('int8').expand_dims(region=[name]) for name, mask in region_masks.items()], dim='region').rename('mask').to_netcdf(output / 'regional_masks.nc')
delta.to_netcdf(output / 'model_changes.nc')
group_fields.to_netcdf(output / 'group_changes.nc')
contrast.to_netcdf(output / 'not_minus_recovered.nc')
loo.to_netcdf(output / 'leave_one_model_out.nc')
pairs.to_netcdf(output / 'pairwise_contrasts.nc')
pair_agreement.to_netcdf(output / 'pairwise_sign_agreement.nc')
loo_sign_stable.to_netcdf(output / 'leave_one_sign_stable.nc')
xr.DataArray(np.stack(normalized_amoc), dims=('model', 'simulation_year'), coords={'model': models, 'simulation_year': np.arange(1, 281)}, name='amoc_over_early').to_netcdf(output / 'amoc_normalized.nc')
membership.to_csv(output / 'membership.csv', index=False)
pd.DataFrame(sensitivity).to_csv(output / 'classification_sensitivity.csv', index=False)
pd.DataFrame(availability).to_csv(output / 'input_availability.csv', index=False)
pd.DataFrame(group_availability).to_csv(output / 'group_availability.csv', index=False)
regional.to_csv(output / 'regional_model_changes.csv', index=False)
observations.to_csv(output / 'observation_reference.csv', index=False)
pd.read_csv(jet / 'n12_endpoint10_temperature_gradient_jet/jet_later10_early10.csv').to_csv(output / 'jet_centroid_changes.csv', index=False)
pd.read_csv(jet / 'n31_low_level_jet_core/core_change_summary.csv').to_csv(output / 'jet_core_changes.csv', index=False)
report = {'method': climate.attrs, 'recovered': recovered, 'not_recovered': not_recovered, 'initial_high_stratification': membership.loc[membership.group == 'high_initial_stratification', 'model'].tolist(), 'observation': 'IAPv4 March2005-2024;mean +/-interannualSD is variability,not measurement uncertainty or confidence of climatology;modern vs idealized initial state', 'initial_screen': 'Provisional 1.5times observed mean,not a statistical observational acceptance test. Sensitivities1.25,1.5,2,2.25. Lower-stratification pool is not called observationally validated.', 'recovery_rule': 'lateAMOC/earlyAMOC >=0.8;not regained fraction of lost strength;not equivalence test. Sensitivities0.7,0.8,0.9.', 'gfdl_warning': 'N2 available March6-15;corrected MLD March6-10 for requested baseline1-10;AMOC/atmosphere/SIC all1-10. Strict ocean initial-window match incomplete;excluding GFDL leaves one recovered model.', 'nor_warning': 'AMOC and LI meanMLD recovery can differ. CDR pr now explicitly reads independently verified daily-derived monthly product;all P/E/PE groups complete2v3.', 'sign_diagnostics': 'Six cross-group model pairs are dependent. All-pair same sign and leave-one sign stability are sensitivity measures,not statistical significance or causal confidence.', 'tipping': 'Weakening and endpoint nonrecovery do not prove a tipping point. No model classified as confirmed tipping.', 'members': members, 'inputs': provenance, 'numpy_composite_error': checks, 'output_dimensions': {name: list(climate[name].dims) for name in climate.data_vars}}
report['user_excluded_models'] = excluded_models
report['selection_caveat'] = climate.attrs['selection_caveat']
report['nor_warning'] = 'AMOC and LI meanMLD recovery can differ. CDR pr reads verified daily-derived monthly product; P/E/PE complete for both two-member groups.'
report['sign_diagnostics'] = 'Four cross-group model pairs are dependent; sign agreement and leave-one stability are membership sensitivities, not statistical significance.'
report['included_model_count'] = len(eligible)
report['cross_group_pair_count'] = len(recovered) * len(not_recovered)
(output / 'qc.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
print(output, flush=True)
