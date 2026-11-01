# %% 1. 导入依赖

from pathlib import Path
import json
import numpy as np
import pandas as pd
import xarray as xr


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
source = root / 'output/s05_atmosphere/p05_amoc_recovery_composites/n13_endpoint10_moisture_budget'
output = root / 'output/s05_atmosphere/p05_amoc_recovery_composites/n14_moisture_group_contrasts'
assert not output.exists()
groups = {'recovered': ['GFDL-ESM4', 'NorESM2-LM'], 'not_recovered': ['CESM2', 'CNRM-ESM2-1']}
models = groups['recovered'] + groups['not_recovered']
# 所有原生输出内部均有插值支撑的固定共同域，不在地图边缘外推。
latitude = np.arange(25., 79., 1.)
longitude = np.arange(-78., 29., 1.)
regions = {'NA_30_60N_60W_0': (30, 60, -60, 0), 'subpolar_50_65N_60W_0': (50, 65, -60, 0), 'Europe_45_60N_10W_20E': (45, 60, -10, 20)}
scientific_variables = ['delta_pr', 'delta_evspsbl', 'delta_pe', 'mfc_total', 'mfc_dynamic', 'mfc_humidity', 'mfc_other', 'mfc_nonlinear', 'mfc_pressure', 'mfc_monthly_covariance', 'delta_storage', 'resolved_pe', 'budget_residual']
member_fields = []
audit = []
surface_comparison = []
with xr.open_dataset(root / 'output/s05_atmosphere/p05_amoc_recovery_composites/n09_member_jet_sections/member_surface_changes.nc') as dataset:
    old_surface = dataset[['pr', 'evspsbl', 'pe']].sel(model=models).load()
old_surface = old_surface.assign_coords(lon=(old_surface.lon + 180) % 360 - 180).sortby('lon').sel(lat=slice(25, 78), lon=slice(-78, 28))


# %% 3. 先核对原生网格恒等式，再统一到同一显示分析网格

for model in models:
    with xr.open_dataset(source / f'{model}.nc') as dataset:
        fields = dataset.load()
    budget_error = float(abs(fields.delta_pe - fields.mfc_total + fields.delta_storage - fields.budget_residual).max())
    decomposition_error = float(abs(fields.mfc_total - fields.mfc_dynamic - fields.mfc_humidity - fields.mfc_other).max())
    assert max(budget_error, decomposition_error) < 1e-9
    # 直接在旧图的原有共同格点比较，不把两次不同网格插值误当成计算误差。
    surface_on_old_grid = fields.interp(lat=old_surface.lat, lon=old_surface.lon)
    for new_name, old_name in [('delta_pr', 'pr'), ('delta_evspsbl', 'evspsbl'), ('delta_pe', 'pe')]:
        maximum = float(abs(surface_on_old_grid[new_name] - old_surface[old_name].sel(model=model)).max())
        surface_comparison.append({'model': model, 'variable': new_name, 'max_difference_from_latest_surface_map_mm_day': maximum})
        assert maximum < 1e-6, (model, new_name, maximum)
    fields = fields.interp(lat=latitude, lon=longitude)
    assert bool(fields.to_array().notnull().all()), (model, '共同区域不可静默填补缺测。')
    member_fields.append(fields.expand_dims(model=[model]))
    audit.append({'model': model, 'native_budget_identity_error': budget_error, 'native_decomposition_error': decomposition_error})
members = xr.concat(member_fields, dim='model')


# %% 4. 两组等权平均与组差，保留四个成员和跨组配对

group_fields = []
for group, names in groups.items():
    group_fields.append(members.sel(model=names).mean('model', skipna=False).expand_dims(group=[group]))
group_fields = xr.concat(group_fields, dim='group')
contrast = group_fields.sel(group='not_recovered', drop=True) - group_fields.sel(group='recovered', drop=True)
pairs = []
for not_recovered in groups['not_recovered']:
    for recovered in groups['recovered']:
        pairs.append((members.sel(model=not_recovered, drop=True) - members.sel(model=recovered, drop=True)).expand_dims(pair=[not_recovered + ' minus ' + recovered]))
pairs = xr.concat(pairs, dim='pair')
sign_agreement = (np.sign(pairs[scientific_variables]) == np.sign(contrast[scientific_variables])).mean('pair')
within_group = []
for group, names in groups.items():
    for variable in scientific_variables:
        region = members[variable].sel(lat=slice(30, 70), lon=slice(-60, 20))
        first = region.sel(model=names[0]).values.ravel()
        second = region.sel(model=names[1]).values.ravel()
        weights = np.broadcast_to(np.cos(np.deg2rad(region.lat.values))[:, None], region.shape[1:]).ravel()
        first_anomaly = first - np.average(first, weights=weights)
        second_anomaly = second - np.average(second, weights=weights)
        correlation = np.sum(weights * first_anomaly * second_anomaly) / np.sqrt(np.sum(weights * first_anomaly**2) * np.sum(weights * second_anomaly**2))
        within_group.append({'group': group, 'variable': variable, 'pattern_r_30_70N_60W_20E': float(correlation)})


# %% 5. 区域均值、空间RMS、残差比例及旧图水文一致性

regional = []
products = {model: members.sel(model=model, drop=True) for model in models}
products.update({group: group_fields.sel(group=group, drop=True) for group in groups})
products['not_minus_recovered'] = contrast
for label, dataset in products.items():
    for region, bounds in regions.items():
        selected = dataset.sel(lat=slice(bounds[0], bounds[1]), lon=slice(bounds[2], bounds[3]))
        weights = np.cos(np.deg2rad(selected.lat)) * xr.ones_like(selected.lon)
        record = {'product': label, 'region': region}
        for variable in scientific_variables:
            record[variable + '_mean'] = float(selected[variable].weighted(weights).mean())
            record[variable + '_rms'] = float(np.sqrt((selected[variable]**2).weighted(weights).mean()))
        record['residual_to_pe_rms'] = record['budget_residual_rms'] / record['delta_pe_rms']
        pe_anomaly = selected.delta_pe - selected.delta_pe.weighted(weights).mean()
        mfc_anomaly = selected.mfc_total - selected.mfc_total.weighted(weights).mean()
        record['mfc_pe_pattern_r'] = float((pe_anomaly*mfc_anomaly).weighted(weights).mean() / np.sqrt((pe_anomaly**2).weighted(weights).mean() * (mfc_anomaly**2).weighted(weights).mean()))
        regional.append(record)


# %% 6. 保存正式组差及可追溯诊断

output.mkdir(parents=True)
members.to_netcdf(output / 'member_fields.nc')
group_fields.to_netcdf(output / 'group_fields.nc')
contrast.to_netcdf(output / 'not_minus_recovered.nc')
pairs.to_netcdf(output / 'pair_contrasts.nc')
sign_agreement.to_netcdf(output / 'pair_sign_agreement.nc')
pd.DataFrame(regional).to_csv(output / 'regional_summary.csv', index=False)
pd.DataFrame(within_group).to_csv(output / 'member_pattern_consistency.csv', index=False)
metadata = json.loads((source / 'metadata.json').read_text(encoding='utf-8'))
report = {'groups': groups, 'window': metadata['window'], 'input': str(source), 'dimensions': {'members': ['model', 'lat', 'lon'], 'groups': ['group', 'lat', 'lon'], 'contrast': ['lat', 'lon']}, 'estimator': 'Each model later10-early10, then equal model group mean; contrast = not recovered minus recovered.', 'horizontal_method': 'Compute column integrals and spherical divergence on each native grid, then bilinear interpolate scalar/vector products onto 1-degree common grid. Interpolation is not conservative; no smoothing or infilling.', 'vector_interpretation': 'Anomaly of column-integrated monthly-mean q times monthly-mean wind, kg m-1 s-1. Not complete submonthly covariance flux.', 'pair_sign_agreement': 'Four cross-group differences, descriptive consistency not statistical significance.', 'audits': audit, 'surface_vs_latest': surface_comparison, 'limitations': metadata['limitations']}
(output / 'qc.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
print(pd.DataFrame(regional)[pd.DataFrame(regional).region == 'NA_30_60N_60W_0'][['product', 'delta_pr_mean', 'delta_evspsbl_mean', 'delta_pe_rms', 'budget_residual_rms', 'residual_to_pe_rms']].to_string(index=False), flush=True)
