# %% 1. 导入依赖

from pathlib import Path
import json
import numpy as np
import pandas as pd
import xarray as xr


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
base = root / 'output/s05_atmosphere/p05_amoc_recovery_composites'
output = base / 'n16_member_mechanism_comparison'
assert not output.exists(), '正式结果存在，先归档。'
models = ['GFDL-ESM4', 'NorESM2-LM', 'CESM2', 'CNRM-ESM2-1']
groups = {'recovered': models[:2], 'not_recovered': models[2:]}
regions = ['LI_coast', 'SPNA_ocean_45_65N_60W_0', 'NA_ocean_30_60N_60W_0']


# %% 3. 读取现有同窗成员结果，保留两种网格

with xr.open_dataset(base / 'n09_member_jet_sections/member_surface_changes.nc') as source:
    context = source[['tas', 'u850', 'gradient850']].sel(model=models).load()
with xr.open_dataset(base / 'n01_endpoint_groups/model_changes.nc') as source:
    context['siconc'] = source.siconc.sel(model=models).load()
with xr.open_dataset(base / 'n14_moisture_group_contrasts/member_fields.nc') as source:
    moisture = source.sel(model=models).load()
with xr.open_dataset(base / 'n01_endpoint_groups/regional_masks.nc') as source:
    masks = source.load()
membership = pd.read_csv(base / 'n01_endpoint_groups/membership.csv').set_index('model').loc[models]
jet = pd.read_csv(base / 'n09_member_jet_sections/na_jet_centroids.csv').set_index('model').loc[models]
assert all(field.dims == ('model', 'lat', 'lon') for field in context.data_vars.values())
assert all(field.dims == ('model', 'lat', 'lon') for field in moisture.data_vars.values())
# 两种网格不混并、不插值；海冰为浓度百分点，温度梯度为同层-dT/dy。
context = context.assign_coords(lon=((context.lon + 180) % 360) - 180).sortby('lon')
masks = masks.assign_coords(lon=((masks.lon + 180) % 360) - 180).sortby('lon')
context = context.sel(lon=slice(-80, 30), lat=slice(20, 80))
masks = masks.sel(lon=context.lon, lat=context.lat)


# %% 4. 等权组均、组差和组内差，缺失不改变成员数量

products = {}
errors = {}
for name, members in [('context', context), ('moisture', moisture)]:
    means = xr.concat([members.sel(model=group_models).mean('model', skipna=False)
                       for group_models in groups.values()], dim=pd.Index(list(groups), name='group'))
    contrast = means.sel(group='not_recovered', drop=True) - means.sel(group='recovered', drop=True)
    within = xr.concat([members.sel(model=pair[1], drop=True) - members.sel(model=pair[0], drop=True)
                        for pair in groups.values()], dim=pd.Index(['Nor_minus_GFDL', 'CNRM_minus_CESM'], name='pair'))
    maximum_error = 0.
    for variable in members.data_vars:
        manual = (members[variable].sel(model=models[2], drop=True) + members[variable].sel(model=models[3], drop=True)
                  - members[variable].sel(model=models[0], drop=True) - members[variable].sel(model=models[1], drop=True)) / 2
        maximum_error = max(maximum_error, float(abs(manual - contrast[variable]).max()))
    assert maximum_error < 1.e-10
    errors[name] = maximum_error
    products[name] = (members, means, contrast, within)


# %% 5. 共同有效格点区域均值和成员场型一致性

regional_rows = []
for region in regions:
    for variable in context.data_vars:
        valid = masks['mask'].sel(region=region).astype(bool) & np.isfinite(context[variable]).all('model')
        weights = np.cos(np.deg2rad(context.lat)).broadcast_like(valid).where(valid, 0.)
        means = context[variable].where(valid).weighted(weights).mean(('lat', 'lon'))
        for model in models:
            regional_rows.append({'region': region, 'variable': variable, 'model': model,
                                  'change': float(means.sel(model=model)), 'common_valid_cells': int(valid.sum())})
regional = pd.DataFrame(regional_rows)
consistency_rows = []
for family, dataset in [('context', context), ('moisture', moisture)]:
    variables = list(context.data_vars) if family == 'context' else ['delta_pr', 'delta_evspsbl', 'delta_pe', 'mfc_total', 'mfc_dynamic', 'mfc_humidity', 'budget_residual']
    domain = (dataset.lat >= 30) & (dataset.lat <= 70) & (dataset.lon >= -60) & (dataset.lon <= 20)
    for group, pair in groups.items():
        for variable in variables:
            first = dataset[variable].sel(model=pair[0], drop=True)
            second = dataset[variable].sel(model=pair[1], drop=True)
            valid = domain & np.isfinite(first) & np.isfinite(second)
            weights = np.cos(np.deg2rad(dataset.lat)).broadcast_like(first).where(valid, 0.)
            first_mean = first.where(valid).weighted(weights).mean()
            second_mean = second.where(valid).weighted(weights).mean()
            first_centered = first - first_mean
            second_centered = second - second_mean
            denominator = np.sqrt((first_centered**2).where(valid).weighted(weights).mean()
                                  * (second_centered**2).where(valid).weighted(weights).mean())
            correlation = (first_centered * second_centered).where(valid).weighted(weights).mean() / denominator
            active = valid & ((abs(first) > 1.e-12) | (abs(second) > 1.e-12))
            active_weights = weights.where(active, 0.)
            agreement = (np.sign(first) == np.sign(second)).astype(float).where(active).weighted(active_weights).mean()
            consistency_rows.append({'family': family, 'group': group, 'variable': variable,
                                     'weighted_pattern_r': float(correlation), 'same_sign_fraction_active': float(agreement),
                                     'valid_cells': int(valid.sum()), 'domain': '30-70N,60W-20E; land and ocean'})


# %% 6. 保存结果与分组限制

output.mkdir(parents=True)
for name, datasets in products.items():
    for suffix, dataset in zip(['members', 'groups', 'contrast', 'within'], datasets):
        dataset.to_netcdf(output / f'{name}_{suffix}.nc')
regional.to_csv(output / 'context_regional_common_support.csv', index=False)
pd.DataFrame(consistency_rows).to_csv(output / 'member_consistency.csv', index=False)
summary = membership[['amoc_early_sv', 'amoc_later_sv', 'amoc_later_over_early', 'li_mld_ratio', 'li_deep_area_ratio', 'nordic_mld_ratio']].join(jet[['early_jet_latitude', 'later_jet_latitude', 'jet_shift_degree']])
summary.to_csv(output / 'recovery_and_jet.csv')
report = {'models': models, 'groups': groups, 'group_label': 'AMOC endpoint recovery, not complete regional deep-convection recovery',
          'window': 'DJF later CDR131-140 minus early1pct1-10;30 calendar months each,day weighted',
          'context_source': 'n09 member_surface_changes, siconc from n01 model_changes',
          'moisture_source': 'n14; Nor daily-derived monthly precipitation verified',
          'dimensions': {'context': 'model/group/pair,lat,lon on existing coarse grid', 'moisture': 'model/group/pair,lat,lon on existing 1-degree grid'},
          'mapping': 'No new mapping; grids retained separately; strict two-model means skipna=False',
          'qc_max_group_identity_error': errors,
          'limitations': ['No smoothing or infilling', 'Member pattern correlation is descriptive, no significance or causal attribution',
                          'Sea-ice and cooling may share ocean-circulation drivers', 'Monthly moisture residual is retained; not an identified eddy term',
                          'Nor AMOC88.1%,LI meanMLD72.9%,deep-mixing area98.8%; recovery metrics are not interchangeable']}
(output / 'qc.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
print(summary.to_string())
print(pd.DataFrame(consistency_rows).to_string(index=False))
