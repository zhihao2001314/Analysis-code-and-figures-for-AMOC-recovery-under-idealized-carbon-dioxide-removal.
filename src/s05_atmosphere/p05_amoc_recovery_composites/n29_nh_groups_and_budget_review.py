# %% 1. 导入依赖

from pathlib import Path
import json
import numpy as np
import pandas as pd
import xarray as xr


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
base = root / 'output/s05_atmosphere/p05_amoc_recovery_composites'
output = base / 'n29_nh_groups_and_budget_review'
assert not output.exists(), '正式结果需先归档。'
output.mkdir(parents=True)
models = ['GFDL-ESM4', 'NorESM2-LM', 'CESM2', 'CNRM-ESM2-1']
groups = {'recovered': models[:2], 'not_recovered': models[2:]}
climates = []


# %% 3. 由每窗十个年度样本重建北半球气候态及组差

for model in models:
    with xr.open_dataset(base / 'n28_nh_endpoint_samples' / (model + '.nc')) as source:
        data = source.load()
    climate = data[['pr', 'evspsbl', 'pe', 'u850']].weighted(data.day_weight).mean('year_in_window', skipna=False)
    climates.append(climate.expand_dims(model=[model]))
climates = xr.concat(climates, dim='model')
changes = climates.sel(period='later', drop=True) - climates.sel(period='early', drop=True)
group_fields = xr.concat([changes.sel(model=names).mean('model', skipna=False).expand_dims(group=[name]) for name, names in groups.items()], dim='group')
contrast = group_fields.sel(group='not_recovered', drop=True) - group_fields.sel(group='recovered', drop=True)
assert float(abs(changes.pe - changes.pr + changes.evspsbl).max()) < 1e-10
climates.to_netcdf(output / 'nh_climatology.nc')
changes.to_netcdf(output / 'nh_members.nc')
group_fields.to_netcdf(output / 'nh_groups.nc')
contrast.to_netcdf(output / 'nh_contrast.nc')


# %% 4. 整理北大西洋收支，显式保留储存与闭合残差

with xr.open_dataset(base / 'n14_moisture_group_contrasts/group_fields.nc') as source:
    moisture_groups = source.load()
with xr.open_dataset(base / 'n14_moisture_group_contrasts/not_minus_recovered.nc') as source:
    moisture_contrast = source.load()
review = xr.concat([moisture_groups.rename({'group': 'product'}), moisture_contrast.expand_dims(product=['not_minus_recovered'])], dim='product')
review = review.assign_coords(product=np.array(['recovered', 'not_recovered', 'not_minus_recovered'], dtype=object))
review['product'].encoding = {}
review['negative_storage'] = -review.delta_storage
review['mass_and_monthly_covariance'] = review.mfc_pressure + review.mfc_monthly_covariance
budget_error = float(abs(review.delta_pe - review.resolved_pe - review.budget_residual).max())
split_error = float(abs(review.mfc_total - review.mfc_dynamic - review.mfc_humidity - review.mfc_nonlinear - review.mass_and_monthly_covariance).max())
assert max(budget_error, split_error) < 1e-9
review.to_netcdf(output / 'na_moisture_review.nc')
with xr.open_dataset(output / 'na_moisture_review.nc') as saved:
    assert saved['product'].values.tolist() == ['recovered', 'not_recovered', 'not_minus_recovered']
rows = []
for product in review.product.values:
    field = review.sel(product=product, lat=slice(30, 60), lon=slice(-60, 0))
    weights = np.cos(np.deg2rad(field.lat)) * xr.ones_like(field.lon)
    row = {'product': str(product)}
    for variable in ['delta_pr', 'delta_evspsbl', 'delta_pe', 'mfc_total', 'negative_storage', 'budget_residual']:
        row[variable + '_mean'] = float(field[variable].weighted(weights).mean())
        row[variable + '_rms'] = float(np.sqrt((field[variable] ** 2).weighted(weights).mean()))
    row['residual_to_pe_rms'] = row['budget_residual_rms'] / row['delta_pe_rms']
    rows.append(row)
pd.DataFrame(rows).to_csv(output / 'closure_review.csv', index=False)


# %% 5. 保存方法、覆盖和缺失资料说明

metadata = {'groups': groups, 'period': 'DJF; CDR131-140 minus1pctCO2 1-10; calendar-year Jan/Feb/Dec,30monthly samples per endpoint, dayweighted.',
    'domain': {'latitude_range': [float(changes.lat.min()), float(changes.lat.max())], 'longitude': 'Full cyclic longitude0-360; no polar/equatorial extrapolation.'},
    'dimensions': {'nh_climatology': ['model','period','lat','lon'], 'nh_members': ['model','lat','lon'], 'nh_groups': ['group','lat','lon'], 'moisture_review': ['product','lat','lon']},
    'nh_grid': 'Original common Gaussian grid retained north of20.93N; positive CNRM precipitation latitudes added below it; cyclic longitude interpolation.',
    'moisture_domain': 'Existing25-78N,78W-28E products only. Moisture is not extended to Northern Hemisphere in this product.',
    'changes': 'North-hemisphere wind and hydrology maps; review existing exact flux-form decomposition, explicitly separate nonlinear and mass/monthly-covariance components.',
    'missing_data': 'No wap_*.nc found in F:/sea filename audit on2026-10-19. Horizontal ua/va, hus,ps do not substitute for wap in direct vertical-advection calculation.',
    'image_equation': 'Five-term expression is a linearized, quasi-steady advection diagnostic; storage, nonlinear, boundary and unresolved submonthly covariance must be retained outside the five terms.',
    'significance': 'No hemispheric significance test performed. Previous North Atlantic BH masks are not transplanted to this enlarged domain.',
    'constraints': 'Four fixed models, UK excluded; no drift correction, no smoothing, no new group selection or jet-centroid extrapolation.',
    'budget_identity_max_error': budget_error, 'decomposition_identity_max_error': split_error}
(output / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(metadata['domain'], indent=2))
print(pd.DataFrame(rows)[['product','residual_to_pe_rms']].to_string(index=False))
