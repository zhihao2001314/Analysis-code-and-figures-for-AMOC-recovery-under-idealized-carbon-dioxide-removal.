# %% 1. 导入依赖

from pathlib import Path
import json
import numpy as np
import pandas as pd
import xarray as xr
from scipy.signal import find_peaks


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
base = root / 'output/s05_atmosphere/p05_amoc_recovery_composites'
jet_base = root / 'output/s05_atmosphere/p02_ua'
output = base / 'n22_pairwise_jet_budget'
assert not output.exists(), '正式产品存在，先归档。'
groups = {'recovered': ['GFDL-ESM4', 'NorESM2-LM'], 'not_recovered': ['CESM2', 'CNRM-ESM2-1']}
models = groups['recovered'] + groups['not_recovered']
variables = ['siconc', 'tas', 'gradient850', 'u850', 'pr', 'evspsbl', 'pe']
regions = ['LI_coast', 'SPNA_ocean_45_65N_60W_0', 'Europe_land_35_65N_10W_30E']
minimum_speed = 3.
peak_gap = 2.
minimum_peak_separation = 5.


# %% 3. 同窗四成员及四跨组配对，共同有效支持不随配对变化

with xr.open_dataset(base / 'n01_endpoint_groups/model_changes.nc') as source:
    members = source[variables].sel(model=models).load()
with xr.open_dataset(base / 'n01_endpoint_groups/regional_masks.nc') as source:
    masks = source.load()
members = members.assign_coords(lon=(members.lon + 180) % 360 - 180).sortby('lon')
masks = masks.assign_coords(lon=(masks.lon + 180) % 360 - 180).sortby('lon')
common = members.notnull().all('model')
pair_fields = []
for model_n in groups['not_recovered']:
    for model_r in groups['recovered']:
        pair_fields.append((members.sel(model=model_n, drop=True) - members.sel(model=model_r, drop=True)).where(common).expand_dims(pair=[model_n + '_minus_' + model_r]))
pairs = xr.concat(pair_fields, 'pair')
contrast = members.sel(model=groups['not_recovered']).mean('model', skipna=False) - members.sel(model=groups['recovered']).mean('model', skipna=False)
contrast = contrast.where(common)
maximum_identity_error = 0.
for variable in variables:
    maximum_identity_error = max(maximum_identity_error, float(abs(pairs[variable].mean('pair', skipna=False) - contrast[variable]).max()))
assert maximum_identity_error < 1.e-10
agreement = ((np.sign(pairs) == np.sign(contrast)) & (abs(pairs) > 1.e-12) & (abs(contrast) > 1.e-12)).sum('pair').where(common)
regional_rows = []
regional_summary = []
for region in regions:
    for variable in variables:
        mask = masks['mask'].sel(region=region).astype(bool) & common[variable]
        if not bool(mask.any()):
            continue
        weights = np.cos(np.deg2rad(members.lat)).broadcast_like(mask).where(mask, 0.)
        values = pairs[variable].where(mask).weighted(weights).mean(('lat', 'lon'))
        group_value = float(contrast[variable].where(mask).weighted(weights).mean(('lat', 'lon')))
        for pair in pairs.pair.values:
            regional_rows.append({'region': region, 'variable': variable, 'pair': str(pair), 'change_difference': float(values.sel(pair=pair))})
        regional_summary.append({'region': region, 'variable': variable, 'group_difference': group_value,
                                 'pair_min': float(values.min()), 'pair_max': float(values.max()),
                                 'same_sign_count_of4': int(((np.sign(values) == np.sign(group_value)) & (abs(values) > 1.e-12)).sum())})


# %% 4. 区域风带位置、峰值、宽度和东西分区

with xr.open_dataset(jet_base / 'n12_endpoint10_temperature_gradient_jet/temperature_gradient_jet.nc') as source:
    wind = source.u850.sel(model=models).load()
wind = wind.assign_coords(lon=(wind.lon + 180) % 360 - 180).sortby('lon').sel(lat=slice(30, 70), lon=slice(-60, 60))
assert wind.dims == ('model', 'period', 'lat', 'lon')
sectors = {'full_60W_60E': wind.lon >= -60, 'west_60W_0': wind.lon < 0, 'east_0_60E': wind.lon >= 0}
profile_fields = []
shape_rows = []
for sector, mask in sectors.items():
    profile = wind.where(mask, drop=True).mean('lon')
    profile_fields.append(profile.expand_dims(sector=[sector]))
    for model in models:
        for period in ['early', 'later']:
            field = profile.sel(model=model, period=period)
            values = field.values
            latitude = field.lat.values
            assert np.isfinite(values).all()
            centroid = np.sum(values**2 * latitude) / np.sum(values**2)
            width = np.sqrt(np.sum(values**2 * (latitude - centroid)**2) / np.sum(values**2))
            shape_rows.append({'model': model, 'period': period, 'sector': sector, 'centroid_degree': centroid,
                               'peak_grid_latitude': latitude[np.argmax(values)], 'peak_speed_m_s': np.max(values), 'squared_wind_rms_width_degree': width})
profiles = xr.concat(profile_fields, 'sector').rename('u850_profile')
shapes = pd.DataFrame(shape_rows)
shape_changes = []
for (model, sector), frame in shapes.groupby(['model', 'sector']):
    initial = frame.set_index('period').loc['early']
    later = frame.set_index('period').loc['later']
    shape_changes.append({'model': model, 'sector': sector, 'centroid_shift_degree': later.centroid_degree - initial.centroid_degree,
                          'peak_grid_shift_degree': later.peak_grid_latitude - initial.peak_grid_latitude,
                          'peak_speed_change_m_s': later.peak_speed_m_s - initial.peak_speed_m_s,
                          'width_change_degree': later.squared_wind_rms_width_degree - initial.squared_wind_rms_width_degree})


# %% 5. 逐经度最大西风候选与多峰标记，不拟合或强行连接轴线

axis_rows = []
for model in models:
    for period in ['early', 'later']:
        field = wind.sel(model=model, period=period)
        for longitude in field.lon.values:
            profile = field.sel(lon=longitude).values
            latitude = field.lat.values
            valid = np.isfinite(profile)
            assert valid.any()
            maximum_index = int(np.nanargmax(profile))
            maxima, _ = find_peaks(np.where(valid, profile, -np.inf))
            separated = maxima[(abs(latitude[maxima] - latitude[maximum_index]) >= minimum_peak_separation) & (profile[maxima] > 0)]
            second_index = int(separated[np.argmax(profile[separated])]) if len(separated) else None
            gap = float(profile[maximum_index] - profile[second_index]) if second_index is not None else np.inf
            boundary = maximum_index in [0, len(latitude) - 1]
            neighbors_valid = not boundary and valid[maximum_index - 1] and valid[maximum_index + 1]
            status = 'ok'
            if valid.mean() < .8:
                status = 'coverage'
            elif boundary:
                status = 'boundary'
            elif not neighbors_valid:
                status = 'missing_neighbor'
            elif profile[maximum_index] < minimum_speed:
                status = 'weak'
            elif gap < peak_gap:
                status = 'multipeak'
            axis_rows.append({'model': model, 'period': period, 'lon': longitude, 'candidate_latitude': latitude[maximum_index],
                              'candidate_speed': profile[maximum_index], 'second_latitude': latitude[second_index] if second_index is not None else np.nan,
                              'peak_gap_m_s': gap, 'valid_latitude_fraction': valid.mean(), 'interior_neighbors_valid': bool(neighbors_valid),
                              'status': status})
axes = pd.DataFrame(axis_rows)
axis_summary = []
sensitivity = []
for model in models:
    for sector, longitude_mask in [('full_60W_60E', axes.lon >= -60), ('west_60W_0', axes.lon < 0), ('east_0_60E', axes.lon >= 0)]:
        selected = axes[(axes.model == model) & longitude_mask]
        early = selected[selected.period == 'early'].set_index('lon')
        later = selected[selected.period == 'later'].set_index('lon')
        paired = (early.status == 'ok') & (later.status == 'ok')
        shift = later.loc[paired, 'candidate_latitude'] - early.loc[paired, 'candidate_latitude']
        axis_summary.append({'model': model, 'sector': sector, 'total_longitudes': len(early), 'paired_reliable_longitudes': int(paired.sum()),
                             'median_paired_shift_degree': shift.median(), 'mean_paired_shift_degree': shift.mean(),
                             'paired_north_count': int((shift > 0).sum()), 'paired_south_count': int((shift < 0).sum()), 'paired_unchanged_count': int((shift == 0).sum())})
        for speed in [2., 3., 4.]:
            for gap_threshold in [1., 2., 3.]:
                usable = (early.valid_latitude_fraction >= .8) & (later.valid_latitude_fraction >= .8)
                usable &= early.interior_neighbors_valid & later.interior_neighbors_valid
                usable &= (early.candidate_speed >= speed) & (later.candidate_speed >= speed)
                usable &= (early.peak_gap_m_s >= gap_threshold) & (later.peak_gap_m_s >= gap_threshold)
                differences = later.loc[usable, 'candidate_latitude'] - early.loc[usable, 'candidate_latitude']
                sensitivity.append({'model': model, 'sector': sector, 'min_speed': speed, 'min_peak_gap': gap_threshold,
                                    'paired_longitudes': int(usable.sum()), 'mean_shift': differences.mean(), 'median_shift': differences.median()})


# %% 6. 保留所有成员和配对的区域降水收支及残差

with xr.open_dataset(base / 'n14_moisture_group_contrasts/member_fields.nc') as source:
    moisture = source.sel(model=models).load()
budget_products = [moisture.sel(model=model, drop=True).expand_dims(product=[model]) for model in models]
for model_n in groups['not_recovered']:
    for model_r in groups['recovered']:
        budget_products.append((moisture.sel(model=model_n, drop=True) - moisture.sel(model=model_r, drop=True)).expand_dims(product=[model_n + '_minus_' + model_r]))
budget_products.append((moisture.sel(model=groups['not_recovered']).mean('model', skipna=False) - moisture.sel(model=groups['recovered']).mean('model', skipna=False)).expand_dims(product=['group_difference']))
budgets = xr.concat(budget_products, 'product')
budget_rows = []
for region, (west, east, south, north) in {'subpolar_50_65N_allcells': (-60, 0, 50, 65), 'Europe_45_60N_allcells': (-10, 20, 45, 60)}.items():
    mask = (budgets.lon >= west) & (budgets.lon <= east) & (budgets.lat >= south) & (budgets.lat <= north)
    weights = np.cos(np.deg2rad(budgets.lat)).broadcast_like(mask).where(mask, 0.)
    for product in budgets['product'].values:
        fields = budgets.sel(product=product)
        record = {'region': region, 'product': str(product)}
        for variable in ['delta_pr', 'delta_evspsbl', 'delta_pe', 'mfc_total', 'delta_storage', 'budget_residual', 'mfc_dynamic', 'mfc_humidity', 'mfc_other']:
            record[variable] = float(fields[variable].where(mask).weighted(weights).mean(('lat', 'lon')))
        record['negative_storage'] = -record['delta_storage']
        record['closure_error'] = record['delta_pr'] - (record['delta_evspsbl'] + record['mfc_total'] + record['negative_storage'] + record['budget_residual'])
        record['residual_to_pe_rms'] = float(np.sqrt((fields.budget_residual**2).where(mask).weighted(weights).mean()) / np.sqrt((fields.delta_pe**2).where(mask).weighted(weights).mean()))
        assert abs(record['closure_error']) < 1.e-10
        budget_rows.append(record)


# %% 7. 保存计算结果和方法边界

output.mkdir(parents=True)
members.to_netcdf(output / 'members.nc')
pairs.to_netcdf(output / 'pair_fields.nc')
contrast.to_netcdf(output / 'contrast.nc')
agreement.to_netcdf(output / 'pair_sign_count.nc')
wind.to_dataset(name='u850').to_netcdf(output / 'wind_climatology.nc')
profiles.to_dataset().to_netcdf(output / 'sector_profiles.nc')
pd.DataFrame(regional_rows).to_csv(output / 'regional_pairs.csv', index=False)
pd.DataFrame(regional_summary).to_csv(output / 'regional_pair_summary.csv', index=False)
shapes.to_csv(output / 'wind_shape_states.csv', index=False)
pd.DataFrame(shape_changes).to_csv(output / 'wind_shape_changes.csv', index=False)
axes.to_csv(output / 'local_axis_candidates.csv', index=False)
pd.DataFrame(axis_summary).to_csv(output / 'local_axis_summary.csv', index=False)
pd.DataFrame(sensitivity).to_csv(output / 'local_axis_threshold_sensitivity.csv', index=False)
pd.DataFrame(budget_rows).to_csv(output / 'regional_budget.csv', index=False)
report = {'groups': groups, 'window': 'DJF laterCDR131-140 minus early1pct1-10;dayweighted30 calendar months each',
          'sources': ['n01 model_changes/regional_masks', 'p02_ua/n12 regional u850 climatology', 'n14 member moisture fields'],
          'dimensions': {'climate': '(model/pair,lat,lon)', 'wind': '(model,period,lat,lon)', 'sector_profiles': '(sector,model,period,lat)'},
          'pair_mean_identity_max_error': maximum_identity_error,
          'axis_definition': 'Exploratory unsmoothed longitude-wise maximum of seasonal-mean u850 over30-70N; not a validated eddy-driven daily jet axis',
          'axis_flags': 'coverage<80%; search-edge max; missing adjacent latitude; max<3m/s; secondary positive local peak >=5deg away within2m/s of max',
          'axis_sensitivity': 'minimum speed2/3/4 and competing-peak gap1/2/3m/s; no optimization to desired sign',
          'width': 'RMS latitude spread of squared regional-mean wind about centroid, not FWHM',
          'budget_regions': '50-65N/60W-0 and45-60N/10W-20E include land and ocean; not the ocean-only/land-only climate masks',
          'limits': ['UK excluded per user, ACCESS not added', 'No smoothing/infilling of wind maps, no forced continuous axis or subgrid peak interpolation',
                     'Four cross-model pairs not independent; sign agreement not significance', 'Fixed coarse grid~2.79deg limits small peak changes',
                     'Residual retained; arithmetic closure not physical budget closure or eddy identification', 'No change to original centroid definition or model selection']}
(output / 'qc.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
print(pd.DataFrame(regional_summary).to_string(index=False))
print(pd.DataFrame(shape_changes).to_string(index=False))
print(pd.DataFrame(axis_summary).to_string(index=False))
