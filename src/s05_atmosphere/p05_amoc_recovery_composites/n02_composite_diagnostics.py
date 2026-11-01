# %% 1. 导入依赖

from pathlib import Path
import json
import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import pearsonr, spearmanr


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
base = root / 'output/s05_atmosphere/p05_amoc_recovery_composites/n01_endpoint_groups'
output = base.parent / 'n02_composite_diagnostics'
assert not output.exists()
membership = pd.read_csv(base / 'membership.csv')
recovered = membership.loc[membership.group == 'amoc_recovered', 'model'].tolist()
not_recovered = membership.loc[membership.group == 'amoc_not_recovered', 'model'].tolist()
eligible = recovered + not_recovered
regional = pd.read_csv(base / 'regional_model_changes.csv')
with xr.open_dataset(base / 'not_minus_recovered.nc') as source:
    contrast = source.load()
with xr.open_dataset(base / 'model_changes.nc') as source:
    delta = source.load()
with xr.open_dataset(base / 'endpoint_climatology.nc') as source:
    climate = source.load()
with xr.open_dataset(base / 'leave_one_model_out.nc') as source:
    loo = source.load()
with xr.open_dataset(base / 'pairwise_sign_agreement.nc') as source:
    agreement = source.load()
with xr.open_dataset(base / 'leave_one_sign_stable.nc') as source:
    sign_stable = source.load()


# %% 3. 区域组间差异及逐一剔除范围，不把模式当独立随机试验

records = []
for region in regional.region.unique():
    for variable in regional.variable.unique():
        selected = regional[(regional.region == region) & (regional.variable == variable)].set_index('model')
        values_r = selected.loc[recovered, 'delta_on_common_model_support'].to_numpy()
        values_n = selected.loc[not_recovered, 'delta_on_common_model_support'].to_numpy()
        complete = np.isfinite(values_r).all() and np.isfinite(values_n).all()
        if not complete:
            records.append({'region': region, 'variable': variable, 'complete': False})
            continue
        difference = float(values_n.mean() - values_r.mean())
        pair_differences = (values_n[:, None] - values_r[None, :]).ravel()
        leave_one = []
        for omitted in eligible:
            group_r = [model for model in recovered if model != omitted]
            group_n = [model for model in not_recovered if model != omitted]
            leave_one.append(float(selected.loc[group_n, 'delta_on_common_model_support'].mean() - selected.loc[group_r, 'delta_on_common_model_support'].mean()))
        records.append({'region': region, 'variable': variable, 'complete': True,
                        'recovered_mean': float(values_r.mean()), 'not_recovered_mean': float(values_n.mean()), 'not_minus_recovered': difference,
                        'pair_min': float(pair_differences.min()), 'pair_max': float(pair_differences.max()),
                        'pair_same_sign_count': int(np.sum(np.sign(pair_differences) == np.sign(difference))), 'pair_count': len(pair_differences),
                        'leave_one_min': min(leave_one), 'leave_one_max': max(leave_one), 'leave_one_all_same_sign': all(np.sign(value) == np.sign(difference) for value in leave_one),
                        'unit': selected.unit.iloc[0], 'n_recovered': len(recovered), 'n_not_recovered': len(not_recovered)})
regional_summary = pd.DataFrame(records)


# %% 4. 两种急流指标保留各模式，分别计算组间差异

centroid = pd.read_csv(base / 'jet_centroid_changes.csv')
core = pd.read_csv(base / 'jet_core_changes.csv')
core = core[(core.kind == '850hPa') & (core.search_range == '15_75N')]
jet_records = []
jet_leave_one = []
for definition, table, column in [('squared_wind_centroid_30_70N', centroid, 'jet_shift_degree'), ('850_westerly_core_15_75N', core, 'climatology_core_shift_degree')]:
    for region in table.region.unique():
        selected = table[table.region == region].set_index('model')
        values_r = selected.loc[recovered, column]
        values_n = selected.loc[not_recovered, column]
        difference = float(values_n.mean() - values_r.mean())
        leave_one = []
        for omitted in eligible:
            names_r = [model for model in recovered if model != omitted]
            names_n = [model for model in not_recovered if model != omitted]
            leave_one.append(float(selected.loc[names_n, column].mean() - selected.loc[names_r, column].mean()))
            jet_leave_one.append({'definition': definition, 'region': region, 'omitted_model': omitted, 'not_minus_recovered_degree': leave_one[-1], 'n_recovered': len(names_r), 'n_not_recovered': len(names_n)})
        jet_records.append({'definition': definition, 'region': region, 'recovered_mean_degree': values_r.mean(), 'not_recovered_mean_degree': values_n.mean(), 'not_minus_recovered_degree': difference, 'leave_one_min_degree': min(leave_one), 'leave_one_max_degree': max(leave_one), 'leave_one_all_same_sign': all(np.sign(value) == np.sign(difference) for value in leave_one)})
jet_summary = pd.DataFrame(jet_records)

# 风场先在入组模式共同有效格点合成，再作海盆经度平均；不把风场均值的质心冒充成员急流纬度均值。
common_wind = climate.u850.sel(model=eligible).notnull().all('model').all('period')
group_wind_fields = []
for group, names in [('recovered', recovered), ('not_recovered', not_recovered)]:
    group_wind_fields.append(climate.u850.sel(model=names).where(common_wind).mean('model', skipna=False).expand_dims(group=[group]))
group_wind_fields = xr.concat(group_wind_fields, 'group')
regional_wind_fields = []
for region in ['North Atlantic-Europe', 'North Pacific']:
    longitude_mask = ((climate.lon >= 300) | (climate.lon <= 60)) if region == 'North Atlantic-Europe' else ((climate.lon >= 140) & (climate.lon <= 240))
    regional_wind_fields.append(group_wind_fields.where(longitude_mask, drop=True).mean('lon', skipna=True).expand_dims(region=[region]))
regional_wind = xr.concat(regional_wind_fields, 'region').transpose('region', 'group', 'period', 'lat')
regional_wind_change = regional_wind.sel(period='later', drop=True) - regional_wind.sel(period='early', drop=True)
regional_wind_contrast = regional_wind_change.sel(group='not_recovered', drop=True) - regional_wind_change.sel(group='recovered', drop=True)


# %% 5. 空间敏感性及分组阈值，不重新选择最优分组

longitude_signed = (delta.lon + 180) % 360 - 180
map_domain = (delta.lat >= 25) & (delta.lat <= 80) & (longitude_signed >= -80) & (longitude_signed <= 30)
weights = np.cos(np.deg2rad(delta.lat)) * xr.ones_like(delta.lon)
map_summary = []
for variable in contrast.data_vars:
    valid = contrast[variable].notnull() & map_domain
    area_weight = weights.where(valid)
    total = float(area_weight.sum())
    if total == 0:
        continue
    same_pairs = float(area_weight.where(agreement[variable] == 1.).sum() / total)
    same_loo = float(area_weight.where(sign_stable[variable] == 1.).sum() / total)
    map_summary.append({'variable': variable, 'all_cross_group_pairs_sign_agree_area_fraction': same_pairs, 'all_leave_one_sign_agree_area_fraction': same_loo,
                        'contrast_min': float(contrast[variable].where(valid).min()), 'contrast_max': float(contrast[variable].where(valid).max())})
threshold_records = []
climate_threshold_records = []
sensitivity = pd.read_csv(base / 'classification_sensitivity.csv')
for _, setting in sensitivity.iterrows():
    names_r = setting.recovered.split(';') if isinstance(setting.recovered, str) else []
    names_n = setting.not_recovered.split(';') if isinstance(setting.not_recovered, str) else []
    if len(names_r) == 0 or len(names_n) == 0:
        continue
    for region in ['North Atlantic-Europe', 'North Pacific']:
        selected = centroid[centroid.region == region].set_index('model')
        threshold_records.append({'n2_ratio_limit': setting.n2_ratio_limit, 'amoc_ratio_limit': setting.amoc_ratio_limit, 'region': region, 'n_recovered': len(names_r), 'n_not_recovered': len(names_n), 'centroid_contrast_degree': float(selected.loc[names_n, 'jet_shift_degree'].mean() - selected.loc[names_r, 'jet_shift_degree'].mean())})
    for region in regional.region.unique():
        for variable in regional.variable.unique():
            selected = regional[(regional.region == region) & (regional.variable == variable)].set_index('model')
            values_r = selected.loc[names_r, 'delta_on_common_model_support'].to_numpy()
            values_n = selected.loc[names_n, 'delta_on_common_model_support'].to_numpy()
            complete = np.isfinite(values_r).all() and np.isfinite(values_n).all()
            value = float(values_n.mean() - values_r.mean()) if complete else np.nan
            climate_threshold_records.append({'n2_ratio_limit': setting.n2_ratio_limit, 'amoc_ratio_limit': setting.amoc_ratio_limit, 'region': region, 'variable': variable, 'n_recovered': len(names_r), 'n_not_recovered': len(names_n), 'at_least_two_each': len(names_r) >= 2 and len(names_n) >= 2, 'not_minus_recovered': value})
ice_initial = regional[(regional.region == 'LI_coast') & (regional.variable == 'siconc')][['model', 'early_mean']].rename(columns={'early_mean': 'early_li_djf_sic_percent'})
ice_check = membership[['model', 'initial_n2_1e5_s2', 'group']].merge(ice_initial, on='model', validate='one_to_one')
ice_pearson = pearsonr(ice_check.initial_n2_1e5_s2, ice_check.early_li_djf_sic_percent)
ice_spearman = spearmanr(ice_check.initial_n2_1e5_s2, ice_check.early_li_djf_sic_percent)


# %% 6. 独立核对全部逐格差值及线性可交换性

audit = []
for variable in climate.data_vars:
    values = climate[variable].transpose('model', 'period', 'lat', 'lon').values
    expected_delta = values[:, 1] - values[:, 0]
    saved_delta = delta[variable].transpose('model', 'lat', 'lon').values
    finite = np.isfinite(saved_delta)
    delta_error = float(np.max(abs(expected_delta[finite] - saved_delta[finite]))) if finite.any() else 0.
    assert np.array_equal(np.isfinite(expected_delta), np.isfinite(saved_delta))
    r_indices = [list(climate.model.values).index(model) for model in recovered]
    n_indices = [list(climate.model.values).index(model) for model in not_recovered]
    difference_of_changes = np.mean(values[n_indices, 1] - values[n_indices, 0], axis=0) - np.mean(values[r_indices, 1] - values[r_indices, 0], axis=0)
    change_of_group_differences = (np.mean(values[n_indices, 1], axis=0) - np.mean(values[r_indices, 1], axis=0)) - (np.mean(values[n_indices, 0], axis=0) - np.mean(values[r_indices, 0], axis=0))
    valid = np.isfinite(contrast[variable].values)
    error = float(np.max(abs(difference_of_changes[valid] - contrast[variable].values[valid]))) if valid.any() else 0.
    linear_error = float(np.max(abs(difference_of_changes[valid] - change_of_group_differences[valid]))) if valid.any() else 0.
    assert max(delta_error, error, linear_error) < 1.e-10
    audit.append({'variable': variable, 'delta_error': delta_error, 'group_error': error, 'linearity_error': linear_error})
output.mkdir(parents=True)
regional_summary.to_csv(output / 'regional_group_summary.csv', index=False)
jet_summary.to_csv(output / 'jet_group_summary.csv', index=False)
pd.DataFrame(jet_leave_one).to_csv(output / 'jet_leave_one_model_out.csv', index=False)
xr.Dataset({'wind_climatology': regional_wind, 'wind_change': regional_wind_change, 'wind_contrast': regional_wind_contrast}).to_netcdf(output / 'regional_wind_profiles.nc')
pd.DataFrame(map_summary).to_csv(output / 'spatial_sign_sensitivity.csv', index=False)
pd.DataFrame(threshold_records).to_csv(output / 'jet_threshold_sensitivity.csv', index=False)
pd.DataFrame(climate_threshold_records).to_csv(output / 'climate_threshold_sensitivity.csv', index=False)
ice_check.to_csv(output / 'initial_ice_stratification_check.csv', index=False)
report = {'independent_arithmetic': audit, 'initial_ice_pearson_r': float(ice_pearson.statistic), 'initial_ice_spearman_rho': float(ice_spearman.statistic), 'initial_ice_interpretation': 'Eight-model DJF LI SIC vs March10-1000m N2 association;not a bias assessment without observation SIC,not same-season causal test;coarse mapped SIC not native-area integral.', 'uncertainty': 'Leave-one ranges are membership sensitivity,not confidence intervals;pairwise contrasts dependent;no significance testing or causal inference.', 'plotting': 'General exploratory figure,not journal-submission certification;all model data preserved.'}
(output / 'independent_qc.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
print(regional_summary[regional_summary.region == 'SPNA_ocean_45_65N_60W_0'].to_string(index=False), flush=True)
print(jet_summary.to_string(index=False), flush=True)
print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
