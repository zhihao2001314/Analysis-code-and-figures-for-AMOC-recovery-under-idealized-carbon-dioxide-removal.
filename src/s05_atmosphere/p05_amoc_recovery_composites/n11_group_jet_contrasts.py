# %% 1. 导入依赖

from pathlib import Path
import json
import numpy as np
import pandas as pd
import xarray as xr


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
base = root / 'output/s05_atmosphere/p05_amoc_recovery_composites'
output = base / 'n11_group_jet_contrasts'
assert not output.exists(), '正式产品存在，需先归档。'
groups = {'recovered': ['GFDL-ESM4', 'NorESM2-LM'], 'not_recovered': ['CESM2', 'CNRM-ESM2-1']}
models = groups['recovered'] + groups['not_recovered']


# %% 3. 读取四成员正式产品，保留既有窗口与维度

with xr.open_dataset(base / 'n09_member_jet_sections/na_member_sections.nc') as source:
    section_climate = source[['wind_profile', 'gradient_profile']].sel(model=models).load()
    longitude_fraction = source.valid_longitude_fraction.sel(model=models).load()
with xr.open_dataset(base / 'n09_member_jet_sections/pressure_change_maps.nc') as source:
    pressure_changes = source.sel(model=models).load()
with xr.open_dataset(base / 'n09_member_jet_sections/member_surface_changes.nc') as source:
    surface_changes = source.sel(model=models).load()
section_changes = section_climate.sel(period='later', drop=True) - section_climate.sel(period='early', drop=True)
assert section_climate.wind_profile.dims == ('model', 'period', 'plev', 'lat')
assert pressure_changes.delta_u.dims == ('model', 'plev', 'lat', 'lon')
assert surface_changes.u850.dims == ('model', 'lat', 'lon')
jet = pd.read_csv(base / 'n09_member_jet_sections/na_jet_centroids.csv').set_index('model')
metadata = json.loads((base / 'n09_member_jet_sections/qc.json').read_text(encoding='utf-8'))


# %% 4. 四成员共同有效支持上等权合成，再作未恢复减恢复

products = {'sections': section_changes, 'pressure_maps': pressure_changes, 'surface_maps': surface_changes}
saved = {}
audit = []
for kind, changes in products.items():
    common = changes.notnull().all('model')
    # 两个剖面变量使用相同的四成员共同剖面支持。
    if kind == 'sections':
        joint = common.wind_profile & common.gradient_profile
        common = xr.Dataset({name: joint for name in changes.data_vars})
    group_fields = []
    for group, members in groups.items():
        group_fields.append(changes.sel(model=members).mean('model', skipna=False).where(common).expand_dims(group=[group]))
    group_fields = xr.concat(group_fields, 'group')
    contrast = group_fields.sel(group='not_recovered', drop=True) - group_fields.sel(group='recovered', drop=True)
    pairs = []
    for model_n in groups['not_recovered']:
        for model_r in groups['recovered']:
            pair = changes.sel(model=model_n, drop=True) - changes.sel(model=model_r, drop=True)
            pairs.append(pair.where(common).expand_dims(pair=[model_n + ' minus ' + model_r]))
    pairs = xr.concat(pairs, 'pair')
    agreement = (np.sign(pairs) == np.sign(contrast)).mean('pair').where(common)
    leave_one = []
    for omitted in models:
        recovered = [model for model in groups['recovered'] if model != omitted]
        not_recovered = [model for model in groups['not_recovered'] if model != omitted]
        value = changes.sel(model=not_recovered).mean('model', skipna=False) - changes.sel(model=recovered).mean('model', skipna=False)
        leave_one.append(value.where(common).expand_dims(omitted_model=[omitted]))
    leave_one = xr.concat(leave_one, 'omitted_model')
    for variable in changes.data_vars:
        first = changes[variable].sel(model=groups['not_recovered']).values
        second = changes[variable].sel(model=groups['recovered']).values
        manual = (first[0] + first[1]) / 2. - (second[0] + second[1]) / 2.
        actual = contrast[variable].values
        finite = np.isfinite(actual)
        error = float(np.max(abs(manual[finite] - actual[finite])))
        assert error < 1.e-10
        audit.append({'product': kind, 'variable': variable, 'numpy_group_contrast_error': error, 'finite_cell_count': int(finite.sum())})
    saved[kind + '_groups'] = group_fields
    saved[kind + '_contrast'] = contrast
    saved[kind + '_pair_sign_agreement'] = agreement
    saved[kind + '_leave_one_out'] = leave_one
    if kind == 'sections':
        climates = []
        for group, members in groups.items():
            climates.append(section_climate.sel(model=members).mean('model', skipna=False).where(common).expand_dims(group=[group]))
        climates = xr.concat(climates, 'group')
        xr.testing.assert_allclose(climates.sel(period='later', drop=True) - climates.sel(period='early', drop=True), group_fields)
        saved['sections_group_climatology'] = climates


# %% 5. 核对既有表面组差与收支恒等式，保存指标成员平均

with xr.open_dataset(base / 'n01_endpoint_groups/not_minus_recovered.nc') as source:
    for variable in ['u850', 'pr', 'evspsbl', 'pe', 'tas', 'gradient850']:
        xr.testing.assert_allclose(saved['surface_maps_contrast'][variable], source[variable])
pe_error = float(abs(saved['surface_maps_contrast'].pe - saved['surface_maps_contrast'].pr + saved['surface_maps_contrast'].evspsbl).max())
assert pe_error < 1.e-10
jet_records = []
for group, members in groups.items():
    selected = jet.loc[members]
    jet_records.append({'group': group, 'early_centroid_degree': selected.early_jet_latitude.mean(), 'later_centroid_degree': selected.later_jet_latitude.mean(), 'mean_member_shift_degree': selected.jet_shift_degree.mean(), 'n_models': len(members)})
jet_summary = pd.DataFrame(jet_records)
shift_contrast = float(jet_summary.set_index('group').loc['not_recovered', 'mean_member_shift_degree'] - jet_summary.set_index('group').loc['recovered', 'mean_member_shift_degree'])


# %% 6. 保存正式组差与支持范围说明

output.mkdir(parents=True)
for name, product in saved.items():
    product.to_netcdf(output / (name + '.nc'))
longitude_fraction.to_netcdf(output / 'member_section_longitude_fraction.nc')
jet_summary.to_csv(output / 'group_jet_centroids.csv', index=False)
pd.DataFrame(audit).to_csv(output / 'arithmetic_qc.csv', index=False)
report = {'groups': groups, 'window': metadata['period'], 'estimator': 'Each member later10-early10, then equal model means, then not_recovered minus recovered; no drift correction.',
          'dimensions': {'section_groups': ['group', 'plev', 'lat'], 'section_contrast': ['plev', 'lat'], 'pressure_groups': ['group', 'plev', 'lat', 'lon'], 'pressure_contrast': ['plev', 'lat', 'lon'], 'surface_groups': ['group', 'lat', 'lon'], 'surface_contrast': ['lat', 'lon']},
          'section_support': metadata['section_support'] + ' Group comparison additionally uses common four-member valid latitude-pressure cells for both wind and gradient. This does not make the original longitude support identical across models.',
          'map_support': 'Common four-member valid grid cells per variable; no infilling. Surface contrasts exactly reproduce n01 products.',
          'sensitivity': 'Four dependent cross-group pairs and four leave-one contrasts retained; not independent replicates or confidence intervals. Two models per group only.',
          'jet_shift_not_minus_recovered_degree': shift_contrast, 'jet_averaging': 'Mean of saved member centroid shifts, not centroid of a composite wind profile.',
          'P_minus_E_contrast_identity_error': pe_error, 'selection': metadata['selection'],
          'limits': 'Member spatial heterogeneity remains. Composite contrast is descriptive, not a common response at every location, statistical significance, or causal identification. Temperature gradient constrains shear, not absolute jet latitude.'}
(output / 'qc.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
print(pd.DataFrame(audit).to_string(index=False), flush=True)
print(jet_summary.to_string(index=False), flush=True)
