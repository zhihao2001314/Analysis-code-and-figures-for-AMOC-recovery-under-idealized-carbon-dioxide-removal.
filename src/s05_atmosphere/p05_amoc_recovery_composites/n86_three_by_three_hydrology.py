# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import false_discovery_control
from threadpoolctl import threadpool_limits


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
base = root / 'output/s05_atmosphere/p05_amoc_recovery_composites'
source, bootstrap_source = base / 'n81_six_model_net_precipitation', base / 'n73_six_model_precipitation_significance'
destination = base / 'n86_three_by_three_hydrology'
assert not destination.exists(), '正式结果存在，须先归档。'
coast = root / 'output/s05_atmosphere/p07_toa_energy_balance/n04_toa_spatial_patterns/coastline_global_110m.geojson'
variables, products = ['pr', 'evspsbl', 'pe'], ['recovered', 'not_recovered', 'contrast']
replicates, lengths, alpha = 9999, [2, 3, 4], .05
coefficients = np.array([-1/3, -1/3, -1/3, 1/3, 1/3, 1/3])


# %% 3. 读取完整三对三年度样本和已有P与P-E检验

with xr.open_dataset(source / 'annual_hydrology_samples.nc') as dataset:
    samples = dataset.load()
with xr.open_dataset(source / 'group_changes.nc') as dataset:
    changes = dataset.sel(product=products).load()
with xr.open_dataset(source / 'member_changes.nc') as dataset:
    individual = dataset.load()
with xr.open_dataset(bootstrap_source / 'pr_significance.nc') as dataset:
    pr_statistics = dataset.load()
with xr.open_dataset(source / 'pe_display_criteria.nc') as dataset:
    pe_statistics = dataset.load()
assert list(samples.model.values) == ['GFDL-ESM4', 'NorESM2-LM', 'CAS-ESM2-0', 'CESM2', 'CNRM-ESM2-1', 'UKESM1-0-LL']
np.testing.assert_allclose(changes.pe, changes.pr - changes.evspsbl, atol=1e-12, rtol=0)
np.testing.assert_allclose(pr_statistics.delta, changes.pr, atol=1e-12, rtol=0)
np.testing.assert_allclose(pe_statistics.delta, changes.pe, atol=1e-12, rtol=0)
weights = samples.day_weight.values
values = samples.evspsbl.values.reshape(6, 2, 10, -1)
observed = changes.evspsbl.sel(product='contrast').values.ravel()
indices_archive = np.load(bootstrap_source / 'bootstrap_year_indices.npz')
probabilities = np.empty((3, values.shape[-1]))


# %% 4. 对蒸发组间差值补算同口径循环块bootstrap

for length_index, length in enumerate(lengths):
    indices = indices_archive['L' + str(length)]
    assert indices.shape == (6, 2, replicates, 10)
    resampling = np.empty((6, 2, replicates, 10))
    for model_index in range(6):
        for period_index in range(2):
            counts = np.eye(10)[indices[model_index, period_index]].sum(axis=1)
            weighted = counts * weights[model_index, period_index][None]
            resampling[model_index, period_index] = weighted / weighted.sum(axis=1, keepdims=True)
    for offset in range(0, values.shape[-1], 256):
        last = min(offset + 256, values.shape[-1])
        with threadpool_limits(limits=4):
            boot = np.zeros((replicates, last - offset))
            for model_index, coefficient in enumerate(coefficients):
                boot += coefficient * (resampling[model_index, 1] @ values[model_index, 1, :, offset:last]
                    - resampling[model_index, 0] @ values[model_index, 0, :, offset:last])
            errors = boot - observed[None, offset:last]
            probabilities[length_index, offset:last] = (1 + np.sum(abs(errors) >= abs(observed[None, offset:last]), axis=0)) / (replicates + 1)
    print('Evaporation contrast bootstrap L=', length, 'complete', flush=True)
p_max = probabilities.max(axis=0).reshape(89, 180)
q_evaporation = false_discovery_control(p_max.ravel(), method='bh').reshape(89, 180)


# %% 5. 合并三种变量的组内一致性与组差检验

delta = np.stack([changes[variable].values for variable in variables])
passed = np.zeros_like(delta, dtype=bool)
p_values = np.stack([pr_statistics.p_max.sel(product='contrast').values, p_max, pe_statistics.p_max.values])
q_values = np.stack([pr_statistics.q.sel(product='contrast').values, q_evaporation, pe_statistics.q.values])
for variable_index, variable in enumerate(variables):
    for group_index, model_indices in enumerate([slice(0, 3), slice(3, 6)]):
        group_fields = individual[variable].isel(model=model_indices).values
        passed[variable_index, group_index] = np.all(group_fields > 0, axis=0) | np.all(group_fields < 0, axis=0)
    passed[variable_index, 2] = q_values[variable_index] < alpha
display = xr.Dataset({'delta': (['variable', 'product', 'lat', 'lon'], delta),
    'passed': (['variable', 'product', 'lat', 'lon'], passed),
    'p_contrast': (['variable', 'lat', 'lon'], p_values),
    'q_contrast': (['variable', 'lat', 'lon'], q_values)},
    coords={'variable': variables, 'product': products, 'lat': changes.lat, 'lon': changes.lon})
caption = ('Annual precipitation, evaporation and net precipitation responses across AMOC recovery groups. '
    'Columns show changes in precipitation(P), upward-positive evaporation(E, including sublimation and evapotranspiration), '
    'and net precipitation(P−E), in mmday−1. Rows show recovery(GFDL-ESM4, NorESM2-LM, CAS-ESM2-0), '
    'non-recovery(CESM2, CNRM-ESM2-1, UKESM1-0-LL), and non-recovery minus recovery. '
    'Each group is an equal-weight mean of three models. Changes are experiment years271–280 minus1–10, '
    'with CF time-duration weighting and periodic bilinear interpolation to2degrees. '
    'Dots indicate unanimous three-model sign agreement in group responses, or group differences passing '
    'the centred two-sided circular block bootstrap(9999resamples; largest p over2,3,4year blocks) '
    'and fieldwise Benjamini–Hochberg FDR atq<0.05. Cross-hatching indicates failure of the applicable criterion. '
    'No significance interpretation of sign agreement as a p-value test is implied.')


# %% 6. 保存正式产品与检验来源

destination.mkdir(parents=True)
display.to_netcdf(destination / 'hydrology_display.nc')
xr.Dataset({'p_by_block': (['block_length', 'lat', 'lon'], probabilities.reshape(3, 89, 180))},
    coords={'block_length': lengths, 'lat': changes.lat, 'lon': changes.lon}).to_netcdf(destination / 'evaporation_bootstrap.nc')
shutil.copy2(coast, destination / 'coastline.geojson')
metadata = {'variables': variables, 'products': products, 'groups': {'recovered': list(samples.model.values[:3]),
    'not_recovered': list(samples.model.values[3:])}, 'units': 'mm day-1', 'caption': caption,
    'group_criterion': 'All three model changes strictly same sign; distinct from statistical significance.',
    'contrast_criterion': 'Centered circular block bootstrap9999, maxp L2/3/4, per-variable spatial BHq<0.05.',
    'uncertainty_scope': 'Temporal sampling conditional on six fixed models; model weights/membership unchanged.',
    'raw_dimensions': 'Inherited source evspsbl/pr(time,lat,lon). No raw fields read here.',
    'saved_dimensions': 'variable3,product3,lat89,lon180', 'piControl_used': False,
    'sources': [{'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in
        [source / 'annual_hydrology_samples.nc', source / 'member_changes.nc', source / 'group_changes.nc',
         source / 'pe_display_criteria.nc', bootstrap_source / 'pr_significance.nc',
         bootstrap_source / 'bootstrap_year_indices.npz', coast]]}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print(destination.resolve(), flush=True)
