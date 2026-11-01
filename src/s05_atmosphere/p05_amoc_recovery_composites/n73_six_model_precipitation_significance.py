# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import false_discovery_control
from threadpoolctl import threadpool_limits


# %% 2. 参数与已保存的六模式逐年降水

root = Path('F:/sea/Antic sea ice constraint')
base = root / 'output/s05_atmosphere/p05_amoc_recovery_composites'
source = base / 'n72_six_model_annual_circulation_precipitation_samples'
destination = base / 'n73_six_model_precipitation_significance'
assert not destination.exists(), '正式结果存在，须先归档。'
models = ['GFDL-ESM4', 'NorESM2-LM', 'CAS-ESM2-0', 'CESM2', 'CNRM-ESM2-1', 'UKESM1-0-LL']
products = ['recovered', 'not_recovered', 'contrast']
replicates, seed, block_lengths, alpha = 9999, 20261030, [2, 3, 4], .05
coefficients = np.array([[1/3, 1/3, 1/3, 0, 0, 0], [0, 0, 0, 1/3, 1/3, 1/3],
                         [-1/3, -1/3, -1/3, 1/3, 1/3, 1/3]])
with xr.open_dataset(source / 'annual_precipitation_samples.nc') as dataset:
    annual = dataset.sel(model=models).load()
with xr.open_dataarray(source / 'group_precipitation_changes.nc') as dataset:
    reference = dataset.sel(product=products).load()
field = annual.pr.transpose('model', 'period', 'year_in_window', 'lat', 'lon')
days = annual.day_weight.transpose('model', 'period', 'year_in_window').values
values = field.values.reshape(6, 2, 10, -1)
assert np.isfinite(values).all() and days.shape == (6, 2, 10)
states = np.sum(values * (days / days.sum(axis=2, keepdims=True))[:, :, :, None], axis=2)
observed = coefficients @ (states[:, 1] - states[:, 0])
np.testing.assert_allclose(observed.reshape(reference.shape), reference, atol=1e-10, rtol=0)


# %% 3. 每模式每时期独立循环年块，共享全空间抽样索引

random = np.random.default_rng(seed)
resampling, index_archive = {}, {}
for length in block_lengths:
    weights = np.empty((6, 2, replicates, 10))
    all_indices = np.empty((6, 2, replicates, 10), dtype=np.int16)
    for model_index in range(6):
        for period_index in range(2):
            starts = random.integers(0, 10, size=(replicates, int(np.ceil(10 / length))))
            indices = ((starts[:, :, None] + np.arange(length)) % 10).reshape(replicates, -1)[:, :10]
            all_indices[model_index, period_index] = indices
            counts = np.eye(10)[indices].sum(axis=1) * days[model_index, period_index][None]
            weights[model_index, period_index] = counts / counts.sum(axis=1, keepdims=True)
    resampling[length] = weights
    index_archive['L' + str(length)] = all_indices


# %% 4. 三个面板的中心化双侧检验与分别BH校正

probabilities = np.empty((3, 3, values.shape[-1]))
for length_index, length in enumerate(block_lengths):
    weights = resampling[length]
    for offset in range(0, values.shape[-1], 256):
        last = min(offset + 256, values.shape[-1])
        with threadpool_limits(limits=4):
            member_bootstrap = np.stack([weights[index, 1] @ values[index, 1, :, offset:last]
                                         - weights[index, 0] @ values[index, 0, :, offset:last]
                                         for index in range(6)])
            for product_index, row in enumerate(coefficients):
                boot = np.einsum('m,mbc->bc', row, member_bootstrap, optimize=True)
                errors = boot - observed[product_index, offset:last][None]
                probabilities[length_index, product_index, offset:last] = (1 + np.sum(
                    abs(errors) >= abs(observed[product_index, offset:last])[None], axis=0)) / (replicates + 1)
    print('Six-model precipitation bootstrap block', length, 'complete', flush=True)
p_max = probabilities.max(axis=0)
q_values = np.stack([false_discovery_control(row, method='bh') for row in p_max])
shape = (3, field.sizes['lat'], field.sizes['lon'])
coordinates = {'product': products, 'lat': field.lat, 'lon': field.lon}
results = xr.Dataset()
for name, array in [('delta', observed), ('p_max', p_max), ('q', q_values),
                    ('tested', np.ones_like(q_values, dtype=bool)), ('significant', q_values < alpha),
                    ('nonsignificant', q_values >= alpha)]:
    results[name] = xr.DataArray(array.reshape(shape), dims=['product', 'lat', 'lon'], coords=coordinates)
results['p_by_block'] = xr.DataArray(probabilities.reshape((3,) + shape), dims=['block_length', 'product', 'lat', 'lon'],
                                     coords={'block_length': block_lengths, **coordinates})
results.attrs = {'inference': 'Temporal sampling conditional on six fixed models; not population or causal inference.',
                 'test': 'Centered two-sided circular year blocks; max p over L2/3/4; BH separately per displayed product.',
                 'units': 'mm day-1'}
destination.mkdir(parents=True)
results.to_netcdf(destination / 'pr_significance.nc')
np.savez_compressed(destination / 'bootstrap_year_indices.npz', **index_archive)
rows = []
for index, product in enumerate(products):
    rows.append({'product': product, 'tested_cells': int(q_values.shape[1]),
                 'significant_cells': int(np.sum(q_values[index] < alpha)), 'min_q': float(q_values[index].min())})
pd.DataFrame(rows).to_csv(destination / 'test_summary.csv', index=False)
metadata = {'status': 'passed', 'models': models, 'groups': {'recovered': models[:3], 'not_recovered': models[3:]},
            'products': products, 'replicates': replicates, 'seed': seed, 'block_lengths': block_lengths,
            'test': 'Independent model/period circular-year blocks; CF day-weighted window estimator; centered two-sided error bootstrap.',
            'FDR': 'max p over L2/3/4 then BH separately for3global2degree panels; q<0.05.',
            'limits': '10annual samples per endpoint; conditional on six fixed models. Rainband peak displacement not tested.',
            'source_files': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in
                             [source / 'annual_precipitation_samples.nc', source / 'group_precipitation_changes.nc', source / 'metadata.json']]}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print(pd.DataFrame(rows).to_string(index=False), flush=True)
print(destination.resolve(), flush=True)
