# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import false_discovery_control
from threadpoolctl import threadpool_limits


# %% 2. 参数、强度门槛与固定成员

root = Path('F:/sea/Antic sea ice constraint')
base = root / 'output/s05_atmosphere/p05_amoc_recovery_composites'
destination = base / 'n71_circulation_precipitation_significance'
assert not destination.exists(), '正式结果存在，须先归档。'
gate_path = base / 'check_circulation_precipitation_strength/metadata.json'
gate = json.loads(gate_path.read_text(encoding='utf-8'))
assert gate['status'] == 'passed'
models = ['GFDL-ESM4', 'NorESM2-LM', 'CESM2', 'CNRM-ESM2-1']
members = ['r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f2']
starts = [[1, 271], [1, 271], [1, 131], [1850, 2120]]
periods = ['early', 'later']
latitude, longitude = np.arange(-88., 90., 2.), np.arange(-180., 180., 2.)
replicates, seed, block_lengths, alpha = 9999, 20261027, [2, 3, 4], .05
input_registry = base / 'n01_endpoint_groups/qc.json'
registry = json.loads(input_registry.read_text(encoding='utf-8'))['inputs']
with xr.open_dataarray(base / 'n65_global_precipitation_rainband/member_precipitation_maps.nc') as dataset:
    reference_maps = dataset.sel(season='Annual', drop=True).load()
with xr.open_dataarray(base / 'n65_global_precipitation_rainband/group_precipitation_changes.nc') as dataset:
    reference_pr_changes = dataset.sel(season='Annual', drop=True).load()
with xr.open_dataset(base / 'n35_global_meridional_samples/group_changes.nc') as dataset:
    reference_psi_changes = dataset.psi.sel(season='Annual', plev=slice(10000, 100000), drop=True).load()
all_pr, all_psi, all_days, source_audit = [], [], [], []


# %% 3. 重建全球降水逐年样本并逐成员复现原气候态

for model_index, model in enumerate(models):
    phases = []
    with xr.open_dataset(base / 'n35_global_meridional_samples' / (model + '_samples.nc')) as dataset:
        sample = dataset.sel(model=model, season='Annual', drop=True).load()
    all_psi.append(sample.psi.sel(plev=slice(10000, 100000)).expand_dims(model=[model]))
    all_days.append(sample.day_weight.expand_dims(model=[model]))
    for period_index, period in enumerate(periods):
        paths = sorted(set(row['file'] for row in registry if row['model'] == model and
                           row['period'] == period and row['variable'] == 'pr'))
        annual_fields = []
        first = starts[model_index][period_index]
        for year_index, year in enumerate(range(first, first + 10)):
            pieces = []
            for filename in paths:
                with xr.open_dataset(filename, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as dataset:
                    assert dataset.pr.dims == ('time', 'lat', 'lon') and dataset.pr.attrs['units'] == 'kg m-2 s-1'
                    assert dataset.attrs['source_id'] == model and dataset.attrs['variant_label'] == members[model_index]
                    take = dataset.time.dt.year == year
                    if not bool(take.any()):
                        continue
                    field = dataset.pr.isel(time=take).load().astype(float) * 86400
                    pieces.append(field)
                    source_audit.append({'model': model, 'period': period, 'year': year, 'path': filename,
                                         'months': int(take.sum()), 'dimensions': 'time,lat,lon',
                                         'selected_values_sha256': hashlib.sha256(field.values.tobytes()).hexdigest()})
            monthly = xr.concat(pieces, 'time').sortby('time')
            np.testing.assert_array_equal(monthly.time.dt.month, np.arange(1, 13))
            assert monthly.time.dt.year.values.tolist() == [year] * 12 and np.isfinite(monthly).all()
            days = monthly.time.dt.days_in_month.astype(float)
            assert float(days.sum()) == float(sample.day_weight.sel(period=period, year_in_window=year_index + 1))
            native = (monthly * days).sum('time', skipna=False) / float(days.sum())
            native = native.assign_coords(lon=((native.lon + 180) % 360) - 180).sortby('lon').sortby('lat')
            periodic = xr.concat([native.isel(lon=[-1]).assign_coords(lon=native.lon.values[[-1]] - 360),
                                  native, native.isel(lon=[0]).assign_coords(lon=native.lon.values[[0]] + 360)], 'lon')
            common = periodic.interp(lat=latitude, lon=longitude).transpose('lat', 'lon')
            assert np.isfinite(common).all()
            annual_fields.append(common.expand_dims(year_in_window=[year_index + 1]))
        phase = xr.concat(annual_fields, 'year_in_window')
        reconstructed = phase.weighted(sample.day_weight.sel(period=period)).mean('year_in_window', skipna=False)
        np.testing.assert_allclose(reconstructed, reference_maps.sel(model=model, period=period), atol=1e-10, rtol=0)
        phases.append(phase.expand_dims(period=[period]))
        print(model + ' ' + period + ': annual precipitation samples verified', flush=True)
    all_pr.append(xr.concat(phases, 'period').expand_dims(model=[model]))
precipitation = xr.concat(all_pr, 'model').transpose('model', 'period', 'year_in_window', 'lat', 'lon').rename('pr')
psi = xr.concat(all_psi, 'model').transpose('model', 'period', 'year_in_window', 'plev', 'lat')
day_weight = xr.concat(all_days, 'model').transpose('model', 'period', 'year_in_window')
days = day_weight.values
assert days.shape == (4, 2, 10)
destination.mkdir(parents=True)
xr.Dataset({'pr': precipitation, 'day_weight': day_weight}).to_netcdf(destination / 'annual_precipitation_samples.nc')
pd.DataFrame(source_audit).to_csv(destination / 'annual_precipitation_source_audit.csv', index=False)


# %% 4. 同一分块样本索引用于空间场和环流

random = np.random.default_rng(seed)
resampling, index_archive = {}, {}
for length in block_lengths:
    weights = np.empty((4, 2, replicates, 10))
    indices_all = np.empty((4, 2, replicates, 10), dtype=np.int16)
    for model_index in range(4):
        for period_index in range(2):
            starts_random = random.integers(0, 10, size=(replicates, int(np.ceil(10 / length))))
            indices = ((starts_random[:, :, None] + np.arange(length)) % 10).reshape(replicates, -1)[:, :10]
            indices_all[model_index, period_index] = indices
            counts = np.eye(10)[indices].sum(axis=1) * days[model_index, period_index][None]
            weights[model_index, period_index] = counts / counts.sum(axis=1, keepdims=True)
    resampling[length] = weights
    index_archive['L' + str(length)] = indices_all
np.savez_compressed(destination / 'bootstrap_year_indices.npz', **index_archive)
test_summaries = []


# %% 5. 中心化双侧误差检验，块长敏感性后逐面板BH校正

for variable, field in [('psi', psi), ('pr', precipitation)]:
    products = ['contrast'] if variable == 'psi' else ['recovered', 'not_recovered', 'contrast']
    coefficients = np.array([[-.5, -.5, .5, .5]]) if variable == 'psi' else np.array([[.5, .5, 0, 0], [0, 0, .5, .5], [-.5, -.5, .5, .5]])
    spatial_dimensions = list(field.dims[3:])
    values = field.values.reshape(4, 2, 10, -1)
    complete_member = np.isfinite(values).all(axis=(1, 2))
    safe = np.where(np.isfinite(values), values, 0.)
    states = np.sum(safe * (days / days.sum(axis=2, keepdims=True))[:, :, :, None], axis=2)
    observed = coefficients @ (states[:, 1] - states[:, 0])
    variance = np.var(safe[:, 0], axis=1) + np.var(safe[:, 1], axis=1)
    tested = np.stack([complete_member[np.flatnonzero(row)].all(axis=0) &
                       np.any(variance[np.flatnonzero(row)] > 1e-16, axis=0) for row in coefficients])
    selected = np.flatnonzero(tested.any(axis=0))
    probabilities = np.full((3, len(products), values.shape[-1]), np.nan)
    for length_index, length in enumerate(block_lengths):
        weights = resampling[length]
        for offset in range(0, len(selected), 256):
            cells = selected[offset:offset + 256]
            with threadpool_limits(limits=4):
                member_bootstrap = np.stack([weights[m, 1] @ safe[m, 1][:, cells] -
                                              weights[m, 0] @ safe[m, 0][:, cells] for m in range(4)])
                for product_index, row in enumerate(coefficients):
                    bootstrap = np.einsum('m,mbc->bc', row, member_bootstrap, optimize=True)
                    error = bootstrap - observed[product_index, cells][None]
                    probabilities[length_index, product_index, cells] = (1 + np.sum(abs(error) >= abs(observed[product_index, cells])[None], axis=0)) / (replicates + 1)
        print(variable + ': block length ' + str(length) + ' completed', flush=True)
    p_max = np.max(probabilities, axis=0)
    q = np.full_like(p_max, np.nan)
    for product_index, product in enumerate(products):
        valid = tested[product_index]
        q[product_index, valid] = false_discovery_control(p_max[product_index, valid], method='bh')
        reference = reference_psi_changes if variable == 'psi' else reference_pr_changes
        np.testing.assert_allclose(observed[product_index].reshape(reference.sel(product=product).shape),
                                   reference.sel(product=product), atol=1e-10, rtol=0)
        test_summaries.append({'variable': variable, 'product': product, 'n_tested': int(valid.sum()),
                               'n_q_below_005': int(np.sum(q[product_index] < alpha)),
                               'min_p_max': float(np.min(p_max[product_index, valid])),
                               'min_q': float(np.min(q[product_index, valid]))})
    dimensions = ['product'] + spatial_dimensions
    coordinates = {'product': products, **{dimension: field[dimension] for dimension in spatial_dimensions}}
    shape = (len(products),) + tuple(field.sizes[dimension] for dimension in spatial_dimensions)
    results = xr.Dataset()
    results['delta'] = xr.DataArray(observed.reshape(shape), dims=dimensions, coords=coordinates)
    results['p_by_block'] = xr.DataArray(probabilities.reshape((3,) + shape), dims=['block_length'] + dimensions,
                                        coords={'block_length': block_lengths, **coordinates})
    for name, array in [('p_max', p_max), ('q', q), ('tested', tested),
                        ('significant', tested & (q < alpha)), ('nonsignificant', tested & (q >= alpha))]:
        results[name] = xr.DataArray(array.reshape(shape), dims=dimensions, coords=coordinates)
    results.attrs = {'inference': 'Temporal sampling conditional on four fixed models; not model population or causal inference.',
                     'test': 'Independent circular blocks per model/period; max p across L2/3/4, then per-panel BH q<0.05.',
                     'null': 'Zero late-minus-early response, or zero contrast of responses; no pairing across experiments.',
                     'units': '1e9 kg s-1' if variable == 'psi' else 'mm day-1'}
    results.to_netcdf(destination / (variable + '_significance.nc'))


# %% 6. 保存方法、门槛和限制

pd.DataFrame(test_summaries).to_csv(destination / 'test_summary.csv', index=False)
source_paths = [gate_path, input_registry, base / 'n35_global_meridional_samples/group_changes.nc',
                base / 'n65_global_precipitation_rainband/member_precipitation_maps.nc',
                base / 'n65_global_precipitation_rainband/group_precipitation_changes.nc']
source_paths.extend([base / 'n35_global_meridional_samples' / (model + '_samples.nc') for model in models])
report = {'status': 'passed', 'strength_gate': gate, 'models': models, 'starts': starts,
          'replicates': replicates, 'seed': seed, 'block_lengths': block_lengths, 'alpha': alpha,
          'test': 'Centered two-sided circular-year-block bootstrap error test; p=(1+count(abs(delta_star-delta_hat)>=abs(delta_hat)))/(B+1).',
          'FDR': 'Maximum p over block lengths2/3/4, then BH independently in4 displayed families: contrast Psi100-1000hPa at88S-88N; three global precipitation panels88S-88N,180longitudes.',
          'resampling_unit': 'Year blocks within each fixed model and period; same draws across all cells and both fields. No resampling of model identity.',
          'effect': 'Existing dayweighted10-year endpoint climates and fixed equal-weight group effects unchanged.',
          'limitations': ['Only10years per endpoint; low-frequency variability not well sampled; circular wrap and approximate within-window stationarity are assumptions.',
                          'Conditional fixed-model temporal inference, not inference across recovery-group model populations or AMOC causation.',
                          'Non-significance is not equivalence or recovery; peak/ITCZ displacement is not tested.',
                          'Global2degree maps use the original nonconservative interpolation; no extrapolation to polar caps; no piControl correction added.'],
          'source_files': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in source_paths],
          'code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(destination / 'metadata.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
print(pd.DataFrame(test_summaries).to_string(index=False), flush=True)
