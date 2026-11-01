# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import netCDF4 as nc
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import false_discovery_control
from threadpoolctl import threadpool_limits


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
source = root / 'output/s05_atmosphere/p07_toa_energy_balance/n18_five_model_figure3'
temperature_source = root / 'output/s03_tem/n13_six_model_global_temperature_groups'
destination = root / 'output/s05_atmosphere/p07_toa_energy_balance/n22_figure3_temperature_significance'
assert not destination.exists(), '正式结果存在，须先归档。'
models = ['CESM2', 'CNRM-ESM2-1', 'UKESM1-0-LL', 'CAS-ESM2-0', 'GFDL-ESM4', 'NorESM2-LM']
period_years = [[[1, 10], [131, 140]], [[1850, 1859], [2120, 2129]],
                [[1850, 1859], [2120, 2129]], [[200, 209], [470, 479]],
                [[6, 15], [271, 280]], [[1, 10], [271, 280]]]
members = ['r1i1p1f1', 'r1i1p1f2', 'r1i1p1f2', 'r1i1p1f1', 'r1i1p1f1', 'r1i1p1f1']
replicates, seed, block_lengths, alpha = 9999, 20261030, [2, 3, 4], .05
cell_batch = 512
annual_fields = np.empty((6, 2, 10, 180, 360))
inventory, map_errors, count_rows = [], [], []
with np.load(source / 'temperature_panels.npz') as saved:
    latitude, longitude = saved['latitude'].copy(), saved['longitude'].copy()
    reference_panels = saved['temperature_panels'].copy()
with np.load(temperature_source / 'global_temperature_groups.npz') as saved:
    reference_models = saved['model_fields'].copy()
temperature_metadata = json.loads((temperature_source / 'metadata.json').read_text(encoding='utf-8'))
assert temperature_metadata['models'] == models
target_lon, target_lat = np.meshgrid(longitude, latitude)
target_lat_rad, target_lon_rad = np.deg2rad(target_lat.ravel()), np.deg2rad(target_lon.ravel())
target_xyz = np.column_stack([np.cos(target_lat_rad) * np.cos(target_lon_rad),
                              np.cos(target_lat_rad) * np.sin(target_lon_rad), np.sin(target_lat_rad)])


# %% 3. 原始月资料计算逐年场，保持原图等月权与相同1度IDW

for model_index, model in enumerate(models):
    for period_index, (experiment, folder) in enumerate([('1pctCO2', 'cmip6tasad'),
                                                        ('1pctCO2-cdr', 'cmip6tas')]):
        first_year, last_year = period_years[model_index][period_index]
        field_parts, month_parts = [], []
        for path in sorted(Path('F:/sea', folder).rglob(
                f'tas_Amon_{model}_{experiment}_{members[model_index]}_*.nc')):
            with nc.Dataset(path) as dataset:
                assert dataset.source_id == model and dataset.experiment_id == experiment
                assert dataset.variant_label == members[model_index]
                field, time = dataset['tas'], dataset['time']
                assert field.dimensions == ('time', 'lat', 'lon') and field.units == 'K'
                dates = nc.num2date(time[:], time.units, calendar=time.calendar)
                selected = np.flatnonzero([first_year <= date.year <= last_year for date in dates])
                if not len(selected):
                    continue
                values = np.asarray(np.ma.filled(field[selected], np.nan), dtype=float)
                assert np.isfinite(values).all()
                field_parts.append(values)
                month_parts.extend([dates[index].year * 12 + dates[index].month - 1 for index in selected])
                native_latitude = np.asarray(dataset['lat'][:], dtype=float)
                native_longitude = np.asarray(dataset['lon'][:], dtype=float)
                inventory.append({'model': model, 'period': ['initial', 'final'][period_index],
                                  'experiment': experiment, 'member': members[model_index], 'path': str(path),
                                  'dimensions': str(field.dimensions), 'shape': str(field.shape),
                                  'calendar': time.calendar, 'time_units': time.units, 'units': field.units,
                                  'first_selected': str(dates[selected[0]]), 'last_selected': str(dates[selected[-1]]),
                                  'months': len(selected), 'file_bytes': path.stat().st_size,
                                  'mtime_ns': path.stat().st_mtime_ns})
        order = np.argsort(month_parts)
        np.testing.assert_array_equal(np.asarray(month_parts)[order],
                                       np.arange(first_year * 12, (last_year + 1) * 12))
        native_annual = np.concatenate(field_parts)[order].reshape(10, 12, len(native_latitude),
                                                                 len(native_longitude)).mean(axis=1) - 273.15
        native_lon, native_lat = np.meshgrid(native_longitude, native_latitude)
        # 复现现有缓存的经度表达，避免等距第4邻居因浮点舍入选择不同格点。
        if model not in ['GFDL-ESM4', 'UKESM1-0-LL']:
            native_lon = (native_lon + 180.) % 360. - 180.
        native_lat_rad, native_lon_rad = np.deg2rad(native_lat.ravel()), np.deg2rad(native_lon.ravel())
        native_xyz = np.column_stack([np.cos(native_lat_rad) * np.cos(native_lon_rad),
                                      np.cos(native_lat_rad) * np.sin(native_lon_rad), np.sin(native_lat_rad)])
        chord, neighbours = cKDTree(native_xyz).query(target_xyz, k=4)
        distance = 6371.0088 * 2 * np.arcsin(np.clip(chord / 2, 0, 1))
        distance_floor = 1e-6 if model in ['GFDL-ESM4', 'UKESM1-0-LL'] else 1.
        neighbour_weights = 1 / np.maximum(distance, distance_floor) ** 2
        neighbour_weights /= neighbour_weights.sum(axis=1, keepdims=True)
        mapped = np.sum(native_annual.reshape(10, -1)[:, neighbours] * neighbour_weights[None], axis=2)
        if model in ['GFDL-ESM4', 'UKESM1-0-LL']:
            exact = distance[:, 0] <= 1e-6
            mapped[:, exact] = native_annual.reshape(10, -1)[:, neighbours[exact, 0]]
        assert np.max(distance[:, 0]) <= 450.
        annual_fields[model_index, period_index] = mapped.reshape(10, 180, 360)
    reconstructed = annual_fields[model_index, 1].mean(axis=0) - annual_fields[model_index, 0].mean(axis=0)
    error = float(np.max(np.abs(reconstructed - reference_models[model_index])))
    map_errors.append({'model': model, 'max_reconstructed_residual_error_degC': error})
    np.testing.assert_allclose(reconstructed, reference_models[model_index], atol=2e-5, rtol=0)
    print(model, 'annual windows complete; max residual error', error, flush=True)


# %% 4. 三模式符号一致性与固定成员的组平均差

residuals = annual_fields[:, 1].mean(axis=1) - annual_fields[:, 0].mean(axis=1)
sign_agreement = np.stack([np.all(reference_models[:3] > 0, axis=0) | np.all(reference_models[:3] < 0, axis=0),
                           np.all(reference_models[3:] > 0, axis=0) | np.all(reference_models[3:] < 0, axis=0)])
observed = residuals[:3].mean(axis=0) - residuals[3:].mean(axis=0)
np.testing.assert_allclose(observed, reference_panels[2], atol=2e-5, rtol=0)
assert np.isfinite(annual_fields).all()
destination.mkdir(parents=True)
np.savez_compressed(destination / 'annual_temperature_samples.npz', latitude=latitude, longitude=longitude,
                    models=np.asarray(models), annual_tas_degC=annual_fields,
                    source_years=np.asarray(period_years), month_weighting=np.asarray('equal calendar months'))
pd.DataFrame(inventory).to_csv(destination / 'annual_source_inventory.csv', index=False)
pd.DataFrame(map_errors).to_csv(destination / 'reconstruction_qc.csv', index=False)


# %% 5. 独立模式和时期循环分块重采样，保持每次组差的三对三等权

random = np.random.default_rng(seed)
resampling, index_archive = {}, {}
coefficients = np.array([1 / 3] * 3 + [-1 / 3] * 3)
for length in block_lengths:
    weights = np.empty((replicates, 6, 2, 10))
    indices_all = np.empty((6, 2, replicates, 10), dtype=np.int16)
    for model_index in range(6):
        for period_index in range(2):
            starts = random.integers(0, 10, size=(replicates, int(np.ceil(10 / length))))
            indices = ((starts[:, :, None] + np.arange(length)) % 10).reshape(replicates, -1)[:, :10]
            indices_all[model_index, period_index] = indices
            weights[:, model_index, period_index] = np.eye(10)[indices].mean(axis=1)
    weights *= coefficients[None, :, None, None] * np.array([-1., 1.])[None, None, :, None]
    resampling[length] = weights.reshape(replicates, 120)
    index_archive['L' + str(length)] = indices_all
np.savez_compressed(destination / 'bootstrap_year_indices.npz', **index_archive)


# %% 6. 双侧中心化bootstrap，块长2/3/4取最大p，再做全球BH校正

centered_fields = annual_fields - annual_fields.mean(axis=2, keepdims=True)
centered_matrix = centered_fields.reshape(120, -1)
observed_flat = observed.ravel()
probabilities = np.empty((3, observed_flat.size))
for length_index, length in enumerate(block_lengths):
    for offset in range(0, observed_flat.size, cell_batch):
        last = min(offset + cell_batch, observed_flat.size)
        with threadpool_limits(limits=4):
            error_samples = resampling[length] @ centered_matrix[:, offset:last]
        probabilities[length_index, offset:last] = (1 + np.sum(
            np.abs(error_samples) >= np.abs(observed_flat[offset:last])[None], axis=0)) / (replicates + 1)
    print('Contrast bootstrap block', length, 'complete', flush=True)
p_max = probabilities.max(axis=0)
q_values = false_discovery_control(p_max, method='bh')
contrast_significant = (q_values < alpha).reshape(180, 360)
significant = np.concatenate([sign_agreement, contrast_significant[None]], axis=0)
np.savez_compressed(destination / 'temperature_significance.npz', latitude=latitude, longitude=longitude,
                    significant=significant, nonsignificant=~significant,
                    model_residuals=reference_models, temperature_panels=reference_panels,
                    contrast_reconstructed=observed, p_by_block=probabilities.reshape(3, 180, 360),
                    p_max=p_max.reshape(180, 360), q=q_values.reshape(180, 360))
area = np.broadcast_to(np.cos(np.deg2rad(latitude))[:, None], (180, 360)).copy()
area /= area.sum()
for index, name in enumerate(['non_recovery_sign_agreement', 'recovery_sign_agreement', 'contrast_bootstrap_fdr']):
    count_rows.append({'product': name, 'valid_cells': 64800, 'passed_cells': int(significant[index].sum()),
                       'passed_area_fraction': float(np.sum(area * significant[index])),
                       'crosshatched_cells': int((~significant[index]).sum())})
pd.DataFrame(count_rows).to_csv(destination / 'significance_summary.csv', index=False)
caption = ('Figure 3. Global radiative and near-surface air-temperature responses. '
           '(a, b) piControl-relative net radiation (black), absorbed shortwave contribution (orange), '
           'and negative outgoing-longwave contribution (blue) for the non-recovery and recovery groups. '
           'Thin lines show individual models and thick lines show equal-weight group means; all components '
           'are plotted as solid lines. (c) Non-recovery minus Recovery group means. Radiation uses two '
           'non-recovery models (CESM2 and CNRM-ESM2-1) and three recovery models (CAS-ESM2-0, GFDL-ESM4, '
           'and NorESM2-LM). Lines show saved 11-year centred means, computed within each forcing phase; '
           'the vertical dotted line marks the forcing reversal. (d, e) Group-mean temperature changes '
           'between experiment years 271–280 and 1–10, with UKESM1-0-LL also included in the non-recovery '
           'temperature group (three models in each group). Cross-hatching marks grid cells where the '
           'three model changes do not all have the same sign. (f) Difference between the two group-mean '
           'temperature changes (non-recovery minus recovery); cross-hatching marks differences that do '
           'not pass a two-sided circular year-block bootstrap test (9,999 resamples; maximum p across '
           'block lengths of 2, 3, and 4 years), followed by global Benjamini–Hochberg FDR correction '
           'at q < 0.05. The test estimates temporal sampling uncertainty conditional on the six fixed '
           'models. Temperature maps use the same colour scale.')
metadata = {'status': 'passed', 'models': models, 'period_years': period_years,
            'annual_dimensions': 'model6,period2,year_in_window10,lat180,lon360',
            'raw_dimensions': 'All six raw tas(time,lat,lon), units K; 1D native lat/lon.',
            'temporal_mean': 'Equal calendar-month means, matching existing displayed residual maps.',
            'spatial_method': 'Native annual means mapped to 1 degree using 4-neighbour spherical IDW power2, max450km; no field smoothing.',
            'd_e_criterion': 'All three individual-model residuals strictly positive or all strictly negative. Sign agreement, not a p-value test.',
            'f_observed': 'Mean of three non-recovery residual maps minus mean of three recovery residual maps.',
            'f_test': 'Two-sided centered circular-block bootstrap; independent within model/period, shared indices across global cells; fixed equal model weights.',
            'bootstrap_replicates': replicates, 'seed': seed, 'block_lengths': block_lengths,
            'multiple_testing': 'Maximum p over L2/3/4, then BH across all64800global grid cells; q<0.05.',
            'inference_limit': 'Temporal sampling conditional on six fixed model identities, not model-population or causal inference. Only ten annual samples per endpoint.',
            'temperature_fields_unchanged': True, 'reconstruction_errors': map_errors,
            'CAS_limitation': 'Existing provisional piControl mapping affects radiation; inconsistent CAS branch metadata also preclude asserting verified restart continuity. No offsets applied to temperature.',
            'caption': caption,
            'source_files': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
                             for path in [source / 'temperature_panels.npz', temperature_source / 'global_temperature_groups.npz',
                                          temperature_source / 'metadata.json']]}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
(destination / 'caption.txt').write_text(caption + '\n', encoding='utf-8')
print(pd.DataFrame(count_rows).to_string(index=False), flush=True)
print(destination.resolve(), flush=True)
