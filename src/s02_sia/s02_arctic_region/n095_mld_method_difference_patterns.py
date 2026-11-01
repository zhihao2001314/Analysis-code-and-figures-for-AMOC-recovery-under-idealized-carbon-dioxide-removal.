# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import shutil
import netCDF4 as nc
import numpy as np
import pandas as pd


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
base = root / 'output/s02_sia/s02_arctic_region/p24_mld_method_comparison'
paired_source = base / 'n079_mld_method_comparison'
pattern_source = base / 'n081_mld_method_pattern_inputs'
destination = base / 'n095_mld_method_difference_patterns'
models = ['ACCESS-ESM1-5', 'CanESM5', 'CESM2', 'CNRM-ESM2-1', 'GFDL-ESM4', 'NorESM2-LM', 'UKESM1-0-LL']
all_models = models + ['MIROC-ES2L', 'CAS-ESM2-0']
quantities = ['initial', 'late', 'change']
experiments = ['1pctCO2', '1pctCO2-cdr']
members = {'ACCESS-ESM1-5': 'r1i1p1f1', 'CanESM5': 'r1i1p2f1', 'CESM2': 'r1i1p1f1',
           'CNRM-ESM2-1': 'r1i1p1f2', 'GFDL-ESM4': 'r1i1p1f1', 'NorESM2-LM': 'r1i1p1f1',
           'UKESM1-0-LL': 'r1i1p1f2', 'MIROC-ES2L': 'r1i1p1f2', 'CAS-ESM2-0': 'r1i1p1f1'}
difference_limit_m = 600.
color_levels = np.arange(-difference_limit_m, difference_limit_m + 25., 25.)
color_ticks = np.arange(-difference_limit_m, difference_limit_m + 200., 200.)
availability_rows, inventory, native_statistics, display_statistics = [], [], [], []
assert not destination.exists(), '正式结果已存在，须先归档。'


# %% 3. 核对已有配对产品和当前原始月份覆盖

pattern_metadata = json.loads((pattern_source / 'metadata.json').read_text(encoding='utf-8'))
assert pattern_metadata['models'] == models and pattern_metadata['methods'] == ['native', 'derived']
for record in pattern_metadata['sources']:
    path = Path(record['path'])
    assert hashlib.sha256(path.read_bytes()).hexdigest() == record['sha256'], path
matched_months = pd.read_csv(paired_source / 'matched_months.csv')
for model in all_models:
    for experiment_index, experiment in enumerate(experiments):
        phase = ['initial', 'late'][experiment_index]
        folder = Path('F:/sea/cmip6mlotstad' if experiment_index == 0 else 'F:/sea/cmip6mlotst')
        paths = sorted(folder.rglob(f'mlotst_Omon_{model}_{experiment}_{members[model]}_*.nc'))
        months = []
        for path in paths:
            with nc.Dataset(path) as dataset:
                assert dataset.source_id == model and dataset.experiment_id == experiment
                assert dataset.variant_label == members[model] and dataset['mlotst'].units == 'm'
                assert dataset['mlotst'].dimensions[0] == 'time'
                dates = nc.num2date(dataset['time'][:], dataset['time'].units, calendar=dataset['time'].calendar)
                months.extend([date.year * 100 + date.month for date in dates])
                inventory.append({'model': model, 'experiment': experiment, 'file': str(path),
                                  'dimensions': ','.join(dataset['mlotst'].dimensions),
                                  'shape': str(dataset['mlotst'].shape), 'units': dataset['mlotst'].units,
                                  'standard_name': dataset['mlotst'].standard_name, 'calendar': dataset['time'].calendar,
                                  'first_month': str(dates[0]), 'last_month': str(dates[-1]),
                                  'file_bytes': path.stat().st_size, 'mtime_ns': path.stat().st_mtime_ns})
        requested = matched_months.loc[matched_months.model.eq(model) & matched_months.period.eq(phase), 'source_year_month'].to_numpy()
        missing = sorted(set(requested) - set(months))
        row = {'model': model, 'experiment': experiment, 'file_count': len(paths),
               'unique_months': len(set(months)), 'paired_march_count': len(requested),
               'paired_march_missing': len(missing), 'paired_endpoint_complete': bool(len(requested) == 10 and not missing),
               'first_available_year_month': min(months) if months else '',
               'last_available_year_month': max(months) if months else ''}
        availability_rows.append(row)
        if model in models:
            assert row['paired_endpoint_complete'], (model, phase, missing)
availability = pd.DataFrame(availability_rows)
assert availability.groupby('model').paired_endpoint_complete.all().sum() == 7


# %% 4. 同网格场直接相减，不重算温盐诊断或插值

with np.load(pattern_source / 'paired_mld_pattern_fields.npz') as saved:
    np.testing.assert_array_equal(saved['models'], models)
    np.testing.assert_array_equal(saved['methods'], ['native', 'derived'])
    longitude, latitude = saved['longitude'].copy(), saved['latitude'].copy()
    initial, late, change = saved['initial'].copy(), saved['late'].copy(), saved['change'].copy()
    land_mask = saved['land_mask'].copy()
assert initial.shape == late.shape == change.shape == (7, 2, 141, 341)
assert np.array_equal(np.isfinite(initial[:, 0]), np.isfinite(initial[:, 1]))
np.testing.assert_allclose(late - initial, change, equal_nan=True, atol=1e-10, rtol=0)
difference = np.stack([initial[:, 1] - initial[:, 0], late[:, 1] - late[:, 0],
                       change[:, 1] - change[:, 0]], axis=1)
np.testing.assert_allclose(difference[:, 2], difference[:, 1] - difference[:, 0], equal_nan=True, atol=1e-10, rtol=0)
paired_cells = pd.read_csv(paired_source / 'paired_native_cells.csv')
reference_statistics = pd.read_csv(paired_source / 'paired_statistics.csv')
for model_index, model in enumerate(models):
    paired = paired_cells.loc[paired_cells.model.eq(model)]
    for quantity_index, quantity in enumerate(quantities):
        source_quantity = 'delta' if quantity == 'change' else quantity
        error = paired['derived_' + source_quantity + '_m'] - paired['native_' + source_quantity + '_m']
        for region in ['LI', 'NS', 'NA_domain']:
            selected = paired['in_' + region].to_numpy(dtype=bool)
            weights = paired.loc[selected, 'area_m2'].to_numpy()
            values = error.to_numpy()[selected]
            mean = float(np.average(values, weights=weights))
            reference = reference_statistics.loc[reference_statistics.model.eq(model) & reference_statistics.region.eq(region)
                                                  & reference_statistics.quantity.eq(source_quantity)].iloc[0]
            assert abs(mean - reference.bias_derived_minus_native_m) < 1e-9
            native_statistics.append({'model': model, 'quantity': quantity, 'region': region,
                                      'paired_cells': len(values), 'area_weighted_difference_m': mean,
                                      'area_weighted_mae_m': float(np.average(abs(values), weights=weights)),
                                      'area_weighted_rmse_m': float(np.sqrt(np.average(values ** 2, weights=weights))),
                                      'minimum_difference_m': float(values.min()), 'maximum_difference_m': float(values.max())})
        display = difference[model_index, quantity_index]
        display_statistics.append({'model': model, 'quantity': quantity, 'finite_cells': int(np.isfinite(display).sum()),
                                   'minimum_m': float(np.nanmin(display)), 'maximum_m': float(np.nanmax(display)),
                                   'median_absolute_m': float(np.nanmedian(abs(display))),
                                   'abs_95th_percentile_m': float(np.nanpercentile(abs(display), 95)),
                                   'outside_color_limit_cells': int(np.sum(abs(display) > difference_limit_m))})


# %% 5. 保存差值场、资料覆盖核对和来源

destination.mkdir(parents=True)
np.savez_compressed(destination / 'mld_method_difference_patterns.npz', models=np.asarray(models), quantities=np.asarray(quantities),
                    longitude=longitude, latitude=latitude, difference=difference, land_mask=land_mask,
                    color_levels=color_levels, color_ticks=color_ticks)
availability.to_csv(destination / 'current_mlotst_availability.csv', index=False)
pd.DataFrame(inventory).to_csv(destination / 'current_mlotst_inventory.csv', index=False)
pd.DataFrame(native_statistics).to_csv(destination / 'native_difference_statistics.csv', index=False)
pd.DataFrame(display_statistics).to_csv(destination / 'display_difference_statistics.csv', index=False)
matched_months.to_csv(destination / 'matched_march_dates.csv', index=False)
for name in ['land_50m.geojson', 'coastline_50m.geojson']:
    shutil.copy2(pattern_source / name, destination / name)
source_paths = [pattern_source / 'paired_mld_pattern_fields.npz', pattern_source / 'metadata.json',
                paired_source / 'paired_native_cells.csv', paired_source / 'paired_statistics.csv',
                paired_source / 'matched_months.csv']
captions = {}
for quantity in quantities:
    explanation = {'initial': 'initial March MLD', 'late': 'final March MLD',
                   'change': 'final-minus-initial March MLD response'}[quantity]
    captions[quantity] = ('Differences in ' + explanation + ' between temperature-salinity-derived mixed-layer depth and native CMIP6 mlotst for seven models. '
                         'Shading is temperature-salinity diagnostics minus mlotst, in metres. Both methods use identical ten-March samples, '
                         'ensemble members, and paired native-grid support. GFDL-ESM4 initial samples are years6-15. '
                         'Initial and final differences compare absolute MLD climatologies; response differences compare the two late-minus-initial fields. '
                         'For response differences, positive values mean a more positive MLD response in the derived product, not necessarily absolute deepening. '
                         'Data are directly subtracted on the existing common0.25degree display grid, with no new interpolation or smoothing. '
                         'The grid spacing is for display, not model resolution. Grey areas represent land or missing paired support. '
                         'Colourbar extensions retain differences beyond600m; no significance inference is added.')
metadata = {'status': 'passed', 'models': models, 'excluded': {'CAS-ESM2-0': 'No local1pctCO2-cdr mlotst files',
            'MIROC-ES2L': 'No local forced mlotst files in either phase'}, 'quantities': quantities,
            'dimensions': 'difference(model7,quantity3,latitude141,longitude341)',
            'difference_sign': 'Temperature-salinity-derived minus native mlotst',
            'source_derived_definition': 'TEOS10 sigma0; density reference interpolated to10m;0.03kg/m3 increase; no crossing uses deepest valid sample.',
            'source_pairing': 'n079 exact native coordinate matching with common valid support and native areacello weights.',
            'time_support': 'Identical10March dates per period and model from matched_march_dates.csv; GFDL initial6-15.',
            'raw_dimension_inventory': 'Actual NetCDF dimensions and calendars reread in current_mlotst_inventory.csv; no raw MLD values recomputed.',
            'operations': 'Read already paired common-grid initial/late/change; directly subtract derived-native; no additional interpolation or smoothing.',
            'extent': pattern_metadata['extent'], 'color_limits_m': [-difference_limit_m, difference_limit_m],
            'uncertainty': 'Method agreement diagnostic, not a statistical significance test; mlotst is not treated as observational truth.',
            'pi_control_correction': False, 'captions': captions,
            'source_files': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in source_paths]}
(destination / 'metadata.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
print(availability.to_string(index=False), flush=True)
print(pd.DataFrame(display_statistics).to_string(index=False), flush=True)
print(destination.resolve(), flush=True)
