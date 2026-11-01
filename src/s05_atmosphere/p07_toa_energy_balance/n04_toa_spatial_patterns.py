# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import netCDF4 as nc
import numpy as np
import pandas as pd


# %% 2. 参数、已核验来源与统一网格

root = Path('F:/sea/Antic sea ice constraint')
source = root / 'output/s05_atmosphere/p07_toa_energy_balance/n01_toa_energy_balance'
destination = root / 'output/s05_atmosphere/p07_toa_energy_balance/n04_toa_spatial_patterns'
assert not destination.exists(), '正式结果已存在，须先归档。'
destination.mkdir(parents=True)
metadata = json.loads((source / 'metadata.json').read_text(encoding='utf-8'))
registry = json.loads((root / 'output/s05_atmosphere/p05_amoc_recovery_composites/n49_group_naming/group_labels.json').read_text(encoding='utf-8'))
monthly = pd.read_csv(source / 'monthly_global_components_long.csv')
endpoint = pd.read_csv(source / 'endpoint_summary.csv')
models = metadata['models']
variables, periods = ['rsdt', 'rsut', 'rlut'], ['initial', 'late']
field_names = ['eei', 'asr', 'negative_olr']
latitude_edges = np.arange(-90, 91, 3, dtype=float)
longitude_edges = np.arange(0, 361, 3, dtype=float)
latitude = (latitude_edges[:-1] + latitude_edges[1:]) / 2
longitude = (longitude_edges[:-1] + longitude_edges[1:]) / 2
target_weights = np.diff(np.sin(np.deg2rad(latitude_edges)))[:, None] * np.diff(np.deg2rad(longitude_edges))[None, :]
target_weights /= target_weights.sum()
remapped = np.full((8, 2, 3, len(latitude), len(longitude)), np.nan)
global_checks, file_records, regional_records, range_records = [], [], [], []


# %% 3. 重读原生月场端点，先按日数求年均，再10年等权

for model_index, model in enumerate(models):
    with np.load(source / ('weights_' + model + '.npz')) as saved:
        native_latitude, native_longitude, native_weights = saved['lat'], saved['lon'], saved['weight']
    native_fields = np.full((2, 3, len(native_latitude), len(native_longitude)), np.nan)
    native_latitude_edges = None
    for period_index, period in enumerate(periods):
        phase_index = 0 if period == 'initial' else 1
        phase = ['increase', 'removal'][phase_index]
        first_year = metadata['starts_increase_removal'][model][phase_index] + (0 if period == 'initial' else 130)
        for variable_index, variable in enumerate(variables):
            selected = monthly[(monthly.model == model) & (monthly.phase == phase) & (monthly.variable == variable)
                               & monthly.source_year.between(first_year, first_year + 9)].sort_values(['source_year', 'month'])
            if selected.empty:
                assert model == 'UKESM1-0-LL' and period == 'late' and variable == 'rsdt'
                continue
            np.testing.assert_array_equal(selected.source_year.to_numpy() * 12 + selected.month.to_numpy(),
                                          np.arange(first_year * 12 + 1, (first_year + 10) * 12 + 1))
            annual_fields = np.full((10, len(native_latitude), len(native_longitude)), np.nan)
            for path_string, file_part in selected.groupby('path', sort=False):
                path = Path(path_string)
                with nc.Dataset(path) as dataset:
                    assert dataset.source_id == model and dataset.experiment_id == ['1pctCO2', '1pctCO2-cdr'][phase_index]
                    assert dataset.variant_label == next(row.member for row in pd.read_csv(source / 'source_inventory.csv').query('model == @model').itertuples())
                    field = dataset.variables[variable]
                    assert field.dimensions == ('time', 'lat', 'lon') and field.units == 'W m-2'
                    np.testing.assert_array_equal(dataset.variables['lat'][:], native_latitude)
                    np.testing.assert_array_equal(dataset.variables['lon'][:], native_longitude)
                    if native_latitude_edges is None:
                        if model == 'CNRM-ESM2-1':
                            nodes, gaussian_weights = np.polynomial.legendre.leggauss(len(native_latitude))
                            # 高斯格点用积分权重构造等效面积边界，不猜测纬度中点。
                            sine_edges = np.concatenate(([-1], -1 + np.cumsum(gaussian_weights)))
                            sine_edges[-1] = 1
                            native_latitude_edges = np.rad2deg(np.arcsin(sine_edges))
                        else:
                            bounds = np.asarray(dataset.variables[dataset.variables['lat'].bounds][:], dtype=float)
                            native_latitude_edges = np.concatenate((bounds[:, 0], bounds[-1:, 1]))
                    time = dataset.variables['time']
                    for year, year_rows in file_part.groupby('source_year'):
                        indices = year_rows.file_index.to_numpy(dtype=int)
                        assert len(indices) == 12
                        dates = nc.num2date(time[indices], time.units, calendar=time.calendar)
                        np.testing.assert_array_equal([date.month for date in dates], np.arange(1, 13))
                        assert all(date.year == year for date in dates)
                        time_bounds = np.asarray(dataset.variables[time.bounds][indices], dtype=float)
                        pairs = nc.num2date(time_bounds.reshape(-1), time.units, calendar=time.calendar).reshape(12, 2)
                        days = np.array([(pair[1] - pair[0]).total_seconds() / 86400 for pair in pairs])
                        np.testing.assert_array_equal(days, year_rows.days)
                        values = np.ma.asarray(field[indices, :, :], dtype=np.float64)
                        assert not np.ma.getmaskarray(values).any() and np.isfinite(values).all()
                        annual_fields[int(year) - first_year] = np.einsum('tyx,t->yx', np.asarray(values), days / days.sum())
                file_records.append({'model': model, 'period': period, 'variable': variable, 'path': str(path),
                                     'selected_months': len(file_part), 'dimensions': 'time,lat,lon', 'units': 'W m-2'})
            assert np.isfinite(annual_fields).all()
            native_fields[period_index, variable_index] = annual_fields.mean(axis=0)
    np.savez_compressed(destination / ('native_endpoints_' + model + '.npz'), model=model,
                        lat=native_latitude, lon=native_longitude, lat_edges=native_latitude_edges,
                        fields=native_fields, variable=np.array(variables), period=np.array(periods))


    # %% 原生网格与3度共同网格的球面面积交叠，严格守恒，不外推

    source_sine = np.sin(np.deg2rad(native_latitude_edges))
    target_sine = np.sin(np.deg2rad(latitude_edges))
    lat_overlap = np.maximum(0, np.minimum(target_sine[1:, None], source_sine[None, 1:])
                             - np.maximum(target_sine[:-1, None], source_sine[None, :-1]))
    lat_transfer = lat_overlap / np.diff(target_sine)[:, None]
    spacing = 360 / len(native_longitude)
    lon_lower, lon_upper = native_longitude - spacing / 2, native_longitude + spacing / 2
    lon_overlap = np.zeros((len(longitude), len(native_longitude)))
    for shift in [-360, 0, 360]:
        lon_overlap += np.maximum(0, np.minimum(longitude_edges[1:, None], lon_upper[None, :] + shift)
                                 - np.maximum(longitude_edges[:-1, None], lon_lower[None, :] + shift))
    lon_transfer = lon_overlap / np.diff(longitude_edges)[:, None]
    np.testing.assert_allclose(lat_transfer.sum(axis=1), 1, atol=1e-12)
    np.testing.assert_allclose(lon_transfer.sum(axis=1), 1, atol=1e-12)
    for period_index, period in enumerate(periods):
        for variable_index, variable in enumerate(variables):
            field = native_fields[period_index, variable_index]
            if np.isnan(field).all():
                continue
            mapped = lat_transfer @ field @ lon_transfer.T
            remapped[model_index, period_index, variable_index] = mapped
            native_mean = float(np.sum(field * native_weights))
            mapped_mean = float(np.sum(mapped * target_weights))
            saved_mean = float(endpoint[(endpoint.model == model) & (endpoint.window == period)].iloc[0][variable])
            assert abs(native_mean - saved_mean) < 5e-10
            assert abs(mapped_mean - native_mean) < 5e-10
            global_checks.append({'model': model, 'period': period, 'variable': variable,
                                  'native_global': native_mean, 'remapped_global': mapped_mean,
                                  'timeseries_endpoint': saved_mean, 'conservation_error': mapped_mean - native_mean,
                                  'timeseries_error': native_mean - saved_mean})
    print(model, 'native endpoints and conservative maps complete', flush=True)


# %% 4. 净辐射和短长波贡献，分组严格说明可用成员

component_climate = np.stack((remapped[:, :, 0] - remapped[:, :, 1] - remapped[:, :, 2],
                              remapped[:, :, 0] - remapped[:, :, 1], -remapped[:, :, 2]), axis=2)
component_change = component_climate[:, 1] - component_climate[:, 0]
np.testing.assert_allclose(component_change[:, 0], component_change[:, 1] + component_change[:, 2], atol=1e-10, equal_nan=True)
group_arrays, membership = [], []
for group, original_members in registry['members'].items():
    available = [model for model in original_members if model != 'UKESM1-0-LL']
    selected_indices = [models.index(model) for model in available]
    group_arrays.append(component_change[selected_indices].mean(axis=0))
    membership.append({'group': group, 'original_members': original_members, 'available_members': available,
                       'expected_n': len(original_members), 'actual_n': len(available),
                       'status': 'complete' if len(available) == len(original_members) else 'explicit_available_subset_not_full_group'})
group_change = np.stack(group_arrays)
core_difference = group_change[0] - group_change[1]


# %% 5. 地区平均、符号覆盖与一致性诊断，禁止由局地辐射直接推断局地升降温

regions = {'global': [-90, 90, 0, 360], 'north_atlantic': [0, 70, 280, 380],
           'subpolar_atlantic': [45, 70, 280, 360], 'tropical_atlantic': [-20, 20, 290, 380],
           'north_pacific': [20, 65, 120, 240], 'southern_ocean': [-90, -45, 0, 360]}
for region, (south, north, west, east) in regions.items():
    latitude_width = np.maximum(0, np.sin(np.deg2rad(np.minimum(latitude_edges[1:], north)))
                                - np.sin(np.deg2rad(np.maximum(latitude_edges[:-1], south))))
    longitude_width = np.zeros(len(longitude))
    for shift in [-360, 0, 360]:
        longitude_width += np.maximum(0, np.minimum(longitude_edges[1:] + shift, east)
                                     - np.maximum(longitude_edges[:-1] + shift, west))
    weights = latitude_width[:, None] * longitude_width[None, :]
    weights /= weights.sum()
    for group_index, member_record in enumerate(membership):
        for field_index, field in enumerate(field_names):
            values = group_change[group_index, field_index]
            regional_records.append({'region': region, 'group': member_record['group'], 'actual_n': member_record['actual_n'],
                                     'field': field, 'mean_Wm2': float(np.sum(values * weights)),
                                     'positive_area_fraction': float(np.sum((values > 0) * weights))})
for group_index, member_record in enumerate(membership):
    for field_index, field in enumerate(field_names):
        values = group_change[group_index, field_index]
        range_records.append({'group': member_record['group'], 'field': field, 'minimum': float(values.min()),
                              'maximum': float(values.max()), 'global_mean': float(np.sum(values * target_weights))})


# %% 6. 保存绘图专用数组与科学审计，不改变已有时间序列

np.savez_compressed(destination / 'toa_pattern_fields.npz', model=np.array(models), field=np.array(field_names),
                    period=np.array(periods), group=np.array(list(registry['members'])),
                    latitude=latitude, longitude=longitude, latitude_edges=latitude_edges, longitude_edges=longitude_edges,
                    target_weights=target_weights, individual_climate=component_climate, individual_change=component_change,
                    group_change=group_change, core_difference=core_difference)
pd.DataFrame(file_records).to_csv(destination / 'source_windows.csv', index=False)
pd.DataFrame(global_checks).to_csv(destination / 'global_conservation_qc.csv', index=False)
pd.DataFrame(regional_records).to_csv(destination / 'regional_summary.csv', index=False)
pd.DataFrame(range_records).to_csv(destination / 'group_ranges.csv', index=False)
(destination / 'membership.json').write_text(json.dumps(membership, ensure_ascii=False, indent=2), encoding='utf-8')
result = {'status': 'passed_with_explicit_UKESM_missing', 'model_count': 8, 'complete_change_models': 7,
          'component_sign': 'All three plotted fields positive downward: EEI, ASR, negativeOLR',
          'endpoint': 'Annual late271-280 minus initial1-10; CF duration per annual mean, ten years equal',
          'source_windows': metadata['starts_increase_removal'],
          'regrid': 'First-order spherical area-overlap conservative; CNRM Gaussian quadrature equivalent-area latitude bounds',
          'target_grid': '3degree regular global60x120, no temporal interpolation, no drift correction',
          'group_composition': membership, 'longwave_same_subset': True,
          'qc_global_rows': len(global_checks),
          'max_conservation_error_Wm2': float(pd.DataFrame(global_checks).conservation_error.abs().max()),
          'max_timeseries_error_Wm2': float(pd.DataFrame(global_checks).timeseries_error.abs().max()),
          'code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(destination / 'metadata.json').write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(result, ensure_ascii=False, indent=2), flush=True)
