# %% 1. 导入依赖

from pathlib import Path
import json
import re
import hashlib
import numpy as np
import pandas as pd
import netCDF4
from scipy.spatial import cKDTree
from cartopy.io import shapereader
from shapely import contains_xy
from shapely.geometry import mapping


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
destination = root / 'output/s03_tem/n12_global_temperature_ice_patterns'
assert not destination.exists(), '正式结果存在，须先归档。'
destination.mkdir(parents=True)
models = ['CESM2', 'CNRM-ESM2-1', 'NorESM2-LM', 'GFDL-ESM4']
members = ['r1i1p1f1', 'r1i1p1f2', 'r1i1p1f1', 'r1i1p1f1']
calendar_starts = [[1, 1], [1850, 1990], [1, 141], [1, 141]]
experiments = ['1pctCO2', '1pctCO2-cdr']
periods = ['initial', 'later']
ice_dimensions = [('time', 'nj', 'ni'), ('time', 'y', 'x'),
                  ('time', 'j', 'i'), ('time', 'y', 'x')]
ice_coordinates = [('lat', 'lon'), ('lat', 'lon'),
                   ('latitude', 'longitude'), ('lat', 'lon')]
latitude = np.arange(-89.5, 90, 1.)
longitude = np.arange(-179.5, 180, 1.)
longitude_grid, latitude_grid = np.meshgrid(longitude, latitude)
target_lat_rad, target_lon_rad = np.deg2rad(latitude_grid.ravel()), np.deg2rad(longitude_grid.ravel())
target_xyz = np.column_stack([np.cos(target_lat_rad) * np.cos(target_lon_rad),
                              np.cos(target_lat_rad) * np.sin(target_lon_rad), np.sin(target_lat_rad)])
radius_km = 6371.0088
threshold = .15
model_fields = np.full((4, 2, 2, 180, 360), np.nan)
inventory, diagnostics = [], []
cartography_root = Path('C:/Users/Lenovo/.local/share/cartopy/shapefiles/natural_earth/physical')
north_path = root / 'output/s05_atmosphere/p05_amoc_recovery_composites/n41_observed_annual_ice_edge/observed_annual_sic_1988_2020.npz'
south_path = Path('F:/sea/seaice-obs-south/sic_pss25_197811-202607_v06r00.nc')


# %% 3. 重读四模式全球气温及海冰的两个完整十年

for model_index, (model, member) in enumerate(zip(models, members)):
    for variable_index, variable in enumerate(['tas', 'siconc']):
        for period_index, (period, experiment) in enumerate(zip(periods, experiments)):
            first_year = calendar_starts[model_index][period_index] + (130 if period_index == 1 else 0)
            last_year = first_year + 9
            suffix = 'ad' if period_index == 0 else ''
            data_root = Path('F:/sea') / (f'cmip6tas{suffix}' if variable == 'tas' else f'cmipseaice{suffix}')
            candidates = sorted(data_root.rglob(f'{variable}_*_{model}_{experiment}_{member}_*.nc'))
            pieces, month_ids = [], []
            for path in candidates:
                dates_in_name = re.search(r'_(\d{6})-(\d{6})', path.name)
                assert dates_in_name is not None
                if int(dates_in_name[1][:4]) > last_year or int(dates_in_name[2][:4]) < first_year:
                    continue
                with netCDF4.Dataset(path) as source:
                    assert source.source_id == model and source.experiment_id == experiment
                    assert source.variant_label == member
                    field = source[variable]
                    expected_dimensions = ('time', 'lat', 'lon') if variable == 'tas' else ice_dimensions[model_index]
                    assert field.dimensions == expected_dimensions
                    assert field.units == ('K' if variable == 'tas' else '%')
                    time = source['time']
                    dates = netCDF4.num2date(time[:], time.units, calendar=time.calendar)
                    selected = np.flatnonzero(np.array([first_year <= date.year <= last_year for date in dates]))
                    assert len(selected) > 0
                    values = np.asarray(np.ma.filled(field[selected], np.nan), dtype=float)
                    if variable == 'siconc':
                        values /= 100.
                        values[(values < 0) | (values > 1)] = np.nan
                    else:
                        values -= 273.15
                    pieces.append(values)
                    month_ids.extend([dates[index].year * 12 + dates[index].month for index in selected])
                    lat_name, lon_name = ('lat', 'lon') if variable == 'tas' else ice_coordinates[model_index]
                    native_latitude = source[lat_name][:].filled(np.nan)
                    native_longitude = source[lon_name][:].filled(np.nan)
                    if variable == 'tas':
                        native_longitude, native_latitude = np.meshgrid(native_longitude, native_latitude)
                    inventory.append({'model': model, 'variable': variable, 'period': period,
                                      'file': str(path), 'member': member, 'raw_dimensions': ','.join(field.dimensions),
                                      'raw_shape': str(field.shape), 'first_selected': str(dates[selected[0]]),
                                      'last_selected': str(dates[selected[-1]]), 'months': len(selected),
                                      'file_bytes': path.stat().st_size, 'mtime_ns': path.stat().st_mtime_ns,
                                      'tracking_id': getattr(source, 'tracking_id', '')})
            monthly = np.concatenate(pieces)
            order = np.argsort(month_ids)
            monthly = monthly[order]
            np.testing.assert_array_equal(np.asarray(month_ids)[order], np.arange(first_year * 12 + 1, (last_year + 1) * 12 + 1))
            assert monthly.shape[0] == 120 and monthly.shape[1:] == native_latitude.shape
            # 沿用之前全年图：每年十二个月等权，十年等权；全窗缺测保持缺测。
            native_mean = np.mean(monthly, axis=0)
            by_year = monthly.reshape((10, 12) + monthly.shape[1:]).mean(axis=1).mean(axis=0)
            mean_error = float(np.nanmax(np.abs(native_mean - by_year)))
            assert mean_error < 1e-10
            valid = np.isfinite(native_mean) & np.isfinite(native_latitude) & np.isfinite(native_longitude)
            source_lat_rad = np.deg2rad(native_latitude[valid])
            source_lon_rad = np.deg2rad(native_longitude[valid])
            xyz = np.column_stack([np.cos(source_lat_rad) * np.cos(source_lon_rad),
                                   np.cos(source_lat_rad) * np.sin(source_lon_rad), np.sin(source_lat_rad)])
            neighbor_count = 4 if variable == 'tas' else 8
            maximum_distance = 450. if variable == 'tas' else 300.
            chord, indices = cKDTree(xyz).query(target_xyz, k=neighbor_count)
            distance = radius_km * 2 * np.arcsin(np.clip(chord / 2, 0, 1))
            weights = 1 / np.maximum(distance, 1e-6) ** 2
            mapped = np.sum(weights * native_mean[valid][indices], axis=1) / weights.sum(axis=1)
            exact = distance[:, 0] <= 1e-6
            mapped[exact] = native_mean[valid][indices[exact, 0]]
            mapped[distance[:, 0] > maximum_distance] = np.nan
            model_fields[model_index, period_index, variable_index] = mapped.reshape(180, 360)
            diagnostics.append({'model': model, 'variable': variable, 'period': period,
                                'first_year': first_year, 'last_year': last_year, 'months': 120,
                                'annual_mean_order_error': mean_error, 'valid_native_cells': int(valid.sum()),
                                'valid_mapped_cells': int(np.isfinite(mapped).sum())})
            print(model, variable, period, first_year, last_year, '120 months complete', flush=True)


# %% 4. 固定两成员组均值与温度组差，陆地仅屏蔽海冰轮廓

group_fields = np.stack([model_fields[:2].mean(axis=0), model_fields[2:].mean(axis=0)])
group_tas_change = group_fields[:, 1, 0] - group_fields[:, 0, 0]
temperature_panels = np.concatenate([group_tas_change, (group_tas_change[0] - group_tas_change[1])[None]])
land_mask = np.zeros((180, 360), dtype=bool)
land_path = cartography_root / 'ne_50m_land.shp'
coast_path = cartography_root / 'ne_50m_coastline.shp'
for geometry in shapereader.Reader(land_path).geometries():
    land_mask |= contains_xy(geometry, longitude_grid, latitude_grid)
group_ice = group_fields[:, :, 1]
paired_support = np.isfinite(group_ice).all(axis=1)
ice_for_contour = np.where(paired_support[:, None] & ~land_mask[None, None], group_ice, np.nan)
coast_features = [{'type': 'Feature', 'properties': {}, 'geometry': mapping(geometry)}
                  for geometry in shapereader.Reader(coast_path).geometries()]
(destination / 'coastline_global_50m.geojson').write_text(json.dumps({'type': 'FeatureCollection', 'features': coast_features}), encoding='utf-8')


# %% 5. 保存南北极同观测时期的原生网格海冰场

with np.load(north_path) as north:
    np.savez_compressed(destination / 'observed_north.npz',
                        x=north['projected_x_m'], y=north['projected_y_m'],
                        sic=north['sic_climatology_fraction'], epsg=np.array(3411))
with netCDF4.Dataset(south_path) as source:
    assert source['cdr_seaice_conc_monthly'].dimensions == ('time', 'y', 'x')
    assert source['cdr_seaice_conc_monthly'].units == '1' and source.geospatial_bounds_crs == 'EPSG:3412'
    time = source['time']
    dates = netCDF4.num2date(time[:], time.units, calendar=getattr(time, 'calendar', 'standard'))
    selected = np.flatnonzero(np.array([1988 <= date.year <= 2020 for date in dates]))
    np.testing.assert_array_equal([dates[index].year * 12 + dates[index].month for index in selected], np.arange(1988 * 12 + 1, 2021 * 12 + 1))
    observed_monthly = np.asarray(np.ma.filled(source['cdr_seaice_conc_monthly'][selected], np.nan), dtype=float)
    observed_monthly[(observed_monthly < 0) | (observed_monthly > 1)] = np.nan
    south_shape = observed_monthly.shape[1:]
    observed_annual = observed_monthly.reshape((33, 12) + south_shape).mean(axis=1)
    valid_years = np.isfinite(observed_annual).sum(axis=0)
    observed_climate = np.divide(np.nansum(observed_annual, axis=0), valid_years,
                                 out=np.full(south_shape, np.nan), where=valid_years >= 27)
    observed_x, observed_y = np.meshgrid(source['x'][:], source['y'][:])
    np.savez_compressed(destination / 'observed_south.npz', x=observed_x, y=observed_y,
                        sic=observed_climate, valid_years=valid_years, epsg=np.array(3412))


# %% 6. 保存全球场、显示数组、输入清单与科学口径

np.savez_compressed(destination / 'global_tas_ice_fields.npz', latitude=latitude, longitude=longitude,
                    models=np.array(models), groups=np.array(['incomplete', 'substantial']),
                    model_fields=model_fields, group_fields=group_fields,
                    temperature_panels=temperature_panels, group_ice=group_ice,
                    ice_for_contour=ice_for_contour, land_mask=land_mask,
                    longitude_cyclic=np.append(longitude, longitude[0] + 360),
                    temperature_cyclic=np.concatenate([temperature_panels, temperature_panels[:, :, :1]], axis=2),
                    ice_cyclic=np.concatenate([ice_for_contour, ice_for_contour[:, :, :, :1]], axis=3))
pd.DataFrame(inventory).to_csv(destination / 'raw_inventory.csv', index=False)
pd.DataFrame(diagnostics).to_csv(destination / 'endpoint_qc.csv', index=False)
metadata = {'status': 'passed', 'groups': {'incomplete': models[:2], 'substantial': models[2:]},
            'variable': 'Near-surface air temperature tas, not sea-surface temperature',
            'periods': 'Initial=up years1-10; later=CDR years131-140; nominal experiment years271-280.',
            'forcing_processing': 'Raw forced fields, no piControl correction.',
            'temporal_mean': '120 unique months, equally weighted; retain missing values; no extra smoothing.',
            'group_mean': 'Two fixed members with equal weight and complete-member support.',
            'temperature': 'Later minus initial in degree C; panel c=incomplete change minus substantial change.',
            'ice': '15% contour of group-mean annual SIC; no subtraction of ice-edge coordinates.',
            'panel_c_ice': 'Incomplete-group later, substantial-group later, and the same modern observed reference.',
            'observations': '1988-2020 annual mean; 12 valid months per year, >=27/33 valid years per cell; no smoothing.',
            'observation_interpretation': 'Modern observational reference, not contemporaneous validation of preindustrial initial or idealized CDR endpoint.',
            'regrid': 'Global 1deg IDW on unit sphere: TAS k4 max450km; SIC k8 max300km. No filling beyond distance threshold.',
            'dimensions': 'model_fields(model4,period2,variable2,latitude180,longitude360); group_fields(group2,period2,variable2,latitude180,longitude360); variables=tas,siconc; temperature_panels(panel3,latitude180,longitude360).',
            'raw_dimensions': {'tas': 'time,lat,lon', 'siconc': dict(zip(models, [','.join(value) for value in ice_dimensions]))},
            'observed_dimensions': 'North (y448,x304), EPSG3411; south(y332,x316), EPSG3412; projected x/y in meters.',
            'tas_ranges': [[float(np.nanmin(field)), float(np.nanmax(field))] for field in temperature_panels],
            'source_files': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
                             for path in [north_path, root / 'output/s05_atmosphere/p05_amoc_recovery_composites/n49_group_naming/group_labels.json']]}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print('OUTPUT', destination.resolve(), metadata['tas_ranges'], flush=True)
