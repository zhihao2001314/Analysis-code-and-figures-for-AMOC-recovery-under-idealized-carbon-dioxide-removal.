# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import pandas as pd


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
source = root / 'output/s02sia/s02_arctic_region/p04_pattern_unadjusted/n02_temperature_pattern'
observation = root / 'output/s02sia/s02_arctic_region/p02_pattern_region/n02/observed_sic_1988_2020_climatology_edge.npz'
coast = root / 'output/s05_atmosphere/p07_toa_energy_balance/n04_toa_spatial_patterns/coastline_global_110m.geojson'
destination = root / 'output/s03_tem/n19_arctic_partition_tas_ice'
assert not destination.exists(), '正式结果存在，须先归档。'
models = ['ACCESS-ESM1-5', 'CanESM5', 'CESM2', 'CNRM-ESM2-1', 'GFDL-ESM4',
          'MIROC-ES2L', 'NorESM2-LM', 'UKESM1-0-LL', 'CAS-ESM2-0']
latitude, longitude = np.arange(45.5, 90., 1.), np.arange(-179.5, 180., 1.)
regions = {'CA': [70., 90., -180., 180.], 'NA': [45., 70., -60., 30.], 'NP': [45., 70., 120., -150.]}


# %% 3. 读取已更新CAS的九模式初末期诊断

cells = pd.read_csv(source / 'arctic_temperature_pattern_cells.csv', float_precision='round_trip')
summary = pd.read_csv(source / 'arctic_temperature_region_summary.csv', float_precision='round_trip')
source_metadata = json.loads((source / 'arctic_temperature_metadata.json').read_text(encoding='utf-8'))
fields = np.empty((9, 3, 45, 360))
for index, model in enumerate(models):
    selected = cells.loc[cells.model.eq(model)].sort_values(['latitude_deg_north', 'longitude_deg_east'])
    assert len(selected) == 45 * 360
    np.testing.assert_array_equal(selected.latitude_deg_north.unique(), latitude)
    np.testing.assert_array_equal(np.sort(selected.longitude_deg_east.unique()), longitude)
    fields[index] = selected[['raw_delta_tas_degC', 'edge_raw_co2_sic_fraction',
                              'edge_raw_cdr_sic_fraction']].to_numpy().T.reshape(3, 45, 360)
assert np.isfinite(fields[:, 0]).all()
caption = ('Arctic near-surface air-temperature residuals and sea-ice edges for nine models. '
    'Shading shows final minus initial annual-mean tas (degreesC). Initial and final are experiment '
    'years1–10 and271–280, using the existing endpoint diagnostics. Purple and green dashed lines '
    'show the initial and final15% sea-ice concentration contours, respectively; magenta lines show '
    'the NOAA/NSIDC1988–2020 observed annual-mean15% edge. Coloured region boundaries identify '
    'CA(70–90N, all longitudes), NA(45–70N,60W–30E), and NP(45–70N,120E–150W). '
    'Regional annotations give existing area-weighted tas changes over fixed geographical regions '
    'on common nine-model valid support, with no sea-ice mask. No significance test is displayed.')


# %% 4. 保存绘图数组及原数据来源

destination.mkdir(parents=True)
np.savez_compressed(destination / 'arctic_partition_tas_ice.npz', models=np.asarray(models),
    latitude=latitude, longitude=longitude, tas_change=fields[:, 0], initial_sic=fields[:, 1], final_sic=fields[:, 2])
summary[['model', 'region', 'region_abbreviation', 'raw_delta_tas_degC']].to_csv(destination / 'regional_tas_changes.csv', index=False)
shutil.copy2(observation, destination / observation.name)
shutil.copy2(coast, destination / 'coastline.geojson')
metadata = {'models': models, 'regions': regions, 'dimensions': 'fields(model9,latitude45,longitude360); no raw NetCDF read.',
    'variable': 'tas: near-surface air temperature over both land and ocean; not SST or SIC filling.',
    'definition': source_metadata['definition'], 'regrid': source_metadata['regrid'],
    'piControl_used': False, 'new_calculation': 'None; exact saved temperature and ice fields, no smoothing or rescaling.',
    'observed_edge': 'Existing NOAA/NSIDC1988–2020 climatology in EPSG3411; existing smoothed-for-edge field retained.',
    'color_limits_degC': [-8., 8.], 'outside_color_limits_cells': int((abs(fields[:, 0]) > 8).sum()),
    'caption': caption, 'sources': [{'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in
        [source / 'arctic_temperature_pattern_cells.csv', source / 'arctic_temperature_region_summary.csv',
         source / 'arctic_temperature_metadata.json', observation, coast]]}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print(destination.resolve(), flush=True)
