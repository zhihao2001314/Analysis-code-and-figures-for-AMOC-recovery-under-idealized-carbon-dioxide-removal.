# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import shutil
import pandas as pd


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
coastal_path = root / 'output/s02sia/s02_arctic_region/p14_mld_fixed_regions_raw/mld_coast_regions_20261012.geojson'
rectangular_directory = root / 'output/s02_sia/s02_arctic_region/p21_grouped_mld/n062_legacy_mld_regions'
destination = root / 'output/s02_sia/s02_arctic_region/p23_stratification_mld_composite/n076_region_locator_inputs'
assert not destination.exists(), '正式分区输入存在，重跑前先归档。'
source_files = [coastal_path, rectangular_directory / 'vertices.csv', rectangular_directory / 'regions.csv']


# %% 3. 核对两种既有区域定义，不重画边界

coastal = json.loads(coastal_path.read_text(encoding='utf-8'))
assert [feature['properties']['region'] for feature in coastal['features']] == ['Labrador-Irminger', 'Nordic Seas']
assert all(feature['geometry']['type'] == 'Polygon' for feature in coastal['features'])
vertices = pd.read_csv(source_files[1])
regions = pd.read_csv(source_files[2])
assert len(regions) == 2 and len(vertices) == 810
assert vertices.groupby('region').size().to_dict() == {'Labrador-Irminger': 405, 'Nordic Seas': 405}
assert regions.loc[regions.region == 'Labrador-Irminger', ['lat_min', 'lat_max', 'lon_min', 'lon_max']].to_numpy().tolist() == [[50., 65., -65., -20.]]
assert regions.loc[regions.region == 'Nordic Seas', ['lat_min', 'lat_max', 'lon_min', 'lon_max']].to_numpy().tolist() == [[65., 80., -30., 20.]]


# %% 4. 原样复制分区与来源，不改动科学输入

destination.mkdir(parents=True)
target_names = ['n2_coastal_regions.geojson', 'mld_rectangular_vertices.csv', 'mld_rectangular_regions.csv']
records = []
for source, name in zip(source_files, target_names):
    target = destination / name
    shutil.copy2(source, target)
    source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    assert hashlib.sha256(target.read_bytes()).hexdigest() == source_hash
    records.append({'source': str(source), 'target': str(target), 'sha256': source_hash})
metadata = {'status': 'passed', 'purpose': 'Existing geographical-boundary inputs for a-panel locator inset; no geometry modification or MLD/N2 recomputation.',
            'coastal': 'Actual n039/n043/n059 N2 coastal polygons for Labrador-Irminger and Nordic-Iceland sectors; not legacy MLD rectangles.',
            'rectangles': 'Actual n061 MLD-series boxes: Labrador50<=lat<65N,65-20W; Nordic65<=lat<=80N,30W-20E. Ocean cells only; not N2 coastal masks.',
            'dimensions': {'coastal': 'FeatureCollection with2Polygon features, coordinates longitude,latitude in degreesE,degreesN',
                           'rectangular_vertices': 'region2 x vertex_order405; longitude,latitude', 'rectangular_regions': 'region2'},
            'selection_policy': 'Plot selection must explicitly name the chosen boundary set; never imply that the two definitions match.',
            'sources': records}
(destination / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print(destination, flush=True)
