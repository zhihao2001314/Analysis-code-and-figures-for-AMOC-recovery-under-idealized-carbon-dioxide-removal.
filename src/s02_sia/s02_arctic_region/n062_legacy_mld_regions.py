# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd


# %% 2. 参数与路径设置

root = Path(r"F:\sea\Antic sea ice constraint")
source = root / "output/s02sia/s02_arctic_region/p06_mld/n04_na_region_timeseries/na_march_mld_region_metadata.json"
destination = root / "output/s02_sia/s02_arctic_region/p21_grouped_mld/n062_legacy_mld_regions"
assert not destination.exists(), "输出已存在，须先归档。"
keys = ["labrador_irminger_box", "nordic_seas_box"]
names = ["Labrador-Irminger", "Nordic Seas"]
labels = ["1  Labrador–Irminger | 50–65°N, 65–20°W", "2  Nordic Seas | 65–80°N, 30°W–20°E"]
colors = ["#D55E00", "#0072B2"]
label_longitudes = [-42, 0]
label_latitudes = [57, 73]
vertex_rows, region_rows = [], []


# %% 3. 读取旧版实际边界并准备示意图顶点

metadata = json.loads(source.read_text(encoding="utf-8"))
assert metadata["labrador_irminger_box"] == {"lat_min": 50.0, "lat_max": 65.0, "lon_min": -65.0, "lon_max": -20.0}
assert metadata["nordic_seas_box"] == {"lat_min": 65.0, "lat_max": 80.0, "lon_min": -30.0, "lon_max": 20.0}
for number, (key, name, label, color, label_lon, label_lat) in enumerate(
        zip(keys, names, labels, colors, label_longitudes, label_latitudes), start=1):
    box = metadata[key]
    corners = np.array([[box["lon_min"], box["lat_min"]], [box["lon_max"], box["lat_min"]],
                        [box["lon_max"], box["lat_max"]], [box["lon_min"], box["lat_max"]],
                        [box["lon_min"], box["lat_min"]]])
    vertices = []
    # 加密仅使投影后的等经纬边界正确弯曲，不改变区域或科学结果。
    for edge in range(4):
        vertices.extend(np.linspace(corners[edge], corners[edge + 1], 101, endpoint=False))
    vertices.append(corners[0])
    for order, (longitude, latitude) in enumerate(vertices):
        vertex_rows.append({"region": name, "vertex_order": order, "longitude": longitude, "latitude": latitude})
    region_rows.append({"region": name, "number": number, "label": label, "color": color,
                        "line_style": "--" if number == 1 else "-", "label_longitude": label_lon,
                        "label_latitude": label_lat, **box,
                        "north_edge_inclusive": number == 2})


# %% 4. 保存分区示意图数据与来源

destination.mkdir(parents=True)
pd.DataFrame(vertex_rows).to_csv(destination / "vertices.csv", index=False)
pd.DataFrame(region_rows).to_csv(destination / "regions.csv", index=False)
report = {"source": str(source), "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
          "purpose": "New schematic faithfully reconstructed from legacy rectangular-region metadata;not a retrieved historical map",
          "dimensions": "vertices CSV:region2 x vertex_order405;regions CSV:2 rows;lon/lat in degrees east/north",
          "selection": "Original native ocean grid-cell centers;positive valid areacello;Labrador50<=lat<65,Nordic65<=lat<=80;both longitude edges inclusive;no MLD-depth threshold",
          "drawing": "Dense latitude-longitude rectangle edges;Natural Earth110m land visually masks fill only;no land-clipped polygon replaces the rectangle",
          "changes_to_existing_series": "None;no raw temperature/salinity read,no MLD recomputation"}
(destination / "metadata.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
print(destination)
