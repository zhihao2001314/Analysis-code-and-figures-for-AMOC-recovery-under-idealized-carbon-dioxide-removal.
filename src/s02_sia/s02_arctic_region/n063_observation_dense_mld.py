# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy.ndimage import label
from cartopy.io.shapereader import Reader
from shapely.geometry import mapping


# %% 2. 参数与路径设置

root = Path(r"F:\sea\Antic sea ice constraint")
destination = root / "output/s02_sia/s02_arctic_region/p21_grouped_mld/n063_observation_dense_mld"
grid_path = root / "output/s02_sia/s02_arctic_region/p21_grouped_mld/n060_mldobs_product_audit/march_reconstructed_grid.npz"
climate_path = Path(r"F:\sea\mldobs\mld_obs\ELG\clim\MLD_clim_month_03.mat")
legacy_path = root / "output/s02_sia/s02_arctic_region/p21_grouped_mld/n061_grouped_mld_timeseries/observation_reference.csv"
land_path = Path(r"C:\Users\Lenovo\.local\share\cartopy\shapefiles\natural_earth\physical\ne_110m_land.shp")
assert not destination.exists(), "已有正式输出须先归档。"
regions = ["Labrador-Irminger", "Nordic Seas"]
minimum_years = 11
thresholds = [8, 11, 14]
# 仅用于定位两个不同海域；实际平均范围是覆盖筛选后的最大8邻域连通片。
locating_windows = [[45, 65, -70, -10], [65, 85, -40, 30]]
annual_rows, summary_rows, cell_rows = [], [], []
primary_masks = []


# %% 3. 读取三月观测重构网格与原始气候态支持字段

climate_data = loadmat(climate_path, squeeze_me=True)
assert climate_data["dimension_order"] == "longitude,latitude"
assert climate_data["units"] == "m" and climate_data["mld_obs"].shape == (360, 180)
climate = climate_data["mld_obs"].astype(float)
with np.load(grid_path) as saved:
    assert saved["product"][0] == "ELG_100m"
    years = saved["year"][45:66]
    reconstructed = saved["reconstructed_mld_m"][0, 45:66].astype(float)
    longitude, latitude = saved["longitude"], saved["latitude"]
assert years.tolist() == list(range(2005, 2026))
assert np.array_equal(longitude, climate_data["longitude"])
assert np.array_equal(latitude, climate_data["latitude"])
longitude_grid, latitude_grid = np.meshgrid((longitude + 180) % 360 - 180, latitude, indexing="ij")
area_weights = np.cos(np.deg2rad(latitude_grid)).astype(float)
valid_year_count = np.isfinite(reconstructed).sum(axis=0)
supported = (np.isfinite(climate) & (climate_data["ocean_mask"] == 1)
             & (climate_data["profile_count"] > 0) & (climate_data["first_guess_only"] == 0))


# %% 4. 按覆盖率划分固定连通范围并计算逐年区域均值

for region, window in zip(regions, locating_windows):
    lower_lat, upper_lat, lower_lon, upper_lon = window
    locating = ((latitude_grid >= lower_lat) & (latitude_grid < upper_lat)
                & (longitude_grid >= lower_lon) & (longitude_grid <= upper_lon))
    for threshold in thresholds:
        candidate = supported & locating & (valid_year_count >= threshold)
        components, component_count = label(candidate, structure=np.ones((3, 3), dtype=int))
        component_sizes = np.bincount(components.ravel())
        component_sizes[0] = 0
        assert component_count > 0
        chosen = components == int(np.argmax(component_sizes))
        # 片区按格点数量选取，与MLD深浅无关；不补洞、不填NaN、不截负值。
        chosen_area = area_weights[chosen].sum()
        for year_index, year in enumerate(years):
            field = reconstructed[year_index]
            valid = chosen & np.isfinite(field)
            assert valid.any()
            annual_rows.append({
                "region": region, "minimum_valid_years": threshold, "year": int(year),
                "mean_m": np.average(field[valid], weights=area_weights[valid]),
                "climate_same_support_m": np.average(climate[valid], weights=area_weights[valid]),
                "valid_cells": int(valid.sum()), "footprint_cells": int(chosen.sum()),
                "valid_area_fraction": area_weights[valid].sum() / chosen_area,
                "negative_area_fraction": np.average(field[valid] < 0, weights=area_weights[valid])})
        annual_table = pd.DataFrame(annual_rows)
        series = annual_table.loc[(annual_table.region == region)
                                  & (annual_table.minimum_valid_years == threshold)]
        summary_rows.append({
            "region": region, "method": "ELG", "minimum_valid_years": threshold,
            "reference_first": 2005, "reference_last": 2025, "n_marches": 21,
            "mean_m": series.mean_m.mean(), "interannual_sd_m": series.mean_m.std(ddof=1),
            "coverage_mean": series.valid_area_fraction.mean(),
            "coverage_min": series.valid_area_fraction.min(), "coverage_max": series.valid_area_fraction.max(),
            "negative_area_fraction_mean": series.negative_area_fraction.mean(),
            "footprint_cells": int(chosen.sum()), "candidate_cells": int(candidate.sum()),
            "component_count": int(component_count),
            "longitude_min": float(longitude_grid[chosen].min()), "longitude_max": float(longitude_grid[chosen].max()),
            "latitude_min": float(latitude_grid[chosen].min()), "latitude_max": float(latitude_grid[chosen].max()),
            "climatology_on_fixed_footprint_m": np.average(climate[chosen], weights=area_weights[chosen]),
            "reference_kind": "ELG climate-plus-anomaly reconstruction;observation-defined connected footprint;not raw absolute MLD"})
        if threshold == minimum_years:
            primary_masks.append(chosen)
            for i, j in np.argwhere(chosen):
                cell_rows.append({"region": region, "longitude_index": int(i), "latitude_index": int(j),
                                  "longitude": float(longitude_grid[i, j]), "latitude": float(latitude_grid[i, j]),
                                  "valid_march_years": int(valid_year_count[i, j]),
                                  "profile_count": int(climate_data["profile_count"][i, j]),
                                  "climate_m": float(climate[i, j]), "area_weight": float(area_weights[i, j])})


# %% 5. 保存主参照、门槛敏感性、地图数组和来源

destination.mkdir(parents=True)
summary = pd.DataFrame(summary_rows)
primary = summary.loc[summary.minimum_valid_years == minimum_years].copy()
annual = pd.DataFrame(annual_rows)
primary.to_csv(destination / "observation_reference.csv", index=False)
summary.to_csv(destination / "threshold_sensitivity.csv", index=False)
annual.to_csv(destination / "annual_reconstruction.csv", index=False)
pd.DataFrame(cell_rows).to_csv(destination / "selected_cells.csv", index=False)
legacy = pd.read_csv(legacy_path)
comparison = primary[["region", "mean_m", "interannual_sd_m"]].merge(
    legacy.loc[legacy.method == "ELG", ["region", "mean_m", "interannual_sd_m"]],
    on="region", suffixes=("_dense_footprint", "_legacy_rectangle"))
comparison.to_csv(destination / "comparison_with_legacy_reference.csv", index=False)
np.savez_compressed(destination / "coverage_and_footprints.npz", longitude=longitude, latitude=latitude,
                    valid_year_count=valid_year_count, supported=supported, masks=np.stack(primary_masks),
                    regions=np.array(regions), locating_windows=np.array(locating_windows))
land = {"type": "FeatureCollection", "features": []}
for geometry in Reader(land_path).geometries():
    land["features"].append({"type": "Feature", "properties": {}, "geometry": mapping(geometry)})
(destination / "map_land.geojson").write_text(json.dumps(land), encoding="utf-8")
metadata = {
    "method": "ELG only;existing reconstruction unchanged", "reference_years": [2005, 2025],
    "regions": regions, "locating_windows_latmin_latmaxexclusive_lonmin_lonmax": locating_windows,
    "selection": "At least11 of21 March reconstructions finite;finite climate,ocean_mask1,profile_count>0,first_guess_only0;largest8-neighbor component by cell count in each broad locating window;no user polygon/legacy rectangle read and no MLD-depth selection",
    "minimum_valid_years": minimum_years, "sensitivity_thresholds": thresholds,
    "mean": "Each year:cos(latitude) weighted over available cells of the fixed selected footprint;then21 annual values weighted equally",
    "spread": "Sample SD of21 processed annual regional reconstructions;not CI or measurement uncertainty;annual available support still varies",
    "limits": ["Not raw absolute MLD;upstream reconstruction caveats remain",
               "Observation footprints differ from unchanged legacy model rectangles;not a matched-support skill comparison",
               "Locating windows separate broad sea vicinities and are analyst choices,not data-derived sea boundaries",
               "11-year coverage is a transparent working rule,not an established validated optimum;8 and14 sensitivity retained",
               "Valid-year count concerns processed gridded output,not21 independent direct profiles;profile_count only screens climate support",
               "Main western footprint may extend beyond the strict Labrador and Irminger seas;map and cell list are authoritative"],
    "dimensions": {"source_reconstruction": "product4,March_year66,longitude360,latitude180;ELG_100m index0,years45:66",
                   "MAT_climate_and_support": "longitude360,latitude180;unitsm for climate", "masks": "region2,longitude360,latitude180"},
    "changes_to_model_series": "None",
    "sources": [{"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                for path in [grid_path, climate_path, legacy_path, land_path]]}
(destination / "metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
print(primary[["region", "footprint_cells", "mean_m", "interannual_sd_m", "coverage_mean"]].round(3).to_string(index=False))
print(summary[["region", "minimum_valid_years", "footprint_cells", "mean_m", "interannual_sd_m"]].round(3).to_string(index=False))
