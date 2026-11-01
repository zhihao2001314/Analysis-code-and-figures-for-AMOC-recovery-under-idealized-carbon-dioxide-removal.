# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import pandas as pd
import xarray as xr
from scipy.io import loadmat


# %% 2. 参数与路径设置

root = Path(r"F:\sea\Antic sea ice constraint")
destination = root / "output/s02_sia/s02_arctic_region/p21_grouped_mld/n065_observed_march_mld_patterns"
seanoe_path = Path(r"F:\sea\mldobs\SEANOE-obs\mld_dr003_ref10m_v2023.nc")
reconstruction_path = root / "output/s02_sia/s02_arctic_region/p21_grouped_mld/n060_mldobs_product_audit/march_reconstructed_grid.npz"
climate_path = Path(r"F:\sea\mldobs\mld_obs\ELG\clim\MLD_clim_month_03.mat")
land_path = root / "output/s02_sia/s02_arctic_region/p21_grouped_mld/n063_observation_dense_mld/map_land.geojson"
assert not destination.exists(), "已有正式输出须先归档。"
map_extent = [-75, 30, 40, 85]
color_limits_m = [0, 1400]
minimum_supported_years_for_display_flag = 11
summary_rows = []


# %% 3. 明确维序读取SEANOE三月气候态

with xr.open_dataset(seanoe_path) as dataset:
    assert dataset["mld_dr003"].dims == ("time", "lat", "lon")
    assert dataset["mask"].dims == ("lat", "lon")
    assert dataset["n_obs"].dims == ("time", "lat", "lon")
    assert dataset["mld_dr003"].attrs["unit"] == "m"
    assert dataset.time.values.tolist() == [14, 45, 73, 104, 134, 165, 195, 226, 257, 287, 318, 348]
    # time为12个月的月中心日序；索引2明确对应3月，不是逐年时间轴。
    longitude = dataset.lon.values.astype(float)
    latitude = dataset.lat.values.astype(float)
    seanoe = dataset.mld_dr003.isel(time=2).transpose("lat", "lon").values.astype(float)
    seanoe_profiles = dataset.n_obs.isel(time=2).transpose("lat", "lon").values
    seanoe_ocean = dataset["mask"].transpose("lat", "lon").values == 1
    source_metadata = {"global": dict(dataset.attrs), "mld_long_name": dataset.mld_dr003.attrs["long_name"],
                       "time_climatology": dataset.time.attrs["climatology"],
                       "fill_value": float(dataset.mld_dr003.encoding["_FillValue"])}
seanoe_masked_land_count = int(np.sum(np.isfinite(seanoe) & ~seanoe_ocean))
seanoe = np.where(seanoe_ocean, seanoe, np.nan)


# %% 4. 重构观测逐格点平均，不使用任何区域掩膜

climate_data = loadmat(climate_path, squeeze_me=True)
assert climate_data["dimension_order"] == "longitude,latitude"
assert climate_data["mld_obs"].shape == (360, 180) and climate_data["units"] == "m"
with np.load(reconstruction_path) as saved:
    assert saved["product"][0] == "ELG_100m"
    assert saved["year"][45:66].tolist() == list(range(2005, 2026))
    values = saved["reconstructed_mld_m"][0, 45:66].astype(float)
    original_longitude, original_latitude = saved["longitude"], saved["latitude"]
assert np.array_equal(original_longitude, climate_data["longitude"])
assert np.array_equal(original_latitude, climate_data["latitude"])
order = np.argsort((original_longitude + 180) % 360 - 180)
assert np.array_equal(((original_longitude + 180) % 360 - 180)[order], longitude)
assert np.array_equal(original_latitude, latitude)
year_count_original = np.isfinite(values).sum(axis=0)
# 分母仅是该格点实际有效年份数；缺测不作0，负值不截零，不做11年门槛筛选。
mean_original = np.divide(np.nansum(values, axis=0), year_count_original,
                          out=np.full((360, 180), np.nan), where=year_count_original > 0)
reconstructed = mean_original[order].T
year_count = year_count_original[order].T
climate = climate_data["mld_obs"].astype(float)[order].T
reconstruction_ocean = climate_data["ocean_mask"][order].T == 1
reconstructed = np.where(reconstruction_ocean, reconstructed, np.nan)
assert reconstructed.shape == seanoe.shape == (180, 360)


# %% 5. 保存覆盖标记和独立量级诊断

longitude_grid, latitude_grid = np.meshgrid(longitude, latitude)
domain = ((longitude_grid >= map_extent[0]) & (longitude_grid <= map_extent[1])
          & (latitude_grid >= map_extent[2]) & (latitude_grid <= map_extent[3]))
fields = np.stack([seanoe, reconstructed])
limited_support = np.stack([np.isfinite(seanoe) & (seanoe_profiles == 0),
                            np.isfinite(reconstructed) & (year_count < minimum_supported_years_for_display_flag)])
names = ["SEANOE_v2023", "ELG_reconstructed_2005_2025"]
for name, field, flag in zip(names, fields, limited_support):
    valid = domain & np.isfinite(field)
    rows, columns = np.where(valid)
    maximum_index = int(np.argmax(field[valid]))
    summary_rows.append({
        "product": name, "valid_cells_in_map_domain": int(valid.sum()),
        "minimum_m": float(field[valid].min()), "maximum_m": float(field[valid].max()),
        "max_longitude": float(longitude[columns[maximum_index]]),
        "max_latitude": float(latitude[rows[maximum_index]]),
        "negative_cells": int(np.sum(field[valid] < 0)), "zero_cells": int(np.sum(field[valid] == 0)),
        "cells_above_500m": int(np.sum(field[valid] > 500)), "cells_above_1000m": int(np.sum(field[valid] > 1000)),
        "limited_support_cells": int(np.sum(valid & flag))})
valid_reconstruction = domain & np.isfinite(reconstructed) & np.isfinite(climate)
climate_identity_error = float(np.max(np.abs(reconstructed[valid_reconstruction] - climate[valid_reconstruction])))
common = domain & np.isfinite(seanoe) & np.isfinite(reconstructed)
common_weights = np.cos(np.deg2rad(latitude_grid[common]))
common_summary = {"common_cells": int(common.sum()),
                  "seanoe_common_area_weighted_mean_m": float(np.average(seanoe[common], weights=common_weights)),
                  "reconstruction_common_area_weighted_mean_m": float(np.average(reconstructed[common], weights=common_weights)),
                  "caution": "Different observation periods and MLD algorithms;not a model bias or proof of reconstruction error"}
destination.mkdir(parents=True)
np.savez_compressed(destination / "march_patterns.npz", longitude=longitude, latitude=latitude,
                    mld_m=fields, limited_support=limited_support, product=np.array(names),
                    seanoe_profile_count=seanoe_profiles, reconstruction_year_count=year_count,
                    original_ELG_climate_m=climate, seanoe_ocean=seanoe_ocean, reconstruction_ocean=reconstruction_ocean)
pd.DataFrame(summary_rows).to_csv(destination / "pattern_summary.csv", index=False)
shutil.copyfile(land_path, destination / "map_land.geojson")
metadata = {
    "products": names, "map_extent_lonmin_lonmax_latmin_latmax": map_extent,
    "color_limits_m": color_limits_m, "month": 3,
    "seanoe": {"period": "1970-01-01 to2021-12-31 pooled monthly climatology", "variable": "mld_dr003",
               "definition": "Density increment0.03kg/m3 relative10m;monthly mapped observational profile MLD",
               "support_flag": "Black dots:finite mapped MLD but n_obs==0 in that1degree March meshbox;not zero MLD",
               "decoded_fill_value_to_NaN": True, "finite_values_masked_on_declared_land": seanoe_masked_land_count,
               "attributes": source_metadata},
    "reconstruction": {"period": "2005-2025 March;equal-year mean per cell over finite values only",
                       "definition": "Existing ELG100m climate-plus-processed-anomaly reconstruction,not raw absolute MLD",
                       "support_flag": "Black dots:1-10 of21 March values available;all such values still plotted,not excluded",
                       "negative_values": "Retained;no negative gridcell means occur in current map domain",
                       "mean_vs_original_climate_max_abs_difference_m": climate_identity_error,
                       "baseline_identity": "Near equality is expected from upstream2005-2025 calendar-month residual centering;mean map mainly reflects input climate"},
    "common_support_diagnostic": common_summary,
    "dimensions": {"SEANOE_mld_and_n_obs": "time12,lat180,lon360;mask lat180,lon360;unitm for MLD",
                   "MAT_climate": "longitude360,latitude180;unitm;explicit transpose and longitude sort",
                   "input_NPZ": "product4,year66,longitude360,latitude180;useproduct0,years45:66",
                   "output_mld_m_and_flags": "product2,latitude180,longitude360;longitude-179.5..179.5,latitude-89.5..89.5"},
    "limits": ["No old user rectangles or observation-dense masks are applied to pattern maps",
               "Each product retains its own available coverage;black-dot meanings differ and are defined separately",
               "Same display scale does not make methods,source periods or spatial support identical",
               "Not an annual time series or interannual uncertainty map;no significance testing",
               "No existing model series or observation reference figures changed"],
    "sources": [{"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                for path in [seanoe_path, reconstruction_path, climate_path, land_path]]}
(destination / "metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
print(pd.DataFrame(summary_rows).round(3).to_string(index=False))
print(json.dumps({"reconstruction_mean_vs_climate_error_m": climate_identity_error,
                  "common_support": common_summary}, ensure_ascii=False, indent=2))
