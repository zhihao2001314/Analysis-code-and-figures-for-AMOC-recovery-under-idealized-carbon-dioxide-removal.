# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.io import loadmat
from shapely.geometry import shape
from shapely import contains_xy


# %% 2. 参数与路径设置

root = Path(r"F:\sea\Antic sea ice constraint")
input_directory = Path(r"F:\sea\mldobs\mld_obs")
output_directory = root / "output/s02_sia/s02_arctic_region/p21_grouped_mld/n060_mldobs_product_audit"
geometry_path = root / "output/s02sia/s02_arctic_region/p14_mld_fixed_regions_raw/mld_coast_regions_20261012.geojson"
assert not output_directory.exists(), "正式输出已存在，须先归档。"
methods = ["ELG", "second_derivative"]
anchors = [100, 300]
baseline_first, baseline_last = 2005, 2025
expected_dates = pd.date_range("1960-01-01", "2025-12-01", freq="MS")
reference = loadmat(input_directory / "ELG/clim/MLD_clim_month_03.mat", squeeze_me=True)
longitude = reference["longitude"]
latitude = reference["latitude"]
assert reference["dimension_order"] == "longitude,latitude"
assert reference["mld_obs"].shape == (360, 180)
longitude_grid, latitude_grid = np.meshgrid((longitude + 180) % 360 - 180, latitude, indexing="ij")
area_weights = np.cos(np.deg2rad(latitude_grid))
features = json.loads(geometry_path.read_text(encoding="utf-8"))["features"]
region_masks = {}
for feature in features:
    region_masks[feature["properties"]["region"]] = contains_xy(
        shape(feature["geometry"]), longitude_grid, latitude_grid)
climate_rows, file_rows, regional_rows, baseline_rows, comparison_rows = [], [], [], [], []
metadata_samples = {}
march_arrays = []
march_names = []
march_years = np.arange(1960, 2026)
output_directory.mkdir(parents=True)


# %% 3. 核查气候态及其区域支持

climate_fields = {}
for method in methods:
    climate_fields[method] = []
    for month in range(1, 13):
        path = input_directory / method / "clim" / f"MLD_clim_month_{month:02d}.mat"
        data = loadmat(path, squeeze_me=True)
        climate = data["mld_obs"].astype(float)
        assert data["dimension_order"] == "longitude,latitude" and data["units"] == "m"
        assert climate.shape == (360, 180)
        assert np.array_equal(data["longitude"], longitude)
        assert np.array_equal(data["latitude"], latitude)
        assert np.array_equal(data["climatology_years"], [baseline_first, baseline_last])
        assert np.array_equal(climate, data["MLD"], equal_nan=True)
        climate_fields[method].append(climate)
        for region, mask in region_masks.items():
            valid = mask & np.isfinite(climate) & (data["ocean_mask"] == 1)
            weights = area_weights[valid]
            climate_rows.append({
                "method": method, "month": month, "region": region,
                "mean_climate_m": np.average(climate[valid], weights=weights),
                "valid_cells": int(valid.sum()),
                "first_guess_only_area_fraction": np.average(data["first_guess_only"][valid] == 1, weights=weights),
                "zero_profile_count_area_fraction": np.average(data["profile_count"][valid] == 0, weights=weights),
                "nearest_observation_distance_mean_km": np.average(data["nearest_observation_distance_km"][valid], weights=weights)})
        if month == 3:
            scalar_fields = {}
            for key, value in data.items():
                if key.startswith("__") or np.size(value) > 30:
                    continue
                scalar_fields[key] = value.tolist() if isinstance(value, np.ndarray) else value
            metadata_samples[method + "_climatology"] = scalar_fields


# %% 4. 核查全部月文件，并以同月份气候态重构

for method in methods:
    for anchor in anchors:
        paths = sorted((input_directory / method / f"depth{anchor}m").glob("*.mat"))
        assert len(paths) == len(expected_dates) == 792
        baseline_sum = np.zeros((12, 360, 180), dtype=float)
        baseline_count = np.zeros((12, 360, 180), dtype=np.int16)
        march_reconstructed = np.full((66, 360, 180), np.nan, dtype=np.float32)
        for path, date in zip(paths, expected_dates):
            data = loadmat(path, squeeze_me=True)
            assert (int(data["year"]), int(data["month"])) == (date.year, date.month)
            assert data["dimension_order"] == "longitude,latitude" and data["units"] == "m"
            assert data["mld_obs"].shape == (360, 180)
            assert np.array_equal(data["longitude"], longitude)
            assert np.array_equal(data["latitude"], latitude)
            assert np.array_equal(data["climatology_years"], [baseline_first, baseline_last])
            assert data["field_not_original_absolute_MLD"] == 1
            assert data["month_centering_applied"] == 1
            assert int(data["anchor_depth_m"]) == anchor
            assert np.array_equal(data["mld_obs"], data["MLD"], equal_nan=True)
            anomaly = data["mld_obs"].astype(float)
            climate = climate_fields[method][date.month - 1]
            # 仅重建基准，不声称等于原始绝对MLD；负值保留并计数，不截零、不补缺测。
            reconstructed = climate + anomaly
            finite = np.isfinite(reconstructed)
            values = reconstructed[finite]
            file_rows.append({
                "method": method, "anchor_depth_m": anchor, "year": date.year, "month": date.month,
                "time_index": int(data["time_index"]), "path": str(path),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "anomaly_valid_cells": int(np.isfinite(anomaly).sum()),
                "anomaly_negative_cells": int((anomaly < 0).sum()),
                "anomaly_zero_cells": int((anomaly == 0).sum()),
                "reconstructed_valid_cells": int(finite.sum()),
                "reconstructed_negative_cells": int((reconstructed < 0).sum()),
                "reconstructed_min_m": float(values.min()) if len(values) else np.nan,
                "reconstructed_max_m": float(values.max()) if len(values) else np.nan})
            if baseline_first <= date.year <= baseline_last:
                baseline_sum[date.month - 1] += np.where(np.isfinite(anomaly), anomaly, 0)
                baseline_count[date.month - 1] += np.isfinite(anomaly)
            if date.month == 3:
                march_reconstructed[date.year - 1960] = reconstructed
                for region, mask in region_masks.items():
                    support = mask & np.isfinite(climate)
                    valid = support & finite
                    weights = area_weights[valid]
                    regional_rows.append({
                        "method": method, "anchor_depth_m": anchor, "year": date.year, "month": 3,
                        "region": region, "valid_cells": int(valid.sum()),
                        "climate_supported_cells": int(support.sum()),
                        "valid_area_fraction": float(weights.sum() / area_weights[support].sum()),
                        "reconstructed_mean_m": np.average(reconstructed[valid], weights=weights) if valid.any() else np.nan,
                        "anomaly_mean_m": np.average(anomaly[valid], weights=weights) if valid.any() else np.nan,
                        "climate_on_valid_cells_mean_m": np.average(climate[valid], weights=weights) if valid.any() else np.nan,
                        "negative_cells": int((valid & (reconstructed < 0)).sum()),
                        "negative_area_fraction": np.average(reconstructed[valid] < 0, weights=weights) if valid.any() else np.nan})
            if date == expected_dates[0]:
                scalar_fields = {}
                for key, value in data.items():
                    if key.startswith("__") or np.size(value) > 30:
                        continue
                    scalar_fields[key] = value.tolist() if isinstance(value, np.ndarray) else value
                metadata_samples[f"{method}_{anchor}m"] = scalar_fields
        baseline_mean = np.divide(baseline_sum, baseline_count,
                                  out=np.full_like(baseline_sum, np.nan), where=baseline_count > 0)
        for month in range(1, 13):
            selected = baseline_count[month - 1] > 0
            baseline_rows.append({"method": method, "anchor_depth_m": anchor, "month": month,
                                  "valid_cells": int(selected.sum()),
                                  "max_abs_baseline_anomaly_mean_m": float(np.max(np.abs(baseline_mean[month - 1][selected]))),
                                  "cells_with_all_21_years": int((baseline_count[month - 1] == 21).sum())})
        march_arrays.append(march_reconstructed)
        march_names.append(f"{method}_{anchor}m")
        print(f"Audited {method}, anchor={anchor}m, 1960-01 to2025-12", flush=True)


# %% 5. 比较两种深度参数，并保存诊断

march_arrays = np.stack(march_arrays)
for method_index, method in enumerate(methods):
    shallow = march_arrays[method_index * 2]
    deeper = march_arrays[method_index * 2 + 1]
    for year_index, year in enumerate(march_years):
        for region, mask in region_masks.items():
            common = mask & np.isfinite(shallow[year_index]) & np.isfinite(deeper[year_index])
            difference = deeper[year_index][common].astype(float) - shallow[year_index][common].astype(float)
            comparison_rows.append({
                "method": method, "year": int(year), "region": region, "common_cells": int(common.sum()),
                "anchor300_minus_anchor100_mean_m": np.average(difference, weights=area_weights[common]) if common.any() else np.nan,
                "max_abs_difference_m": float(np.max(np.abs(difference))) if common.any() else np.nan})

inventory = pd.DataFrame(file_rows)
regional = pd.DataFrame(regional_rows)
inventory.to_csv(output_directory / "monthly_file_audit.csv", index=False)
pd.DataFrame(climate_rows).to_csv(output_directory / "regional_climatology.csv", index=False)
regional.to_csv(output_directory / "march_regional_reconstruction.csv", index=False)
pd.DataFrame(baseline_rows).to_csv(output_directory / "baseline_centering_check.csv", index=False)
pd.DataFrame(comparison_rows).to_csv(output_directory / "anchor_comparison_common_cells.csv", index=False)
np.savez_compressed(output_directory / "march_reconstructed_grid.npz", reconstructed_mld_m=march_arrays,
                    product=np.array(march_names), year=march_years, longitude=longitude, latitude=latitude)
report = {
    "status": "audit_completed_not_approved_as_raw_absolute_MLD",
    "time": "1960-01 to2025-12 inclusive;792 monthly files per method/anchor;confirmed by user and all MAT year/month fields",
    "dimensions": {"MAT_field": ["longitude:360", "latitude:180"],
                   "longitude": "0.5 to359.5degreesE", "latitude": "-89.5 to89.5degreesN",
                   "NPZ_reconstructed_mld_m": ["product:4", "year:66 March values", "longitude:360", "latitude:180"]},
    "formula": "reconstructed_MLD(year,month,lon,lat)=clim(method,month,lon,lat)+anomaly(method,anchor,year,month,lon,lat)",
    "important": ["Reconstruction is not recovery of raw unsmoothed absolute MLD; residual month centering was applied upstream.",
                  "100m/300m are anchor_depth_m parameters, not added constants or limits on output MLD.",
                  "Negative anomalies and zeros are valid. Negative reconstructed absolute depths are flagged, not clipped to zero.",
                  "NaN retained. Spatial coverage varies with date; regional mean climate on available cells can vary even with fixed climate.",
                  "Methods ELG and second_derivative are distinct from the model sigma0 increase0.03kg/m3 relative10m definition.",
                  "Upstream implementation, QC thresholds, exact ELG algorithm, and No04 time-window code not supplied in this folder; metadata is evidence, not executable verification of upstream processing.",
                  "Attached time-window table is contextual only; its units and exact match to this V6 product require upstream source confirmation."],
    "region_geometry": str(geometry_path), "geometry_sha256": hashlib.sha256(geometry_path.read_bytes()).hexdigest(),
    "region_weighting": "same user-drawn coastal polygons;1degree cell centers;cos(latitude);valid reconstruction cells only;no depth threshold;coverage and negative fractions saved",
    "metadata_samples": metadata_samples,
    "summary": inventory.groupby(["method", "anchor_depth_m"])[["reconstructed_valid_cells", "reconstructed_negative_cells"]].sum().reset_index().to_dict(orient="records")}
(output_directory / "metadata.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps(report["summary"], ensure_ascii=False, indent=2))
