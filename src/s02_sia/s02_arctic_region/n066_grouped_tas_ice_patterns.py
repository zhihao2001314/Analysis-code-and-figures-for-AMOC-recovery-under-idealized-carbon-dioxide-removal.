# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import pandas as pd
from shapely import contains_xy
from shapely.geometry import shape


# %% 2. 参数与路径设置

root = Path(r"F:\sea\Antic sea ice constraint")
destination = root / "output/s02_sia/s02_arctic_region/p22_grouped_temperature_ice/n066_grouped_tas_ice_patterns"
membership_path = root / "output/s02_sia/s02_arctic_region/p20_grouped_initial_stratification/n059_grouped_initial_n2_profiles/membership.csv"
ice_path = root / "data/s02_sia/s02_arctic_region/arctic_sia_pattern_values_1deg.csv"
temperature_directory = root / "output/s03tem_n01"
land_path = root / "output/s02_sia/s02_arctic_region/p21_grouped_mld/n063_observation_dense_mld/map_land.geojson"
assert not destination.exists() or not list(destination.iterdir()), "正式输出已存在，须先归档。"
group_order = ["better_nonrecovery", "better_recovery", "larger_bias"]
group_titles = ["Better-simulated stratification\nNon-recovery (n = 2)",
                "Better-simulated stratification\nRecovery (n = 2)",
                "Larger stratification bias\n(n = 4)"]
latitude = np.arange(45.5, 90.0, 1.0)
longitude = np.arange(-179.5, 180.0, 1.0)
longitude_grid, latitude_grid = np.meshgrid(longitude, latitude)
grid_index = pd.MultiIndex.from_product([latitude, longitude],
                                      names=["latitude_deg_north", "longitude_deg_east"])
ice_threshold = 0.15
temperature_color_limit = 7.0
source_paths = [membership_path, ice_path, land_path,
                temperature_directory / "arctic_tas_pattern_summary.csv",
                temperature_directory / "arctic_tas_pattern_qc.csv"]


# %% 3. 读取固定分组及既有端点场

original_membership = pd.read_csv(membership_path)
membership = original_membership[["model", "group", "panel_index", "member_order"]].copy()
membership = membership.sort_values(["panel_index", "member_order"]).reset_index(drop=True)
assert membership.groupby("group").size().reindex(group_order).tolist() == [2, 2, 4]
assert "CAS-ESM2-0" not in membership.model.tolist()
models = membership.model.to_numpy(dtype=str)
temperature_summary = pd.read_csv(temperature_directory / "arctic_tas_pattern_summary.csv").set_index("model")
temperature_qc = pd.read_csv(temperature_directory / "arctic_tas_pattern_qc.csv").set_index("model")
ice_table = pd.read_csv(ice_path, usecols=["model", "latitude_deg_north", "longitude_deg_east",
                                        "edge_raw_co2_sic_fraction", "edge_raw_cdr_sic_fraction"])
ice_table["longitude_deg_east"] = (ice_table.longitude_deg_east + 180) % 360 - 180
member_temperature = []
member_initial_ice = []
member_later_ice = []
member_rows = []
for model in models:
    path = temperature_directory / f"arctic_tas_{model}_pattern_1deg.npz"
    source_paths.append(path)
    with np.load(path) as saved:
        assert saved["lat"].shape == saved["lon"].shape == (45, 360)
        assert np.array_equal(saved["lat"], latitude_grid)
        assert np.array_equal(saved["lon"], longitude_grid)
        temperature = saved["raw_difference"].astype(float)
        np.testing.assert_allclose(temperature, saved["cdr_final"] - saved["co2_initial"],
                                   atol=1e-10, rtol=0, equal_nan=True)
    assert temperature_summary.loc[model, "co2_window_year_index"] == "0:10"
    assert temperature_summary.loc[model, "cdr_window_year_index"] == "130:140"
    assert temperature_summary.loc[model, "tas_vertical_definition"] == "standard_tas"
    assert temperature_qc.loc[model, "co2_start_step"] == 0
    assert temperature_qc.loc[model, "co2_stop_step_exclusive"] == 120
    assert temperature_qc.loc[model, "cdr_start_step"] == 1560
    assert temperature_qc.loc[model, "cdr_stop_step_exclusive"] == 1680
    model_ice = ice_table.loc[ice_table.model.eq(model)].set_index(
        ["latitude_deg_north", "longitude_deg_east"])
    assert model_ice.index.is_unique
    model_ice = model_ice.reindex(grid_index)
    # 恢复原CSV的实际支持；缺失格点不沿用旧绘图中填0的处理。
    initial = model_ice.edge_raw_co2_sic_fraction.to_numpy(float).reshape(45, 360)
    later = model_ice.edge_raw_cdr_sic_fraction.to_numpy(float).reshape(45, 360)
    for field in [initial, later]:
        assert np.nanmin(field) >= 0 and np.nanmax(field) <= 1
    member_temperature.append(temperature)
    member_initial_ice.append(initial)
    member_later_ice.append(later)
    member_rows.append({"model": model, "tas_finite_cells": int(np.isfinite(temperature).sum()),
                        "initial_sic_finite_cells": int(np.isfinite(initial).sum()),
                        "later_sic_finite_cells": int(np.isfinite(later).sum()),
                        "tas_minimum_degC": float(np.nanmin(temperature)),
                        "tas_maximum_degC": float(np.nanmax(temperature))})

member_temperature = np.stack(member_temperature)
member_initial_ice = np.stack(member_initial_ice)
member_later_ice = np.stack(member_later_ice)
land_json = json.loads(land_path.read_text(encoding="utf-8"))
land_mask = np.zeros((45, 360), dtype=bool)
for feature in land_json["features"]:
    land_mask |= contains_xy(shape(feature["geometry"]), longitude_grid, latitude_grid)


# %% 4. 模式等权复合，保持固定成员支持

group_temperature = []
group_initial_ice = []
group_later_ice = []
diagnostics = []
for group in group_order:
    selected = membership.group.eq(group).to_numpy()
    # 普通mean保留任一成员的NaN，不用nanmean改变逐格点入组数。
    temperature = np.mean(member_temperature[selected], axis=0)
    initial = np.mean(member_initial_ice[selected], axis=0)
    later = np.mean(member_later_ice[selected], axis=0)
    paired_support = np.isfinite(initial) & np.isfinite(later)
    initial = np.where(paired_support, initial, np.nan)
    later = np.where(paired_support, later, np.nan)
    group_temperature.append(temperature)
    group_initial_ice.append(initial)
    group_later_ice.append(later)
    diagnostics.append({"group": group, "models": ";".join(models[selected]),
                        "n_models": int(selected.sum()),
                        "tas_finite_cells": int(np.isfinite(temperature).sum()),
                        "tas_minimum_degC": float(np.nanmin(temperature)),
                        "tas_maximum_degC": float(np.nanmax(temperature)),
                        "tas_outside_color_limits_cells": int(np.sum(np.abs(temperature) > temperature_color_limit)),
                        "sic_paired_support_cells": int(paired_support.sum()),
                        "initial_sic_at_least_15pct_cells": int(np.sum(initial >= ice_threshold)),
                        "later_sic_at_least_15pct_cells": int(np.sum(later >= ice_threshold))})

group_temperature = np.stack(group_temperature)
group_initial_ice = np.stack(group_initial_ice)
group_later_ice = np.stack(group_later_ice)


# %% 5. 保存计算结果与来源记录

destination.mkdir(parents=True, exist_ok=True)
np.savez_compressed(destination / "grouped_tas_ice_fields.npz", longitude=longitude, latitude=latitude,
                    model=models, group=np.array(group_order), group_title=np.array(group_titles),
                    member_tas_change_degC=member_temperature,
                    member_initial_sic_fraction=member_initial_ice,
                    member_later_sic_fraction=member_later_ice,
                    group_tas_change_degC=group_temperature,
                    group_initial_sic_fraction=group_initial_ice,
                    group_later_sic_fraction=group_later_ice,
                    initial_sic_for_contour=np.where(land_mask[None], np.nan, group_initial_ice),
                    later_sic_for_contour=np.where(land_mask[None], np.nan, group_later_ice),
                    display_land_mask=land_mask)
membership.to_csv(destination / "membership.csv", index=False)
pd.DataFrame(member_rows).to_csv(destination / "member_diagnostics.csv", index=False)
pd.DataFrame(diagnostics).to_csv(destination / "group_diagnostics.csv", index=False)
shutil.copyfile(land_path, destination / "map_land.geojson")
metadata = {
    "definition": "Annual mean;1pctCO2-cdr years131-140 minus1pctCO2 years1-10;nominal concatenated years271-280 minus1-10",
    "temperature": "Standard near-surface air temperature tas;raw_difference,not corrected_difference;land and ocean both displayed",
    "ice": "Existing raw initial/later annual-mean SIC fields from original three-pattern source CSV;15% contour after equal-model group averaging",
    "grouping": "Fixed user-approved2/2/4 membership copied fromn059;not a new skill or recovery test;third group not asserted universally recovered",
    "aggregation": "Equal-model arithmetic mean;all fixed members must be finite;initial/later SIC further restricted to paired support",
    "missing_data": "Source missing SIC retained asNaN;no zero filling;missing TAS distinct from zero change",
    "display_land_mask": "Natural Earth110m land polygons inherited fromn063;only contour display fields masked;TAS remains visible overland",
    "smoothing": "No new temporal or spatial smoothing;no new interpolation;contours connect existing1degree fields",
    "source_grid_history": "Retains original temperature and SIC remapping from earlier pattern products;not raw native-grid contours",
    "adjustment": "none;raw forced differences used;upstream raw TAS validity mask also required piControl support,not drift subtraction",
    "color_limits_degC": [-temperature_color_limit, temperature_color_limit],
    "ice_threshold_fraction": ice_threshold,
    "dimensions": {"source_TAS_NPZ_fields": "latitude45,longitude360;lat/lon coordinates2D45x360",
                   "source_SIC_CSV": "model,latitude_deg_north,longitude_deg_east;explicit coordinate reindex;no inferred array axes",
                   "output_member_fields": "model8,latitude45,longitude360",
                   "output_group_fields": "group3,latitude45,longitude360",
                   "coordinates": "lat45.5..89.5N,lon-179.5..179.5E;1degree"},
    "limits": ["Annual15%SIC contour is not March ice edge,maximum extent,or a sea-ice-area change map",
               "Contour of group-mean SIC is not mean contour position or a contour shared by all individual models",
               "No observation contour requested or added;no significance test or inferential uncertainty shown",
               "Three small groups are descriptive composites;spatial juxtaposition alone does not establish causality",
               "GFDL6-15 exception belonged to earlierN2 profiles,not reused for these temperature/ice endpoints"],
    "sources": [{"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()} for path in source_paths]}
(destination / "metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
print(pd.DataFrame(diagnostics).to_string(index=False))
