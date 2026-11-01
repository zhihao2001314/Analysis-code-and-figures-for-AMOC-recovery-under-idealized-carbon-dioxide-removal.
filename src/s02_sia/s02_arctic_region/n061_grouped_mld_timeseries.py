# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.io import loadmat


# %% 2. 参数与路径设置

root = Path(r"F:\sea\Antic sea ice constraint")
ocean_directory = root / "output/s02sia/s02_arctic_region/p06_mld/n05_na_region_absolute_timeseries"
membership_path = root / "output/s02_sia/s02_arctic_region/p20_grouped_initial_stratification/n059_grouped_initial_n2_profiles/membership.csv"
observation_directory = root / "output/s02_sia/s02_arctic_region/p21_grouped_mld/n060_mldobs_product_audit"
output_directory = root / "output/s02_sia/s02_arctic_region/p21_grouped_mld/n061_grouped_mld_timeseries"
assert not output_directory.exists(), "正式输出已存在，须先归档。"
model_path = ocean_directory / "na_march_mld_region_absolute_timeseries.csv"
observation_path = observation_directory / "march_reconstructed_grid.npz"
climate_directory = Path(r"F:\sea\mldobs\mld_obs")
regions = ["Labrador-Irminger", "Nordic Seas"]
region_titles = ["Labrador–Irminger", "Nordic Seas"]
groups = ["better_nonrecovery", "better_recovery", "larger_bias"]
running_window = 11
first_reference_year, last_reference_year = 2005, 2025
member_rows, mean_rows, reference_rows, observation_rows = [], [], [], []
region_boxes = {"Labrador-Irminger": {"lat_min": 50, "lat_max_exclusive": 65, "lon_min": -65, "lon_max": -20},
                "Nordic Seas": {"lat_min": 65, "lat_max_inclusive": 80, "lon_min": -30, "lon_max": 20}}
region_columns = {"Labrador-Irminger": "li", "Nordic Seas": "nordic"}


# %% 3. 读取固定分组、旧矩形序列及观测重构网格

membership = pd.read_csv(membership_path)
models = pd.read_csv(model_path)
with np.load(observation_path) as saved:
    reconstructed = saved["reconstructed_mld_m"].astype(float)
    longitude = saved["longitude"]
    latitude = saved["latitude"]
    observation_years = saved["year"]
longitude_grid, latitude_grid = np.meshgrid((longitude + 180) % 360 - 180, latitude, indexing="ij")
weights = np.cos(np.deg2rad(latitude_grid))
region_masks = {"Labrador-Irminger": (latitude_grid >= 50) & (latitude_grid < 65) & (longitude_grid >= -65) & (longitude_grid <= -20),
                "Nordic Seas": (latitude_grid >= 65) & (latitude_grid <= 80) & (longitude_grid >= -30) & (longitude_grid <= 20)}
climate_paths = []
assert len(membership) == membership.model.nunique() == 8
assert membership.groupby("group").size().to_dict() == {"better_nonrecovery": 2, "better_recovery": 2, "larger_bias": 4}


# %% 4. 整理逐年原值与完整窗口平滑，计算固定成员均值

for member in membership.itertuples(index=False):
    for region in regions:
        subset = models.loc[models.model == member.model].sort_values("simulation_year")
        prefix = region_columns[region]
        # 用户要求复用旧曲线；不把沿岸多边形产品混入矩形序列，也不静默修正旧GFDL年份。
        series = subset.set_index("simulation_year")[prefix + "_march_mld_m"].reindex(range(1, 281))
        assert len(series) == 280
        assert series.notna().all()
        smoothed = subset.set_index("simulation_year")[prefix + "_march_mld_absolute_11yr_m"]
        member_rows.append(pd.DataFrame({"model": member.model, "group": member.group, "region": region,
                                         "simulation_year": series.index, "march_mld_m": series.to_numpy(),
                                         "march_mld_11yr_m": smoothed.to_numpy()}))
members = pd.concat(member_rows, ignore_index=True)
for group in groups:
    names = membership.loc[membership.group == group, "model"].tolist()
    for region in regions:
        subset = members.loc[(members.group == group) & (members.region == region)]
        matrix = subset.pivot(index="simulation_year", columns="model", values="march_mld_m").loc[:, names]
        # 固定2/2/4成员；有成员缺失时保留NaN，不以单个模式代替两模式组均值。
        mean = np.mean(matrix.to_numpy(), axis=1)
        smooth = pd.Series(mean, index=matrix.index).rolling(running_window, center=True, min_periods=running_window).mean()
        mean_rows.append(pd.DataFrame({"group": group, "region": region, "simulation_year": matrix.index,
                                       "n_members": len(names), "n_valid_members": matrix.notna().sum(axis=1),
                                       "group_mean_m": mean, "group_mean_11yr_m": smooth.to_numpy()}))
means = pd.concat(mean_rows, ignore_index=True)


# %% 5. 汇总观测重构参考，不把自然历与实验年对齐

for region in regions:
    mask = region_masks[region]
    for method_index, method in enumerate(["ELG", "second_derivative"]):
        climate_path = climate_directory / method / "clim/MLD_clim_month_03.mat"
        climate_paths.append(climate_path)
        climate_data = loadmat(climate_path, squeeze_me=True)
        assert climate_data["dimension_order"] == "longitude,latitude"
        assert climate_data["units"] == "m" and climate_data["mld_obs"].shape == (360, 180)
        assert np.array_equal(climate_data["longitude"], longitude)
        assert np.array_equal(climate_data["latitude"], latitude)
        climate = climate_data["mld_obs"].astype(float)
        support = mask & np.isfinite(climate)
        for year_index, year in enumerate(observation_years):
            field = reconstructed[method_index * 2, year_index]
            valid = support & np.isfinite(field)
            # 负重构值原样保留；无额外MLD深度筛选，不补缺测。
            observation_rows.append({"region": region, "method": method, "year": int(year),
                                     "reconstructed_mean_m": np.average(field[valid], weights=weights[valid]),
                                     "valid_cells": int(valid.sum()),
                                     "valid_area_fraction": weights[valid].sum() / weights[support].sum(),
                                     "negative_area_fraction": np.average(field[valid] < 0, weights=weights[valid])})
            if first_reference_year <= year <= last_reference_year:
                assert np.array_equal(reconstructed[method_index * 2, year_index][support],
                                      reconstructed[method_index * 2 + 1, year_index][support], equal_nan=True)
        subset = pd.DataFrame(observation_rows)
        modern = subset.loc[(subset.region == region) & (subset.method == method) & subset.year.between(first_reference_year, last_reference_year)]
        assert len(modern) == 21
        reference_rows.append({
            "region": region, "method": method, "reference_first": first_reference_year,
            "reference_last": last_reference_year, "n_marches": len(modern),
            "mean_m": modern.reconstructed_mean_m.mean(),
            "interannual_sd_m": modern.reconstructed_mean_m.std(ddof=1),
            "coverage_mean": modern.valid_area_fraction.mean(),
            "negative_area_fraction_mean": modern.negative_area_fraction.mean(),
            "climatology_full_region_m": np.average(climate[support], weights=weights[support]),
            "anchor_choice": "100m;300m yields identical modern regional March values",
            "reference_kind": "climatology_plus_anomaly_reconstruction;valid-cell regional mean;not raw absolute observation"})
references = pd.DataFrame(reference_rows)
axis_limits = {}
for region in regions:
    plotted_values = members.loc[members.region == region, "march_mld_11yr_m"]
    reference = references.loc[references.region == region]
    largest = max(float(plotted_values.max()), float((reference.mean_m + reference.interannual_sd_m).max()))
    upper = float(np.ceil(largest / 100) * 100)
    axis_limits[region] = {"y_limits_m": [0, upper], "y_ticks_m": np.arange(0, upper + 1, 100).tolist()}


# %% 6. 保存正式计算与来源记录

output_directory.mkdir(parents=True)
membership.to_csv(output_directory / "membership.csv", index=False)
members.to_csv(output_directory / "member_timeseries.csv", index=False)
means.to_csv(output_directory / "group_mean_timeseries.csv", index=False)
references.to_csv(output_directory / "observation_reference.csv", index=False)
pd.DataFrame(observation_rows).to_csv(output_directory / "observation_box_march_annual.csv", index=False)
metadata = {
    "regions": regions, "region_titles": region_titles, "groups": groups,
    "model_definition": "User-selected legacy rectangular-region series:March native-area-weighted mean,resolved sigma0+0.03kg/m3 crossings;nearest valid level to10m reference,not exact10m interpolation;no anomaly or drift correction",
    "region_boxes": region_boxes,
    "series": "Legacy experiment indices1-280 as saved;nominal ramp-up1-140 and CDR141-280;reuse saved11year centered means,min_periods11;no endpoint padding",
    "group_mean": "Equal model weights within fixed2/2/4 membership;no skipna or time-varying membership",
    "gfdl": "Known legacy timing limitation retained for exact historical-curve comparison:old index-based alignment does not establish correct physical years;not corrected with n055 coastal products because masks differ. No new quantitative timing inference is validated here.",
    "observation": "Re-average n060 climate-plus-anomaly March grids within the SAME legacy rectangles;2005-2025 means plus/minus sample interannual SD;modern anchor100/300 equal;one reference per method",
    "observation_limits": "Only contextual reference:coverage roughly28%Nordic and44%Labrador-Irminger,negative reconstructed cells retained,ELG/second_derivative differ from model density-threshold method;SD not measurement uncertainty or CI;no calendar-to-experiment alignment",
    "reference_placement": "Observation means and SD annotated on shared physical depth scale,not tied to experiment time;current display placement recorded in figure_manifest.json",
    "axis_limits": axis_limits, "dimensions": {"members": "model8 x region2 x simulation_year280", "means": "group3 x region2 x simulation_year280", "references": "method2 x region2"},
    "sources": [{"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                for path in [membership_path, model_path, observation_path] + list(dict.fromkeys(climate_paths))]}
(output_directory / "metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
print(references[["region", "method", "mean_m", "interannual_sd_m"]].round(2).to_string(index=False))
print(axis_limits)
