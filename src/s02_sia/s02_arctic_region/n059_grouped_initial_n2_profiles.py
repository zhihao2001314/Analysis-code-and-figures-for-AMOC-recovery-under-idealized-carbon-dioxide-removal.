# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd


# %% 2. 参数与路径设置

root = Path(r"F:\sea\Antic sea ice constraint")
source_directory = root / "output/s02sia/s02_arctic_region/p14_mld_fixed_regions_raw"
output_directory = root / "output/s02_sia/s02_arctic_region/p20_grouped_initial_stratification/n059_grouped_initial_n2_profiles"
assert not output_directory.exists(), "正式输出已存在；重新运行前须先归档。"
model_path = source_directory / "coast_region_20261012_initial_n2_profiles.csv"
observation_path = source_directory / "coast_region_20261012_iapv4_profiles.csv"
model_metadata_path = source_directory / "coast_region_20261012_initial_n2_metadata.json"
observation_metadata_path = source_directory / "coast_region_20261012_iapv4_metadata.json"
geometry_path = source_directory / "mld_coast_regions_20261012.geojson"
regions = ["Labrador-Irminger Seas", "Nordic Seas"]
group_keys = ["better_nonrecovery", "better_recovery", "larger_bias"]
group_titles = ["Better-simulated stratification\nNon-recovery",
                "Better-simulated stratification\nRecovery", "Larger stratification bias"]
group_members = [["CESM2", "CNRM-ESM2-1"], ["NorESM2-LM", "GFDL-ESM4"],
                 ["ACCESS-ESM1-5", "CanESM5", "MIROC-ES2L", "UKESM1-0-LL"]]
models = [model for group in group_members for model in group]
depths = np.arange(5.0, 1000.0, 10.0)
member_colors = {"CESM2": "#348F46", "CNRM-ESM2-1": "#CC4C49",
                 "NorESM2-LM": "#9D755D", "GFDL-ESM4": "#16857A",
                 "ACCESS-ESM1-5": "#4C78A8", "CanESM5": "#C97B00",
                 "MIROC-ES2L": "#9467BD", "UKESM1-0-LL": "#C0567D"}
member_styles = ["solid", "dashed", "dashdot", "dotted"]
membership_rows, mean_rows, metric_rows = [], [], []


# %% 3. 读取既有初始三月剖面

source_models = pd.read_csv(model_path)
source_observations = pd.read_csv(observation_path)
model_metadata = json.loads(model_metadata_path.read_text(encoding="utf-8"))
observation_metadata = json.loads(observation_metadata_path.read_text(encoding="utf-8"))
profiles = source_models.loc[source_models.model.isin(models) & source_models.region.isin(regions)].copy()
observations = source_observations.loc[source_observations.region.isin(regions)].copy()
assert len(models) == len(set(models)) == 8
assert not profiles.duplicated(["model", "region", "depth_m"]).any()
assert not observations.duplicated(["region", "depth_m"]).any()
assert np.allclose(profiles.initial_march_n2_s2 * 1e5, profiles.initial_march_n2_1e5_s2, equal_nan=True)
assert np.allclose(observations.iapv4_march_mean_n2_s2 * 1e5,
                   observations.iapv4_march_mean_n2_1e5_s2, equal_nan=True)


# %% 4. 固定成员并逐层等权平均

for group_index, (group_key, names) in enumerate(zip(group_keys, group_members)):
    for member_index, model in enumerate(names):
        source_info = model_metadata["models_metadata"][model]
        membership_rows.append({"model": model, "group": group_key, "panel_index": group_index,
                                "panel_title": group_titles[group_index], "group_size": len(names),
                                "member_order": member_index, "member": source_info["member"],
                                "color": member_colors[model], "line_style": member_styles[member_index],
                                "calendar_first": source_info["initial_first_year"],
                                "calendar_last": source_info["initial_last_year"],
                                "experiment_year_first": 6 if model == "GFDL-ESM4" else 1,
                                "experiment_year_last": 15 if model == "GFDL-ESM4" else 10,
                                "n_marches": 10,
                                "baseline_exception": "User approved reuse of years6-15; local first5 years missing" if model == "GFDL-ESM4" else "",
                                "grouping_basis": "User-specified members, including UKESM in larger-bias group; no new recovery test"})

    for region in regions:
        values = []
        for model in names:
            profile = profiles.loc[(profiles.model == model) & (profiles.region == region)].sort_values("depth_m")
            assert np.array_equal(profile.depth_m.to_numpy(), depths)
            assert profile.initial_march_n2_s2.iloc[1:].notna().all()
            assert profile.valid_march_count.iloc[1:].eq(10).all()
            assert pd.isna(profile.initial_march_n2_s2.iloc[0])
            values.append(profile.initial_march_n2_s2.to_numpy())
        values = np.asarray(values)
        # 不用nanmean替换缺失成员；每一有效深度保持2、2或4个固定成员。
        complete = np.isfinite(values).all(axis=0)
        group_mean = np.full(len(depths), np.nan)
        group_mean[complete] = values[:, complete].mean(axis=0)
        for index, depth in enumerate(depths):
            mean_rows.append({"group": group_key, "panel_index": group_index,
                              "region": region, "depth_m": depth, "n_models_expected": len(names),
                              "n_models_valid": int(np.isfinite(values[:, index]).sum()),
                              "group_mean_n2_s2": group_mean[index],
                              "group_mean_n2_1e5_s2": group_mean[index] * 1e5})

        observed = observations.loc[observations.region == region].sort_values("depth_m")
        assert np.array_equal(observed.depth_m.to_numpy(), depths)
        assert observed.iapv4_march_mean_n2_s2.notna().all()
        observed_values = observed.iapv4_march_mean_n2_s2.to_numpy()
        common = complete & np.isfinite(observed_values)
        metric_rows.append({"group": group_key, "region": region, "n_models": len(names),
                            "n_shared_depths": int(common.sum()),
                            "shared_depth_min_m": depths[common].min(),
                            "shared_depth_max_m": depths[common].max(),
                            "group_depth_mean_n2_1e5_s2": group_mean[common].mean() * 1e5,
                            "observed_depth_mean_n2_1e5_s2": observed_values[common].mean() * 1e5,
                            "profile_rmse_1e5_s2": np.sqrt(np.mean((group_mean[common] - observed_values[common]) ** 2)) * 1e5})


# %% 5. 保存正式绘图输入与口径

membership = pd.DataFrame(membership_rows)
means = pd.DataFrame(mean_rows)
profiles = profiles.merge(membership[["model", "group", "panel_index"]], on="model", validate="many_to_one")
output_directory.mkdir(parents=True)
membership.to_csv(output_directory / "membership.csv", index=False, encoding="utf-8-sig")
profiles.to_csv(output_directory / "member_profiles.csv", index=False, encoding="utf-8-sig")
means.to_csv(output_directory / "group_mean_profiles.csv", index=False, encoding="utf-8-sig")
observations.to_csv(output_directory / "observation_profiles.csv", index=False, encoding="utf-8-sig")
pd.DataFrame(metric_rows).to_csv(output_directory / "group_profile_diagnostics.csv", index=False, encoding="utf-8-sig")
sources = [model_path, observation_path, model_metadata_path, observation_metadata_path, geometry_path]
metadata = {"purpose": "Initial March N2 profiles grouped by the user's revised eight-model membership, with equal-model mean and IAPv4 observation in each panel",
            "group_sizes": [2, 2, 4], "groups": dict(zip(group_keys, group_members)),
            "group_label_caution": "Recovery/non-recovery titles refer to user-selected groups, not a new regional deep-convection recovery diagnosis; larger-bias group is not labeled recovered. Stratification skill labels are user-selected, not a statistical pass/fail test; Nordic panel uses same fixed memberships without asserting a new Nordic-skill ranking",
            "source_files": [{"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()} for path in sources],
            "dimensions": {"member_profiles": "model x region x depth_m, 8x2x100", "group_mean_profiles": "group x region x depth_m, 3x2x100", "observation_profiles": "region x depth_m, 2x100", "membership": "model,8 rows"},
            "original_raw_dimensions": "No new NetCDF/MAT read. Model thetao dimensions retained in n039 metadata; IAPv4 source N2(depth_mid,lat,lon), original10^-4 s^-2",
            "model_method": "Reuse unchanged n039 TEOS-10 N2; T/S interpolated to common geometric-depth levels every10m, without extrapolation; regional native areacello-weighted N2 for each March then equal mean of10 March years",
            "model_baseline": "Experiment years1-10 March except GFDL-ESM4 years6-15 as explicitly approved by user; ten March samples each",
            "observation": "IAPv4 March2005-2024, 20 years; same saved coastal polygons, cos(latitude) weighted regional profiles then20-year mean; native vertical interpolation already performed in n043, no new interpolation",
            "regions": {"Labrador-Irminger Seas": "Labrador-Irminger polygon in mld_coast_regions_20261012.geojson", "Nordic Seas": "Nordic Seas polygon in mld_coast_regions_20261012.geojson; Nordic-Iceland sector"},
            "group_mean": "Arithmetic equal-model mean at each saved depth, only where every fixed member is valid; mean of N2 profiles, not N2 recomputed from mean T/S, not time-sample or grid-cell pooled mean",
            "missing_depth": "All models lack N2 at5m, retainedNaN; means valid at15..995m,99 depths. Observation valid5..995m,100 depths. No 0-10m extrapolation or surface fill",
            "units": "N2 in s^-2 with separate display column scaled by1e5; depth positive downward in m",
            "drift_correction": "none", "smoothing": "none; exact source profiles retained",
            "uncertainty": "No confidence bands or significance claims. IAP temporal interannual SD retained in CSV but not plotted; it is not observational error or a CI",
            "comparison_limit": "Modern observations are a structural climatology reference for initial model states, not contemporaneous bias or a causal attribution",
            "plot_x_limits": observation_metadata["plot_x_limits_n2_1e5_s2"],
            "plot_y_limits": observation_metadata["plot_y_limits_depth_m"],
            "previous_products": "Original per-model N2 profiles, membership tables and prior atmosphere composites not edited"}
(output_directory / "metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
print(membership[["model", "group", "experiment_year_first", "experiment_year_last"]].to_string(index=False))
print(pd.DataFrame(metric_rows).round(4).to_string(index=False))
