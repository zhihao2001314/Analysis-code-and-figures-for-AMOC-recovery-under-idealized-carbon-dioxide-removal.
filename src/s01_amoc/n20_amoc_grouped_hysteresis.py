# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd


# %% 2. 参数与路径设置

root = Path(r"F:\sea\Antic sea ice constraint")
source_path = root / "output/amoc_s01_n01/l_40/CMIP6_AMOC_9models_annual_series.csv"
source_metadata_path = root / "output/amoc_s01_n01/l_40/CMIP6_AMOC_9models_processing_metadata.csv"
membership_path = root / "output/s05_atmosphere/p05_amoc_recovery_composites/n01_endpoint_groups/membership.csv"
output_directory = root / "output/s01_amoc/n20_amoc_grouped_hysteresis"
assert not output_directory.exists(), "正式结果已存在；重新运行前须先归档。"

models = ["ACCESS-ESM1-5", "CanESM5", "MIROC-ES2L", "GFDL-ESM4",
          "NorESM2-LM", "CESM2", "CNRM-ESM2-1", "UKESM1-0-LL"]
group_names = ["high_initial_stratification", "amoc_recovered", "amoc_not_recovered"]
group_titles = ["High initial stratification", "AMOC recovered", "AMOC not recovered"]
group_descriptions = ["Initial stratification > 1.5 times the observational reference",
                      "Endpoint return: late / initial AMOC >= 0.8",
                      "Persistent deficit: late / initial AMOC < 0.8"]
first_year, peak_year, final_year = 1, 140, 280
smoothing_windows = [1, 5, 11, 21]
display_window = 11
initial_years = np.arange(1, 11)
late_years = np.arange(271, 281)
arrow_year_pairs = [(45, 55), (100, 110), (165, 175), (220, 230)]


# %% 3. 读取年序列与既有分组

source = pd.read_csv(source_path)
source_metadata = pd.read_csv(source_metadata_path).set_index("Model")
membership = pd.read_csv(membership_path).set_index("model").loc[models]
# UKESM仅为这次八模式展示恢复原类别，不修改四模式大气复合成员表。
membership["figure_group"] = membership["group_before_user_exclusion"]
years = np.arange(first_year, final_year + 1)
strength = source.pivot(index="Simulation year", columns="Model",
                        values="Forced AMOC strength (Sv)").loc[years, models]
co2_table = source.pivot(index="Simulation year", columns="Model",
                         values="Atmospheric CO2 (ppm)").loc[years, models]
co2 = co2_table[models[0]]
assert np.isfinite(strength.to_numpy()).all()
assert np.allclose(co2_table.to_numpy(), co2.to_numpy()[:, None], rtol=0, atol=1e-7)
assert co2.loc[1:140].is_monotonic_increasing
assert co2.loc[140:280].is_monotonic_decreasing


# %% 4. 计算端点与平滑序列

initial_mean = strength.loc[initial_years].mean()
late_mean = strength.loc[late_years].mean()
assert np.allclose(initial_mean, membership.amoc_early_sv, atol=1e-7)
assert np.allclose(late_mean, membership.amoc_later_sv, atol=1e-7)
plot_rows, paired_rows, metric_rows, arrow_rows, area_rows = [], [], [], [], []
group_rows = []
panel_letter = 0

for group_index, group_name in enumerate(group_names):
    group_models = membership.index[membership.figure_group == group_name]
    group_rows.append({"group": group_name, "column": group_index,
                       "title": group_titles[group_index],
                       "description": group_descriptions[group_index],
                       "n_models": len(group_models)})

    for row_index, model in enumerate(group_models):
        panel_letter += 1
        raw = strength[model]
        smooth = raw.rolling(display_window, center=True, min_periods=1).mean()
        model_initial = initial_mean[model]
        for year in years:
            plot_rows.append({"model": model, "group": group_name,
                              "panel": chr(96 + panel_letter),
                              "row": row_index, "column": group_index,
                              "simulation_year": year, "co2_ppm": co2.loc[year],
                              "phase": "increase" if year <= peak_year else "decrease",
                              "amoc_raw_sv": raw.loc[year], "amoc_smooth_sv": smooth.loc[year],
                              "amoc_smooth_percent_initial": 100 * smooth.loc[year] / model_initial})

        for arrow_number, (start_year, end_year) in enumerate(arrow_year_pairs):
            arrow_rows.append({"model": model, "arrow": arrow_number,
                               "phase": "increase" if end_year <= peak_year else "decrease",
                               "start_year": start_year, "end_year": end_year,
                               "start_co2_ppm": co2.loc[start_year],
                               "end_co2_ppm": co2.loc[end_year],
                               "start_amoc_sv": smooth.loc[start_year],
                               "end_amoc_sv": smooth.loc[end_year]})

        # 相同CO₂配对：升高第k年对应降低阶段合并第280-k年。
        # 第140年为共享转折点；没有补造升高第0年的AMOC观测。
        up_years = np.arange(1, 141)
        down_years = final_year - up_years
        paired_co2 = co2.loc[up_years].to_numpy()
        assert np.allclose(paired_co2, co2.loc[down_years].to_numpy(), atol=1e-7)
        for year_up, year_down, concentration in zip(up_years, down_years, paired_co2):
            paired_rows.append({"model": model, "up_year": year_up, "down_year": year_down,
                                "co2_ppm": concentration,
                                "up_raw_sv": raw.loc[year_up], "down_raw_sv": raw.loc[year_down],
                                "up_smooth_sv": smooth.loc[year_up],
                                "down_smooth_sv": smooth.loc[year_down],
                                "gap_smooth_sv": smooth.loc[year_down] - smooth.loc[year_up]})

        for window in smoothing_windows:
            filtered = raw.rolling(window, center=True, min_periods=1).mean()
            up_values = filtered.loc[up_years].to_numpy()
            down_values = filtered.loc[down_years].to_numpy()
            gaps = down_values - up_values
            # 插入两段折线的交点后，对|差值|做梯形积分，避免跨零区间面积高估。
            crossing_indices = np.flatnonzero(gaps[:-1] * gaps[1:] < 0)
            crossing_fraction = -gaps[crossing_indices] / (gaps[crossing_indices + 1] - gaps[crossing_indices])
            crossing_co2 = paired_co2[crossing_indices] + crossing_fraction * np.diff(paired_co2)[crossing_indices]
            nodes = np.sort(np.concatenate([paired_co2, crossing_co2]))
            up_nodes = np.interp(nodes, paired_co2, up_values)
            down_nodes = np.interp(nodes, paired_co2, down_values)
            area = np.trapezoid(np.abs(down_nodes - up_nodes), nodes)
            concentration_range = nodes[-1] - nodes[0]
            area_rows.append({"model": model, "group": group_name, "smoothing_years": window,
                              "common_co2_min_ppm": nodes[0], "common_co2_max_ppm": nodes[-1],
                              "hysteresis_area_sv_ppm": area,
                              "mean_absolute_branch_gap_sv": area / concentration_range,
                              "normalized_area": area / (model_initial * concentration_range),
                              "n_interior_crossings": len(crossing_indices)})
            if window == display_window:
                for concentration, up_value, down_value in zip(nodes, up_nodes, down_nodes):
                    metric_rows.append({"model": model, "co2_ppm": concentration,
                                        "up_smooth_sv": up_value, "down_smooth_sv": down_value})

endpoint = pd.DataFrame({"initial_mean_sv": initial_mean, "late_mean_sv": late_mean,
                         "late_minus_initial_sv": late_mean - initial_mean,
                         "late_over_initial": late_mean / initial_mean,
                         "late_percent_initial": 100 * late_mean / initial_mean,
                         "group": membership.figure_group,
                         "member": membership.member,
                         "actual_latitude_degN": source_metadata.loc[models, "Actual latitude (degN)"],
                         "source_variable": source_metadata.loc[models, "Variable"]})
endpoint.index.name = "model"
area_table = pd.DataFrame(area_rows)
endpoint = endpoint.join(area_table.loc[area_table.smoothing_years == display_window].set_index("model")
                         [["hysteresis_area_sv_ppm", "mean_absolute_branch_gap_sv", "normalized_area",
                           "n_interior_crossings"]])


# %% 5. 保存正式结果与元数据

output_directory.mkdir(parents=True)
pd.DataFrame(plot_rows).to_csv(output_directory / "annual_plot_data.csv", index=False, encoding="utf-8-sig")
pd.DataFrame(paired_rows).to_csv(output_directory / "same_co2_pairs.csv", index=False, encoding="utf-8-sig")
pd.DataFrame(metric_rows).to_csv(output_directory / "area_fill_nodes.csv", index=False, encoding="utf-8-sig")
pd.DataFrame(arrow_rows).to_csv(output_directory / "direction_arrows.csv", index=False, encoding="utf-8-sig")
pd.DataFrame(group_rows).to_csv(output_directory / "groups.csv", index=False, encoding="utf-8-sig")
area_table.to_csv(output_directory / "area_smoothing_sensitivity.csv", index=False, encoding="utf-8-sig")
endpoint.to_csv(output_directory / "endpoint_and_hysteresis_metrics.csv", encoding="utf-8-sig")
metadata = {
    "source": str(source_path),
    "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
    "source_metadata": str(source_metadata_path), "membership_source": str(membership_path),
    "source_value_column": "Forced AMOC strength (Sv)", "drift_correction": "none",
    "definition": "Maximum Atlantic-Arctic overturning at the latitude nearest 40N, depth >= 500 m; not a fixed-depth measurement",
    "source_dimensions": "CSV Model x Simulation year; pivot Simulation year x Model",
    "saved_dimensions": "annual_plot_data: model x simulation_year; same_co2_pairs: model x up_year; area_fill_nodes: model x irregular co2 node",
    "models": models, "excluded_model": "CAS-ESM2-0",
    "groups": "Existing hierarchical grouping: initial stratification > 1.5 observational mean first; other models split by late/initial AMOC >= 0.8",
    "UKESM_scope": "Original amoc_not_recovered category restored only for this requested eight-model display; previous four-model composites unchanged",
    "years": [first_year, final_year], "turning_point_year": peak_year,
    "co2_C0_ppm": float(co2.loc[final_year]),
    "co2_sequence": "C0*1.01**t for t=1..140; C0*1.01**(280-t) for t=141..280",
    "smoothing": "11-year centered mean on years 1..280 only; min_periods=1 (six years at endpoints); raw and 1/5/11/21-year area sensitivity retained",
    "endpoint_statistics": "Unsmoothed years 1-10 and 271-280 means, not endpoints of the smoothed curve",
    "shared_peak": "Decrease branch begins at the same measured year140 turning point; not a second independent observation",
    "same_co2_pairing": "Increasing year k=1..140 paired with combined year280-k=279..140; actual final year280 retained for curve endpoint only",
    "hysteresis_area": "Integral over common linear CO2 domain of absolute (decrease minus increase) AMOC difference; exact for piecewise linear branches after inserting crossings; units Sv ppm",
    "normalized_area": "A / [(initial years1-10 AMOC mean) * (common Cmax-Cmin)]; dimensionless mean absolute relative path separation",
    "uncertainty": "One existing member per model; no CI, significance, equilibrium irreversibility, or equivalence-to-initial test claimed",
    "classification": "Recovered is the existing >=0.8 ratio rule, not an exact Closed-loop criterion. No arbitrary large/small area thresholds applied",
    "reference": "Kim et al. (2022), Nature Climate Change, https://doi.org/10.1038/s41558-022-01452-z",
}
(output_directory / "metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
print(endpoint.round(4).to_string())
