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
output_directory = root / "output/s01_amoc/n22_amoc_recovery_definition_audit"
assert not output_directory.exists(), "正式输出已存在；重新执行前须先归档。"
models = ["ACCESS-ESM1-5", "CanESM5", "CAS-ESM2-0", "GFDL-ESM4", "NorESM2-LM",
          "UKESM1-0-LL", "MIROC-ES2L", "CNRM-ESM2-1", "CESM2"]
short_names = ["ACCESS", "CanESM5", "CAS", "GFDL", "Nor", "UKESM", "MIROC", "CNRM", "CESM2"]
years = np.arange(1, 281)
peak_year, final_year = 140, 280
initial_window = 10
late_windows = [10, 20, 30]
smoothing_windows = [1, 5, 11, 21]
relative_thresholds = [0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00]
absolute_tolerances_sv = [1.0, 2.0, 3.0, 4.0, 5.0]
reference_relative_threshold = 0.80
endpoint_rows, window_rows, smoothing_rows = [], [], []
relative_rows, absolute_rows, return_rows, lag_threshold_rows = [], [], [], []
normalized_curves = []


# %% 3. 读取九模式原值

source = pd.read_csv(source_path)
source_metadata = pd.read_csv(source_metadata_path).set_index("Model")
strength = source.pivot(index="Simulation year", columns="Model",
                        values="Forced AMOC strength (Sv)").loc[years, models]
forcing = source.pivot(index="Simulation year", columns="Model",
                       values="Atmospheric CO2 (ppm)").loc[years, models]
assert strength.notna().all().all()
assert np.allclose(forcing.to_numpy(), forcing[models[0]].to_numpy()[:, None], rtol=0, atol=1e-7)
co2 = forcing[models[0]]
initial_means = strength.loc[1:10].mean()
up_years = np.arange(1, 141)
down_years = 280 - up_years
co2_nodes = co2.loc[up_years].to_numpy()
assert np.allclose(co2_nodes, co2.loc[down_years], atol=1e-7)


# %% 4. 计算端点、最低点与路径指标

for model_index, model in enumerate(models):
    raw = strength[model]
    initial_mean = initial_means[model]
    late_mean = raw.loc[271:280].mean()
    ratio = late_mean / initial_mean
    residual = late_mean - initial_mean
    smoothed = raw.rolling(11, center=True, min_periods=1).mean()
    minimum_year = int(smoothed.idxmin())
    minimum_mean = smoothed.loc[minimum_year]
    maximum_loss = initial_mean - minimum_mean
    # R为剩余/恢复到的初始强度比例；F才是最大损失被补回的比例，均不截断到100%。
    regained_fraction = (late_mean - minimum_mean) / maximum_loss

    endpoint = {"model": model, "short_name": short_names[model_index],
                "initial_mean_sv": initial_mean, "late_mean_sv": late_mean,
                "late_minus_initial_sv": residual, "absolute_endpoint_difference_sv": abs(residual),
                "endpoint_ratio": ratio, "endpoint_percent_initial": 100 * ratio,
                "regained_fraction_of_maximum_loss": regained_fraction,
                "minimum_11yr_sv": minimum_mean, "minimum_11yr_year": minimum_year,
                "minimum_lag_from_co2_peak_years": minimum_year - peak_year,
                "initial_annual_sd_sv": raw.loc[1:10].std(ddof=1),
                "late_annual_sd_sv": raw.loc[271:280].std(ddof=1),
                "last_decade_minus_previous_decade_sv": late_mean - raw.loc[261:270].mean(),
                "candidate_recovered_at_80pct": ratio >= reference_relative_threshold,
                "included_in_previous_eight_model_hysteresis": model != "CAS-ESM2-0",
                "member": source_metadata.loc[model, "Selected member"],
                "actual_latitude_degN": source_metadata.loc[model, "Actual latitude (degN)"]}
    for window in [10, 20, 30]:
        trend_years = np.arange(final_year - window + 1, final_year + 1)
        # 末期仍受变动CO₂强迫；OLS斜率仅描述轨迹，不作显著性或稳态判定。
        slope = np.polyfit(trend_years, raw.loc[trend_years].to_numpy(), 1)[0]
        endpoint[f"late_{window}yr_slope_sv_per_decade"] = 10 * slope

    for year in years:
        normalized_curves.append({"model": model, "short_name": short_names[model_index],
                                  "simulation_year": year, "co2_ppm": co2.loc[year],
                                  "raw_sv": raw.loc[year], "smooth_11yr_sv": smoothed.loc[year],
                                  "raw_percent_initial": 100 * raw.loc[year] / initial_mean,
                                  "smooth_11yr_percent_initial": 100 * smoothed.loc[year] / initial_mean,
                                  "candidate_recovered_at_80pct": ratio >= reference_relative_threshold})

    for smoothing_window in smoothing_windows:
        filtered = raw.rolling(smoothing_window, center=True, min_periods=1).mean()
        min_year = int(filtered.idxmin())
        increase = filtered.loc[up_years].to_numpy()
        decrease = filtered.loc[down_years].to_numpy()
        difference = decrease - increase
        dx = np.diff(co2_nodes)
        left, right = np.abs(difference[:-1]), np.abs(difference[1:])
        interval_area = dx * (left + right) / 2
        crossing = difference[:-1] * difference[1:] < 0
        # 跨零时使用两个三角形的总面积，避免绝对差梯形高估。
        interval_area[crossing] = dx[crossing] * (left[crossing] ** 2 + right[crossing] ** 2) / (2 * (left[crossing] + right[crossing]))
        area = interval_area.sum()
        normalized_area = area / (initial_mean * (co2_nodes[-1] - co2_nodes[0]))
        smoothing_rows.append({"model": model, "smoothing_years": smoothing_window,
                               "minimum_year": min_year, "minimum_lag_years": min_year - peak_year,
                               "minimum_sv": filtered.loc[min_year],
                               "hysteresis_area_sv_ppm": area, "normalized_hysteresis_area": normalized_area,
                               "interior_crossings": int(crossing.sum())})
        if smoothing_window == 11:
            endpoint["hysteresis_area_sv_ppm"] = area
            endpoint["normalized_hysteresis_area"] = normalized_area

        if smoothing_window in [5, 11, 21]:
            for lag_threshold in [5, 10, 15]:
                lag_threshold_rows.append({"model": model, "smoothing_years": smoothing_window,
                                           "lag_threshold_years": lag_threshold,
                                           "minimum_lag_years": min_year - peak_year,
                                           "delayed_minimum_candidate": min_year - peak_year > lag_threshold})

    for late_window in late_windows:
        first_late_year = final_year - late_window + 1
        late_window_mean = raw.loc[first_late_year:final_year].mean()
        matched_early_mean = raw.loc[1:late_window].mean()
        for baseline_type in ["fixed_initial_10yr", "matched_initial_window"]:
            baseline_mean = initial_mean if baseline_type == "fixed_initial_10yr" else matched_early_mean
            window_ratio = late_window_mean / baseline_mean
            window_rows.append({"model": model, "baseline_type": baseline_type,
                                "late_window_years": late_window, "initial_mean_sv": baseline_mean,
                                "late_mean_sv": late_window_mean, "endpoint_ratio": window_ratio,
                                "late_minus_initial_sv": late_window_mean - baseline_mean,
                                "initial_mean_co2_ppm": co2.loc[1:10 if baseline_type == "fixed_initial_10yr" else late_window].mean(),
                                "late_mean_co2_ppm": co2.loc[first_late_year:280].mean()})
            for threshold in relative_thresholds:
                relative_rows.append({"model": model, "baseline_type": baseline_type,
                                      "late_window_years": late_window, "threshold": threshold,
                                      "endpoint_ratio": window_ratio,
                                      "recovered_candidate": window_ratio >= threshold})

    for tolerance in absolute_tolerances_sv:
        absolute_rows.append({"model": model, "tolerance_sv": tolerance,
                              "late_minus_initial_sv": residual,
                              "shortfall_within_tolerance": residual >= -tolerance,
                              "within_symmetric_band": abs(residual) <= tolerance})

    # 最低点之后，首个十年滑动均值越过阈值且之后所有完整十年均值均维持阈值以上。
    trailing_mean = raw.rolling(10, min_periods=10).mean()
    eligible = trailing_mean.loc[minimum_year + 10:280] / initial_mean
    for threshold in relative_thresholds:
        above = eligible.to_numpy() >= threshold
        remains_above = np.logical_and.accumulate(above[::-1])[::-1]
        candidate_years = eligible.index[remains_above]
        return_year = float(candidate_years[0]) if len(candidate_years) else np.nan
        return_rows.append({"model": model, "threshold": threshold,
                            "first_sustained_trailing10_window_end_year": return_year,
                            "years_after_co2_peak": return_year - peak_year,
                            "minimum_reference_year": minimum_year})

    endpoint_rows.append(endpoint)


# %% 5. 核对与保存诊断

endpoints = pd.DataFrame(endpoint_rows)
windows = pd.DataFrame(window_rows)
relative = pd.DataFrame(relative_rows)
smoothing = pd.DataFrame(smoothing_rows)
previous = pd.read_csv(root / "output/s01_amoc/n20_amoc_grouped_hysteresis/endpoint_and_hysteresis_metrics.csv").set_index("model")
current = endpoints.set_index("model").loc[previous.index]
assert np.allclose(current.endpoint_ratio, previous.late_over_initial, atol=1e-12)
assert np.allclose(current.hysteresis_area_sv_ppm, previous.hysteresis_area_sv_ppm, atol=1e-8)

# 先检验70%至85%与80%是否同分组；不把这一样本空档解释为物理临界点。
comparison = relative.loc[relative.threshold.between(0.70, 0.85)].copy()
reference = endpoints.set_index("model").candidate_recovered_at_80pct
comparison["same_as_reference_80pct_10yr"] = comparison.recovered_candidate.to_numpy() == reference.loc[comparison.model].to_numpy()
stability = comparison.groupby(["baseline_type", "late_window_years", "threshold"], as_index=False).same_as_reference_80pct_10yr.all()

output_directory.mkdir(parents=True)
endpoints.to_csv(output_directory / "nine_model_endpoint_diagnostics.csv", index=False, encoding="utf-8-sig")
pd.DataFrame(normalized_curves).to_csv(output_directory / "nine_model_normalized_series.csv", index=False, encoding="utf-8-sig")
windows.to_csv(output_directory / "endpoint_window_sensitivity.csv", index=False, encoding="utf-8-sig")
relative.to_csv(output_directory / "relative_threshold_sensitivity.csv", index=False, encoding="utf-8-sig")
pd.DataFrame(absolute_rows).to_csv(output_directory / "absolute_tolerance_sensitivity.csv", index=False, encoding="utf-8-sig")
smoothing.to_csv(output_directory / "area_and_minimum_smoothing_sensitivity.csv", index=False, encoding="utf-8-sig")
pd.DataFrame(lag_threshold_rows).to_csv(output_directory / "minimum_lag_threshold_sensitivity.csv", index=False, encoding="utf-8-sig")
pd.DataFrame(return_rows).to_csv(output_directory / "sustained_threshold_return_years.csv", index=False, encoding="utf-8-sig")
stability.to_csv(output_directory / "classification_stability.csv", index=False, encoding="utf-8-sig")
metadata = {"purpose": "Exploratory audit of recovery definitions before manuscript classification; not a new adopted grouping",
            "source": str(source_path), "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
            "value_column": "Forced AMOC strength (Sv)", "drift_correction": "none",
            "definition": "Atlantic-Arctic overturning maximum at nearest 40N latitude, depth >=500m",
            "dimensions": "Source CSV model x simulation_year; pivot simulation_year x model; output nine_model_normalized_series model x simulation_year",
            "models": models, "scope": "All nine models including CAS; previous eight-model exclusion and four-model atmosphere composites unchanged",
            "years": [1, 280], "initial_reference": "Raw forced years1-10 mean; not an equilibrium piControl mean",
            "R": "Mean years271-280 / mean years1-10: fraction of initial strength attained, not fraction of maximum lost strength regained",
            "F": "(late10 - minimum centered11) / (initial10 - minimum centered11); exploratory, selected minimum and smoothing affect this metric; not clipped at1",
            "relative_thresholds": relative_thresholds, "absolute_tolerances_sv": absolute_tolerances_sv,
            "window_sensitivity": "Last10/20/30 years compared with fixed first10 baseline and, separately, matched first10/20/30 windows. All windows sample different CO2 along a transient forcing trajectory; not uncertainty intervals or equilibrium estimates",
            "area": "Absolute branch separation integrated over common linear CO2 domain, increasing years1-140 paired to combined years279-140; crossings treated exactly for piecewise linear branches",
            "Astar": "Area / (initial10 AMOC * common CO2 range); dimensionless mean relative path separation; distinct from endpoint ratio",
            "lag": "Year of smoothed AMOC minimum minus CO2 peak year140; not the time to full recovery; sensitivity5/11/21 years",
            "sustained_return_year": "Earliest end of a full trailing10 mean window after minimum+10, exceeding threshold with all later complete trailing10 means still above until year280; descriptive and right-censored if not attained",
            "trend": "OLS slope on raw last10/20/30 years in Sv per decade; descriptive only, no significance/steady-state inference",
            "uncertainty": "No equivalence test, confidence interval or claim of permanent irreversibility; short forced windows are not long stationary internal-variability records",
            "selection_caution": "Threshold sensitivity diagnoses current-sample separation, not a universal or objectively estimated physical recovery threshold; do not select cutoffs to maximize downstream significance",
            "adopted_changes": "None: previous memberships, schematic and atmosphere results untouched"}
(output_directory / "metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
qc = {"all_nine_models_present": len(endpoints) == 9, "source_values_complete": True,
      "reproduces_previous_eight_endpoint_ratios_and_areas": True,
      "n_series_rows": len(normalized_curves), "n_window_rows": len(windows),
      "n_relative_threshold_rows": len(relative),
      "all_baseline_types_late_windows_stable_at_80pct": bool(stability.loc[stability.threshold == 0.8, "same_as_reference_80pct_10yr"].all())}
(output_directory / "qc_report.json").write_text(json.dumps(qc, indent=2), encoding="utf-8")
columns = ["model", "endpoint_percent_initial", "late_minus_initial_sv",
           "regained_fraction_of_maximum_loss", "minimum_lag_from_co2_peak_years",
           "normalized_hysteresis_area", "late_20yr_slope_sv_per_decade"]
print(endpoints[columns].round(4).to_string(index=False))
print(stability.to_string(index=False))
