# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd


# %% 2. 参数与路径设置

root = Path(r"F:\sea\Antic sea ice constraint")
directory = root / "output/s01_amoc/n23_amoc_two_class_schematic"
audit_path = root / "output/s01_amoc/n22_amoc_recovery_definition_audit/nine_model_endpoint_diagnostics.csv"
assert not directory.exists(), "正式输出已存在；重新执行前须先归档。"
titles = ["Substantial recovery", "Incomplete recovery"]
recovery_threshold = 0.80
final_coordinates = [0.92, 0.43]
ellipse_angle = np.linspace(0, np.pi, 601)
concentration = (1 - np.cos(ellipse_angle)) / 2
ellipse_half_width = 0.22
arrow_indices = {"increase": (300, 355), "decrease": (300, 245)}
coordinate_rows, arrow_rows = [], []


# %% 3. 保存九模式端点二分类

# 用户选择两类后，仅在本阶段记录纯AMOC端点判据；不覆盖旧层结/大气复合成员表。
audit = pd.read_csv(audit_path)
membership = audit[["model", "member", "initial_mean_sv", "late_mean_sv",
                    "late_minus_initial_sv", "endpoint_ratio",
                    "included_in_previous_eight_model_hysteresis"]].copy()
membership["endpoint_class"] = titles[1]
membership.loc[membership.endpoint_ratio >= recovery_threshold, "endpoint_class"] = titles[0]
membership["criterion"] = "Final years271-280 / initial years1-10 mean >=0.8 versus <0.8"
assert membership.endpoint_class.value_counts().to_dict() == {titles[0]: 5, titles[1]: 4}


# %% 4. 构造两类概念路径

# 上弧取椭圆四分弧，起点为最高点，避免倾斜椭圆引入虚假的初始增强。
ellipse_centerline = 1 - 0.72 * concentration
increase = 0.28 + 0.72 * np.sqrt(1 - concentration ** 2)
assert increase.max() <= 1
assert np.all(np.diff(increase) < 0)
for panel, title in enumerate(titles):
    endpoint_gap = 1 - final_coordinates[panel]
    decrease = ellipse_centerline - ellipse_half_width * np.sin(ellipse_angle) - endpoint_gap * (1 - concentration)
    decrease[0] = final_coordinates[panel]
    decrease[-1] = increase[-1]
    assert np.all(decrease <= increase + 1e-12)
    assert np.all(decrease > 0)
    assert (decrease[0] >= recovery_threshold) == (panel == 0)

    for index, value in enumerate(concentration):
        coordinate_rows.append({"panel": panel, "title": title,
                                "co2_coordinate": value,
                                "increase_amoc_coordinate": increase[index],
                                "decrease_amoc_coordinate": decrease[index]})
    for phase, (start, end) in arrow_indices.items():
        curve = increase if phase == "increase" else decrease
        arrow_rows.append({"panel": panel, "phase": phase,
                           "start_x": concentration[start], "start_y": curve[start],
                           "end_x": concentration[end], "end_y": curve[end]})


# %% 5. 保存坐标、分组和元数据

directory.mkdir(parents=True)
membership.to_csv(directory / "amoc_endpoint_two_class_membership.csv", index=False, encoding="utf-8-sig")
pd.DataFrame(coordinate_rows).to_csv(directory / "schematic_coordinates.csv", index=False, encoding="utf-8-sig")
pd.DataFrame(arrow_rows).to_csv(directory / "schematic_arrows.csv", index=False, encoding="utf-8-sig")
metadata = {"type": "Two-class conceptual schematic, not model trajectories or a composite mean",
            "user_decision": "Two classes selected after nine-model recovery audit; user requested tilted elliptical schematic arcs for more intuitive loops",
            "titles": titles, "recovery_threshold": recovery_threshold,
            "R_definition": "Unsmoothed forced AMOC mean years271-280 divided by mean years1-10, nearest40N Atlantic-Arctic overturning maximum at depths>=500m",
            "membership_source": str(audit_path), "membership_source_sha256": hashlib.sha256(audit_path.read_bytes()).hexdigest(),
            "membership_scope": "All nine source models including CAS; no change to previous stratification-screened or four-model atmospheric composites",
            "data_treatment": "Source audit uses forced raw AMOC, no piControl drift correction; this schematic reads only saved audit outputs",
            "dimensions": "Coordinates: panel x drawing node; arrows: panel x phase; membership: model rows. No new NetCDF/MAT read",
            "schematic_endpoints": final_coordinates,
            "geometry": "Ellipse-derived rounded loop; upper branch is a monotone ellipse quadrant starting at its maximum1, lower branch is an opened tilted-ellipse arc; no forced closure",
            "endpoint_caution": "0.92 and0.43 are chosen illustrative drawing coordinates, not estimated group means; actual recovery criterion uses ten-year means",
            "formulas": {"parameter": "theta from0 to pi; x=(1-cos(theta))/2; ellipse_half_width=0.22",
                         "centerline": "1-0.72*x",
                         "increase": "0.28+0.72*sqrt(1-x**2), an ellipse quadrant; max1 at x0 and strictly decreasing",
                         "substantial_recovery_decrease": "centerline-0.22*sin(theta)-0.08*(1-x)",
                         "incomplete_recovery_decrease": "centerline-0.22*sin(theta)-0.57*(1-x)"},
            "geometry_caution": "The upper arc is constructed analytically to remain <=1 and decrease strictly; no clipping or flat cap is applied. Curve shape remains conceptual, not a model fit. Common endpoints and recovery definition remain unchanged",
            "interpretation": "Both classes retain hysteresis. Substantial recovery need not return exactly to initial strength; incomplete recovery means below adopted endpoint threshold, not permanent irreversibility or absence of rebound",
            "illustrative_scope": "No model labels, fitted timescales, estimated loop areas, uncertainty bands or significance claims in schematic",
            "style": "Restrained two-color vector drawing with common axes; generic journal-style design, no journal-specific compliance claim",
            "reference_framework": "Kim et al.2022, https://doi.org/10.1038/s41558-022-01452-z; endpoint return distinct from path separation; 80% is this study's operational threshold, not a threshold from that paper"}
(directory / "metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
print(membership[["model", "endpoint_ratio", "endpoint_class"]].to_string(index=False))
