# %% 1. 导入依赖

from pathlib import Path
import ast
import difflib
import hashlib
import io
import json
import re
import tokenize


# %% 2. 参数与路径设置

root = Path("F:/sea/Antic sea ice constraint")
destination = root / "output/s05_atmosphere/p05_amoc_recovery_composites/n49_group_naming"
assert not destination.exists(), "正式命名清单已存在，须先归档。"
atmosphere = "plot/s05_atmosphere/p05_amoc_recovery_composites/"
scripts = ["plot/s01_amoc/n20_amoc_grouped_hysteresis.py",
           "plot/S02_sia/s02_arctic_region/n059_grouped_initial_n2_profiles.py",
           "plot/S02_sia/s02_arctic_region/n061_grouped_mld_timeseries.py",
           "plot/S02_sia/s02_arctic_region/n066_grouped_tas_ice_patterns.py"]
scripts += [atmosphere + name for name in ["n03_endpoint_groups.py", "n10_member_jet_hydrology.py", "n12_group_jet_contrasts.py",
            "n15_moisture_group_contrasts.py", "n17_member_mechanism_comparison.py", "n23_pairwise_jet_budget.py",
            "n30_nh_maps_and_budget_review.py", "n31_nh_land_precipitation.py", "n33_nh_circulation_sections.py",
            "n37_global_meridional_circulation.py", "n38_grouped_climate_precipitation_jet.py",
            "n40_net_precipitation_circulation_patterns.py", "n46_jet_position_strength.py",
            "n47_hydrology_rainbelt_comparison.py", "n48_climate_patterns_djf_u250.py"]]
registry = {"basis": "User-approved original two-class AMOC hysteresis terminology;2026-10-22",
            "recovery_only": {"amoc_recovered": "Substantial recovery", "amoc_not_recovered": "Incomplete recovery",
                              "recovered": "Substantial recovery", "not_recovered": "Incomplete recovery"},
            "three_group_titles": {"better_nonrecovery": "Better-simulated\nIncomplete recovery",
                                   "better_recovery": "Better-simulated\nSubstantial recovery", "larger_bias": "Larger-bias models"},
            "core_names": {"recovered": "Better-simulated\nSubstantial recovery", "not_recovered": "Better-simulated\nIncomplete recovery",
                           "contrast": "Incomplete − substantial recovery"},
            "criterion": "R=mean(AMOC271-280)/mean(AMOC1-10);Substantial recovery R>=0.8,Incomplete recovery R<0.8",
            "members": {"better_nonrecovery": ["CESM2", "CNRM-ESM2-1"], "better_recovery": ["NorESM2-LM", "GFDL-ESM4"],
                        "larger_bias": ["ACCESS-ESM1-5", "CanESM5", "MIROC-ES2L", "UKESM1-0-LL"]},
            "scope": "Better-simulated describes study-region sea-ice/upper-ocean diagnostics,not global skill orcompleted strictMLDvalidation",
            "legacy": "Earlier alternative screening/member counts remain unchanged;only AMOC recovery terminology is standardized there",
            "identifiers": "Internal recovered/not_recovered/better_nonrecovery keys and scientific files stay unchanged for downstream compatibility"}
replacements = [("AMOC not recovered", "Incomplete recovery"), ("AMOC recovered", "Substantial recovery"),
                ("Not recovered", "Incomplete recovery"), ("not recovered", "incomplete recovery"),
                ("Non-recovery", "Incomplete recovery"), ("non-recovery", "incomplete recovery"),
                ("Better-simulated stratification", "Better-simulated"),
                ("Larger stratification bias", "Larger-bias models"), ("Larger bias", "Larger-bias models")]


# %% 3. 仅生成显示文字的修改补丁，不修改科学成员或文件键

patches, change_records, figure_directories = [], [], []
for relative_path in scripts:
    path = root / relative_path
    original = path.read_text(encoding="utf-8")
    revised_tokens = []
    for token in tokenize.generate_tokens(io.StringIO(original).readline):
        if token.type == tokenize.STRING:
            value = token.string
            for previous, current in replacements:
                value = value.replace(previous, current)
            value = re.sub(r"\bRecovered\b", "Substantial recovery", value)
            value = re.sub(r"\bRecovery\b", "Substantial recovery", value)
            # 含下划线的机器键与独立小写recovered字符串保持原样。
            if " " in value:
                value = re.sub(r"\brecovered\b", "substantial recovery", value)
            token = token._replace(string=value)
        revised_tokens.append(token)
    revised = tokenize.untokenize(revised_tokens)
    name = path.name
    needs_registry = name in ["n20_amoc_grouped_hysteresis.py", "n059_grouped_initial_n2_profiles.py", "n061_grouped_mld_timeseries.py",
                              "n066_grouped_tas_ice_patterns.py", "n37_global_meridional_circulation.py", "n38_grouped_climate_precipitation_jet.py",
                              "n40_net_precipitation_circulation_patterns.py", "n46_jet_position_strength.py", "n47_hydrology_rainbelt_comparison.py",
                              "n48_climate_patterns_djf_u250.py"]
    if needs_registry:
        if not re.search(r"^import json$", revised, flags=re.MULTILINE):
            revised = revised.replace("from pathlib import Path\n", "from pathlib import Path\nimport json\n", 1)
        revised = re.sub(r"^(root = Path\(.+\))$", r'\1\nnaming = json.loads((root / "output/s05_atmosphere/p05_amoc_recovery_composites/n49_group_naming/group_labels.json").read_text(encoding="utf-8"))', revised, count=1, flags=re.MULTILINE)
    if name == "n059_grouped_initial_n2_profiles.py":
        revised = revised.replace("panel_title = members.panel_title.iloc[0]", 'panel_title = naming["three_group_titles"][group_key]')
    elif name == "n061_grouped_mld_timeseries.py":
        revised = revised.replace("ax.set_title(selection.panel_title.iloc[0],", 'ax.set_title(naming["three_group_titles"][group],')
    elif name == "n20_amoc_grouped_hysteresis.py":
        revised = revised.replace("1.18, group.title,", '1.18, naming["recovery_only"].get(group.group, group.title),')
    elif name == "n37_global_meridional_circulation.py":
        revised = re.sub(r"^names = .+$", 'names = naming["core_names"]', revised, count=1, flags=re.MULTILINE)
        revised = revised.replace("fontsize=11.5)", "fontsize=10.5)")
        revised = revised.replace("fontsize=12)\n        axis.grid", "fontsize=10.5)\n        axis.grid")
    elif name in ["n066_grouped_tas_ice_patterns.py", "n38_grouped_climate_precipitation_jet.py"]:
        variable = "main_titles" if name.startswith("n066") else "titles"
        revised = re.sub(r"^" + variable + r" = .+$", variable + ' = list(naming["three_group_titles"].values())', revised, count=1, flags=re.MULTILINE)
        if name.startswith("n066"):
            revised = revised.replace('fontsize=12, color="#20282D"', 'fontsize=10.2, color="#20282D"')
    elif name in ["n40_net_precipitation_circulation_patterns.py", "n46_jet_position_strength.py", "n47_hydrology_rainbelt_comparison.py", "n48_climate_patterns_djf_u250.py"]:
        revised = re.sub(r"^group_titles = .+$", 'group_titles = [naming["three_group_titles"][key] + " (n = " + str(count) + ")" for key, count in zip(naming["members"], [2, 2, 4])]', revised, count=1, flags=re.MULTILINE)
    # 将实际图名注册到新manifest；旧输入metadata中的显示字段只作历史记录，不覆盖。
    if needs_registry and 'manifest = {' in revised:
        revised = revised.replace('manifest = {', 'manifest = {"group_naming_registry": str(root / "output/s05_atmosphere/p05_amoc_recovery_composites/n49_group_naming/group_labels.json"), "canonical_group_titles": naming["three_group_titles"], ', 1)
    ast.parse(revised)
    assert not any(isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef, ast.Lambda)) for node in ast.walk(ast.parse(revised)))
    if revised != original:
        diff = list(difflib.unified_diff(original.splitlines(keepends=True), revised.splitlines(keepends=True), n=3))
        patches.append("*** Update File: " + str(path).replace("\\", "/") + "\n" + "".join(diff[2:]))
        change_records.append({"path": relative_path, "before_sha256": hashlib.sha256(original.encode()).hexdigest(),
                               "after_sha256": hashlib.sha256(revised.encode()).hexdigest(), "registry": needs_registry})
    # 对脚本中明确写出的正式图目录建立归档清单；n31为stage变量单独处理。
    matches = re.findall(r"['\"](figure/[^'\"]+)['\"]", original)
    if name == "n31_nh_land_precipitation.py":
        matches = ["figure/s05_atmosphere/p05_amoc_recovery_composites/n31_nh_land_precipitation"]
    for match in matches:
        folder = root / match
        if folder.is_dir() and match not in figure_directories:
            figure_directories.append(match)


# %% 4. 记录未改动科学输入的哈希和当前图件目录

science_directories = [root / "output/s01_amoc/n20_amoc_grouped_hysteresis", root / "output/s01_amoc/n23_amoc_two_class_schematic",
                      root / "output/s02_sia/s02_arctic_region/p20_grouped_initial_stratification/n059_grouped_initial_n2_profiles",
                      root / "output/s02_sia/s02_arctic_region/p21_grouped_mld/n061_grouped_mld_timeseries",
                      root / "output/s02_sia/s02_arctic_region/p21_grouped_mld/n063_observation_dense_mld",
                      root / "output/s02_sia/s02_arctic_region/p22_grouped_temperature_ice/n066_grouped_tas_ice_patterns"]
stage_output = root / "output/s05_atmosphere/p05_amoc_recovery_composites"
science_directories += [folder for folder in stage_output.iterdir() if folder.is_dir() and folder.name.startswith(
    ("n01_", "n02_", "n09_", "n11_", "n14_", "n16_", "n22_", "n28_", "n29_", "n32_", "n35_", "n36_", "n37_", "n39_", "n41_", "n44_", "n45_"))]
science_hashes = []
for folder in science_directories:
    for path in sorted(folder.iterdir()):
        if path.is_file() and path.suffix.lower() in [".csv", ".nc", ".npz", ".json"]:
            digest = hashlib.sha256()
            with path.open("rb") as stream:
                for chunk in iter(stream.read, b""):
                    # 文件较小，按一整块读取；不修改文件。
                    digest.update(chunk)
            science_hashes.append({"path": str(path), "sha256": digest.hexdigest(), "bytes": path.stat().st_size})
destination.mkdir(parents=True)
(destination / "group_labels.json").write_text(json.dumps(registry, ensure_ascii=False, indent=2), encoding="utf-8")
(destination / "source_label_changes.patch").write_text("*** Begin Patch\n" + "".join(patches) + "*** End Patch\n", encoding="utf-8")
(destination / "rename_plan.json").write_text(json.dumps({"script_changes": change_records, "figure_directories": figure_directories,
                                                          "science_hashes": science_hashes}, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps({"scripts": len(change_records), "active_figure_directories": figure_directories,
                  "science_inputs": len(science_hashes)}, ensure_ascii=False, indent=2))
