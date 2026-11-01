# -*- coding: utf-8 -*-
"""
从原始CMIP6 thetao/tas重新计算当前北极未订正温度endpoint pattern与CA/NA/NP区域值。

【是否需要运行】
需要。当前正式SIA pattern完成后运行本脚本。

【输入文件】
F:/sea/cmipthetaoad
F:/sea/cmipthetao
F:/sea/cmiptasad 或 F:/sea/cmip6tasad
F:/sea/cmiptas 或 F:/sea/cmip6tas
output/s02sia/s02_arctic_region/p04_pattern_unadjusted/n01_arctic_sia_pattern/arctic_sia_pattern_cells.csv

【输出文件】
output/s02sia/s02_arctic_region/p04_pattern_unadjusted/n02_temperature_pattern/
    arctic_temperature_pattern_cells.csv
    arctic_temperature_region_summary.csv
    arctic_temperature_qa.csv
    arctic_temperature_failures.csv
    arctic_temperature_metadata.json

【什么时候需要重新运行】
原始thetao/tas、成员、CDR时间规则、表层选择、重网格规则或区域定义变化后。

科学定义
1. 不做piControl订正。
2. raw ΔT = mean(CDR experiment years 131–140) - mean(1pctCO2 experiment years 1–10)。
3. 普通CDR-only文件固定取实验年索引[130:140]，绝不取文件最后10年。
4. thetao取最浅有效海洋层。
5. CAS thetao/tas只读标准原生月资料，排除旧年度和重网格文件。
6. CAS降CO2窗口为月索引[1560:1680]；thetao为0131–0140年，2米tas为0470–0479年。
7. NorESM2-LM thetao/tas的CO2和CDR均强制使用对齐成员r1i1p1f1，禁止回退旧p4数据。
8. thetao和tas均使用9模式，包括GFDL-ESM4。
9. 区域温度均值使用固定地理区域 + 各变量九模式共同有效1°格点，不使用模式自身SIC掩膜。
10. CA=70–90°N；NA=45–70°N, 60°W–30°E；NP=45–70°N, 120°E–150°W。
11. --models CAS-ESM2-0仅重算CAS，其余二维场从已存结果读取；正式输出先归档再更新。
"""

# %% 1. 导入依赖
from pathlib import Path
from datetime import datetime
import argparse
import hashlib
import shutil
import json
import os
import re
import warnings

import numpy as np
import pandas as pd
import xarray as xr
from scipy.spatial import cKDTree

warnings.filterwarnings("ignore", category=RuntimeWarning)

# %% 2. 路径与参数
ROOT = Path(r"F:\sea\Antic sea ice constraint")
SEA_ROOT = Path(r"F:\sea")

THETAO_CO2_ROOT = SEA_ROOT / "cmipthetaoad"
THETAO_CDR_ROOT = SEA_ROOT / "cmipthetao"

TAS_CO2_CANDIDATES = [SEA_ROOT / "cmiptasad", SEA_ROOT / "cmip6tasad"]
TAS_CDR_CANDIDATES = [SEA_ROOT / "cmiptas", SEA_ROOT / "cmip6tas"]

TAS_CO2_ROOT = next((path for path in TAS_CO2_CANDIDATES if path.is_dir()), TAS_CO2_CANDIDATES[0])
TAS_CDR_ROOT = next((path for path in TAS_CDR_CANDIDATES if path.is_dir()), TAS_CDR_CANDIDATES[0])

SIA_DIR = ROOT / "output" / "s02sia" / "s02_arctic_region" / "p04_pattern_unadjusted" / "n01_arctic_sia_pattern"
SIA_CELL_FILE = SIA_DIR / "arctic_sia_pattern_cells.csv"

OUT = ROOT / "output" / "s02sia" / "s02_arctic_region" / "p04_pattern_unadjusted" / "n02_temperature_pattern"

CELL_FILE = OUT / "arctic_temperature_pattern_cells.csv"
SUMMARY_FILE = OUT / "arctic_temperature_region_summary.csv"
QA_FILE = OUT / "arctic_temperature_qa.csv"
FAIL_FILE = OUT / "arctic_temperature_failures.csv"
META_FILE = OUT / "arctic_temperature_metadata.json"

for path in [THETAO_CO2_ROOT, THETAO_CDR_ROOT, TAS_CO2_ROOT, TAS_CDR_ROOT]:
    if not path.is_dir():
        raise FileNotFoundError(path)

if not SIA_CELL_FILE.is_file():
    raise FileNotFoundError(f"先运行当前SIA pattern src，缺少：{SIA_CELL_FILE}")

OUT.mkdir(parents=True, exist_ok=True)

MODELS = [
    "ACCESS-ESM1-5",
    "CanESM5",
    "CESM2",
    "CNRM-ESM2-1",
    "GFDL-ESM4",
    "MIROC-ES2L",
    "UKESM1-0-LL",
    "CAS-ESM2-0",
    "NorESM2-LM",
]

THETAO_MODELS = MODELS.copy()
TAS_MODELS = MODELS.copy()

ALIASES = {
    "ACCESS-ESM1-5": ("ACCESS-ESM1-5", "ACCESSESM5", "ACCESS"),
    "CanESM5": ("CanESM5", "CANESM5"),
    "CESM2": ("CESM2",),
    "CNRM-ESM2-1": ("CNRM-ESM2-1", "CNRM-ESM2"),
    "GFDL-ESM4": ("GFDL-ESM4", "GFDL"),
    "MIROC-ES2L": ("MIROC-ES2L", "MIROC"),
    "UKESM1-0-LL": ("UKESM1-0-LL", "UKESM1"),
    "CAS-ESM2-0": ("CAS-ESM2-0", "CAS-ESM2", "CAS"),
    "NorESM2-LM": ("NorESM2-LM", "NORESM2-LM", "NorESM2-LM2", "NORESM2-LM2"),
}

MEMBER_PRIORITY = ("r1i1p1f1", "r1i1p1f2", "r1i1p1f3")
NOR_REQUIRED_MEMBER = "r1i1p1f1"
MEMBER_RE = re.compile(r"r\d+i\d+p\d+f\d+", re.IGNORECASE)

WINDOW_YEARS = 10
CO2_START_YEAR_INDEX = 0
CDR_START_YEAR_INDEX = 130
parser = argparse.ArgumentParser()
parser.add_argument('--models', nargs='+', choices=MODELS, default=MODELS)
RECALCULATE_MODELS = parser.parse_args().models

EARTH_RADIUS_KM = 6371.0088
MAX_REGRID_DISTANCE_KM = 450.0
IDW_NEIGHBOURS = 4
IDW_POWER = 2.0

TARGET_LAT = np.arange(45.5, 90.0, 1.0)
TARGET_LON = np.arange(-179.5, 180.0, 1.0)
TARGET_LON_GRID, TARGET_LAT_GRID = np.meshgrid(TARGET_LON, TARGET_LAT)

target_latr = np.deg2rad(TARGET_LAT_GRID.ravel())
target_lonr = np.deg2rad(TARGET_LON_GRID.ravel())
TARGET_XYZ = np.column_stack([np.cos(target_latr) * np.cos(target_lonr), np.cos(target_latr) * np.sin(target_lonr), np.sin(target_latr)])

TIME_CODER = xr.coders.CFDatetimeCoder(use_cftime=True)

# %% 3. 读取当前SIA格点，只用于共享1°网格与海冰边缘绘图字段
sia_cells = pd.read_csv(SIA_CELL_FILE)
sia_required = {"model", "latitude_deg_north", "longitude_deg_east", "edge_raw_co2_sic_fraction", "edge_raw_cdr_sic_fraction"}

if not sia_required.issubset(sia_cells.columns):
    raise KeyError(f"当前SIA cells缺少字段：{sorted(sia_required - set(sia_cells.columns))}")

sia_cells = sia_cells[["model", "latitude_deg_north", "longitude_deg_east", "edge_raw_co2_sic_fraction", "edge_raw_cdr_sic_fraction"]].copy()
sia_cells = sia_cells.sort_values(["model", "latitude_deg_north", "longitude_deg_east"]).reset_index(drop=True)

# %% 4. 从原始NC重新计算thetao和tas
# 注意：本段不使用SIC决定温度采样区域。
temperature_fields = {"thetao": {}, "tas": {}}
qa_rows = []
failure_rows = []

# 部分更新只读已有二维场，其他模式不再读取原始数据。
if set(RECALCULATE_MODELS) != set(MODELS):
    saved_cells = pd.read_csv(CELL_FILE, float_precision='round_trip')
    saved_qa = pd.read_csv(QA_FILE)
    qa_rows = saved_qa.loc[~saved_qa['model'].isin(RECALCULATE_MODELS)].to_dict('records')
    for variable in ['thetao', 'tas']:
        for model in MODELS:
            if model not in RECALCULATE_MODELS:
                saved_model = saved_cells.loc[saved_cells['model'] == model].sort_values(
                    ['latitude_deg_north', 'longitude_deg_east'])
                temperature_fields[variable][model] = saved_model[
                    f'raw_delta_{variable}_degC'].to_numpy().reshape(TARGET_LAT_GRID.shape)

for variable in ["thetao", "tas"]:
    models = THETAO_MODELS if variable == "thetao" else TAS_MODELS
    co2_root = THETAO_CO2_ROOT if variable == "thetao" else TAS_CO2_ROOT
    cdr_root = THETAO_CDR_ROOT if variable == "thetao" else TAS_CDR_ROOT

    print("\n" + "#" * 110)
    print(f"{variable} | CO2={co2_root} | CDR={cdr_root}")
    print("#" * 110)

    for model in models:
        if model not in RECALCULATE_MODELS:
            continue
        print("\n" + "=" * 110)
        print(f"[{variable}] {model}")
        print("=" * 110)

        try:
            stage_files = {}
            stage_members = {}

            # %% 4.1 查找CO2/CDR目录与文件
            for stage, base_root in [("co2", co2_root), ("cdr", cdr_root)]:
                matched_dirs = {}

                for folder in sorted(path for path in base_root.iterdir() if path.is_dir()):
                    normalized_folder = re.sub(r"[^A-Z0-9]", "", folder.name.upper())

                    if any(normalized_folder == re.sub(r"[^A-Z0-9]", "", alias.upper()) for alias in ALIASES[model]):
                        key = os.path.normcase(os.path.realpath(str(folder)))
                        matched_dirs.setdefault(key, folder)

                if not matched_dirs:
                    for folder in sorted(path for path in base_root.iterdir() if path.is_dir()):
                        normalized_folder = re.sub(r"[^A-Z0-9]", "", folder.name.upper())

                        if any(
                            re.sub(r"[^A-Z0-9]", "", alias.upper()) in normalized_folder
                            or normalized_folder in re.sub(r"[^A-Z0-9]", "", alias.upper())
                            for alias in ALIASES[model]
                        ):
                            key = os.path.normcase(os.path.realpath(str(folder)))
                            matched_dirs.setdefault(key, folder)

                if len(matched_dirs) == 0:
                    raise FileNotFoundError(f"{model} {variable} {stage}未找到模式目录：{base_root}")

                if len(matched_dirs) > 1:
                    raise RuntimeError(f"{model} {variable} {stage}命中多个真实不同目录：{list(matched_dirs.values())}")

                folder = list(matched_dirs.values())[0]
                files = [path for path in sorted(folder.rglob("*.nc")) if variable in path.name.lower()]

                if not files:
                    raise FileNotFoundError(f"{model} {variable} {stage}目录无变量文件：{folder}")

                # CAS只用标准月文件，排除老师的年度及重网格派生文件。
                if model == 'CAS-ESM2-0':
                    experiment = '1pctCO2' if stage == 'co2' else '1pctCO2-cdr'
                    table = 'Omon' if variable == 'thetao' else 'Amon'
                    pattern = rf'{variable}_{table}_CAS-ESM2-0_{experiment}_r1i1p1f1_gn_\d{{6}}-\d{{6}}\.nc'
                    files = [path for path in files if re.fullmatch(pattern, path.name)]
                    if not files:
                        raise FileNotFoundError(f'CAS缺少标准月文件：{variable} {stage}')

                stage_files[stage] = files
                print(f"[{stage}] folder={folder} | raw files={len(files)}")

            # %% 4.2 member选择；Nor严格p1/p1，不允许fallback
            available = {}

            for stage in ["co2", "cdr"]:
                members = set()

                for path in stage_files[stage]:
                    match = MEMBER_RE.search(str(path))

                    if match:
                        members.add(match.group(0).lower())

                available[stage] = members

            # Nor：thetao/tas均必须使用新下载的r1i1p1f1
            if model == "NorESM2-LM":
                for stage in ["co2", "cdr"]:
                    if NOR_REQUIRED_MEMBER not in available[stage]:
                        raise FileNotFoundError(
                            f"NorESM2-LM {variable} {stage}必须使用对齐成员{NOR_REQUIRED_MEMBER}；"
                            f"当前可用member={sorted(available[stage])}，禁止回退旧成员"
                        )

                    stage_files[stage] = [
                        path
                        for path in stage_files[stage]
                        if (match := MEMBER_RE.search(str(path))) and match.group(0).lower() == NOR_REQUIRED_MEMBER
                    ]
                    stage_members[stage] = NOR_REQUIRED_MEMBER

            # 其他模式要求CO2/CDR使用共同member
            else:
                if available["co2"] and available["cdr"]:
                    common_members = available["co2"] & available["cdr"]
                    selected_member = next(
                        (candidate for candidate in MEMBER_PRIORITY if candidate in common_members),
                        sorted(common_members)[0] if common_members else None,
                    )

                    if selected_member is None:
                        raise ValueError(
                            f"{model} {variable} CO2/CDR无共同member："
                            f"CO2={sorted(available['co2'])}, CDR={sorted(available['cdr'])}"
                        )

                    for stage in ["co2", "cdr"]:
                        stage_files[stage] = [
                            path
                            for path in stage_files[stage]
                            if (match := MEMBER_RE.search(str(path))) and match.group(0).lower() == selected_member
                        ]
                        stage_members[stage] = selected_member

                elif not available["co2"] and not available["cdr"]:
                    stage_members["co2"] = "unlabeled"
                    stage_members["cdr"] = "unlabeled"

                else:
                    raise ValueError(f"{model} {variable}只有一侧文件带member标签，禁止含糊拼接")

            print(f"[member] CO2={stage_members['co2']} | CDR={stage_members['cdr']}")

            # %% 4.3 读取固定端点并重网格
            stage_result = {}

            for stage in ["co2", "cdr"]:
                files = stage_files[stage]
                ds = None
                open_error = None

                for open_options in [
                    {"combine": "by_coords", "decode_times": TIME_CODER},
                    {"combine": "by_coords", "decode_times": False},
                    {"combine": "nested", "concat_dim": "time", "decode_times": False},
                ]:
                    try:
                        ds = xr.open_mfdataset(
                            [str(path) for path in files],
                            chunks=({'time': 12, **({'lev': 1} if variable == 'thetao' else {})}
                                    if model == 'CAS-ESM2-0' else None),
                            parallel=False,
                            data_vars="minimal",
                            coords="minimal",
                            compat="override",
                            mask_and_scale=True,
                            **open_options,
                        )
                        break

                    except Exception as error:
                        open_error = error

                if ds is None:
                    raise RuntimeError(f"{model} {variable} {stage} NetCDF读取失败：{open_error}")

                try:
                    # %% 4.3.1 变量识别
                    if variable in ds.data_vars:
                        raw = ds[variable]

                    else:
                        target_standard = "sea_water_potential_temperature" if variable == "thetao" else "air_temperature"
                        candidates = []

                        for name, data_array in ds.data_vars.items():
                            score = 12 * (str(data_array.attrs.get("standard_name", "")).lower() == target_standard)
                            score += 8 * (variable in name.lower())

                            if score > 0:
                                candidates.append((score, name))

                        if not candidates:
                            raise KeyError(f"{model} {stage}未找到{variable}")

                        candidates.sort(reverse=True)
                        raw = ds[candidates[0][1]]

                    time_dim = next((dim for dim in raw.dims if dim.lower() == "time"), None)

                    if time_dim is None:
                        raise KeyError(f"{model} {variable} {stage}无time维：{raw.dims}")

                    field = raw

                    # 删除非time单例辅助维
                    for dim in list(field.dims):
                        if dim != time_dim and field.sizes[dim] == 1:
                            field = field.isel({dim: 0}, drop=True)

                    # %% 4.3.2 thetao取最浅有效层
                    surface_note = "not_applicable"

                    if variable == "thetao":
                        non_time_dims = [dim for dim in field.dims if dim != time_dim]

                        if len(non_time_dims) > 2:
                            vertical_candidates = []

                            for dim in non_time_dims:
                                coordinate = field.coords.get(dim)
                                score = 10 if dim.lower() in {"lev", "olevel", "depth", "deptht", "z", "zt", "z_t"} else 0

                                if coordinate is not None:
                                    score += 8 * (str(coordinate.attrs.get("axis", "")).upper() == "Z")
                                    score += 8 * ("depth" in str(coordinate.attrs.get("standard_name", "")).lower())

                                if score > 0:
                                    vertical_candidates.append((score, dim))

                            if not vertical_candidates:
                                raise ValueError(f"{model} {stage} thetao无法识别垂向维：{field.dims}")

                            vertical_candidates.sort(reverse=True)
                            zdim = vertical_candidates[0][1]
                            zcoord = field.coords.get(zdim)
                            zindex = 0

                            if zcoord is not None and zcoord.ndim == 1:
                                zvalues = np.asarray(zcoord.values, dtype=float)
                                finite_z = np.flatnonzero(np.isfinite(zvalues))

                                if finite_z.size:
                                    zindex = int(finite_z[np.argmin(np.abs(zvalues[finite_z]))])

                                surface_note = f"{zdim}[{zindex}]={zvalues[zindex]:g} {zcoord.attrs.get('units', '')}".strip()

                            else:
                                surface_note = f"{zdim}[0]"

                            field = field.isel({zdim: zindex}, drop=True)

                        else:
                            surface_note = "thetao_no_vertical_dimension"

                    spatial_dims = [dim for dim in field.dims if dim != time_dim]

                    if len(spatial_dims) != 2:
                        raise ValueError(f"{model} {variable} {stage} 表层选择后不是time+二维空间：{field.dims}")

                    field = field.astype("float64")
                    field = field.where((field > -150.0) & (field < 450.0))

                    # %% 4.4 时间频率与固定窗口
                    time_count = int(field.sizes[time_dim])
                    steps_per_year = None
                    time_method = None

                    if model == 'CAS-ESM2-0':
                        steps_per_year = 12
                        time_method = 'decoded_monthly'
                        if ds.attrs['source_id'] != model or ds.attrs['variant_label'] != 'r1i1p1f1':
                            raise ValueError('CAS标准月文件内部模式或成员不匹配')
                        expected_experiment = '1pctCO2' if stage == 'co2' else '1pctCO2-cdr'
                        if ds.attrs['experiment_id'] != expected_experiment:
                            raise ValueError('CAS标准月文件内部实验不匹配')

                    else:
                        try:
                            years = np.asarray(ds[time_dim].dt.year.values, dtype=int)
                            _, counts = np.unique(years, return_counts=True)
                            typical = int(np.median(counts)) if counts.size else 0

                            if typical >= 10:
                                steps_per_year = 12
                                time_method = "decoded_monthly"

                            elif typical == 1:
                                steps_per_year = 1
                                time_method = "decoded_annual"

                        except Exception:
                            pass

                        if steps_per_year is None:
                            metadata_text = " ".join([
                                str(ds.attrs.get("frequency", "")),
                                str(ds.attrs.get("table_id", "")),
                                str(field.attrs.get("frequency", "")),
                                " ".join(path.name for path in files[:5]),
                            ]).lower()

                            if any(token in metadata_text for token in ("oyr", "ayr", "year", "annual")):
                                steps_per_year = 1
                                time_method = "metadata_annual"

                            elif any(token in metadata_text for token in ("omon", "amon", "month", "monthly")):
                                steps_per_year = 12
                                time_method = "metadata_monthly"

                            else:
                                raise ValueError(f"{model} {variable} {stage} 无法可靠确定时间频率")

                    # CO2始终是实验年1–10
                    if stage == "co2":
                        start_year_index = CO2_START_YEAR_INDEX
                        experiment_window = "1pctCO2 years 1-10"

                    # 标准CDR月文件统一使用本实验第131–140年。
                    else:
                        start_year_index = CDR_START_YEAR_INDEX
                        experiment_window = "CDR years 131-140"

                    start_step = start_year_index * steps_per_year
                    stop_step = (start_year_index + WINDOW_YEARS) * steps_per_year

                    if time_count < stop_step:
                        raise ValueError(
                            f"{model} {variable} {stage}仅{time_count}步，固定窗口需要至少{stop_step}步；禁止使用文件尾部"
                        )

                    # %% 4.4.1 连续性检查
                    prefix_time = ds[time_dim].isel({time_dim: slice(0, stop_step)})
                    continuity = None

                    try:
                        years = np.asarray(prefix_time.dt.year.values, dtype=int)

                        if steps_per_year == 12:
                            months = np.asarray(prefix_time.dt.month.values, dtype=int)
                            serial = years.astype(np.int64) * 12 + months.astype(np.int64)

                            if serial.size > 1 and not np.all(np.diff(serial) == 1):
                                raise ValueError(f"{model} {variable} {stage}固定窗口前存在缺月")

                            continuity = "decoded_calendar_monthly"

                        else:
                            if years.size > 1 and not np.all(np.diff(years) == 1):
                                raise ValueError(f"{model} {variable} {stage}固定窗口前存在缺年")

                            continuity = "decoded_calendar_annual"

                    except (AttributeError, TypeError, OverflowError):
                        time_values = np.asarray(prefix_time.values)

                        if not np.issubdtype(time_values.dtype, np.number):
                            raise ValueError(f"{model} {variable} {stage} 时间轴无法验证连续性")

                        differences = np.diff(time_values.astype(float))
                        finite_diff = differences[np.isfinite(differences)]

                        if finite_diff.size != differences.size or finite_diff.size == 0:
                            raise ValueError(f"{model} {variable} {stage} 数值时间轴含NaN或长度不足")

                        median_step = float(np.median(finite_diff))

                        if (
                            median_step <= 0.0
                            or np.any(
                                (differences <= 0.0)
                                | (differences < 0.5 * median_step)
                                | (differences > 1.5 * median_step)
                            )
                        ):
                            raise ValueError(f"{model} {variable} {stage} 数值时间轴不连续")

                        continuity = "numeric_regular_prefix"

                    endpoint = field.isel({time_dim: slice(start_step, stop_step)}).mean(dim=time_dim, skipna=True).load()

                    # %% 4.5 单位统一到degC
                    units = str(raw.attrs.get("units", "")).strip().lower().replace(" ", "")

                    if units in {"k", "kelvin", "degrees_k", "degree_k"}:
                        endpoint = endpoint - 273.15
                        unit_method = "kelvin_to_degC"

                    elif units in {"degc", "degree_celsius", "degrees_celsius", "celsius", "degreec", "degreesc"}:
                        unit_method = "degC"

                    else:
                        raise ValueError(f"{model} {variable} {stage} 温度单位无法确认：{raw.attrs.get('units', '')!r}")

                    values = np.asarray(endpoint.values, dtype=float)
                    shape = endpoint.shape

                    # %% 4.6 经纬度
                    coordinates = {}

                    for kind in ["lat", "lon"]:
                        if kind == "lat":
                            standard_target = "latitude"
                            axis_target = "Y"
                            names = {"lat", "latitude", "nav_lat", "geolat", "tlat", "ulat"}

                        else:
                            standard_target = "longitude"
                            axis_target = "X"
                            names = {"lon", "longitude", "nav_lon", "geolon", "tlon", "ulon"}

                        candidates = []

                        for name in ds.variables:
                            coordinate = ds[name].squeeze(drop=True)

                            if coordinate.ndim not in [1, 2]:
                                continue

                            score = 12 * (str(coordinate.attrs.get("standard_name", "")).lower() == standard_target)
                            score += 5 * (str(coordinate.attrs.get("axis", "")).upper() == axis_target)
                            score += 6 * (name.lower() in names)
                            score += 4 * ("degree" in str(coordinate.attrs.get("units", "")).lower())

                            if score > 0:
                                candidates.append((score, name))

                        candidates.sort(reverse=True)
                        coordinate_2d = None

                        for _, name in candidates:
                            coordinate_values = np.asarray(ds[name].squeeze(drop=True).values, dtype=float)

                            if coordinate_values.ndim == 1 and coordinate_values.size == shape[0]:
                                coordinate_2d = np.repeat(coordinate_values[:, None], shape[1], axis=1)

                            elif coordinate_values.ndim == 1 and coordinate_values.size == shape[1]:
                                coordinate_2d = np.repeat(coordinate_values[None, :], shape[0], axis=0)

                            elif coordinate_values.ndim == 2 and coordinate_values.shape == shape:
                                coordinate_2d = coordinate_values

                            elif coordinate_values.ndim == 2 and coordinate_values.T.shape == shape:
                                coordinate_2d = coordinate_values.T

                            if coordinate_2d is not None:
                                break

                        if coordinate_2d is None:
                            raise KeyError(f"{model} {variable} {stage}无法识别{kind}坐标匹配shape={shape}")

                        coordinates[kind] = coordinate_2d

                    source_lat = coordinates["lat"]
                    source_lon = coordinates["lon"]

                    bad_coordinate = (
                        ~np.isfinite(source_lat)
                        | ~np.isfinite(source_lon)
                        | (source_lat < -90.0)
                        | (source_lat > 90.0)
                        | (np.abs(source_lon) > 1.0e6)
                    )

                    source_lon = np.where(bad_coordinate, np.nan, ((source_lon + 180.0) % 360.0) - 180.0)
                    source_lat = np.where(bad_coordinate, np.nan, source_lat)

                    # %% 4.7 球面IDW到共同1°北极网格
                    geo = (
                        np.isfinite(source_lat)
                        & np.isfinite(source_lon)
                        & (source_lat >= 40.0)
                        & (source_lat <= 90.0)
                    )

                    valid = geo & np.isfinite(values)

                    if int(valid.sum()) < 50:
                        raise ValueError(f"{model} {variable} {stage}可重网格有效点仅{int(valid.sum())}")

                    # 原生网格几何支持，用于防止跨陆地外插
                    geo_latr = np.deg2rad(source_lat[geo])
                    geo_lonr = np.deg2rad(source_lon[geo])

                    geo_xyz = np.column_stack([
                        np.cos(geo_latr) * np.cos(geo_lonr),
                        np.cos(geo_latr) * np.sin(geo_lonr),
                        np.sin(geo_latr),
                    ])

                    support_tree = cKDTree(geo_xyz)
                    _, support_index = support_tree.query(TARGET_XYZ, k=1)
                    support_valid = np.isfinite(values[geo])[support_index]

                    valid_latr = np.deg2rad(source_lat[valid])
                    valid_lonr = np.deg2rad(source_lon[valid])

                    valid_xyz = np.column_stack([
                        np.cos(valid_latr) * np.cos(valid_lonr),
                        np.cos(valid_latr) * np.sin(valid_lonr),
                        np.sin(valid_latr),
                    ])

                    valid_values = values[valid]
                    k = min(IDW_NEIGHBOURS, valid_values.size)
                    tree = cKDTree(valid_xyz)
                    chord, source_index = tree.query(TARGET_XYZ, k=k)

                    if k == 1:
                        chord = chord[:, None]
                        source_index = source_index[:, None]

                    chord = np.clip(np.asarray(chord, dtype=float), 0.0, 2.0)
                    distance = 2.0 * EARTH_RADIUS_KM * np.arcsin(chord / 2.0)
                    neighbour_values = valid_values[source_index]
                    idw_weights = 1.0 / np.maximum(distance, 1.0) ** IDW_POWER
                    regridded = np.sum(idw_weights * neighbour_values, axis=1) / np.sum(idw_weights, axis=1)

                    exact = distance[:, 0] < 1.0e-6
                    regridded[exact] = neighbour_values[exact, 0]
                    regridded[(distance[:, 0] > MAX_REGRID_DISTANCE_KM) | (~support_valid)] = np.nan
                    regridded = regridded.reshape(TARGET_LAT_GRID.shape)

                    stage_result[stage] = {
                        "field": regridded,
                        "member": stage_members[stage],
                        "file_count": len(files),
                        "time_count": time_count,
                        "steps_per_year": steps_per_year,
                        "time_method": time_method,
                        "continuity": continuity,
                        "source_window_start_step": start_step,
                        "source_window_stop_step": stop_step,
                        "experiment_window": experiment_window,
                        "surface": surface_note,
                        "units": str(raw.attrs.get("units", "")),
                        "unit_method": unit_method,
                        "valid_target_cells": int(np.isfinite(regridded).sum()),
                        "nearest_distance_p95_km": float(
                            np.nanpercentile(
                                distance[:, 0][distance[:, 0] <= MAX_REGRID_DISTANCE_KM],
                                95.0,
                            )
                        ),
                        "source_files": ";".join(path.name for path in files),
                    }

                    print(
                        f"[{stage}] steps={time_count} | {steps_per_year}/yr | "
                        f"source index={start_step}:{stop_step} | {experiment_window} | "
                        f"valid target={stage_result[stage]['valid_target_cells']}"
                    )

                finally:
                    ds.close()

            # %% 4.8 raw endpoint difference
            co2_field = stage_result["co2"]["field"]
            cdr_field = stage_result["cdr"]["field"]
            common = np.isfinite(co2_field) & np.isfinite(cdr_field)

            if int(common.sum()) < 300:
                raise ValueError(f"{model} {variable} CO2/CDR共同1°有效格点仅{int(common.sum())}")

            raw_delta = np.where(common, cdr_field - co2_field, np.nan)
            temperature_fields[variable][model] = raw_delta

            domain_weight = np.cos(np.deg2rad(TARGET_LAT_GRID))
            domain_valid = common & np.isfinite(domain_weight) & (domain_weight > 0.0)
            domain_mean = float(
                np.sum(raw_delta[domain_valid] * domain_weight[domain_valid])
                / np.sum(domain_weight[domain_valid])
            )

            qa_rows.append({
                "variable": variable,
                "model": model,
                "co2_member": stage_result["co2"]["member"],
                "cdr_member": stage_result["cdr"]["member"],
                "co2_file_count": stage_result["co2"]["file_count"],
                "cdr_file_count": stage_result["cdr"]["file_count"],
                "co2_source_files": stage_result["co2"]["source_files"],
                "cdr_source_files": stage_result["cdr"]["source_files"],
                "co2_time_steps": stage_result["co2"]["time_count"],
                "cdr_time_steps": stage_result["cdr"]["time_count"],
                "co2_steps_per_year": stage_result["co2"]["steps_per_year"],
                "cdr_steps_per_year": stage_result["cdr"]["steps_per_year"],
                "co2_source_window_start_step_0based": stage_result["co2"]["source_window_start_step"],
                "co2_source_window_stop_step_exclusive_0based": stage_result["co2"]["source_window_stop_step"],
                "cdr_source_window_start_step_0based": stage_result["cdr"]["source_window_start_step"],
                "cdr_source_window_stop_step_exclusive_0based": stage_result["cdr"]["source_window_stop_step"],
                "co2_experiment_window": stage_result["co2"]["experiment_window"],
                "cdr_experiment_window": stage_result["cdr"]["experiment_window"],
                "co2_time_method": stage_result["co2"]["time_method"],
                "cdr_time_method": stage_result["cdr"]["time_method"],
                "co2_continuity": stage_result["co2"]["continuity"],
                "cdr_continuity": stage_result["cdr"]["continuity"],
                "co2_surface": stage_result["co2"]["surface"],
                "cdr_surface": stage_result["cdr"]["surface"],
                "co2_units": stage_result["co2"]["units"],
                "cdr_units": stage_result["cdr"]["units"],
                "co2_valid_target_cells": stage_result["co2"]["valid_target_cells"],
                "cdr_valid_target_cells": stage_result["cdr"]["valid_target_cells"],
                "common_target_cells": int(common.sum()),
                "raw_delta_45_90N_mean_degC": domain_mean,
            })

            print(f"[完成] {variable} {model}: raw ΔT 45–90N={domain_mean:+.3f} °C")

        except Exception as error:
            failure_rows.append({
                "variable": variable,
                "model": model,
                "error_type": type(error).__name__,
                "error": str(error),
            })

            print(f"[失败] {variable} {model}: {type(error).__name__}: {error}")

# %% 5. 正式完整性检查
failures = pd.DataFrame(failure_rows)

if failure_rows:
    raise RuntimeError(f'thetao/tas计算失败，正式结果未改动：{failure_rows}')

for variable in ["thetao", "tas"]:
    missing_models = [model for model in MODELS if model not in temperature_fields[variable]]

    if missing_models:
        raise RuntimeError(f"{variable}缺少正式模式：{missing_models}")

print("[完整性] thetao=9/9，tas=9/9")
print(f"[Nor] thetao/tas CO2/CDR均强制 {NOR_REQUIRED_MEMBER}")

# %% 6. 合并到当前SIA网格
cell_frames = []

for model in MODELS:
    model_sia = sia_cells[sia_cells["model"] == model].sort_values(
        ["latitude_deg_north", "longitude_deg_east"]
    ).reset_index(drop=True)

    if len(model_sia) != TARGET_LAT.size * TARGET_LON.size:
        raise ValueError(
            f"当前SIA网格 {model} 应为{TARGET_LAT.size * TARGET_LON.size}格，"
            f"实际={len(model_sia)}"
        )

    actual_lat = model_sia["latitude_deg_north"].to_numpy(dtype=float)
    actual_lon = model_sia["longitude_deg_east"].to_numpy(dtype=float)

    if (
        not np.allclose(actual_lat, TARGET_LAT_GRID.ravel(), rtol=0.0, atol=1.0e-8)
        or not np.allclose(actual_lon, TARGET_LON_GRID.ravel(), rtol=0.0, atol=1.0e-8)
    ):
        raise ValueError(
            f"{model} 当前SIA格点坐标与温度目标1°网格不一致，禁止按行强行拼接"
        )

    model_frame = model_sia.copy()
    model_frame["raw_delta_thetao_degC"] = temperature_fields["thetao"][model].ravel()
    model_frame["raw_delta_tas_degC"] = temperature_fields["tas"][model].ravel()
    cell_frames.append(model_frame)

temperature_cells = pd.concat(cell_frames, ignore_index=True)

# %% 7. 固定地理区域 + 九模式共同有效格点的温度平均
common_supports = {}

for variable in ["thetao", "tas"]:
    common_support = np.ones(TARGET_LAT_GRID.shape, dtype=bool)

    for model in MODELS:
        common_support &= np.isfinite(temperature_fields[variable][model])

    common_supports[variable] = common_support

    print(
        f"[共同支持] {variable}: "
        f"{int(common_support.sum())}/{common_support.size} target cells"
    )

region_masks = {
    "North Atlantic": (
        (TARGET_LAT_GRID >= 45.0)
        & (TARGET_LAT_GRID < 70.0)
        & (TARGET_LON_GRID >= -60.0)
        & (TARGET_LON_GRID < 30.0)
    ),
    "Central Arctic": (
        (TARGET_LAT_GRID >= 70.0)
        & (TARGET_LAT_GRID <= 90.0)
    ),
    "North Pacific": (
        (TARGET_LAT_GRID >= 45.0)
        & (TARGET_LAT_GRID < 70.0)
        & ((TARGET_LON_GRID >= 120.0) | (TARGET_LON_GRID < -150.0))
    ),
}

abbreviations = {
    "North Atlantic": "NA",
    "Central Arctic": "CA",
    "North Pacific": "NP",
}

grid_weights = np.cos(np.deg2rad(TARGET_LAT_GRID))
summary_rows = []

for region, region_mask in region_masks.items():
    thetao_count = int((region_mask & common_supports["thetao"]).sum())
    tas_count = int((region_mask & common_supports["tas"]).sum())

    if thetao_count < 10:
        raise RuntimeError(f"{region} thetao九模式共同有效格点过少：{thetao_count}")

    if tas_count < 10:
        raise RuntimeError(f"{region} tas九模式共同有效格点过少：{tas_count}")

    for model in MODELS:
        thetao = temperature_fields["thetao"][model]
        tas = temperature_fields["tas"][model]

        thetao_valid = (
            region_mask
            & common_supports["thetao"]
            & np.isfinite(thetao)
            & np.isfinite(grid_weights)
            & (grid_weights > 0.0)
        )

        tas_valid = (
            region_mask
            & common_supports["tas"]
            & np.isfinite(tas)
            & np.isfinite(grid_weights)
            & (grid_weights > 0.0)
        )

        thetao_mean = float(
            np.sum(thetao[thetao_valid] * grid_weights[thetao_valid])
            / np.sum(grid_weights[thetao_valid])
        )

        tas_mean = float(
            np.sum(tas[tas_valid] * grid_weights[tas_valid])
            / np.sum(grid_weights[tas_valid])
        )

        summary_rows.append({
            "model": model,
            "region": region,
            "region_abbreviation": abbreviations[region],
            "raw_delta_thetao_degC": thetao_mean,
            "raw_delta_tas_degC": tas_mean,
            "temperature_sampling_method": "fixed_geographic_region_and_9model_common_valid_support",
            "thetao_common_valid_grid_cell_count": int(thetao_valid.sum()),
            "tas_common_valid_grid_cell_count": int(tas_valid.sum()),
        })

summary = pd.DataFrame(summary_rows)

if len(summary) != 27:
    raise RuntimeError(f"区域summary应为27行，实际={len(summary)}")

print("\n[区域共同支持格点]")

for region, region_mask in region_masks.items():
    print(
        f"{region}: "
        f"thetao={int((region_mask & common_supports['thetao']).sum())}, "
        f"tas={int((region_mask & common_supports['tas']).sum())}"
    )

# %% 8. 覆盖正式输出
archive_dir = ROOT / 'archive' / (datetime.now().strftime('%Y%m%d%H%M%S') + '_cas_arctic_temperature_monthly')
for old_file in [CELL_FILE, SUMMARY_FILE, QA_FILE, FAIL_FILE, META_FILE]:
    if old_file.exists():
        archived_file = archive_dir / old_file.relative_to(ROOT)
        archived_file.parent.mkdir(parents=True, exist_ok=True)
        old_hash = hashlib.sha256(old_file.read_bytes()).hexdigest()
        shutil.move(str(old_file), str(archived_file))
        assert not old_file.exists() and hashlib.sha256(archived_file.read_bytes()).hexdigest() == old_hash

qa = pd.DataFrame(qa_rows)

empty_failures = pd.DataFrame(
    columns=["variable", "model", "error_type", "error"]
)

temperature_cells.to_csv(
    CELL_FILE,
    index=False,
    encoding="utf-8-sig",
    float_format="%.8g",
)

summary.to_csv(
    SUMMARY_FILE,
    index=False,
    encoding="utf-8-sig",
    float_format="%.8f",
)

qa.to_csv(
    QA_FILE,
    index=False,
    encoding="utf-8-sig",
    float_format="%.8g",
)

empty_failures.to_csv(
    FAIL_FILE,
    index=False,
    encoding="utf-8-sig",
)

metadata = {
    "adjustment": "none",
    "definition": "raw delta T = CDR experiment years 131-140 minus 1pctCO2 experiment years 1-10",
    "ordinary_co2_year_index_0based": [0, 10],
    "ordinary_cdr_year_index_0based": [130, 140],
    "use_file_last10": False,
    "thetao_models": THETAO_MODELS,
    "tas_models": TAS_MODELS,
    "tas_excluded_models": [],
    "thetao_surface": "shallowest valid ocean layer",
    "NorESM2-LM": {
        "required_member": NOR_REQUIRED_MEMBER,
        "co2_member": NOR_REQUIRED_MEMBER,
        "cdr_member": NOR_REQUIRED_MEMBER,
        "fallback_allowed": False,
    },
    "CAS_thetao_CDR": {
        "file_rule": "standard native monthly Omon files; exclude annual/regridded derivatives",
        "source_index_0based": [1560, 1680],
        "selected_calendar_years": "0131-0140",
        "experiment_window": "CDR years 131-140",
    },
    "CAS_tas_CDR": {
        "file_rule": "standard 2 m tas_Amon monthly file",
        "source_index_0based": [1560, 1680],
        "selected_calendar_years": "0470-0479",
        "experiment_window": "CDR years 131-140",
    },
    "target_grid": "45.5-89.5N x -179.5-179.5E, 1 degree",
    "regrid": {
        "method": "spherical inverse-distance weighting",
        "neighbours": IDW_NEIGHBOURS,
        "power": IDW_POWER,
        "max_nearest_distance_km": MAX_REGRID_DISTANCE_KM,
    },
    "temperature_region_sampling": (
        "fixed geographic region AND variable-specific 9-model common valid "
        "target-grid support; no SIC mask"
    ),
    "common_valid_target_cells": {
        "thetao": int(common_supports["thetao"].sum()),
        "tas": int(common_supports["tas"].sum()),
    },
    "regions": {
        "Central Arctic": "70-90N, all longitudes",
        "North Atlantic": "45-70N, 60W-30E",
        "North Pacific": "45-70N, 120E-150W across dateline",
    },
    "overwrite_existing_outputs": True,
    "recalculated_models": RECALCULATE_MODELS,
    "previous_outputs_archive": str(archive_dir),
}

META_FILE.write_text(
    json.dumps(metadata, ensure_ascii=False, indent=2),
    encoding="utf-8",
)

# %% 9. 终端汇总
print("\n" + "=" * 110)
print("Arctic thetao / tas fixed-region recalculation complete")
print("=" * 110)

print(
    summary.pivot(
        index="model",
        columns="region",
        values=["raw_delta_thetao_degC", "raw_delta_tas_degC"],
    ).reindex(MODELS).to_string()
)

print(f"\n[覆盖 cells]    {CELL_FILE}")
print(f"[覆盖 summary]  {SUMMARY_FILE}")
print(f"[覆盖 QA]       {QA_FILE}")
print(f"[覆盖 failures] {FAIL_FILE}")
print(f"[覆盖 metadata] {META_FILE}")
