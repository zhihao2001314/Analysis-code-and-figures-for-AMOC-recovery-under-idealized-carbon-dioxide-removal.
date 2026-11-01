# -*- coding: utf-8 -*-
"""
AMOC 9模式重新计算
【是否需要运行】需要
【输入文件】input/cmip6amocad、input/cmip6amoc、input/cmip6amoc-picontrol；CNRM额外读取同模式nav_lat/nav_lon
【输出文件】output/amoc/CMIP6_AMOC_9models_*.csv；output/amoc_figure1/amoc_figure1_9models_*
【什么时候需要重新运行】AMOC原始数据、member、piControl、CNRM网格映射或AMOC定义变化后
定义：Atlantic-Arctic；固定40°N最近模式纬度；depth >=500 m；先年度平均流函数剖面再取深度极值；
piControl分支对齐线性订正；11年中心滑动平均；NorESM2-LM使用p1/p1/p1。
"""


# %% 1. 导入依赖
from datetime import datetime
from operator import itemgetter
from pathlib import Path
import json
import re
import cftime
import numpy as np
import pandas as pd
import xarray as xr


# %% 2. 参数与路径
ROOT = Path(__file__).resolve().parents[2]
UP_ROOT, CDR_ROOT, PI_ROOT = ROOT/"input"/"cmip6amocad", ROOT/"input"/"cmip6amoc", ROOT/"input"/"cmip6amoc-picontrol"
OUT = ROOT / "output" / "amoc_s01_n01"/"l_40"
CNRM_CACHE = OUT/"cnrm_grid_cache"/"CNRM-ESM2-1_j294_latitude_mapping_verified_v4.npz"

MODELS = ["ACCESS-ESM1-5","CanESM5","CESM2","CNRM-ESM2-1","GFDL-ESM4","MIROC-ES2L","NorESM2-LM","UKESM1-0-LL","CAS-ESM2-0"]
DIRS = {"ACCESS-ESM1-5":"ACCESSESM5","CanESM5":"CanESM5","CESM2":"CESM2","CNRM-ESM2-1":"CNRM-ESM2","GFDL-ESM4":"GFDL","MIROC-ES2L":"MIROC-ES2L","NorESM2-LM":"NORESM2-LM2","UKESM1-0-LL":"UKESM1-0-LL","CAS-ESM2-0":"CAS"}
VARS = {"ACCESS-ESM1-5":"msftmz","CanESM5":"msftmz","CESM2":"msftmz","CNRM-ESM2-1":"msftyz","GFDL-ESM4":"msftyz","MIROC-ES2L":"msftmz","NorESM2-LM":"msftmz","UKESM1-0-LL":"msftyz","CAS-ESM2-0":"msftmz"}
STAGES = [("piControl",PI_ROOT),("1pctCO2",UP_ROOT),("1pctCO2-cdr",CDR_ROOT)]

TARGET_LAT, MIN_DEPTH, RAMP_YEARS, MAX_YEARS, SMOOTH = 40, 500.0, 140, 301, 11
CO2_0, START_YEAR = 284.7, 1850
TIME_CODER = xr.coders.CFDatetimeCoder(use_cftime=True)

ANNUAL = OUT/"CMIP6_AMOC_9models_annual_series.csv"
META = OUT/"CMIP6_AMOC_9models_processing_metadata.csv"
FAILED = OUT/"CMIP6_AMOC_9models_failures.csv"
PLOT_CSV = OUT/"amoc_figure1_9models_plot_data.csv"
PLOT_JSON = OUT/"amoc_figure1_9models_metadata.json"

CNRM_GRID_ROOTS = [
    ROOT/"input"/"areacello", ROOT/"input"/"cmipseaicead", ROOT/"input"/"cmipseaice", ROOT/"input"/"cmipseaice-picontrol", UP_ROOT, CDR_ROOT, PI_ROOT,
    Path(r"F:\sea\areacello"), Path(r"F:\sea\cmipseaicead"), Path(r"F:\sea\cmipseaice"), Path(r"F:\sea\cmipseaice-picontrol"),
    Path(r"F:\sea\cmip6amocad"), Path(r"F:\sea\cmip6amoc"), Path(r"F:\sea\cmip6amoc-picontrol"),
]


# %% 3. 归档旧结果
archive = ROOT/"archive"/f"{datetime.now():%Y%m%d%H%M}_amoc_9model_recompute"

for path in (ANNUAL,META,FAILED,PLOT_CSV,PLOT_JSON):
    if path.exists():
        target = archive/path.relative_to(ROOT)
        target.parent.mkdir(parents=True,exist_ok=True)
        path.replace(target)

OUT.mkdir(parents=True,exist_ok=True)

CNRM_CACHE.parent.mkdir(parents=True,exist_ok=True)

results, metadata_rows, failure_rows = {}, [], []


# %% 4. 逐模式重新计算
for model in MODELS:
    print(f"\n{'='*84}\n{model} | fixed {TARGET_LAT:.1f}°N\n{'='*84}")

    try:
        raw_files, groups = {}, {}

        # 三阶段member选择；Nor固定使用补齐后的r1i1p1f1。
        for stage, root in STAGES:
            folder = root/DIRS[model]

            if not folder.is_dir():
                raise FileNotFoundError(f"{model} {stage}目录不存在：{folder}")

            files = sorted(folder.rglob("*.nc"))
            amoc_files = [p for p in files if any(x in p.name.lower() for x in ("msftmz","msftyz","amoc","moc"))]
            files = amoc_files or files

            if not files:
                raise FileNotFoundError(f"{model} {stage}目录没有NetCDF")

            raw_files[stage], groups[stage] = files, {}

            for path in files:
                match = re.search(r"r(\d+)i(\d+)p(\d+)f(\d+)",path.as_posix(),re.I)

                if match:
                    label = tuple(int(x) for x in match.groups())
                    groups[stage].setdefault(label,[]).append(path)

        if model == "NorESM2-LM":
            selected = (1,1,1,1)

            for stage, _ in STAGES:
                if selected not in groups[stage]:
                    available = [f"r{r}i{i}p{p}f{f}" for r,i,p,f in sorted(groups[stage])]
                    raise ValueError(f"NorESM2-LM {stage}必须有r1i1p1f1；可用={available}")

        else:
            labelled = [stage for stage,_ in STAGES if groups[stage]]
            common = set(groups[labelled[0]]) if labelled else set()

            for stage in labelled[1:]:
                common &= set(groups[stage])

            if not common:
                available = {stage:[f"r{r}i{i}p{p}f{f}" for r,i,p,f in sorted(groups[stage])] for stage in labelled}
                raise ValueError(f"{model}三阶段没有共同member：{available}")

            selected = sorted(common)[0]

        member = f"r{selected[0]}i{selected[1]}p{selected[2]}f{selected[3]}"
        stage_files = {stage:(sorted(groups[stage][selected]) if selected in groups[stage] else sorted(raw_files[stage])) for stage,_ in STAGES}

        print(f"member={member}")

        series, stage_meta = {}, {}
        basin_from_pi, cnrm_orientation = None, "max"

        for stage, _ in STAGES:
            files = stage_files[stage]
            default_var = VARS[model]
            priority = ("amoc","msftmz","msftyz","moc","MOC") if model=="CAS-ESM2-0" and stage=="1pctCO2-cdr" else (default_var,"msftmz","msftyz","amoc","moc","MOC")

            print(f"[{stage}] files={len(files)}")

            with xr.open_mfdataset(
                [str(p) for p in files],
                combine="by_coords",
                decode_times=TIME_CODER,
                data_vars="minimal",
                coords="minimal",
                compat="override",
                parallel=False,
                chunks={"time":120},
            ) as ds:

                variable_name = next((name for name in priority if name in ds.data_vars),None)

                if variable_name is None:
                    raise KeyError(f"{model} {stage}未找到AMOC流函数变量")

                psi = ds[variable_name].squeeze().astype("float64")
                psi = psi.where(np.isfinite(psi)&(np.abs(psi)<1e20))

                if int(psi.count().compute().item()) == 0:
                    raise ValueError(f"{model} {stage} {variable_name}全部NaN")

                time_dim = next((d for d in psi.dims if d.lower()=="time"),None)

                if time_dim is None:
                    time_dim = next((d for d in psi.dims if d in psi.coords and str(psi[d].attrs.get("axis","")).upper()=="T"),None)

                if time_dim is None:
                    raise ValueError(f"{model} {stage}无法识别时间维：{psi.dims}")

                variant = str(ds.attrs.get("variant_label","")).strip().lower() or member

                if variant != member.lower() and not (model=="CAS-ESM2-0" and stage=="1pctCO2-cdr"):
                    raise ValueError(f"{model} {stage} metadata member={variant}，选择={member}")

                # Atlantic-Arctic：恢复原代码requested/标签/flag识别，不按幅度猜测。
                basin_dim = next((d for d in psi.dims if "basin" in d.lower() or "region" in d.lower() or "sector" in d.lower() or d.lower() in {"transport_reg","transport_region"}),None)
                basin_index, basin_text = None, "no basin dimension; already basin-selected"

                if basin_dim is not None:
                    nbasin = int(psi.sizes[basin_dim])

                    if stage!="piControl" and basin_from_pi is not None:
                        basin_index, basin_text = basin_from_pi, f"reuse piControl basin[{basin_from_pi}]"

                    elif model=="CESM2" and basin_dim.lower() in {"transport_reg","transport_region"} and nbasin==2:
                        basin_index, basin_text = 1, "CESM2 POP transport_reg[1] Atlantic-Arctic"

                    else:
                        candidates, seen = [], set()

                        for name in ("transport_regions","transport_region","sector","basin","region",basin_dim):
                            if name in ds and ds[name].ndim in (1,2) and ds[name].dims[0]==basin_dim:
                                candidates.append((name,ds[name]))
                                seen.add(name)

                        for name, coord in ds.variables.items():
                            if name not in seen and coord.ndim in (1,2) and coord.dims[0]==basin_dim:
                                candidates.append((name,coord))
                                seen.add(name)

                        for name, coord in candidates:
                            requested = str(coord.attrs.get("requested","")).lower()
                            match = re.search(r"(?:atlantic_arctic_ocean|atlantic_ocean)\s*=\s*(-?\d+)",requested)

                            if match:
                                value = int(match.group(1))
                                flat = np.asarray(coord.values).reshape(-1)

                                if flat.size==nbasin and np.issubdtype(flat.dtype,np.number):
                                    pos = np.where(np.isclose(flat.astype(float),value))[0]
                                    basin_index = int(pos[0]) if pos.size else None

                                elif 0<=value<nbasin:
                                    basin_index = value

                                if basin_index is not None:
                                    basin_text = f"{name}:requested Atlantic-Arctic"
                                    break

                            names = re.findall(r"(?:atlantic_arctic_ocean|atlantic_ocean|indian_pacific_ocean|global_ocean)",requested)
                            hits = [i for i,text in enumerate(names[:nbasin]) if text in {"atlantic_arctic_ocean","atlantic_ocean"}]

                            if hits:
                                basin_index, basin_text = hits[0], f"{name}:requested Atlantic-Arctic"
                                break

                            values = np.asarray(coord.values)
                            labels = []

                            if values.dtype.kind in {"U","S","O"}:
                                if values.ndim==1 and values.size==nbasin:
                                    labels = [
                                        x.decode("utf-8","ignore").strip().lower()
                                        if isinstance(x,(bytes,np.bytes_))
                                        else str(x).strip().lower()
                                        for x in values
                                    ]

                                elif values.ndim==2 and values.shape[0]==nbasin:
                                    labels = [
                                        "".join(
                                            x.decode("utf-8","ignore")
                                            if isinstance(x,(bytes,np.bytes_))
                                            else str(x)
                                            for x in row
                                        ).strip().lower()
                                        for row in values
                                    ]

                            hits = [i for i,text in enumerate(labels) if "atlantic" in text or text.startswith("atl") or "atl_" in text]

                            if hits:
                                basin_index, basin_text = hits[0], f"{name}:{labels[hits[0]]}"
                                break

                            meanings = str(coord.attrs.get("flag_meanings",""))
                            flags = coord.attrs.get("flag_values")

                            if meanings and flags is not None and values.ndim==1:
                                for meaning, flag in zip(meanings.split(),np.atleast_1d(flags)):
                                    pos = np.where(values==flag)[0]

                                    if "atl" in meaning.lower() and pos.size:
                                        basin_index, basin_text = int(pos[0]), f"{name}-flag:{meaning}"
                                        break

                                if basin_index is not None:
                                    break

                        if basin_index is None and model=="CESM2" and nbasin>=3:
                            basin_index, basin_text = 0, "CESM2 verified basin[0] Atlantic-Arctic"

                        if basin_index is None and nbasin==1:
                            basin_index, basin_text = 0, "single basin"

                        if basin_index is None:
                            raise ValueError(f"{model} {stage}无法从metadata识别Atlantic-Arctic basin；禁止按幅度猜测")

                    psi = psi.isel({basin_dim:basin_index})

                    if stage=="piControl":
                        basin_from_pi = basin_index

                # CESM2 POP总residual MOC。
                comp_dim = next((d for d in psi.dims if re.sub(r"[^A-Z0-9]","",d.upper()) in {"MOCCOMP","MOCCOMPONENT","MOCCOMPONENTS"}),None)
                component_text = "none"

                if comp_dim:
                    if model=="CESM2":
                        valid = psi.notnull().any(comp_dim)
                        component_text = f"sum({comp_dim}, n={psi.sizes[comp_dim]})"
                        psi = psi.sum(comp_dim,skipna=True).where(valid)

                    elif psi.sizes[comp_dim]==1:
                        psi = psi.isel({comp_dim:0})
                        component_text = f"{comp_dim}[0]"

                    else:
                        raise ValueError(f"{model} {stage}未验证MOC分量：{comp_dim}={psi.sizes[comp_dim]}")

                # 纬度：恢复原find_latitude；CNRM cache不存在时自动重建。
                lat_candidates = []
                variables = dict(ds.variables)
                variables.update({name:psi.coords[name] for name in psi.coords})

                for name, coord in variables.items():
                    coord = coord.squeeze()

                    score = 10 if str(coord.attrs.get("standard_name","")).lower()=="latitude" else 0
                    score += 6 if str(coord.attrs.get("axis","")).upper()=="Y" else 0
                    score += 5 if name.lower() in {"lat","latitude","nav_lat","geolat","rlat","lat_aux_grid","y"} else 0

                    units = str(coord.attrs.get("units","")).lower()
                    score += 4 if "degree_north" in units or "degrees_north" in units else 0

                    if score>0 and coord.ndim==1 and coord.dims[0] in psi.dims:
                        lat_candidates.append((score,name,coord))

                latitude = None

                if lat_candidates:
                    lat_candidates.sort(key=itemgetter(0),reverse=True)
                    _, lat_name, latitude = lat_candidates[0]
                    lat_method, lat_source = f"original find_latitude: {lat_name}", str(files[0])

                if latitude is None and model=="CNRM-ESM2-1":
                    j_dim = next((d for d in psi.dims if d.lower() in {"j-mean","j_mean","jmean"}),None)

                    if j_dim is None or psi.sizes[j_dim]!=294:
                        raise ValueError(f"CNRM {stage}要求294行j-mean，实际={psi.dims}")

                    row_lat = None

                    # cache只是加速项，不是前置条件。
                    if CNRM_CACHE.is_file():
                        cached = np.asarray(np.load(CNRM_CACHE,allow_pickle=False)["row_latitude"],float).reshape(-1)
                        j40 = int(np.nanargmin(np.abs(cached-40))) if cached.size==294 and np.isfinite(cached).sum()>=20 else -1

                        if j40==207 and abs(float(cached[207])-39.970472)<=0.25:
                            row_lat = cached
                            lat_method = "verified CNRM j294 cache"
                            lat_source = str(CNRM_CACHE)

                    # cache不存在时，从原始同模式海洋网格自动重建。
                    if row_lat is None:
                        grid_files, seen = [], set()

                        for root in CNRM_GRID_ROOTS:
                            if not root.is_dir():
                                continue

                            for alias in ("CNRM-ESM2","CNRM-ESM2-1"):
                                folder = root/alias

                                if folder.is_dir():
                                    for path in sorted(folder.rglob("*.nc")):
                                        key = str(path.resolve()).lower()

                                        if key not in seen:
                                            seen.add(key)
                                            grid_files.append(path)

                        for grid_file in grid_files:
                            try:
                                with xr.open_dataset(grid_file,decode_times=False,cache=False) as grid:
                                    lat2d = next(
                                        (
                                            grid[n].squeeze()
                                            for n in ("nav_lat","latitude","lat","geolat")
                                            if n in grid and grid[n].squeeze().ndim==2
                                        ),
                                        None,
                                    )

                                    lon2d = next(
                                        (
                                            grid[n].squeeze()
                                            for n in ("nav_lon","longitude","lon","geolon")
                                            if n in grid and grid[n].squeeze().ndim==2
                                        ),
                                        None,
                                    )

                                    if lat2d is None or lon2d is None or set(lat2d.dims)!=set(lon2d.dims):
                                        continue

                                    distance = [abs(int(lat2d.sizes[d])-294) for d in lat2d.dims]
                                    grid_j = lat2d.dims[int(np.argmin(distance))]
                                    source_n = int(lat2d.sizes[grid_j])

                                    if abs(source_n-294)>6:
                                        continue

                                    grid_i = next(d for d in lat2d.dims if d!=grid_j)

                                    latv = np.asarray(lat2d.transpose(grid_j,grid_i).values,float)
                                    lonv = np.asarray(lon2d.transpose(grid_j,grid_i).values,float)

                                    lon180 = ((lonv+180)%360)-180
                                    mask = (
                                        np.isfinite(latv)
                                        & np.isfinite(lon180)
                                        & (lon180>=-100)
                                        & (lon180<=20)
                                        & (latv>=-20)
                                        & (latv<=80)
                                    )

                                    source = np.full(source_n,np.nan)

                                    for row in range(source_n):
                                        values = latv[row,mask[row]]

                                        if values.size>=3:
                                            source[row] = float(np.nanmedian(values))

                                    finite = np.isfinite(source)

                                    if finite.sum()<20:
                                        continue

                                    if source_n==294:
                                        trial = source

                                    else:
                                        sx = np.linspace(0,1,source_n)
                                        tx = np.linspace(0,1,294)
                                        trial = np.interp(tx,sx[finite],source[finite],left=np.nan,right=np.nan)

                                    j40 = int(np.nanargmin(np.abs(trial-40)))
                                    actual40 = float(trial[j40])

                                    if j40==207 and abs(actual40-39.970472)<=0.25:
                                        row_lat = trial
                                        lat_method = f"Atlantic row median {source_n}->294"
                                        lat_source = str(grid_file)

                                        np.savez_compressed(CNRM_CACHE,row_latitude=row_lat)

                                        print(f"  CNRM纬度映射已重建：{grid_file}")
                                        print(f"  40°N QA：j={j40}, actual={actual40:.6f}°N")
                                        break

                            except (OSError,ValueError,KeyError):
                                continue

                    if row_lat is None:
                        raise ValueError(
                            f"CNRM j-mean映射失败：cache不存在且网格重建未通过40°N QA；"
                            f"候选文件={len(grid_files)}"
                        )

                    jcoord = np.asarray(psi[j_dim].values) if j_dim in psi.coords else np.arange(294)
                    latitude = xr.DataArray(row_lat,dims=(j_dim,),coords={j_dim:jcoord})

                if latitude is None:
                    raise ValueError(f"{model} {stage}无法按原代码规则识别一维纬度")

                lat_values = np.asarray(latitude.values,float).reshape(-1)
                lat_index = int(np.nanargmin(np.abs(lat_values-TARGET_LAT)))
                actual_lat = float(lat_values[lat_index])

                if abs(actual_lat-TARGET_LAT)>2.5:
                    raise ValueError(f"{model} {stage}最近纬度={actual_lat:.3f}°N，距{TARGET_LAT:g}°N过远")

                lat_dim = latitude.dims[0]
                psi = psi.isel({lat_dim:lat_index})

                print(f"  {variable_name} | basin={basin_text} | lat={lat_method} | actual={actual_lat:.3f}°N")

                # 深度：恢复原CF metadata识别；centimeters正确换算。
                depth_candidates = []
                variables = dict(ds.variables)
                variables.update({name:psi.coords[name] for name in psi.coords})

                for name, coord in variables.items():
                    coord = coord.squeeze()

                    if coord.ndim!=1 or coord.dims[0] not in psi.dims:
                        continue

                    score = 10 if str(coord.attrs.get("standard_name","")).lower() in {"depth","ocean_depth","depth_below_geoid"} else 0
                    score += 6 if str(coord.attrs.get("axis","")).upper()=="Z" else 0
                    score += 2 if str(coord.attrs.get("positive","")).lower() in {"up","down"} else 0
                    score += 5 if name.lower() in {"lev","olevel","depth","depth_coord","moc_z","z_t","z"} else 0

                    if score>0:
                        depth_candidates.append((score,name,coord))

                if not depth_candidates:
                    raise ValueError(f"{model} {stage}无法识别深度")

                depth_candidates.sort(key=itemgetter(0),reverse=True)
                _, depth_name, depth = depth_candidates[0]

                depth_dim = depth.dims[0]
                depth_units = str(depth.attrs.get("units","")).lower().replace(" ","")
                values = np.abs(np.asarray(depth.values,float))

                if depth_units in {"m","meter","meters","metre","metres"}:
                    depth_m = values

                elif depth_units in {"cm","centimeter","centimeters","centimetre","centimetres"}:
                    depth_m = values/100

                elif depth_units in {"km","kilometer","kilometers","kilometre","kilometres"}:
                    depth_m = values*1000

                else:
                    raise ValueError(f"{model} {stage}深度单位无法确认：{depth.attrs.get('units','')!r}")

                depth_index = np.flatnonzero(np.isfinite(depth_m)&(depth_m>=MIN_DEPTH))

                if not depth_index.size:
                    raise ValueError(f"{model} {stage}在500m以下无有效层")

                psi = psi.isel({depth_dim:depth_index})

                extra = [d for d in psi.dims if d not in {time_dim,depth_dim}]

                if extra:
                    raise ValueError(f"{model} {stage}固定盆地/纬度后仍有维度：{extra}")

                # 完整年份；先年度平均剖面再取深度极值。
                years = np.asarray(psi[time_dim].dt.year.values,int)
                unique_years, counts = np.unique(years,return_counts=True)

                steps = int(pd.Series(counts).mode().iloc[0])
                complete = unique_years[counts==steps]
                required = 20 if stage=="piControl" else RAMP_YEARS

                if steps not in {1,12} or complete.size<required:
                    raise ValueError(
                        f"{model} {stage}时间覆盖异常："
                        f"steps/year={steps}, complete={complete.size}"
                    )

                if stage!="piControl" and np.any(np.diff(complete[:required])!=1):
                    raise ValueError(f"{model} {stage}前{required}个完整年份不连续")

                psi = psi.isel({time_dim:np.flatnonzero(np.isin(years,complete))})

                if steps==12:
                    days = psi[time_dim].dt.days_in_month.astype("float64")
                    annual = (
                        (psi*days).groupby(f"{time_dim}.year").sum(time_dim,skipna=True)
                        / days.where(psi.notnull()).groupby(f"{time_dim}.year").sum(time_dim,skipna=True)
                    )

                else:
                    annual = psi

                positive = np.asarray(annual.max(depth_dim,skipna=True).compute().values,float).reshape(-1)
                negative = np.asarray(-annual.min(depth_dim,skipna=True).compute().values,float).reshape(-1)

                units = str(ds[variable_name].attrs.get("units","")).strip()
                key = units.lower().replace(" ","").replace("**","^")

                if "sv" in key or "sverdrup" in key:
                    factor, conversion = 1.0, "already Sv"

                elif "kg" in key and ("s-1" in key or "/s" in key):
                    factor, conversion = 1e-9, "kg s-1 / 1e9"

                elif ("m3" in key or "m^3" in key) and ("s-1" in key or "/s" in key):
                    factor, conversion = 1e-6, "m3 s-1 / 1e6"

                else:
                    raise ValueError(f"{model} {stage}AMOC单位无法转换：{units!r}")

                positive, negative = positive*factor, negative*factor

                if model=="CNRM-ESM2-1" and stage=="piControl":
                    max_mean = float(np.nanmean(positive))
                    neg_mean = float(np.nanmean(negative))

                    cnrm_orientation = (
                        "negative_min"
                        if (not np.isfinite(max_mean) or max_mean<1)
                        and np.isfinite(neg_mean)
                        and 2<=neg_mean<=80
                        else "max"
                    )

                    print(
                        f"  CNRM sign QA: max={max_mean:.3f}, "
                        f"-min={neg_mean:.3f} Sv -> {cnrm_orientation}"
                    )

                orientation = cnrm_orientation if model=="CNRM-ESM2-1" else "max"
                series[stage] = negative if orientation=="negative_min" else positive

                time = ds[time_dim]
                time_units = str(time.encoding.get("units","") or time.attrs.get("units",""))
                calendar = str(time.encoding.get("calendar","") or time.attrs.get("calendar","") or "standard")
                first_time = np.asarray(time.values).reshape(-1)[0]

                time_start = (
                    float(cftime.date2num(first_time,units=time_units,calendar=calendar))
                    if time_units
                    else np.nan
                )

                stage_meta[stage] = {
                    "variant":variant,
                    "parent_exp":str(ds.attrs.get("parent_experiment_id","")).strip(),
                    "parent_variant":str(ds.attrs.get("parent_variant_label","")).strip().lower(),
                    "branch":ds.attrs.get("branch_time_in_parent",np.nan),
                    "time_start":time_start,
                    "time_units":time_units,
                    "calendar":calendar,
                    "variable":variable_name,
                    "basin":basin_text,
                    "basin_index":basin_index,
                    "component":component_text,
                    "lat":actual_lat,
                    "lat_index":lat_index,
                    "lat_dim":lat_dim,
                    "lat_method":lat_method,
                    "lat_source":lat_source,
                    "depth":depth_name,
                    "depth_units":depth_units,
                    "conversion":conversion,
                    "orientation":orientation,
                }

        # 三阶段QA与piControl线性订正。
        variants = [stage_meta[stage]["variant"] for stage,_ in STAGES]

        if variants != [member.lower()]*3:
            raise ValueError(f"{model}三阶段member不一致：{variants}")

        if model=="NorESM2-LM" and variants != ["r1i1p1f1"]*3:
            raise ValueError(f"NorESM2-LM必须p1/p1/p1：{variants}")

        cdr_parent_exp = stage_meta["1pctCO2-cdr"]["parent_exp"]
        cdr_parent = stage_meta["1pctCO2-cdr"]["parent_variant"]

        if cdr_parent_exp and re.sub(r"[^A-Z0-9]","",cdr_parent_exp.upper())!="1PCTCO2":
            raise ValueError(f"{model} CDR parent_experiment={cdr_parent_exp}")

        if cdr_parent and cdr_parent!=member.lower():
            raise ValueError(f"{model} CDR parent_variant={cdr_parent}，选择={member}")

        lats = np.asarray([stage_meta[stage]["lat"] for stage,_ in STAGES],float)

        if np.ptp(lats)>1:
            raise ValueError(f"{model}三阶段实际纬度不一致：{lats.tolist()}")

        pi = np.asarray(series["piControl"],float)
        up = np.asarray(series["1pctCO2"],float)
        cdr = np.asarray(series["1pctCO2-cdr"],float)

        sign_flipped = float(np.nanmean(pi))<0

        if sign_flipped:
            pi, up, cdr = -pi, -up, -cdr

        total = min(MAX_YEARS,RAMP_YEARS+cdr.size)
        combined = np.concatenate([up[:RAMP_YEARS],cdr[:total-RAMP_YEARS]])

        branch_index = 0

        try:
            delta = float(stage_meta["1pctCO2"]["branch"])-float(stage_meta["piControl"]["time_start"])
            pi_units = stage_meta["piControl"]["time_units"].lower()
            pi_calendar = stage_meta["piControl"]["calendar"].lower()

            if np.isfinite(delta) and "day" in pi_units:
                divisor = (
                    360
                    if "360" in pi_calendar
                    else 365.2425
                    if pi_calendar in {"gregorian","standard","proleptic_gregorian"}
                    else 365
                )
                branch_index = int(np.floor(delta/divisor+1e-6))

            elif np.isfinite(delta) and "month" in pi_units:
                branch_index = int(np.floor(delta/12+1e-6))

            elif np.isfinite(delta) and "year" in pi_units:
                branch_index = int(np.floor(delta+1e-6))

        except (TypeError,ValueError):
            branch_index = 0

        branch_index = int(np.clip(branch_index,0,max(0,pi.size-1)))

        pi_segment = pi[branch_index:]
        valid = np.isfinite(pi_segment)

        if valid.sum()<20:
            branch_index = 0
            pi_segment = pi
            valid = np.isfinite(pi)

        if valid.sum()<20:
            raise ValueError(f"{model} piControl有效年数少于20年")

        t = np.arange(pi_segment.size,dtype=float)
        slope, intercept = np.polyfit(t[valid],pi_segment[valid],1)

        pi_mean = float(np.nanmean(pi_segment))
        reference = intercept+slope*np.arange(total,dtype=float)
        anomaly = combined-reference
        mean_anomaly = combined-pi_mean

        smooth = (
            pd.Series(anomaly)
            .rolling(SMOOTH,center=True,min_periods=SMOOTH//2)
            .mean()
            .to_numpy()
        )

        initial = float(np.nanmean(anomaly[:10]))
        jump = float(np.nanmean(cdr[:5])-np.nanmean(up[RAMP_YEARS-5:RAMP_YEARS]))

        fail = []

        if abs(initial)>3:
            fail.append(f"initial10={initial:+.2f} Sv >3")

        if abs(jump)>4:
            fail.append(f"UP→CDR jump={jump:+.2f} Sv >4")

        if not 2<=pi_mean<=50:
            fail.append(f"piControl={pi_mean:.2f} Sv不在2–50")

        if fail:
            raise ValueError("；".join(fail))

        results[model] = {
            "combined":combined,
            "reference":reference,
            "anomaly":anomaly,
            "mean_anomaly":mean_anomaly,
            "smooth":smooth,
            "length":total,
        }

        m = stage_meta["piControl"]

        metadata_rows.append({
            "Model":model,
            "Selected member":member,
            "Variable":m["variable"],
            "Basin selection":m["basin"],
            "Basin index":m["basin_index"],
            "MOC component treatment":m["component"],
            "Actual latitude (degN)":m["lat"],
            "Latitude index":m["lat_index"],
            "Latitude dimension":m["lat_dim"],
            "Latitude mapping method":m["lat_method"],
            "Latitude mapping source":m["lat_source"],
            "Depth coordinate":m["depth"],
            "Depth units":m["depth_units"],
            "AMOC unit conversion":m["conversion"],
            "Streamfunction orientation":m["orientation"],
            "Sign flipped after extraction":sign_flipped,
            "piControl climatology (Sv)":pi_mean,
            "piControl trend (Sv yr-1)":float(slope),
            "Initial 10-yr anomaly (Sv)":initial,
            "1pctCO2-to-CDR jump (Sv)":jump,
        })

        print(
            f"[OK] pi={pi_mean:.2f} Sv | "
            f"trend={slope*100:+.3f} Sv/century | "
            f"initial10={initial:+.2f} | jump={jump:+.2f}"
        )

    except Exception as error:
        failure_rows.append({
            "Model":model,
            "Error":f"{type(error).__name__}: {error}",
        })

        print(f"[FAILED] {model}: {type(error).__name__}: {error}")


# %% 5. 九模式QA与保存
pd.DataFrame(failure_rows,columns=["Model","Error"]).to_csv(
    FAILED,
    index=False,
    encoding="utf-8-sig",
)

success = [model for model in MODELS if model in results]

if success != MODELS:
    raise RuntimeError(f"9模式未全部通过。成功={success}；失败报告={FAILED}")

metadata = pd.DataFrame(metadata_rows)

nor_member = str(
    metadata.loc[
        metadata["Model"]=="NorESM2-LM",
        "Selected member",
    ].iloc[0]
).lower()

if nor_member != "r1i1p1f1":
    raise RuntimeError(f"NorESM2-LM最终member QA失败：{nor_member}")

common = min(results[model]["length"] for model in MODELS)

if common<280:
    raise RuntimeError(f"9模式共同长度={common}年，少于280年")

sim_year = np.arange(1,common+1)
calendar_year = START_YEAR+sim_year-1
index = np.arange(common)

co2_exp = np.where(
    index<RAMP_YEARS,
    index+1,
    np.maximum(2*RAMP_YEARS-1-index,0),
)

co2_ppm = CO2_0*1.01**co2_exp

frames = []

plot = pd.DataFrame({
    "simulation_year":sim_year,
    "calendar_year":calendar_year,
    "co2_ppm":co2_ppm,
})

for model in MODELS:
    r = results[model]

    plot[model] = r["smooth"][:common]

    frames.append(pd.DataFrame({
        "Model":model,
        "Simulation year":sim_year,
        "Calendar year":calendar_year,
        "Atmospheric CO2 (ppm)":co2_ppm,
        "Forced AMOC strength (Sv)":r["combined"][:common],
        "piControl linear reference (Sv)":r["reference"][:common],
        "AMOC anomaly (Sv)":r["anomaly"][:common],
        "AMOC smoothed anomaly (Sv)":r["smooth"][:common],
        "AMOC anomaly using piControl mean (Sv)":r["mean_anomaly"][:common],
    }))

stack = plot[MODELS].to_numpy(float)

if not np.all(np.isfinite(stack)):
    bad = np.argwhere(~np.isfinite(stack))[0]
    raise RuntimeError(
        f"绘图期非有限值：year={sim_year[bad[0]]}, "
        f"model={MODELS[bad[1]]}"
    )

plot["mmm"] = np.mean(stack,axis=1)
plot["std"] = np.std(stack,axis=1,ddof=0)
plot["lower"] = plot["mmm"]-plot["std"]
plot["upper"] = plot["mmm"]+plot["std"]

pd.concat(frames,ignore_index=True).to_csv(
    ANNUAL,
    index=False,
    encoding="utf-8-sig",
    float_format="%.8f",
)

metadata.to_csv(
    META,
    index=False,
    encoding="utf-8-sig",
    float_format="%.8f",
)

plot.to_csv(
    PLOT_CSV,
    index=False,
    encoding="utf-8-sig",
    float_format="%.8f",
)

PLOT_JSON.write_text(
    json.dumps({
        "models":MODELS,
        "model_count":9,
        "nor_variants":{
            "piControl":"r1i1p1f1",
            "1pctCO2":"r1i1p1f1",
            "1pctCO2-cdr":"r1i1p1f1",
        },
        "amoc_definition":f"Atlantic-Arctic, fixed{TARGET_LAT:g}°N nearest latitude, maximum below 500 m",
        "piControl_correction":"branch-aligned linear",
        "smoothing":"11-year centered",
        "spread":"±1 standard deviation across 9 models",
        "common_length_years":int(common),
    },ensure_ascii=False,indent=2),
    encoding="utf-8",
)

print(f"\n[完成] {ANNUAL}")
print(f"[完成] {META}")
print(f"[完成] {PLOT_CSV}")
print("9/9模式通过；CNRM cache不存在时自动重建，不再把cache作为前置条件。")