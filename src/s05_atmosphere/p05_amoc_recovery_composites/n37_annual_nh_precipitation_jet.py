# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import xarray as xr
from scipy.spatial import cKDTree
from cartopy.io import shapereader
from shapely import contains_xy
from shapely.geometry import mapping, shape


# %% 2. 参数与路径设置

root = Path(r"F:\sea\Antic sea ice constraint")
destination = root / "output/s05_atmosphere/p05_amoc_recovery_composites/n37_annual_nh_precipitation_jet"
assert not destination.exists() or not list(destination.iterdir()), "正式计算结果已存在，须先归档。"
upper_directory = root / "output/s02_sia/s02_arctic_region/p22_grouped_temperature_ice/n066_grouped_tas_ice_patterns"
membership_path = upper_directory / "membership.csv"
membership = pd.read_csv(membership_path)
models = membership.model.to_numpy(dtype=str)
groups = ["better_nonrecovery", "better_recovery", "larger_bias"]
members = {"CESM2": "r1i1p1f1", "CNRM-ESM2-1": "r1i1p1f2", "NorESM2-LM": "r1i1p1f1",
           "GFDL-ESM4": "r1i1p1f1", "ACCESS-ESM1-5": "r1i1p1f1", "CanESM5": "r1i1p2f1",
           "MIROC-ES2L": "r1i1p1f2", "UKESM1-0-LL": "r1i1p1f2"}
calendar_starts = {"CESM2": [1, 1], "CNRM-ESM2-1": [1850, 1990], "NorESM2-LM": [1, 141],
                   "GFDL-ESM4": [1, 141], "ACCESS-ESM1-5": [101, 241], "CanESM5": [1850, 1991],
                   "MIROC-ES2L": [1850, 1989], "UKESM1-0-LL": [1850, 1990]}
periods = ["initial", "later"]
experiments = ["1pctCO2", "1pctCO2-cdr"]
target_pressure = 25000.0
wind_contour_levels = [20.0, 30.0]
pr_color_limits = [-1.0, 1.0]
monthly_file = root / "output/s05_atmosphere/p05_amoc_recovery_composites/n05_nor_pr_monthly_from_daily/nor_pr_cdr_027101-028012_monthly_derived.nc"
can_grid = Path(r"F:\sea\cmip6prad\CANSESM5\pr_Amon_CanESM5_1pctCO2_r1i1p2f1_gn_185001-200012.nc")
with xr.open_dataset(can_grid, decode_times=False) as source:
    assert source.pr.dims == ("time", "lat", "lon")
    latitude = source.lat.values[source.lat.values > 0].astype(float)
    longitude = source.lon.values.astype(float)
    reference_latitude_global = source.lat.values.astype(float)
assert latitude.size == 32 and longitude.size == 128
longitude_grid, latitude_grid = np.meshgrid(longitude, latitude)
source_rows, diagnostics, yearly_fields = [], [], []


# %% 3. 按真实历年读取八模式的全年端点样本

for model in models:
    variables = {}
    for variable in ["pr", "ua"]:
        phase_fields = []
        for period_index, period in enumerate(periods):
            first_year = calendar_starts[model][period_index] + (130 if period_index == 1 else 0)
            last_year = first_year + 9
            suffix = "ad" if period_index == 0 else ""
            folder = Path("F:/sea") / f"cmip6{variable}{suffix}"
            table = "AERmon" if model == "NorESM2-LM" and variable == "ua" else "Amon"
            if model == "NorESM2-LM" and variable == "pr" and period_index == 1:
                paths = [monthly_file]
            else:
                candidates = sorted(folder.rglob(f"{variable}_{table}_{model}_{experiments[period_index]}_{members[model]}_*.nc"))
                paths = []
                for path in candidates:
                    date_tokens = path.stem.rsplit("_", 1)[1].split("-")
                    if int(date_tokens[0][:4]) <= last_year and int(date_tokens[1][:4]) >= first_year:
                        paths.append(path)
            assert paths, (model, variable, period)
            pieces = []
            for path in paths:
                with xr.open_dataset(path, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as source:
                    assert source.attrs["source_id"] == model
                    assert source.attrs["experiment_id"] == experiments[period_index]
                    assert source.attrs["variant_label"] == members[model]
                    selected = source.time.dt.year.isin(range(first_year, last_year + 1))
                    assert bool(selected.any())
                    if variable == "pr":
                        assert source.pr.dims == ("time", "lat", "lon")
                        assert source.pr.attrs["units"] == "kg m-2 s-1"
                        field = source.pr.isel(time=selected).sortby("lat").sel(lat=slice(-5, 90)).load().astype(float) * 86400.0
                    elif model != "NorESM2-LM":
                        assert source.ua.dims == ("time", "plev", "lat", "lon")
                        assert source.plev.attrs["units"] == "Pa" and source.ua.attrs["units"] == "m s-1"
                        # ACCESS气压坐标含浮点尾差，按已查明的Pa层值取整后明确取250hPa。
                        field = source.ua.assign_coords(plev=np.rint(source.plev)).sel(plev=target_pressure, drop=True)
                        field = field.isel(time=selected).sortby("lat").sel(lat=slice(-5, 90)).load().astype(float)
                    else:
                        assert source.ua.dims == ("time", "lev", "lat", "lon")
                        assert source.ps.dims == ("time", "lat", "lon")
                        assert source.ua.attrs["units"] == "m s-1" and source.ps.attrs["units"] == "Pa"
                        native = source.ua.isel(time=selected).sortby("lat").sel(lat=slice(-5, 90)).load().astype(float)
                        surface_pressure = source.ps.isel(time=selected).sel(lat=native.lat).load().values
                        # Nor原始层由近地面向高空排列，显式反转为压力递增顺序。
                        offset = (source.a.values * float(source.p0))[::-1]
                        fraction = source.b.values[::-1]
                        wind = native.values[:, ::-1]
                        converted = np.full((native.sizes["time"], native.sizes["lat"], native.sizes["lon"]), np.nan)
                        for month_index in range(native.sizes["time"]):
                            pressure = offset[:, None, None] + fraction[:, None, None] * surface_pressure[month_index][None]
                            assert np.all(np.diff(pressure, axis=0) > 0)
                            lower = np.clip(np.sum(pressure <= target_pressure, axis=0) - 1, 0, len(offset) - 2)
                            low_p = np.take_along_axis(pressure, lower[None], axis=0)[0]
                            high_p = np.take_along_axis(pressure, (lower + 1)[None], axis=0)[0]
                            low_u = np.take_along_axis(wind[month_index], lower[None], axis=0)[0]
                            high_u = np.take_along_axis(wind[month_index], (lower + 1)[None], axis=0)[0]
                            amount = (np.log(target_pressure) - np.log(low_p)) / (np.log(high_p) - np.log(low_p))
                            valid = ((target_pressure >= pressure[0]) & (target_pressure <= pressure[-1])
                                     & (target_pressure <= surface_pressure[month_index]))
                            converted[month_index] = np.where(valid, low_u + amount * (high_u - low_u), np.nan)
                        field = xr.DataArray(converted, dims=("time", "lat", "lon"),
                                             coords={"time": native.time, "lat": native.lat, "lon": native.lon})
                    pieces.append(field.reset_coords(drop=True))
                    source_rows.append({"model": model, "variable": variable, "period": period,
                                        "file": str(path), "raw_dimensions": str(source[variable].dims),
                                        "raw_shape": str(source[variable].shape), "units": source[variable].attrs["units"],
                                        "first_selected_month": str(field.time.values[0]),
                                        "last_selected_month": str(field.time.values[-1]),
                                        "months": int(selected.sum()), "calendar": source.time.encoding["calendar"],
                                        "tracking_id": source.attrs.get("tracking_id", "local_daily_derived_monthly"),
                                        "file_bytes": path.stat().st_size})
            monthly = xr.concat(pieces, "time", join="exact").sortby("time")
            month_ids = monthly.time.dt.year.values * 12 + monthly.time.dt.month.values
            expected_ids = np.arange(first_year * 12 + 1, (last_year + 1) * 12 + 1)
            np.testing.assert_array_equal(month_ids, expected_ids)
            assert monthly.dims == ("time", "lat", "lon") and monthly.sizes["time"] == 120
            # 与上排原有年平均一致：12个月等权，10个年份等权；不混用旧DJF日数权重。
            annual = monthly.coarsen(time=12, boundary="exact").mean(skipna=False)
            annual = annual.rename(time="year_in_window").assign_coords(year_in_window=np.arange(1, 11))
            annual = annual.assign_coords(lon=annual.lon % 360).sortby("lon")
            if model == "MIROC-ES2L":
                # 两套同名64点高斯纬度仅差约1e-13度；明确对齐，防止边界被误判为外推。
                same_gaussian_latitude = reference_latitude_global[reference_latitude_global >= -5]
                np.testing.assert_allclose(annual.lat.values, same_gaussian_latitude, atol=1e-10, rtol=0)
                annual = annual.assign_coords(lat=same_gaussian_latitude)
            assert np.all(np.diff(annual.lat.values) > 0) and np.all(np.diff(annual.lon.values) > 0)
            cyclic = xr.concat([annual.isel(lon=[-1]).assign_coords(lon=[float(annual.lon[-1]) - 360]),
                                annual, annual.isel(lon=[0]).assign_coords(lon=[float(annual.lon[0]) + 360])], "lon")
            mapped = cyclic.interp(lat=latitude, lon=longitude)
            assert mapped.dims == ("year_in_window", "lat", "lon")
            phase_fields.append(mapped.expand_dims(period=[period]))
            diagnostics.append({"model": model, "variable": variable, "period": period,
                                "calendar_first_year": first_year, "calendar_last_year": last_year,
                                "months": 120, "native_monthly_selected_sha256": hashlib.sha256(np.ascontiguousarray(monthly.values).tobytes()).hexdigest(),
                                "mapped_finite_cells_all_years": int(np.isfinite(mapped.values).all(axis=0).sum()),
                                "minimum": float(np.nanmin(mapped)), "maximum": float(np.nanmax(mapped))})
            print(model, variable, period, f"{first_year}-{last_year};120months;complete", flush=True)
        variable_name = "pr" if variable == "pr" else "u250"
        values = xr.concat(phase_fields, "period", join="exact")
        support = values.notnull().all(["period", "year_in_window"])
        variables[variable_name] = values.where(support)
    yearly_fields.append(xr.Dataset(variables).expand_dims(model=[model]))


# %% 4. 核查并修正旧UKESM气温缓存的历年错位

with np.load(upper_directory / "grouped_tas_ice_fields.npz") as saved:
    upper_arrays = {key: saved[key] for key in saved.files}
uk_index = int(np.where(upper_arrays["model"] == "UKESM1-0-LL")[0][0])
target_lon, target_lat = np.meshgrid(upper_arrays["longitude"], upper_arrays["latitude"])
lat_rad, lon_rad = np.deg2rad(target_lat.ravel()), np.deg2rad(target_lon.ravel())
target_xyz = np.column_stack([np.cos(lat_rad) * np.cos(lon_rad), np.cos(lat_rad) * np.sin(lon_rad), np.sin(lat_rad)])
tas_windows = ["initial1850_1859", "correct_later2120_2129", "old_cache_later2180_2189"]
tas_first_years = [1850, 2120, 2180]
uk_tas_endpoints = []
uk_tas_audit = []
for window_index, (window, first_year) in enumerate(zip(tas_windows, tas_first_years)):
    folder = Path("F:/sea/cmip6tasad" if window_index == 0 else "F:/sea/cmip6tas")
    experiment = "1pctCO2" if window_index == 0 else "1pctCO2-cdr"
    paths = sorted(folder.rglob(f"tas_Amon_UKESM1-0-LL_{experiment}_r1i1p1f2_gn_*.nc"))
    parts = []
    for path in paths:
        date_tokens = path.stem.rsplit("_", 1)[1].split("-")
        if int(date_tokens[0][:4]) > first_year + 9 or int(date_tokens[1][:4]) < first_year:
            continue
        with xr.open_dataset(path, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as source:
            assert source.tas.dims == ("time", "lat", "lon") and source.tas.attrs["units"] == "K"
            assert source.attrs["source_id"] == "UKESM1-0-LL" and source.attrs["variant_label"] == "r1i1p1f2"
            selected = source.time.dt.year.isin(range(first_year, first_year + 10))
            part = source.tas.isel(time=selected).load().astype(float) - 273.15
            parts.append(part)
            uk_tas_audit.append({"window": window, "file": str(path), "raw_dimensions": str(source.tas.dims),
                                 "first_month": str(part.time.values[0]), "last_month": str(part.time.values[-1]),
                                 "months": int(selected.sum()), "tracking_id": source.attrs.get("tracking_id", "")})
    monthly = xr.concat(parts, "time", join="exact").sortby("time")
    np.testing.assert_array_equal(monthly.time.dt.year.values * 12 + monthly.time.dt.month.values,
                                  np.arange(first_year * 12 + 1, (first_year + 10) * 12 + 1))
    native_mean = monthly.mean("time", skipna=False)
    native_lon, native_lat = np.meshgrid(native_mean.lon.values, native_mean.lat.values)
    native_lon = (native_lon + 180) % 360 - 180
    selected = (native_lat >= 45) & (native_lat <= 90) & np.isfinite(native_mean.values)
    source_lat_rad, source_lon_rad = np.deg2rad(native_lat[selected]), np.deg2rad(native_lon[selected])
    xyz = np.column_stack([np.cos(source_lat_rad) * np.cos(source_lon_rad),
                           np.cos(source_lat_rad) * np.sin(source_lon_rad), np.sin(source_lat_rad)])
    chord, indices = cKDTree(xyz).query(target_xyz, k=4)
    distance = 6371.0 * 2 * np.arcsin(np.clip(chord / 2, 0, 1))
    neighbors = native_mean.values[selected][indices]
    weights = 1 / np.maximum(distance, 1.0) ** 2
    mapped = np.sum(weights * neighbors, axis=1) / np.sum(weights, axis=1)
    exact = distance[:, 0] < 1e-6
    mapped[exact] = neighbors[exact, 0]
    mapped[distance[:, 0] > 450] = np.nan
    uk_tas_endpoints.append(mapped.reshape(45, 360))
    print("UKESM tas verified", window, flush=True)
with np.load(root / "output/s03tem_n01/arctic_tas_UKESM1-0-LL_pattern_1deg.npz") as old:
    initial_error = float(np.nanmax(np.abs(uk_tas_endpoints[0] - old["co2_initial"])))
    wrong_later_error = float(np.nanmax(np.abs(uk_tas_endpoints[2] - old["cdr_final"])))
assert initial_error < 1e-8 and wrong_later_error < 1e-8
correct_uk_change = uk_tas_endpoints[1] - uk_tas_endpoints[0]
previous_uk_change = upper_arrays["member_tas_change_degC"][uk_index].copy()
upper_arrays["member_tas_change_degC"][uk_index] = correct_uk_change
for group_index, group in enumerate(groups):
    selected = membership.group.eq(group).to_numpy()
    upper_arrays["group_tas_change_degC"][group_index] = np.mean(upper_arrays["member_tas_change_degC"][selected], axis=0)
tas_correction = {"model": "UKESM1-0-LL", "reason": "Old cached TAS used index130:140 after missing initial60CDRyears;actual2180-2189 instead of2120-2129",
                  "initial_reproduction_max_error_degC": initial_error,
                  "old_wrong_later_reproduction_max_error_degC": wrong_later_error,
                  "max_abs_change_from_cached_UK_TAS_degC": float(np.nanmax(np.abs(correct_uk_change - previous_uk_change))),
                  "max_abs_effect_on_four_member_group_degC": float(np.nanmax(np.abs(correct_uk_change - previous_uk_change))) / 4,
                  "old_n066_products": "Retained unchanged for audit;current six-panel product uses corrected calendar window"}


# %% 5. 用原始120个月核验CanESM5初始海冰并更新端点缓存

can_ice_directory = root / "output/s02sia/s02_arctic_region/p04_pattern_unadjusted/n01_arctic_sia_pattern"
can_ice_path = can_ice_directory / "arctic_sia_pattern_cells.csv"
ice_source = pd.read_csv(can_ice_path)
can_ice = ice_source.loc[ice_source.model.eq("CanESM5")].sort_values(["latitude_deg_north", "longitude_deg_east"])
np.testing.assert_array_equal(can_ice.latitude_deg_north.to_numpy().reshape(45, 360), target_lat)
np.testing.assert_array_equal(can_ice.longitude_deg_east.to_numpy().reshape(45, 360), target_lon)
canonical_initial = can_ice.edge_raw_co2_sic_fraction.to_numpy().reshape(45, 360)
canonical_later = can_ice.edge_raw_cdr_sic_fraction.to_numpy().reshape(45, 360)
can_sic_file = Path(r"F:\sea\cmipseaicead\CANESM5\siconc_SImon_CanESM5_1pctCO2_r1i1p2f1_gn_185001-200012.nc")
with xr.open_dataset(can_sic_file, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as source:
    assert source.siconc.dims == ("time", "j", "i") and source.siconc.attrs["units"] == "%"
    assert source.latitude.dims == source.longitude.dims == ("j", "i")
    selected = source.time.dt.year.isin(range(1850, 1860))
    assert int(selected.sum()) == 120
    native = source.siconc.isel(time=selected).astype(float).mean("time", skipna=True).values / 100
    native_lat = source.latitude.values
    native_lon = (source.longitude.values + 180) % 360 - 180
    can_sic_tracking_id = source.attrs.get("tracking_id", "")
valid = np.isfinite(native) & np.isfinite(native_lat) & np.isfinite(native_lon) & (native_lat >= 45) & (native_lat <= 90)
lat_rad, lon_rad = np.deg2rad(native_lat[valid]), np.deg2rad(native_lon[valid])
xyz = np.column_stack([np.cos(lat_rad) * np.cos(lon_rad), np.cos(lat_rad) * np.sin(lon_rad), np.sin(lat_rad)])
chord, indices = cKDTree(xyz).query(target_xyz, k=8)
distance = 6371.0088 * 2 * np.arcsin(np.clip(chord / 2, 0, 1))
weights = 1 / np.maximum(distance, 1e-6) ** 2
checked_initial = np.sum(weights * native[valid][indices], axis=1) / np.sum(weights, axis=1)
exact = distance[:, 0] <= 1e-6
checked_initial[exact] = native[valid][indices[exact, 0]]
checked_initial[distance[:, 0] > 300] = np.nan
checked_initial = checked_initial.reshape(45, 360)
can_source_error = float(np.nanmax(np.abs(checked_initial - canonical_initial)))
assert can_source_error < 1e-6
can_index = int(np.where(upper_arrays["model"] == "CanESM5")[0][0])
old_can_initial = upper_arrays["member_initial_sic_fraction"][can_index].copy()
old_can_later = upper_arrays["member_later_sic_fraction"][can_index].copy()
upper_arrays["member_initial_sic_fraction"][can_index] = canonical_initial
upper_arrays["member_later_sic_fraction"][can_index] = canonical_later
for group_index, group in enumerate(groups):
    selected = membership.group.eq(group).to_numpy()
    initial = np.mean(upper_arrays["member_initial_sic_fraction"][selected], axis=0)
    later = np.mean(upper_arrays["member_later_sic_fraction"][selected], axis=0)
    paired = np.isfinite(initial) & np.isfinite(later)
    upper_arrays["group_initial_sic_fraction"][group_index] = np.where(paired, initial, np.nan)
    upper_arrays["group_later_sic_fraction"][group_index] = np.where(paired, later, np.nan)
upper_arrays["initial_sic_for_contour"] = np.where(upper_arrays["display_land_mask"][None], np.nan, upper_arrays["group_initial_sic_fraction"])
upper_arrays["later_sic_for_contour"] = np.where(upper_arrays["display_land_mask"][None], np.nan, upper_arrays["group_later_sic_fraction"])
can_ice_correction = {"model": "CanESM5", "raw_file": str(can_sic_file), "tracking_id": can_sic_tracking_id,
                      "initial_months": "1850-01 through1859-12;120months", "raw_dimensions": "siconc(time,j,i),latitude/longitude(j,i)",
                      "canonical_initial_vs_raw_regridded_max_error_fraction": can_source_error,
                      "old_initial_vs_correct_max_difference_fraction": float(np.nanmax(np.abs(old_can_initial - canonical_initial))),
                      "later_source_difference_max_fraction": float(np.nanmax(np.abs(old_can_later - canonical_later))),
                      "action": "Use verified current raw initial/later SIC endpoint pair for Can only;other7member SIC fields unchanged;recompute fixed4member group and paired support",
                      "legacy_initial_error_cause": "Not localized in obsolete CSV generation;do not assume it is solely regridding"}
print("CanESM5 initial SIC raw120month verification", can_source_error, flush=True)


# %% 6. 等权分组及北半球大陆区域诊断

samples = xr.concat(yearly_fields, "model", join="exact")
samples["pr"].attrs["units"] = "mm day-1"
samples["u250"].attrs["units"] = "m s-1"
climate = samples.mean("year_in_window", skipna=False)
group_climate = xr.concat([climate.sel(model=membership.loc[membership.group.eq(group), "model"].tolist()).mean("model", skipna=False).expand_dims(group=[group])
                           for group in groups], "group")
changes = group_climate.sel(period="later", drop=True) - group_climate.sel(period="initial", drop=True)
coastline_source = Path(r"C:\Users\Lenovo\.local\share\cartopy\shapefiles\natural_earth\physical\ne_50m_coastline.shp")
coastline_reader = shapereader.Reader(coastline_source)
coastline_features = []
for geometry in coastline_reader.geometries():
    if geometry.bounds[3] >= -2:
        coastline_features.append({"type": "Feature", "properties": {}, "geometry": mapping(geometry)})
coastline_reader.close()
land_json = json.loads((upper_directory / "map_land.geojson").read_text(encoding="utf-8"))
land_mask = np.zeros(longitude_grid.shape, dtype=bool)
signed_longitude = (longitude_grid + 180) % 360 - 180
for feature in land_json["features"]:
    land_mask |= contains_xy(shape(feature["geometry"]), signed_longitude, latitude_grid)
lat_bounds = np.concatenate([[0.0], (latitude[:-1] + latitude[1:]) / 2, [90.0]])
area_weights = np.repeat((np.sin(np.deg2rad(lat_bounds[1:])) - np.sin(np.deg2rad(lat_bounds[:-1])))[:, None], 128, axis=1)
land_boxes = {"North_America": [-170, -50, 15, 75], "Europe": [-10, 45, 35, 72],
              "Asia": [45, 180, 5, 75], "North_Africa": [-20, 45, 0, 35]}
region_rows, group_rows = [], []
for group in groups:
    precipitation = changes.pr.sel(group=group).values
    initial_wind = group_climate.u250.sel(group=group, period="initial").values
    later_wind = group_climate.u250.sel(group=group, period="later").values
    group_rows.append({"group": group, "n_models": int(membership.group.eq(group).sum()),
                       "pr_minimum_mm_day": float(np.nanmin(precipitation)), "pr_maximum_mm_day": float(np.nanmax(precipitation)),
                       "pr_below_color_limit_cells": int(np.sum(precipitation < pr_color_limits[0])),
                       "pr_above_color_limit_cells": int(np.sum(precipitation > pr_color_limits[1])),
                       "initial_u250_maximum_m_s": float(np.nanmax(initial_wind)), "later_u250_maximum_m_s": float(np.nanmax(later_wind))})
    for region, bounds in land_boxes.items():
        valid = (land_mask & (signed_longitude >= bounds[0]) & (signed_longitude < bounds[1])
                 & (latitude_grid >= bounds[2]) & (latitude_grid < bounds[3]) & np.isfinite(precipitation))
        region_rows.append({"group": group, "region": region, "land_cells": int(valid.sum()),
                            "area_weighted_pr_change_mm_day": float(np.average(precipitation[valid], weights=area_weights[valid]))})


# %% 7. 保存计算产品、原始维序与图外说明所需口径

destination.mkdir(parents=True, exist_ok=True)
samples.to_netcdf(destination / "annual_member_samples.nc")
climate.to_netcdf(destination / "member_climatology.nc")
group_climate.to_netcdf(destination / "group_climatology.nc")
changes.to_netcdf(destination / "group_changes.nc")
np.savez_compressed(destination / "upper_tas_ice_fields.npz", **upper_arrays)
np.savez_compressed(destination / "UKESM_tas_window_audit.npz", window=np.array(tas_windows),
                    endpoint_degC=np.stack(uk_tas_endpoints), previous_change_degC=previous_uk_change,
                    corrected_change_degC=correct_uk_change)
np.savez_compressed(destination / "CanESM5_sic_source_audit.npz", previous_initial_fraction=old_can_initial,
                    previous_later_fraction=old_can_later, correct_initial_fraction=canonical_initial,
                    correct_later_fraction=canonical_later, independent_raw_initial_fraction=checked_initial)
np.savez_compressed(destination / "nh_land_support.npz", longitude=longitude, latitude=latitude,
                    land_mask=land_mask, cell_area_weights=area_weights)
pd.DataFrame(source_rows).to_csv(destination / "input_audit.csv", index=False)
pd.DataFrame(uk_tas_audit).to_csv(destination / "UKESM_tas_input_audit.csv", index=False)
pd.DataFrame(diagnostics).to_csv(destination / "endpoint_diagnostics.csv", index=False)
pd.DataFrame(group_rows).to_csv(destination / "group_diagnostics.csv", index=False)
pd.DataFrame(region_rows).to_csv(destination / "land_region_precipitation.csv", index=False)
membership.to_csv(destination / "membership.csv", index=False)
(destination / "coastline_nh_50m.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": coastline_features}), encoding="utf-8")
metadata = {"models": models.tolist(), "groups": groups,
            "season": "Annual;12calendar months equally weighted,10years equally weighted;matches retained top-row convention,not oldDJFdayweighting",
            "period": "initial1pctCO2years1-10;laterCDRyears131-140;nominal concatenatedyears271-280",
            "calendar_starts": calendar_starts, "wind_pressure_Pa": target_pressure,
            "wind_contour_levels_m_s": wind_contour_levels,
            "wind_definition": "250hPa zonal wind u,positive eastward;20and30m/s strong-westerly band contours,not total wind speed or a uniquely tracked jet axis",
            "geographic_extent": {"upper": [-180, 180, 45, 90], "lower": [-180, 180, 0, 90]},
            "precipitation": "pr*86400 to mm/day;land and ocean retained;not P-E or moisture convergence",
            "precipitation_color_limits_mm_day": pr_color_limits,
            "groups_method": "fixed2/2/4members;model endpoints first,then equal-model means;all fixed members required finite",
            "horizontal_mapping": "Native yearly pr/u250 bilinear periodic interpolation to CanESM564Gaussian-grid northern32latitudes and128longitudes;no extrapolation;not an area-conservative water-budget remap",
            "coordinate_precision": "MIROC64Gaussian latitudes checked against Can64to1e-10degree,then assigned same nominal coordinates;difference about1e-13degree is metadata roundoff,not spatial extrapolation",
            "native_coordinate_details": "pr(time,lat,lon);ua for7models(time,plev,lat,lon),Pa rounded for floating tails;Nor(time,lev,lat,lon),ps(time,lat,lon),a/b(lev),p0scalar;monthly log-pressure interpolation using a*p0+b*ps",
            "saved_dimension_order": {"samples": ["model8", "period2", "year_in_window10", "lat32", "lon128"],
                                      "group_climate": ["group3", "period2", "lat32", "lon128"],
                                      "group_change": ["group3", "lat32", "lon128"],
                                      "upper": "same45x360grid asn066;UKESMTAScalendar corrected;Caninitial/laterSICfromverifiedrawendpointcache;correspondingfour-membergroupupdated"},
            "latitude_range": [float(latitude.min()), float(latitude.max())],
            "land_region_boxes": land_boxes, "land_region_weighting": "110mland mask;sin(latitude-boundary) cell areas;not significance-selected",
            "Nor_CDR_pr": "Existing verified duration-weighted monthly product from complete published daily records;not original publishedAmon",
            "UKESM_TAS_correction": tas_correction,
            "CanESM5_SIC_correction": can_ice_correction,
            "no_new_spatial_smoothing": True, "drift_correction": "none",
            "limitations": ["Top and bottom have intentionally different map extents;not identical Arctic zooms",
                            "A30m/s contour is not jet-axis latitude;strength and meridional displacement can both move it",
                            "No significance dots,causal attribution,or statistical wave-train diagnosis in this exploratory six-panel figure",
                            "Regional land boxes are descriptive approximations,not exact political continents",
                            "Upper original standalone n066 retainedforaudit;newcombinedproductcorrectsUKESMwindowandCanSICsource"],
            "small_sources": [{"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                              for path in [membership_path, upper_directory / "grouped_tas_ice_fields.npz", monthly_file,
                                           can_ice_path, can_ice_directory / "arctic_sia_pattern_metadata.json",
                                           can_ice_directory / "arctic_sia_pattern_qa.csv"]],
            "large_raw_sources": "Identified by input_audit.csv with tracking_id,file size,real selected dates,dimensions;selected monthly arrays haveSHA256inendpoint_diagnostics.csv;no full multiGB file hash implied"}
(destination / "metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
print(pd.DataFrame(group_rows).to_string(index=False), flush=True)
print(json.dumps(tas_correction, ensure_ascii=False, indent=2), flush=True)
