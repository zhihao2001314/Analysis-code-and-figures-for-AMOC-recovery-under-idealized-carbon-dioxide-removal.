# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import xarray as xr
from scipy.signal import find_peaks


# %% 2. 参数与路径设置

root = Path("F:/sea/Antic sea ice constraint")
base = root / "output/s05_atmosphere/p05_amoc_recovery_composites"
source_directory = base / "n44_seasonal_jet_hydrology_samples"
destination = base / "n45_jet_hydrology_diagnostics"
assert not destination.exists(), "正式诊断已存在，须先归档。"
ridge_latitude_bounds = [25.0, 65.0]
ridge_minimum_u = 15.0
ridge_ambiguity_fraction = 0.90
ridge_peak_separation_degrees = 8.0
ridge_maximum_neighbor_jump = 8.0
profile_windows = {"primary25_65": [25.0, 65.0], "sensitivity20_70": [20.0, 70.0], "sensitivity30_65": [30.0, 65.0]}
jet_sectors = {"North Atlantic": [-60.0, 0.0], "North Pacific": [140.0, 240.0], "Northern Hemisphere": [0.0, 360.0]}
rain_sectors = {"Global": [0.0, 360.0], "Atlantic sector": [-60.0, 20.0], "Pacific sector": [140.0, 240.0]}
land_boxes = {"North America": [-170, -50, 15, 75], "Europe": [-10, 45, 35, 72],
              "Asia": [45, 180, 5, 75], "North Africa": [-20, 45, 0, 35]}
ridge_sensitivity = {"primary": [15.0, 0.90], "minimum10": [10.0, 0.90],
                     "minimum20": [20.0, 0.90], "ambiguity85": [15.0, 0.85], "ambiguity95": [15.0, 0.95]}


# %% 3. 读取固定成员的季节端点

metadata = json.loads((source_directory / "metadata.json").read_text(encoding="utf-8"))
membership = pd.read_csv(source_directory / "membership.csv")
with xr.open_dataset(source_directory / "member_climatology.nc") as source:
    member_climate = source.load()
with xr.open_dataset(source_directory / "group_climatology.nc") as source:
    group_climate = source.load()
with xr.open_dataset(source_directory / "member_changes.nc") as source:
    member_changes = source.load()
with xr.open_dataset(source_directory / "group_changes.nc") as source:
    group_changes = source.load()
product_climate = xr.concat([member_climate.rename(model="product"), group_climate.rename(group="product")], "product", join="exact")
product_changes = xr.concat([member_changes.rename(model="product"), group_changes.rename(group="product")], "product", join="exact")
latitude, longitude = product_climate.lat.values, product_climate.lon.values
signed_longitude = (longitude + 180) % 360 - 180
products = product_climate.product.values.tolist()
# model改名product后不可继承原NetCDF的短字符串编码，否则长组名在保存时会截断。
product_climate = product_climate.assign_coords(product=np.asarray(products, dtype=str))
product_changes = product_changes.assign_coords(product=np.asarray(products, dtype=str))
periods, seasons = metadata["periods"], metadata["seasons"]
with np.load(source_directory / "nh_land_support.npz") as source:
    land_latitude, land_longitude = source["latitude"], source["longitude"]
    nh_land = source["land_mask"]
    nh_area = source["cell_area_weights"]
np.testing.assert_array_equal(land_latitude, latitude[latitude > 0])
np.testing.assert_array_equal(land_longitude, longitude)
latitude_edges = np.concatenate([[-90.0], (latitude[:-1] + latitude[1:]) / 2, [90.0]])


# %% 4. 沿经度识别局地强风脊线并拒绝双峰和跳跃

ridge_shape = (len(products), len(periods), len(seasons), len(longitude))
candidate_latitude = np.full(ridge_shape, np.nan)
candidate_speed = np.full(ridge_shape, np.nan)
second_peak_ratio = np.full(ridge_shape, np.nan)
base_flags = np.zeros(ridge_shape, dtype=np.int16)
ridge_indices = np.flatnonzero((latitude >= ridge_latitude_bounds[0]) & (latitude <= ridge_latitude_bounds[1]))
ridge_latitude = latitude[ridge_indices]
for product_index, product in enumerate(products):
    for period_index, period in enumerate(periods):
        for season_index, season in enumerate(seasons):
            field = product_climate.u250.sel(product=product, period=period, season=season).values[ridge_indices]
            for column in range(len(longitude)):
                values = field[:, column]
                key = (product_index, period_index, season_index, column)
                if not np.isfinite(values).all():
                    base_flags[key] = 16
                    continue
                peak = int(np.argmax(values))
                if peak == 0 or peak == len(values) - 1:
                    base_flags[key] = 1
                    continue
                coefficients = np.polyfit(ridge_latitude[peak - 1:peak + 2], values[peak - 1:peak + 2], 2)
                vertex = -coefficients[1] / (2 * coefficients[0])
                assert coefficients[0] < 0 and ridge_latitude[peak - 1] <= vertex <= ridge_latitude[peak + 1]
                candidate_latitude[key] = vertex
                candidate_speed[key] = float(np.polyval(coefficients, vertex))
                peaks, properties = find_peaks(values)
                competitors = peaks[np.abs(ridge_latitude[peaks] - ridge_latitude[peak]) >= ridge_peak_separation_degrees]
                second_peak_ratio[key] = float(np.max(values[competitors]) / values[peak]) if len(competitors) else 0.0
ridge_versions, ridge_flags = [], []
for version, (minimum_u, ambiguity_fraction) in ridge_sensitivity.items():
    flags = base_flags.copy()
    flags[(flags == 0) & (candidate_speed < minimum_u)] |= 2
    flags[(base_flags == 0) & (second_peak_ratio >= ambiguity_fraction)] |= 4
    accepted = np.where(flags == 0, candidate_latitude, np.nan)
    # 经度周期邻接；大跳跃两侧都断线，不人为选取能产生预期位移的分支。
    jump = np.abs(accepted - np.roll(accepted, 1, axis=-1)) > ridge_maximum_neighbor_jump
    bad_neighbor = jump | np.roll(jump, -1, axis=-1)
    flags[bad_neighbor] |= 8
    ridge_versions.append(np.where(flags == 0, candidate_latitude, np.nan))
    ridge_flags.append(flags)
ridge = xr.Dataset({"latitude": (("version", "product", "period", "season", "lon"), np.stack(ridge_versions)),
                    "quality_flag": (("version", "product", "period", "season", "lon"), np.stack(ridge_flags)),
                    "candidate_latitude": (("product", "period", "season", "lon"), candidate_latitude),
                    "candidate_speed": (("product", "period", "season", "lon"), candidate_speed),
                    "second_peak_ratio": (("product", "period", "season", "lon"), second_peak_ratio)},
                   coords={"version": list(ridge_sensitivity), "product": products, "period": periods, "season": seasons, "lon": longitude})
ridge.latitude.attrs["units"] = "degree_north"
ridge["paired_latitude_change"] = ridge.latitude.sel(period="later", drop=True) - ridge.latitude.sel(period="initial", drop=True)
ridge.paired_latitude_change.attrs["units"] = "degree_north"
ridge.candidate_speed.attrs["units"] = "m s-1"
ridge.quality_flag.attrs["meaning"] = "bit1 search-boundary;bit2 weak;bit4 comparable separated peak;bit8 adjacent longitude jump;bit16 missing"
ridge_summary = []
for version in ridge_sensitivity:
    for product in products:
        for season in seasons:
            initial = ridge.latitude.sel(version=version, product=product, season=season, period="initial").values
            later = ridge.latitude.sel(version=version, product=product, season=season, period="later").values
            for sector, (west, east) in jet_sectors.items():
                coordinate = signed_longitude if west < 0 else longitude
                chosen = (coordinate >= west) & (coordinate < east)
                paired = chosen & np.isfinite(initial) & np.isfinite(later)
                ridge_summary.append({"product": product, "season": season, "version": version, "sector": sector,
                                      "longitude_count": int(chosen.sum()), "paired_valid_count": int(paired.sum()),
                                      "paired_coverage_fraction": float(paired.sum() / chosen.sum()),
                                      "mean_local_displacement_degree_north": float(np.mean(later[paired] - initial[paired])) if paired.any() else np.nan,
                                      "northward_count": int(np.sum((later - initial)[paired] > 0)),
                                      "southward_count": int(np.sum((later - initial)[paired] < 0))})


# %% 5. 分海盆纬度剖面、风峰和加权中心分开保存

jet_profiles = []
for sector, (west, east) in jet_sectors.items():
    coordinate = signed_longitude if west < 0 else longitude
    chosen = (coordinate >= west) & (coordinate < east)
    jet_profiles.append(product_climate.u250.isel(lon=np.flatnonzero(chosen)).mean("lon", skipna=False).expand_dims(sector=[sector]))
jet_profile = xr.concat(jet_profiles, "sector", join="exact").transpose("product", "period", "season", "sector", "lat")
profile_records = []
for product in products:
    for season in seasons:
        for sector in jet_sectors:
            for window, bounds in profile_windows.items():
                chosen = np.flatnonzero((latitude >= bounds[0]) & (latitude <= bounds[1]))
                selected_latitude = latitude[chosen]
                # 一维纬度中心的权重为正u²乘纬度带宽，不额外混入cos纬度的面积权重。
                edges = np.concatenate([[bounds[0]], (selected_latitude[:-1] + selected_latitude[1:]) / 2, [bounds[1]]])
                latitude_width = np.diff(edges)
                for period in periods:
                    values = jet_profile.sel(product=product, period=period, season=season, sector=sector).values[chosen]
                    assert np.isfinite(values).all()
                    peak = int(np.argmax(values))
                    status, peak_latitude, peak_u = "ok", np.nan, np.nan
                    if peak == 0 or peak == len(values) - 1:
                        status = "search_boundary"
                    elif values[peak] <= 0:
                        status = "no_westerly_peak"
                    else:
                        coefficients = np.polyfit(selected_latitude[peak - 1:peak + 2], values[peak - 1:peak + 2], 2)
                        vertex = -coefficients[1] / (2 * coefficients[0])
                        assert coefficients[0] < 0 and selected_latitude[peak - 1] <= vertex <= selected_latitude[peak + 1]
                        peak_latitude, peak_u = float(vertex), float(np.polyval(coefficients, vertex))
                        peaks, properties = find_peaks(values)
                        competitors = peaks[np.abs(selected_latitude[peaks] - selected_latitude[peak]) >= ridge_peak_separation_degrees]
                        if len(competitors) and np.max(values[competitors]) >= ridge_ambiguity_fraction * values[peak]:
                            status = "comparable_two_peaks"
                    weights = np.maximum(values, 0) ** 2 * latitude_width
                    centroid = float(np.sum(weights * selected_latitude) / weights.sum()) if weights.sum() > 0 else np.nan
                    profile_records.append({"product": product, "season": season, "sector": sector, "window": window, "period": period,
                                            "status": status, "peak_latitude": peak_latitude if status == "ok" else np.nan,
                                            "candidate_peak_latitude": peak_latitude, "peak_speed": peak_u,
                                            "grid_maximum_latitude": float(selected_latitude[peak]), "u2_latitude_centroid": centroid})
profile_metrics = pd.DataFrame(profile_records)
initial_metrics = profile_metrics.loc[profile_metrics.period.eq("initial")].drop(columns="period")
later_metrics = profile_metrics.loc[profile_metrics.period.eq("later")].drop(columns="period")
profile_change = initial_metrics.merge(later_metrics, on=["product", "season", "sector", "window"], suffixes=("_initial", "_later"), validate="one_to_one")
for quantity in ["peak_latitude", "peak_speed", "u2_latitude_centroid"]:
    profile_change["delta_" + quantity] = profile_change[quantity + "_later"] - profile_change[quantity + "_initial"]


# %% 6. 热带降水剖面与雨量中位纬度

rain_profiles = []
for sector, (west, east) in rain_sectors.items():
    coordinate = signed_longitude if west < 0 else longitude
    chosen = (coordinate >= west) & (coordinate < east)
    rain_profiles.append(product_climate.pr.isel(lon=np.flatnonzero(chosen)).mean("lon", skipna=False).expand_dims(sector=[sector]))
rain_profile = xr.concat(rain_profiles, "sector", join="exact").transpose("product", "period", "season", "sector", "lat")
south_edges, north_edges = np.maximum(latitude_edges[:-1], -20), np.minimum(latitude_edges[1:], 20)
tropical = np.flatnonzero(north_edges > south_edges)
tropical_weights = np.sin(np.deg2rad(north_edges[tropical])) - np.sin(np.deg2rad(south_edges[tropical]))
rain_records = []
for product in products:
    for season in seasons:
        for sector in rain_sectors:
            for period in periods:
                values = rain_profile.sel(product=product, period=period, season=season, sector=sector).values[tropical]
                assert np.isfinite(values).all() and np.all(values >= 0)
                mass = values * tropical_weights
                cumulative = np.cumsum(mass)
                median_index = int(np.searchsorted(cumulative, cumulative[-1] / 2))
                previous = 0 if median_index == 0 else cumulative[median_index - 1]
                fraction = (cumulative[-1] / 2 - previous) / mass[median_index]
                median_sin = (np.sin(np.deg2rad(south_edges[tropical][median_index]))
                              + fraction * tropical_weights[median_index])
                median = float(np.rad2deg(np.arcsin(median_sin)))
                centroid = float(np.sum(latitude[tropical] * mass) / mass.sum())
                rain_records.append({"product": product, "season": season, "sector": sector, "period": period,
                                     "precipitation_median_latitude": median, "precipitation_centroid_latitude": centroid,
                                     "grid_precipitation_peak_latitude": float(latitude[tropical[np.argmax(values)]]),
                                     "north_fraction": float(mass[latitude[tropical] > 0].sum() / mass.sum())})
rain_metrics = pd.DataFrame(rain_records)
rain_initial = rain_metrics.loc[rain_metrics.period.eq("initial")].drop(columns="period")
rain_later = rain_metrics.loc[rain_metrics.period.eq("later")].drop(columns="period")
rain_changes = rain_initial.merge(rain_later, on=["product", "season", "sector"], suffixes=("_initial", "_later"), validate="one_to_one")
for quantity in ["precipitation_median_latitude", "precipitation_centroid_latitude", "north_fraction"]:
    rain_changes["delta_" + quantity] = rain_changes[quantity + "_later"] - rain_changes[quantity + "_initial"]


# %% 7. 陆地及副极地海洋水文、组差和符号一致性

regional_records = []
nh_lon, nh_lat = np.meshgrid(signed_longitude, land_latitude)
region_masks = {}
for region, bounds in land_boxes.items():
    region_masks[region] = nh_land & (nh_lon >= bounds[0]) & (nh_lon < bounds[1]) & (nh_lat >= bounds[2]) & (nh_lat < bounds[3])
region_masks["Subpolar Atlantic ocean"] = ~nh_land & (nh_lon >= -60) & (nh_lon < 0) & (nh_lat >= 45) & (nh_lat < 65)
for product in products:
    for season in seasons:
        for region, mask in region_masks.items():
            row = {"product": product, "season": season, "region": region, "grid_cells": int(mask.sum())}
            for variable in ["pr", "evspsbl", "pe"]:
                field = product_changes[variable].sel(product=product, season=season, lat=land_latitude).values
                assert np.isfinite(field[mask]).all()
                row["delta_" + variable + "_mm_day"] = float(np.average(field[mask], weights=nh_area[mask]))
            assert abs(row["delta_pe_mm_day"] - row["delta_pr_mm_day"] + row["delta_evspsbl_mm_day"]) < 1e-12
            regional_records.append(row)
contrast = group_changes.sel(group="better_nonrecovery", drop=True) - group_changes.sel(group="better_recovery", drop=True)
agreement_groups = []
for group in metadata["groups"]:
    selected_models = membership.loc[membership.group.eq(group), "model"].tolist()
    variables = {}
    for variable in ["u250", "pr", "evspsbl", "pe"]:
        member_signs = np.sign(member_changes[variable].sel(model=selected_models))
        same_sign = member_signs == np.sign(group_changes[variable].sel(group=group))
        variables[variable] = same_sign.mean("model").where(group_changes[variable].sel(group=group).notnull())
    agreement_groups.append(xr.Dataset(variables).expand_dims(group=[group]))
agreement = xr.concat(agreement_groups, "group", join="exact")


# %% 8. 保存诊断、阈值与解释限制

destination.mkdir(parents=True)
ridge.to_netcdf(destination / "u250_ridges.nc")
jet_profile.to_dataset(name="u250").to_netcdf(destination / "u250_sector_profiles.nc")
rain_profile.to_dataset(name="pr").to_netcdf(destination / "tropical_precipitation_profiles.nc")
contrast.to_netcdf(destination / "nonrecovery_minus_recovery.nc")
agreement.to_netcdf(destination / "model_sign_agreement.nc")
profile_metrics.to_csv(destination / "jet_endpoint_metrics.csv", index=False)
profile_change.to_csv(destination / "jet_changes_and_window_sensitivity.csv", index=False)
pd.DataFrame(ridge_summary).to_csv(destination / "ridge_coverage_and_sensitivity.csv", index=False)
rain_metrics.to_csv(destination / "rainbelt_endpoint_metrics.csv", index=False)
rain_changes.to_csv(destination / "rainbelt_changes.csv", index=False)
pd.DataFrame(regional_records).to_csv(destination / "regional_hydrology.csv", index=False)
diagnostic_metadata = {"source": str(source_directory), "products": products, "periods": periods, "seasons": seasons,
                       "jet_sectors": jet_sectors, "rain_sectors": rain_sectors,
                       "primary_ridge": {"latitude_window": ridge_latitude_bounds, "minimum_u_m_s": ridge_minimum_u,
                                         "comparable_peak_ratio": ridge_ambiguity_fraction, "peak_separation_degrees": ridge_peak_separation_degrees,
                                         "max_neighbor_jump_degrees": ridge_maximum_neighbor_jump, "sensitivity_versions": ridge_sensitivity},
                       "ridge_definition": "Local maximum of climatological U250 within25-65N at each longitude;3point concave quadratic refinement. Reject boundary/weak/comparable separated peaks;mask both sides of>8degree adjacent longitude jumps. No spatial smoothing,forced continuity or branch selection by desired change direction.",
                       "profile_definition": "Average U250 over longitude first,then locate interior maximum;3point quadratic refinement;reject boundary/comparable peaks. Speed and latitude reported separately. U-squared latitude centroid is supplementary,not a tracked axis.",
                       "profile_search_windows": profile_windows,
                       "rain_definition": "P notP-E;median over20S-20N using clipped latitude-cell area sin(north)-sin(south);P assumed constant within each cell,median interpolated in sin(latitude);also retain centroid,discretepeak andnorthfraction",
                       "land_ocean": "Rain sectors include land and ocean;not an ocean-only Atlantic mask. Regional continent diagnostics use saved110mland mask andfixedboxes;SPNA45-65N60W-0oceanonly.",
                       "dimensions": {"ridge": ["version", "product", "period", "season", "lon"], "profiles": ["product", "period", "season", "sector", "lat"],
                                      "group_contrast_and_agreement": "season,lat,lon or group,season,lat,lon"},
                       "normalization": "physical units;no division by residualGSAT",
                       "limitations": ["U250 ridge is an operational strong-westerly peak,not uniquely isolated eddy-driven or subtropical jet",
                                       "Ridge movement and sector-average peak movement are different estimators;groupfielddiagnostics are notmeanofmembermetrics",
                                       "A few degrees of inferred motion is not subgrid observational precision;native commonlatspacing about2.8degrees",
                                       "2/2/4 fixed model members;sign agreement is not significance;no borrowedt-testdots orbootstrapintervals",
                                       "9completeDJF means differfromold30calendarDJFmonth/dayweightedproducts;Annualreproducesexistingphysicalfields",
                                       "Rainmedian changes can reflect rainband redistribution/doubleITCZ,not necessarilyrigidtranslation",
                                       "Association doesnotisolateAMOCcausality;no imposedwave-trainorcontinentalcausalchain"],
                       "sources": [{"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                                   for path in [source_directory / "metadata.json", source_directory / "member_climatology.nc",
                                                source_directory / "group_climatology.nc", source_directory / "member_changes.nc",
                                                source_directory / "group_changes.nc", source_directory / "membership.csv", source_directory / "nh_land_support.npz"]]}
(destination / "metadata.json").write_text(json.dumps(diagnostic_metadata, ensure_ascii=False, indent=2), encoding="utf-8")
print(profile_change.loc[profile_change["product"].isin(metadata["groups"]) & profile_change.window.eq("primary25_65"),
                         ["product", "season", "sector", "delta_peak_latitude", "delta_peak_speed", "delta_u2_latitude_centroid"]].to_string(index=False))
print(rain_changes.loc[rain_changes["product"].isin(metadata["groups"]), ["product", "season", "sector", "delta_precipitation_median_latitude"]].to_string(index=False))
