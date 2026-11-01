# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import pandas as pd
import xarray as xr
from scipy.ndimage import label
from shapely.geometry import box, mapping
from shapely.ops import unary_union


# %% 2. 参数与路径设置

root = Path(r"F:\sea\Antic sea ice constraint")
source = Path(r"F:\sea\mld-argo\Argo_mixedlayers_all_04142022.nc")
destination = root / "output/s02_sia/s02_arctic_region/p21_grouped_mld/n069_argo_mld_patterns_timeseries"
land_source = root / "output/s02_sia/s02_arctic_region/p21_grouped_mld/n065_observed_march_mld_patterns/map_land.geojson"
assert not destination.exists(), "正式结果重跑前先归档。"
years = np.arange(2005, 2023)
month = 3
map_extent = [-75, 30, 40, 85]
map_grid_degrees = 1
series_grid_degrees = 2
primary_minimum_years = 9
minimum_year_sensitivity = [6, 9]
line_minimum_coverage = .5
low_map_count = 5
region_names = ["Labrador-Irminger", "Nordic Seas"]
# 只用作两个海区的地理定位窗；实际平均是Argo覆盖决定的连通格点集合。
locating_windows = [[45, 65, -70, -10], [65, 85, -40, 30]]
methods = ["density_threshold", "density_hybrid"]
variables = ["dt_mld", "da_mld"]
annual_rows, footprint_rows, support_rows, summary_rows = [], [], [], []
footprint_features = []


# %% 3. 读取逐剖面结果，明确时间、经纬度与算法字段

with xr.open_dataset(source, decode_times=False) as dataset:
    assert dict(dataset.sizes) == {"iNPROF": 2480040}
    fields = ["profiledate", "profilelat", "profilelon", "floatnumber", "profilenumber", "status", "dt_mld", "da_mld"]
    for name in fields:
        assert dataset[name].dims == ("iNPROF",)
    assert dataset.dt_mld.attrs["units"] == dataset.da_mld.attrs["units"] == "meters"
    assert dataset.status.attrs["units"] == "1 corresponds to delayed mode profile"
    date_numbers = dataset.profiledate.values.astype(float)
    # MATLAB datenum(1970-01-01)=719529，不能当作Unix日数或普通Python ordinal。
    dates = pd.to_datetime(date_numbers - 719529, unit="D", origin="unix")
    latitude = dataset.profilelat.values.astype(float)
    longitude = (dataset.profilelon.values.astype(float) + 180) % 360 - 180
    threshold = dataset.dt_mld.values.astype(float)
    hybrid = dataset.da_mld.values.astype(float)
    valid = (np.isfinite(date_numbers) & np.isfinite(latitude) & np.isfinite(longitude)
             & (latitude >= -90) & (latitude <= 90) & np.isfinite(threshold) & (threshold > 0)
             & np.isfinite(hybrid) & (hybrid > 0))
    selected = (valid & (dates.month == month) & np.isin(dates.year, years)
                & (longitude >= map_extent[0]) & (longitude < map_extent[1])
                & (latitude >= map_extent[2]) & (latitude < map_extent[3]))
    selected_indices = np.flatnonzero(selected)
    profiles = pd.DataFrame({
        "source_profile_index": selected_indices, "date": dates[selected], "matlab_datenum": date_numbers[selected],
        "year": dates.year[selected], "month": dates.month[selected], "latitude": latitude[selected],
        "longitude": longitude[selected], "dt_mld": threshold[selected], "da_mld": hybrid[selected],
        "float_number": dataset.floatnumber.values[selected].astype(np.int64),
        "profile_number": dataset.profilenumber.values[selected].astype(np.int64),
        "delayed_mode": dataset.status.values[selected].astype(np.int8)})
    source_attributes = {"global": dict(dataset.attrs),
                         "variables": {name: dict(dataset[name].attrs) for name in fields}}
source_date_minimum, source_date_maximum = str(dates.min()), str(dates.max())
global_profile_year_counts = pd.Series(dates.year).value_counts().sort_index()
selected_status_counts = dict(zip(*np.unique(profiles.delayed_mode.values, return_counts=True)))
raw_selected_profile_count = len(profiles)
profile_keys = ["float_number", "profile_number", "matlab_datenum"]
duplicate_flag = profiles.duplicated(profile_keys, keep=False)
duplicate_key_count = len(profiles.loc[duplicate_flag, profile_keys].drop_duplicates())
# 同一观测优先采用延迟质控版本；同优先级结果不一致时不任意挑一个。
preferred = profiles.loc[profiles.delayed_mode == profiles.groupby(profile_keys).delayed_mode.transform("max")].copy()
preferred_uniqueness = preferred.groupby(profile_keys)[["latitude", "longitude", "dt_mld", "da_mld"]].nunique()
ambiguous_keys = preferred_uniqueness.index[(preferred_uniqueness > 1).any(axis=1)]
ambiguous = preferred.loc[pd.MultiIndex.from_frame(preferred[profile_keys]).isin(ambiguous_keys)].copy()
preferred = preferred.loc[~pd.MultiIndex.from_frame(preferred[profile_keys]).isin(ambiguous_keys)]
retained = preferred.sort_values("source_profile_index").drop_duplicates(profile_keys, keep="first").copy()
duplicate_audit = profiles.loc[duplicate_flag].copy()
duplicate_audit["retained"] = duplicate_audit.source_profile_index.isin(retained.source_profile_index)
duplicate_audit["ambiguous_equal_priority_key"] = pd.MultiIndex.from_frame(duplicate_audit[profile_keys]).isin(ambiguous_keys)
profiles = retained.reset_index(drop=True)
assert not profiles.duplicated(profile_keys).any()
profiles["map_iy"] = np.minimum(np.floor(profiles.latitude + 90).astype(int), 179)
profiles["map_ix"] = np.floor(profiles.longitude + 180).astype(int)
profiles["series_iy"] = np.minimum(np.floor((profiles.latitude + 90) / series_grid_degrees).astype(int), 89)
profiles["series_ix"] = np.floor((profiles.longitude + 180) / series_grid_degrees).astype(int)
print(f"[source] {len(date_numbers)} profiles, {source_date_minimum} to {source_date_maximum}", flush=True)
print(f"[selected] {len(profiles)} paired-method March profiles in map domain", flush=True)


# %% 4. 构建1度三月pattern，不填海洋空白、不平滑

map_longitude, map_latitude = np.arange(-179.5, 180, 1), np.arange(-89.5, 90, 1)
map_groups = profiles.groupby(["map_iy", "map_ix"])
map_table = map_groups.agg(
    profile_count=("dt_mld", "size"), float_count=("float_number", "nunique"),
    year_count=("year", "nunique"), delayed_fraction=("delayed_mode", "mean"),
    dt_mean_m=("dt_mld", "mean"), da_mean_m=("da_mld", "mean"),
    dt_median_m=("dt_mld", "median"), da_median_m=("da_mld", "median"),
    dt_individual_max_m=("dt_mld", "max"), da_individual_max_m=("da_mld", "max")).reset_index()
map_table["longitude"] = map_longitude[map_table.map_ix]
map_table["latitude"] = map_latitude[map_table.map_iy]
map_fields = np.full((2, 180, 360), np.nan)
map_counts = np.zeros((180, 360), dtype=int)
map_year_counts = np.zeros((180, 360), dtype=int)
for index, column in enumerate(["dt_mean_m", "da_mean_m"]):
    map_fields[index, map_table.map_iy, map_table.map_ix] = map_table[column]
map_counts[map_table.map_iy, map_table.map_ix] = map_table.profile_count
map_year_counts[map_table.map_iy, map_table.map_ix] = map_table.year_count
assert np.nanmax(map_fields) < 2000


# %% 5. 先按年和2度格点平均，再按有效观测格点面积加权

series_latitude = np.arange(-89., 90., series_grid_degrees)
series_longitude = np.arange(-179., 180., series_grid_degrees)
series_lon_grid, series_lat_grid = np.meshgrid(series_longitude, series_latitude)
cell_groups = profiles.groupby(["year", "series_iy", "series_ix"])
annual_cells = cell_groups.agg(
    dt_mean_m=("dt_mld", "mean"), da_mean_m=("da_mld", "mean"),
    profile_count=("dt_mld", "size"), float_count=("float_number", "nunique"),
    delayed_fraction=("delayed_mode", "mean")).reset_index()
annual_cells["longitude"] = series_longitude[annual_cells.series_ix]
annual_cells["latitude"] = series_latitude[annual_cells.series_iy]
annual_cells["area_weight"] = np.cos(np.deg2rad(annual_cells.latitude))
cell_support = annual_cells.groupby(["series_iy", "series_ix"]).agg(
    march_year_count=("year", "nunique"), dt_reference_mean_m=("dt_mean_m", "mean"),
    da_reference_mean_m=("da_mean_m", "mean"), total_profiles=("profile_count", "sum")).reset_index()
year_support = np.zeros(series_lat_grid.shape, dtype=int)
year_support[cell_support.series_iy, cell_support.series_ix] = cell_support.march_year_count
cell_reference = np.full((2,) + series_lat_grid.shape, np.nan)
cell_reference[0, cell_support.series_iy, cell_support.series_ix] = cell_support.dt_reference_mean_m
cell_reference[1, cell_support.series_iy, cell_support.series_ix] = cell_support.da_reference_mean_m

for region_index, (region, window) in enumerate(zip(region_names, locating_windows)):
    lat_min, lat_max, lon_min, lon_max = window
    locating = ((series_lat_grid >= lat_min) & (series_lat_grid < lat_max)
                & (series_lon_grid >= lon_min) & (series_lon_grid < lon_max))
    for minimum_years in minimum_year_sensitivity:
        candidate = locating & (year_support >= minimum_years)
        components, component_count = label(candidate, structure=np.ones((3, 3), dtype=int))
        component_sizes = np.bincount(components.ravel())
        component_sizes[0] = 0
        assert component_count > 0
        footprint = components == int(component_sizes.argmax())
        footprint_weights = np.cos(np.deg2rad(series_lat_grid[footprint]))
        footprint_total_weight = footprint_weights.sum()
        if minimum_years == primary_minimum_years:
            cells = []
            for row, column in np.argwhere(footprint):
                lon, lat = series_longitude[column], series_latitude[row]
                footprint_rows.append({"region": region, "series_iy": int(row), "series_ix": int(column),
                                       "longitude": lon, "latitude": lat,
                                       "march_year_count": int(year_support[row, column]),
                                       "area_weight": float(np.cos(np.deg2rad(lat)))})
                cells.append(box(lon - 1, lat - 1, lon + 1, lat + 1))
            footprint_features.append({"type": "Feature", "properties": {"region": region, "number": region_index + 1},
                                       "geometry": mapping(unary_union(cells))})
        support_rows.append({"region": region, "minimum_march_years": minimum_years,
                             "candidate_cells": int(candidate.sum()), "connected_components": int(component_count),
                             "footprint_cells": int(footprint.sum()),
                             "latitude_center_min": float(series_lat_grid[footprint].min()),
                             "latitude_center_max": float(series_lat_grid[footprint].max()),
                             "longitude_center_min": float(series_lon_grid[footprint].min()),
                             "longitude_center_max": float(series_lon_grid[footprint].max())})
        footprint_profile = footprint[profiles.series_iy, profiles.series_ix]
        for year in years:
            region_profiles = profiles.loc[footprint_profile & (profiles.year == year)]
            year_cells = annual_cells.loc[annual_cells.year == year]
            year_cells = year_cells.loc[footprint[year_cells.series_iy, year_cells.series_ix]]
            cell_weights = year_cells.area_weight.to_numpy()
            coverage = float(cell_weights.sum() / footprint_total_weight)
            for method_index, (method, column, variable) in enumerate(zip(methods, ["dt_mean_m", "da_mean_m"], variables)):
                mean = np.nan
                same_support_reference = np.nan
                if len(year_cells):
                    mean = float(np.average(year_cells[column], weights=cell_weights))
                    same_support_reference = float(np.average(
                        cell_reference[method_index, year_cells.series_iy, year_cells.series_ix], weights=cell_weights))
                annual_rows.append({"region": region, "method": method, "minimum_march_years": minimum_years,
                                    "year": int(year), "mean_m": mean,
                                    "same_support_climatology_m": same_support_reference,
                                    "same_support_anomaly_m": mean - same_support_reference,
                                    "profile_unweighted_mean_m": float(region_profiles[variable].mean()),
                                    "profile_median_m": float(region_profiles[variable].median()),
                                    "profile_p90_m": float(region_profiles[variable].quantile(.9)),
                                    "profile_maximum_m": float(region_profiles[variable].max()),
                                    "profile_count": len(region_profiles), "float_count": region_profiles.float_number.nunique(),
                                    "sampled_cells": len(year_cells), "footprint_cells": int(footprint.sum()),
                                    "area_coverage": coverage,
                                    "delayed_mode_fraction": float(region_profiles.delayed_mode.mean()),
                                    "connect_line": bool(coverage >= line_minimum_coverage)})
annual = pd.DataFrame(annual_rows)
primary = annual.loc[annual.minimum_march_years == primary_minimum_years]
for region in region_names:
    for method in methods:
        records = primary.loc[(primary.region == region) & (primary.method == method)]
        adequate = records.loc[records.connect_line & np.isfinite(records.mean_m)]
        peak = adequate.loc[adequate.mean_m.idxmax()]
        minimum = adequate.loc[adequate.mean_m.idxmin()]
        summary_rows.append({"region": region, "method": method, "adequate_coverage_years": len(adequate),
                             "years_total": len(records), "maximum_mean_m": peak.mean_m, "maximum_year": int(peak.year),
                             "minimum_mean_m": minimum.mean_m, "minimum_year": int(minimum.year),
                             "minimum_coverage": float(records.area_coverage.min()),
                             "mean_coverage": float(records.area_coverage.mean()),
                             "maximum_individual_profile_m": float(records.profile_maximum_m.max()),
                             "total_profiles": int(records.profile_count.sum()),
                             "mean_of_annual_sampled_area_means_m": float(records.mean_m.mean()),
                             "mean_of_adequate_coverage_annual_means_m": float(adequate.mean_m.mean()),
                             "profile_count_2015": int(records.loc[records.year == 2015, "profile_count"].iloc[0]),
                             "coverage_2015": float(records.loc[records.year == 2015, "area_coverage"].iloc[0])})


# %% 6. 保存图用数组、年序列、覆盖和源信息

destination.mkdir(parents=True)
profiles.to_csv(destination / "march_profiles_in_map_domain.csv", index=False)
duplicate_audit.to_csv(destination / "duplicate_profile_audit.csv", index=False)
ambiguous.to_csv(destination / "ambiguous_profile_records.csv", index=False)
map_table.to_csv(destination / "march_1degree_map_cells.csv", index=False)
annual_cells.to_csv(destination / "march_2degree_annual_cells.csv", index=False)
cell_support.to_csv(destination / "march_2degree_cell_support.csv", index=False)
annual.to_csv(destination / "regional_march_timeseries.csv", index=False)
pd.DataFrame(footprint_rows).to_csv(destination / "observation_footprint_cells.csv", index=False)
pd.DataFrame(support_rows).to_csv(destination / "footprint_support_summary.csv", index=False)
pd.DataFrame(summary_rows).to_csv(destination / "regional_summary.csv", index=False)
global_profile_year_counts.rename_axis("year").rename("profile_count").to_csv(destination / "global_source_year_counts.csv")
np.savez_compressed(destination / "march_patterns.npz", longitude=map_longitude, latitude=map_latitude,
                    mld_m=map_fields, profile_count=map_counts, march_year_count=map_year_counts,
                    low_profile_support=(map_counts > 0) & (map_counts < low_map_count),
                    method=np.array(methods))
(destination / "observation_footprints.geojson").write_text(
    json.dumps({"type": "FeatureCollection", "features": footprint_features}, ensure_ascii=False), encoding="utf-8")
shutil.copyfile(land_source, destination / "map_land.geojson")
metadata = {
    "status": "calculation_completed_independent_QC_separate", "source": str(source),
    "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(), "source_bytes": source.stat().st_size,
    "source_profiles": int(date_numbers.size), "source_actual_date_minimum": source_date_minimum,
    "source_actual_date_maximum": source_date_maximum,
    "source_dimensions": "All fields are1-D iNPROF2480040; no time/lat/lon grid or vertical T/S profile dimension",
    "source_attributes": source_attributes,
    "date_conversion": "MATLAB datenum minus719529 to Unix days; fractions retained; UTC-naive source dates",
    "years": years.tolist(), "month": month, "selected_map_profiles_before_deduplication": raw_selected_profile_count,
    "selected_map_profiles": len(profiles),
    "selected_profile_status_counts_before_deduplication": {str(k): int(v) for k, v in selected_status_counts.items()},
    "duplicate_profile_keys": duplicate_key_count, "ambiguous_equal_priority_keys_excluded": len(ambiguous_keys),
    "retained_profile_key_duplicates": 0,
    "duplicate_policy": "Key=float number/profile number/exact source date; prefer delayed mode, collapse identical equal-priority records, exclude conflicting equal-priority key from paired-method comparison; complete source-index audit saved",
    "map": {"extent": map_extent, "bin_degrees": map_grid_degrees,
            "estimator": "Arithmetic mean of paired valid individual profile MLDs in each1degree March bin over2005-2022",
            "missing": "No profile means NaN, not zero; no spatial mapping or smoothing",
            "low_support_dots": f"Fewer than{low_map_count} pooled profiles in cell; not significance",
            "method_maxima_m": np.nanmax(map_fields, axis=(1, 2)).tolist(), "color_limits_m": [0, 2000]},
    "series": {"bin_degrees": series_grid_degrees, "primary_minimum_distinct_march_years": primary_minimum_years,
               "minimum_years_sensitivity": minimum_year_sensitivity,
               "regions": region_names, "locating_windows_latmin_latmax_lonmin_lonmax": locating_windows,
               "footprint": "Largest8-neighbor connected set of2degree cells with data in>=9 of18 March years; determined only from sampling, not MLD depth; fixed across methods/years",
               "estimator": "Within each March/year/cell profile mean, then cos(latitude) weighted mean of sampled footprint cells only",
               "denominator": "Only that year's observed cells; absent cells not zero; coverage denominator is the full fixed footprint",
               "line_minimum_area_coverage": line_minimum_coverage,
               "low_coverage_display": "Finite low-coverage means retained as open markers, not joined by lines; missing means remain gaps",
               "sampling_composition": "Same-observed-cell climatology and anomalies saved as diagnostics; main series are actual sampled-area means, not gap-filled estimates"},
    "methods": {"density_threshold": "Author-provided dt_mld", "density_hybrid": "Author-provided da_mld (hybrid density algorithm)",
                "not_recomputed": "No complete T/S-depth profiles in this file; author algorithm not rerun; reference depth/details not inferred from variable name alone"},
    "limitations": ["Argo sampling is uneven; binning reduces dense-track domination but cannot remove all within-bin float/time dependence.",
                    "No bootstrap confidence intervals or trend significance tests; profile counts are not independent replicate counts.",
                    "Real-time and delayed-mode profiles both retained as delivered; fractions saved by year.",
                    "2015 has a marked source/sample reduction; low regional coverage is flagged, not interpreted as a physical event.",
                    "Two-degree fixed sampling cores are not whole-sea or old/user-drawn polygon averages; mapped outlines show their actual support.",
                    "March2022 exists although the source stops inApril2022; this is not a full2022 annual mean.",
                    "Argo typically does not sample beyond roughly2000m; source maxima near2000m do not establish the absence of deeper mixing.",
                    "Profile-pooled map means and equal-year cell-reference means have different temporal weights; they need not be equal."],
    "references": ["https://mixedlayer.ucsd.edu/", "https://doi.org/10.1002/2017GL073426", "https://doi.org/10.1175/2009JTECHO543.1"],
    "software": {"numpy": np.__version__, "pandas": pd.__version__, "xarray": xr.__version__}}
(destination / "metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
print(pd.DataFrame(support_rows).to_string(index=False), flush=True)
print(pd.DataFrame(summary_rows).to_string(index=False), flush=True)
