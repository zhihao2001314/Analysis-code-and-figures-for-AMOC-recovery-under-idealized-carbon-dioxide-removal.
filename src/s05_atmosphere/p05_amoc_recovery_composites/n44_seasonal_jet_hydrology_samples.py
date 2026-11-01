# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import xarray as xr


# %% 2. 参数与路径设置

root = Path("F:/sea/Antic sea ice constraint")
base = root / "output/s05_atmosphere/p05_amoc_recovery_composites"
destination = base / "n44_seasonal_jet_hydrology_samples"
assert not destination.exists(), "正式结果已存在，须先归档。"
annual_directory = base / "n37_annual_nh_precipitation_jet"
netp_directory = base / "n39_annual_net_precipitation_circulation"
membership = pd.read_csv(annual_directory / "membership.csv")
models = membership.model.tolist()
groups = ["better_nonrecovery", "better_recovery", "larger_bias"]
periods, seasons = ["initial", "later"], ["Annual", "DJF"]
sample_count = {"Annual": 10, "DJF": 9}
first_years = {"CESM2": [1, 131], "CNRM-ESM2-1": [1850, 2120], "NorESM2-LM": [1, 271],
               "GFDL-ESM4": [1, 271], "ACCESS-ESM1-5": [101, 371], "CanESM5": [1850, 2121],
               "MIROC-ES2L": [1850, 2119], "UKESM1-0-LL": [1850, 2120]}
members = {"CESM2": "r1i1p1f1", "CNRM-ESM2-1": "r1i1p1f2", "NorESM2-LM": "r1i1p1f1",
           "GFDL-ESM4": "r1i1p1f1", "ACCESS-ESM1-5": "r1i1p1f1", "CanESM5": "r1i1p2f1",
           "MIROC-ES2L": "r1i1p1f2", "UKESM1-0-LL": "r1i1p1f2"}
pressure_Pa = 25000.0
pr_ua_audit = pd.read_csv(annual_directory / "input_audit.csv")
evaporation_audit = pd.read_csv(netp_directory / "input_audit.csv")
can_file = pr_ua_audit.loc[pr_ua_audit.model.eq("CanESM5") & pr_ua_audit.variable.eq("pr") & pr_ua_audit.period.eq("initial"), "file"].iloc[0]
with xr.open_dataset(can_file, decode_times=False) as source:
    latitude, longitude = source.lat.values.astype(float), source.lon.values.astype(float)
assert latitude.size == 64 and longitude.size == 128
wind_latitude = latitude[latitude > 0]
sample_indices = np.arange(1, 11)
model_results, input_rows, endpoint_rows = [], [], []


# %% 3. 读取月度原生场并生成完整冬季样本

for model in models:
    variables = {}
    for variable in ["pr", "evspsbl", "ua"]:
        audit = evaporation_audit if variable == "evspsbl" else pr_ua_audit
        phase_results = []
        for period_index, period in enumerate(periods):
            start_year = first_years[model][period_index]
            chosen = audit.model.eq(model) & audit.variable.eq(variable) & audit.period.eq(period)
            paths = [Path(value) for value in audit.loc[chosen, "file"].drop_duplicates()]
            assert paths
            pieces = []
            for path in paths:
                with xr.open_dataset(path, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as source:
                    derived_nor_pr = model == "NorESM2-LM" and variable == "pr" and period == "later"
                    if not derived_nor_pr:
                        assert source.attrs["source_id"] == model and source.attrs["variant_label"] == members[model]
                        assert source.attrs["experiment_id"] == ("1pctCO2" if period == "initial" else "1pctCO2-cdr")
                    selection = source.time.dt.year.isin(range(start_year, start_year + 10))
                    assert bool(selection.any())
                    if variable in ["pr", "evspsbl"]:
                        assert source[variable].dims == ("time", "lat", "lon")
                        assert source[variable].attrs["units"] == "kg m-2 s-1"
                        field = source[variable].isel(time=selection).sortby("lat").load().astype(float) * 86400
                    elif model != "NorESM2-LM":
                        assert source.ua.dims == ("time", "plev", "lat", "lon")
                        assert source.plev.attrs["units"] == "Pa" and source.ua.attrs["units"] == "m s-1"
                        field = source.ua.assign_coords(plev=np.rint(source.plev)).sel(plev=pressure_Pa, drop=True)
                        field = field.isel(time=selection).sortby("lat").sel(lat=slice(-5, 90)).load().astype(float)
                    else:
                        assert source.ua.dims == ("time", "lev", "lat", "lon")
                        assert source.ps.dims == ("time", "lat", "lon")
                        native = source.ua.isel(time=selection).sortby("lat").sel(lat=slice(-5, 90)).load().astype(float)
                        surface_pressure = source.ps.isel(time=selection).sel(lat=native.lat).load().values
                        offset, fraction = (source.a.values * float(source.p0))[::-1], source.b.values[::-1]
                        wind = native.values[:, ::-1]
                        values = np.full((native.sizes["time"], native.sizes["lat"], native.sizes["lon"]), np.nan)
                        for month in range(native.sizes["time"]):
                            pressure = offset[:, None, None] + fraction[:, None, None] * surface_pressure[month][None]
                            assert np.all(np.diff(pressure, axis=0) > 0)
                            lower = np.clip(np.sum(pressure <= pressure_Pa, axis=0) - 1, 0, len(offset) - 2)
                            low_p = np.take_along_axis(pressure, lower[None], axis=0)[0]
                            high_p = np.take_along_axis(pressure, (lower + 1)[None], axis=0)[0]
                            low_u = np.take_along_axis(wind[month], lower[None], axis=0)[0]
                            high_u = np.take_along_axis(wind[month], (lower + 1)[None], axis=0)[0]
                            amount = (np.log(pressure_Pa) - np.log(low_p)) / (np.log(high_p) - np.log(low_p))
                            valid = (pressure_Pa >= pressure[0]) & (pressure_Pa <= pressure[-1]) & (pressure_Pa <= surface_pressure[month])
                            values[month] = np.where(valid, low_u + amount * (high_u - low_u), np.nan)
                        field = xr.DataArray(values, dims=("time", "lat", "lon"), coords={"time": native.time, "lat": native.lat, "lon": native.lon})
                    field = field.reset_coords(drop=True)
                    pieces.append(field)
                    input_rows.append({"model": model, "variable": variable, "period": period, "path": str(path),
                                       "raw_dimensions": str(source[variable].dims), "raw_shape": str(source[variable].shape),
                                       "units": source[variable].attrs["units"], "months_selected": int(selection.sum()),
                                       "first_month": str(field.time.values[0]), "last_month": str(field.time.values[-1]),
                                       "tracking_id": source.attrs.get("tracking_id", "verified_daily_derived_monthly"), "file_bytes": path.stat().st_size})
            monthly = xr.concat(pieces, "time", join="exact").sortby("time")
            np.testing.assert_array_equal(monthly.time.dt.year.values * 12 + monthly.time.dt.month.values,
                                          np.arange(start_year * 12 + 1, (start_year + 10) * 12 + 1))
            assert monthly.dims == ("time", "lat", "lon") and monthly.sizes["time"] == 120
            native_hash = hashlib.sha256(np.ascontiguousarray(monthly.values).tobytes()).hexdigest()
            season_results = []
            for season in seasons:
                if season == "Annual":
                    seasonal = monthly.coarsen(time=12, boundary="exact").mean(skipna=False)
                    seasonal = seasonal.rename(time="sample_index").assign_coords(sample_index=sample_indices)
                else:
                    winters = []
                    for winter in range(9):
                        indices = [winter * 12 + 11, (winter + 1) * 12, (winter + 1) * 12 + 1]
                        assert monthly.time.dt.month.values[indices].tolist() == [12, 1, 2]
                        winter_mean = monthly.isel(time=indices).mean("time", skipna=False)
                        winters.append(winter_mean.expand_dims(sample_index=[winter + 1]))
                    seasonal = xr.concat(winters, "sample_index", join="exact")
                seasonal = seasonal.assign_coords(lon=seasonal.lon % 360).sortby("lon")
                if model == "MIROC-ES2L":
                    expected_latitude = latitude[latitude >= -5] if variable == "ua" else latitude
                    np.testing.assert_allclose(seasonal.lat.values, expected_latitude, atol=1e-10, rtol=0)
                    seasonal = seasonal.assign_coords(lat=expected_latitude)
                cyclic = xr.concat([seasonal.isel(lon=[-1]).assign_coords(lon=[float(seasonal.lon[-1]) - 360]),
                                    seasonal, seasonal.isel(lon=[0]).assign_coords(lon=[float(seasonal.lon[0]) + 360])], "lon")
                target_latitude = wind_latitude if variable == "ua" else latitude
                mapped = cyclic.interp(lat=target_latitude, lon=longitude)
                season_results.append(mapped.reindex(sample_index=sample_indices, lat=latitude).expand_dims(season=[season]))
                endpoint_rows.append({"model": model, "variable": variable, "period": period, "season": season,
                                      "calendar_first_year": start_year, "calendar_last_year": start_year + 9,
                                      "monthly_window_count": 120, "sample_count": sample_count[season],
                                      "first_winter_end_year": start_year + 1 if season == "DJF" else None,
                                      "last_winter_end_year": start_year + 9 if season == "DJF" else None,
                                      "native_selected_values_sha256": native_hash})
            phase_results.append(xr.concat(season_results, "season", join="exact").expand_dims(period=[period]))
            print(model, variable, period, "Annual10 / complete DJF9", flush=True)
        combined = xr.concat(phase_results, "period", join="exact")
        supported_seasons = []
        for season in seasons:
            selected = combined.sel(season=season, drop=True).isel(sample_index=slice(0, sample_count[season]))
            paired_support = selected.notnull().all(["period", "sample_index"])
            supported_seasons.append(combined.sel(season=season, drop=True).where(paired_support).expand_dims(season=[season]))
        name = "u250" if variable == "ua" else variable
        variables[name] = xr.concat(supported_seasons, "season", join="exact")
    # P和E使用共同有效样本支持后相减，缺测不填零。
    variables["pr"], variables["evspsbl"] = variables["pr"].where(variables["evspsbl"].notnull()), variables["evspsbl"].where(variables["pr"].notnull())
    variables["pe"] = variables["pr"] - variables["evspsbl"]
    model_results.append(xr.Dataset(variables).expand_dims(model=[model]))


# %% 4. 固定成员组平均与既有年平均逐值复核

samples = xr.concat(model_results, "model", join="exact").transpose("model", "period", "season", "sample_index", "lat", "lon")
climate_seasons = []
for season in seasons:
    climate_seasons.append(samples.sel(season=season, drop=True).isel(sample_index=slice(0, sample_count[season]))
                           .mean("sample_index", skipna=False).expand_dims(season=[season]))
member_climate = xr.concat(climate_seasons, "season", join="exact").transpose("model", "period", "season", "lat", "lon")
group_results = []
for group in groups:
    selected_models = membership.loc[membership.group.eq(group), "model"].tolist()
    group_results.append(member_climate.sel(model=selected_models).mean("model", skipna=False).expand_dims(group=[group]))
group_climate = xr.concat(group_results, "group", join="exact")
member_changes = member_climate.sel(period="later", drop=True) - member_climate.sel(period="initial", drop=True)
group_changes = group_climate.sel(period="later", drop=True) - group_climate.sel(period="initial", drop=True)
annual_errors = []
for folder, variables_to_check in [(annual_directory, ["u250"]), (netp_directory, ["pr", "evspsbl", "pe"])]:
    with xr.open_dataset(folder / "annual_member_samples.nc") as source:
        for variable in variables_to_check:
            old = source[variable].transpose("model", "period", "year_in_window", "lat", "lon").values
            new = samples[variable].sel(season="Annual", lat=wind_latitude).values
            np.testing.assert_allclose(new, old, atol=1e-10, rtol=0, equal_nan=True)
            annual_errors.append({"variable": variable, "max_abs_error": float(np.nanmax(np.abs(new - old)))})
for dataset in [samples, member_climate, group_climate, member_changes, group_changes]:
    np.testing.assert_allclose(dataset.pe, dataset.pr - dataset.evspsbl, atol=1e-12, rtol=0, equal_nan=True)
    for variable in dataset.data_vars:
        dataset[variable].attrs["units"] = "m s-1" if variable == "u250" else "mm day-1"


# %% 5. 保存正式样本与来源维序

destination.mkdir(parents=True)
for name, dataset in [("seasonal_member_samples", samples), ("member_climatology", member_climate),
                      ("group_climatology", group_climate), ("member_changes", member_changes), ("group_changes", group_changes)]:
    encoding = {variable: {"zlib": True, "complevel": 3} for variable in dataset.data_vars}
    dataset.to_netcdf(destination / f"{name}.nc", encoding=encoding)
membership.to_csv(destination / "membership.csv", index=False)
pd.DataFrame(input_rows).to_csv(destination / "input_audit.csv", index=False)
pd.DataFrame(endpoint_rows).to_csv(destination / "endpoint_diagnostics.csv", index=False)
for name in ["coastline_nh_50m.geojson", "nh_land_support.npz"]:
    (destination / name).write_bytes((annual_directory / name).read_bytes())
metadata = {"models": models, "groups": groups, "periods": periods, "seasons": seasons, "sample_count": sample_count,
            "first_calendar_years": first_years, "pressure_Pa": pressure_Pa,
            "time": "Initial1pctCO2years1-10;laterCDRyears131-140. Annual10calendar-year means,12months equally weighted. DJF9complete contiguous winters entirely inside each10year window,Dec ofyears1-9 andJanFeb ofyears2-10,3months equally weighted.",
            "padding": "sample_index10 is NaN in DJF;not a missing scientific winter,never included in9winter means or support checks",
            "native_dimensions": {"P_E": ["time", "lat", "lon"], "UA_seven": ["time", "plev", "lat", "lon"],
                                  "UA_Nor": ["time", "lev", "lat", "lon"], "PS_Nor": ["time", "lat", "lon"], "Nor_coefficients": "a/b(lev),p0 scalar"},
            "saved_dimensions": {"samples": ["model", "period", "season", "sample_index", "lat", "lon"],
                                 "member_climate": ["model", "period", "season", "lat", "lon"],
                                 "group_climate": ["group", "period", "season", "lat", "lon"]},
            "grid": "CanESM5global64Gaussian latitudes x128longitudes;u250available onlyNH32latitudes;SHwind deliberatelyNaN",
            "mapping": "Seasonal native means then periodic bilinear mapping;no extrapolation/no smoothing;MIROC nominal latitudes checked to1e-10degree before roundoff alignment",
            "wind": "250hPa eastward component,not total horizontal speed;Nor monthly actualp=a*p0+b*ps,log-pressure interpolation without extrapolation",
            "hydrology": "P and E multiplied86400 fromkg/m2/s;E upwardpositive;PE=P-E oncommon finite support",
            "grouping": "fixed2/2/4 equalmodel;all group members required;one realization permodel;no newgroup classification",
            "annual_reproduction": annual_errors, "drift_correction": "none", "significance": "not computed",
            "references": ["https://doi.org/10.1038/s41467-021-24015-w"],
            "small_source_hashes": [{"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                                    for path in [annual_directory / "membership.csv", annual_directory / "input_audit.csv",
                                                 netp_directory / "input_audit.csv", annual_directory / "annual_member_samples.nc",
                                                 netp_directory / "annual_member_samples.nc"]]}
(destination / "metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps({"annual_reproduction": annual_errors, "destination": str(destination)}, indent=2), flush=True)
