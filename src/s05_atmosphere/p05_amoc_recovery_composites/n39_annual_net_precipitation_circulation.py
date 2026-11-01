# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import xarray as xr


# %% 2. 参数与路径设置

root = Path("F:/sea/Antic sea ice constraint")
previous = root / "output/s05_atmosphere/p05_amoc_recovery_composites/n37_annual_nh_precipitation_jet"
destination = root / "output/s05_atmosphere/p05_amoc_recovery_composites/n39_annual_net_precipitation_circulation"
assert not destination.exists(), "正式结果已存在，须先归档。"
membership_path = previous / "membership.csv"
membership = pd.read_csv(membership_path)
models = membership.model.tolist()
groups = ["better_nonrecovery", "better_recovery", "larger_bias"]
periods = ["initial", "later"]
experiments = ["1pctCO2", "1pctCO2-cdr"]
members = {"CESM2": "r1i1p1f1", "CNRM-ESM2-1": "r1i1p1f2", "NorESM2-LM": "r1i1p1f1",
           "GFDL-ESM4": "r1i1p1f1", "ACCESS-ESM1-5": "r1i1p1f1", "CanESM5": "r1i1p2f1",
           "MIROC-ES2L": "r1i1p1f2", "UKESM1-0-LL": "r1i1p1f2"}
first_years = {"CESM2": [1, 131], "CNRM-ESM2-1": [1850, 2120], "NorESM2-LM": [1, 271],
               "GFDL-ESM4": [1, 271], "ACCESS-ESM1-5": [101, 371], "CanESM5": [1850, 2121],
               "MIROC-ES2L": [1850, 2119], "UKESM1-0-LL": [1850, 2120]}
target_pressure = 85000.0
input_audit, endpoint_qc, member_results = [], [], []
with xr.open_dataset(previous / "annual_member_samples.nc") as source:
    precipitation = source.pr.load()
assert precipitation.dims == ("model", "period", "year_in_window", "lat", "lon")
latitude, longitude = precipitation.lat.values, precipitation.lon.values
wind_sources = pd.read_csv(previous / "input_audit.csv")


# %% 3. 按真实历年读取蒸发、低层西风和海平面气压

for model in models:
    variables = {"pr": precipitation.sel(model=model, drop=True)}
    for variable in ["evspsbl", "ua", "psl"]:
        phase_fields = []
        for period_index, period in enumerate(periods):
            start_year = first_years[model][period_index]
            if variable == "ua":
                selected_rows = wind_sources.model.eq(model) & wind_sources.variable.eq("ua") & wind_sources.period.eq(period)
                paths = [Path(value) for value in wind_sources.loc[selected_rows, "file"]]
            else:
                suffix = "ad" if period_index == 0 else ""
                folder = Path("F:/sea") / f"cmip6{variable}{suffix}"
                candidates = sorted(folder.rglob(f"{variable}_Amon_{model}_{experiments[period_index]}_{members[model]}_*.nc"))
                paths = []
                for path in candidates:
                    tokens = path.stem.rsplit("_", 1)[1].split("-")
                    if int(tokens[0][:4]) <= start_year + 9 and int(tokens[1][:4]) >= start_year:
                        paths.append(path)
            assert paths, (model, variable, period)
            pieces = []
            for path in paths:
                with xr.open_dataset(path, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as source:
                    assert source.attrs["source_id"] == model
                    assert source.attrs["experiment_id"] == experiments[period_index]
                    assert source.attrs["variant_label"] == members[model]
                    selected = source.time.dt.year.isin(range(start_year, start_year + 10))
                    assert bool(selected.any())
                    native_dimensions = source[variable].dims
                    if variable in ["evspsbl", "psl"]:
                        assert native_dimensions == ("time", "lat", "lon")
                        expected_unit = "kg m-2 s-1" if variable == "evspsbl" else "Pa"
                        assert source[variable].attrs["units"] == expected_unit
                        # E保持向上为正；kg/m²/s转mm/day。PSL只转hPa，不使用PS代替。
                        conversion = 86400.0 if variable == "evspsbl" else 0.01
                        field = source[variable].isel(time=selected).sortby("lat").sel(lat=slice(-5, 90)).load().astype(float) * conversion
                    elif model != "NorESM2-LM":
                        assert native_dimensions == ("time", "plev", "lat", "lon")
                        assert source.plev.attrs["units"] == "Pa" and source.ua.attrs["units"] == "m s-1"
                        field = source.ua.assign_coords(plev=np.rint(source.plev)).sel(plev=target_pressure, drop=True)
                        field = field.isel(time=selected).sortby("lat").sel(lat=slice(-5, 90)).load().astype(float)
                    else:
                        assert native_dimensions == ("time", "lev", "lat", "lon")
                        assert source.ps.dims == ("time", "lat", "lon") and source.ps.attrs["units"] == "Pa"
                        native = source.ua.isel(time=selected).sortby("lat").sel(lat=slice(-5, 90)).load().astype(float)
                        surface_pressure = source.ps.isel(time=selected).sel(lat=native.lat).load().values
                        offset = (source.a.values * float(source.p0))[::-1]
                        fraction = source.b.values[::-1]
                        native_wind = native.values[:, ::-1]
                        converted = np.full((native.sizes["time"], native.sizes["lat"], native.sizes["lon"]), np.nan)
                        for month_index in range(native.sizes["time"]):
                            pressure = offset[:, None, None] + fraction[:, None, None] * surface_pressure[month_index][None]
                            assert np.all(np.diff(pressure, axis=0) > 0)
                            lower = np.clip(np.sum(pressure <= target_pressure, axis=0) - 1, 0, len(offset) - 2)
                            low_p = np.take_along_axis(pressure, lower[None], axis=0)[0]
                            high_p = np.take_along_axis(pressure, (lower + 1)[None], axis=0)[0]
                            low_u = np.take_along_axis(native_wind[month_index], lower[None], axis=0)[0]
                            high_u = np.take_along_axis(native_wind[month_index], (lower + 1)[None], axis=0)[0]
                            amount = (np.log(target_pressure) - np.log(low_p)) / (np.log(high_p) - np.log(low_p))
                            valid = ((target_pressure >= pressure[0]) & (target_pressure <= pressure[-1])
                                     & (target_pressure <= surface_pressure[month_index]))
                            converted[month_index] = np.where(valid, low_u + amount * (high_u - low_u), np.nan)
                        field = xr.DataArray(converted, dims=("time", "lat", "lon"),
                                             coords={"time": native.time, "lat": native.lat, "lon": native.lon})
                    field = field.reset_coords(drop=True)
                    pieces.append(field)
                    input_audit.append({"model": model, "variable": variable, "period": period,
                                        "file": str(path), "raw_dimensions": str(native_dimensions),
                                        "raw_shape": str(source[variable].shape), "units": source[variable].attrs["units"],
                                        "standard_name": source[variable].attrs.get("standard_name", ""),
                                        "months": field.sizes["time"], "first_month": str(field.time.values[0]),
                                        "last_month": str(field.time.values[-1]), "calendar": source.time.encoding["calendar"],
                                        "tracking_id": source.attrs.get("tracking_id", ""), "file_bytes": path.stat().st_size})
            monthly = xr.concat(pieces, "time", join="exact").sortby("time")
            np.testing.assert_array_equal(monthly.time.dt.year.values * 12 + monthly.time.dt.month.values,
                                          np.arange(start_year * 12 + 1, (start_year + 10) * 12 + 1))
            assert monthly.dims == ("time", "lat", "lon") and monthly.sizes["time"] == 120
            annual = monthly.coarsen(time=12, boundary="exact").mean(skipna=False)
            annual = annual.rename(time="year_in_window").assign_coords(year_in_window=np.arange(1, 11))
            annual = annual.assign_coords(lon=annual.lon % 360).sortby("lon")
            if model == "MIROC-ES2L":
                nominal_latitude = np.concatenate([-latitude[:2][::-1], latitude])
                np.testing.assert_allclose(annual.lat, nominal_latitude, atol=1e-10, rtol=0)
                annual = annual.assign_coords(lat=nominal_latitude)
            cyclic = xr.concat([annual.isel(lon=[-1]).assign_coords(lon=[float(annual.lon[-1]) - 360]),
                                annual, annual.isel(lon=[0]).assign_coords(lon=[float(annual.lon[0]) + 360])], "lon")
            mapped = cyclic.interp(lat=latitude, lon=longitude)
            assert mapped.dims == ("year_in_window", "lat", "lon")
            phase_fields.append(mapped.expand_dims(period=[period]))
            endpoint_qc.append({"model": model, "variable": variable, "period": period,
                                "first_year": start_year, "last_year": start_year + 9, "months": 120,
                                "monthly_selected_sha256": hashlib.sha256(np.ascontiguousarray(monthly.values).tobytes()).hexdigest(),
                                "all_year_finite_cells": int(np.isfinite(mapped.values).all(axis=0).sum()),
                                "minimum": float(np.nanmin(mapped)), "maximum": float(np.nanmax(mapped))})
            print(model, variable, period, start_year, "120 months", flush=True)
        values = xr.concat(phase_fields, "period", join="exact")
        support = values.notnull().all(["period", "year_in_window"])
        name = "u850" if variable == "ua" else variable
        variables[name] = values.where(support)
    # P和E逐成员在完全相同的经纬度/时期/年度样本上相减。
    hydro_support = (variables["pr"].notnull() & variables["evspsbl"].notnull()).all(["period", "year_in_window"])
    variables["pr"] = variables["pr"].where(hydro_support)
    variables["evspsbl"] = variables["evspsbl"].where(hydro_support)
    variables["pe"] = variables["pr"] - variables["evspsbl"]
    member_results.append(xr.Dataset(variables).expand_dims(model=[model]))


# %% 4. 模式等权分组并核对净降水恒等式

samples = xr.concat(member_results, "model", join="exact")
units = {"pr": "mm day-1", "evspsbl": "mm day-1", "pe": "mm day-1", "u850": "m s-1", "psl": "hPa"}
for variable, unit in units.items():
    samples[variable].attrs["units"] = unit
climate = samples.mean("year_in_window", skipna=False)
group_climate = xr.concat([climate.sel(model=membership.loc[membership.group.eq(group), "model"].tolist())
                           .mean("model", skipna=False).expand_dims(group=[group]) for group in groups], "group")
changes = group_climate.sel(period="later", drop=True) - group_climate.sel(period="initial", drop=True)
member_changes = climate.sel(period="later", drop=True) - climate.sel(period="initial", drop=True)
for dataset in [samples, climate, group_climate, changes, member_changes]:
    error = float(np.nanmax(np.abs((dataset.pe - dataset.pr + dataset.evspsbl).values)))
    assert error < 1e-10
    for variable, unit in units.items():
        dataset[variable].attrs["units"] = unit
group_diagnostics = []
color_limits = {"pe": [-1.2, 1.2], "u850": [-2.4, 2.4], "psl": [-3.0, 3.0]}
for group in groups:
    for variable, limits in color_limits.items():
        values = changes[variable].sel(group=group).values
        group_diagnostics.append({"group": group, "variable": variable, "units": units[variable],
                                  "minimum": float(np.nanmin(values)), "maximum": float(np.nanmax(values)),
                                  "below_color_limit_cells": int(np.sum(values < limits[0])),
                                  "above_color_limit_cells": int(np.sum(values > limits[1])),
                                  "missing_cells": int(np.isnan(values).sum())})


# %% 5. 保存正式结果及实际维序

destination.mkdir(parents=True)
samples.to_netcdf(destination / "annual_member_samples.nc")
climate.to_netcdf(destination / "member_climatology.nc")
member_changes.to_netcdf(destination / "member_changes.nc")
group_climate.to_netcdf(destination / "group_climatology.nc")
changes.to_netcdf(destination / "group_changes.nc")
membership.to_csv(destination / "membership.csv", index=False)
pd.DataFrame(input_audit).to_csv(destination / "input_audit.csv", index=False)
pd.DataFrame(endpoint_qc).to_csv(destination / "endpoint_diagnostics.csv", index=False)
pd.DataFrame(group_diagnostics).to_csv(destination / "group_diagnostics.csv", index=False)
coast_path = previous / "coastline_nh_50m.geojson"
(destination / "coastline_nh_50m.geojson").write_bytes(coast_path.read_bytes())
metadata = {"models": models, "groups": groups, "periods": periods, "first_calendar_years": first_years,
            "time": "Annual,12months equal weight and10years equal weight;1pctCO2years1-10 vsCDRyears131-140 (nominal271-280)",
            "net_precipitation": "pe=pr-evspsbl;E upward positive,86400 converts kg/m2/s to mm/day;land and ocean retained",
            "wind": "850hPa eastward wind,u positive eastward;Nor monthly actual pressure=a*p0+b*ps,log-pressure interpolation without extrapolation",
            "pressure": "psl is sea-level pressure,Pa/100 to hPa;not surface pressure ps",
            "grouping": "fixed2/2/4,first each model endpoint then equal-model average;all fixed members required",
            "saved_dimension_order": {"samples": ["model", "period", "year_in_window", "lat", "lon"],
                                      "climate": ["group", "period", "lat", "lon"], "changes": ["group", "lat", "lon"]},
            "raw_dimension_order": "E/PSL(time,lat,lon);7modelsUA(time,plev,lat,lon);NorUA(time,lev,lat,lon),PS(time,lat,lon),a/b(lev),p0scalar",
            "horizontal_mapping": "CanGaussian32NHlatx128lon,periodicbilinear,annualnative means before remap;noextrapolation;MIROCnominalcoords aligned within1e-10degree",
            "terrain": "Native supplied missing wind retained;Nor additionally requires850hPa at or above surface;no uniform extra ps-based terrain masking for the seven fixed-pressure products",
            "geographic_extent": [-180, 180, 0, 90], "color_limits": color_limits,
            "climatology_contours": {"u850": [5, 10, 15], "psl": [1000, 1005, 1010, 1015, 1020, 1025]},
            "u850_change_contours": [-2, -1, -0.5, 0.5, 1, 2],
            "reference": "Yu et al.2025,GRL,10.1029/2024GL113541,Fig2:green climatology at5m/s and5hPa intervals;ourannual three-group/netP context differs from original winter paper",
            "drift_correction": "none", "smoothing": "none", "significance": "not computed;no borrowed hatching or dots",
            "interpretation": "Contours are wind-speed component isolines,not tracked jet axes;P-E increase is not necessarily P increase;association is not causality",
            "small_sources": [{"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                              for path in [membership_path, previous / "annual_member_samples.nc", previous / "input_audit.csv", coast_path]]}
(destination / "metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
print(pd.DataFrame(group_diagnostics).to_string(index=False), flush=True)
