# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import numpy as np
import xarray as xr


# %% 2. 参数与路径设置

root = Path("F:/sea/Antic sea ice constraint")
source_path = Path("F:/sea/seaice-obs-north-monthly/sic_psn25_197811-202607_v06r00.nc")
destination = root / "output/s05_atmosphere/p05_amoc_recovery_composites/n41_observed_annual_ice_edge"
assert not destination.exists(), "正式观测结果已存在，须先归档。"
first_year, last_year = 1988, 2020
years = np.arange(first_year, last_year + 1)
minimum_valid_years = int(np.ceil(len(years) * 0.8))
threshold = 0.15


# %% 3. 读取原生投影网格及396个唯一月份

with xr.open_dataset(source_path, decode_times=True, mask_and_scale=True) as source:
    assert source.cdr_seaice_conc_monthly.dims == ("time", "y", "x")
    assert source.cdr_seaice_conc_monthly.attrs["units"] == "1"
    assert source.x.dims == ("x",) and source.y.dims == ("y",)
    assert source.x.attrs["units"] == source.y.attrs["units"] == "meters"
    assert source.attrs["geospatial_bounds_crs"] == "EPSG:3411"
    selected = source.time.dt.year.isin(years)
    monthly = source.cdr_seaice_conc_monthly.isel(time=selected).load().astype(float)
    month_ids = monthly.time.dt.year.values * 12 + monthly.time.dt.month.values
    np.testing.assert_array_equal(month_ids, np.arange(first_year * 12 + 1, (last_year + 1) * 12 + 1))
    assert monthly.shape == (396, 448, 304)
    x, y = source.x.values.copy(), source.y.values.copy()
    product_title = source.attrs["title"]
    product_version = source.attrs["product_version"]
    product_doi = source.attrs["id"]
    projection_attributes = dict(source.crs.attrs)


# %% 4. 排除特殊编码，逐年12个月等权再逐年等权

# CF解码已执行scale_factor；特殊状态可能仍>1，不能当浓度或填为0。
monthly_values = monthly.values
valid = np.isfinite(monthly_values) & (monthly_values >= 0) & (monthly_values <= 1)
values_by_year = np.where(valid, monthly_values, np.nan).reshape(33, 12, 448, 304)
valid_months = valid.reshape(33, 12, 448, 304).sum(axis=1)
annual = np.sum(np.where(np.isfinite(values_by_year), values_by_year, 0), axis=1) / 12
annual[valid_months != 12] = np.nan
valid_years = np.isfinite(annual).sum(axis=0)
climatology = np.full((448, 304), np.nan)
accepted = valid_years >= minimum_valid_years
climatology[accepted] = np.nansum(annual, axis=0)[accepted] / valid_years[accepted]
assert np.nanmin(climatology) >= 0 and np.nanmax(climatology) <= 1
assert np.nanmin(climatology) < threshold < np.nanmax(climatology)
projected_x, projected_y = np.meshgrid(x, y)
assert np.all(np.diff(x) == 25000) and np.all(np.diff(y) == -25000)


# %% 5. 保存可直接绘图的观测结果与数据口径

destination.mkdir(parents=True)
np.savez_compressed(destination / "observed_annual_sic_1988_2020.npz", projected_x_m=projected_x,
                    projected_y_m=projected_y, sic_climatology_fraction=climatology,
                    valid_year_count=valid_years.astype(np.int16), valid_month_count=valid_months.astype(np.int8),
                    years=years, threshold_fraction=np.array(threshold))
metadata = {"product": product_title, "version": product_version, "doi": product_doi,
            "source_path": str(source_path), "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
            "source_file_bytes": source_path.stat().st_size, "raw_dimension_order": ["time", "y", "x"],
            "selected_shape": [396, 448, 304], "years": [first_year, last_year], "months": 396,
            "season": "Annual,12months equally weighted per year,then accepted years equally weighted",
            "source_units": "1,CF mask-and-scale decoded;do not divide by100 again",
            "validity": "Only finite0<=SIC<=1;12valid months required per year;at least27of33years required per cell;land/specialflags/missing not zero",
            "minimum_valid_years": minimum_valid_years, "threshold_fraction": threshold,
            "projection": "EPSG:3411", "grid_spacing_m": 25000,
            "projection_a_m": float(projection_attributes["semi_major_axis"]),
            "projection_inverse_flattening": float(projection_attributes["inverse_flattening"]),
            "grid_order": "projected_x_m,projected_y_m,sic_climatology_fraction(y448,x304);valid_month_count(year33,y448,x304)",
            "smoothing": "none", "horizontal_regridding": "none;nativeprojectedcontour transformed byCartopy",
            "valid_climatology_cells": int(accepted.sum()), "valid_years_minimum_on_support": int(valid_years[accepted].min()),
            "range_fraction": [float(np.nanmin(climatology)), float(np.nanmax(climatology))],
            "not_reused_legacy_cache": "Earlierannualedgecache used183samples/year from timestamp-freebinary files and sigma0.8display smoothing;not used for this new overlay",
            "interpretation": "1988-2020modernobservationalreference,not contemporaneous validation ofpreindustrialinitialoridealizedCDRlater;edgeofannualmeanSIC not meanofmonthlyedgepositions"}
(destination / "metadata.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps(metadata, ensure_ascii=False, indent=2))
