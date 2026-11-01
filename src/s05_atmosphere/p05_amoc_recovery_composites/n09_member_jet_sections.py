# %% 1. 导入依赖

from pathlib import Path
import json
import numpy as np
import pandas as pd
import xarray as xr


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
base = root / 'output/s05_atmosphere/p05_amoc_recovery_composites'
jet_base = root / 'output/s05_atmosphere/p02_ua'
output = base / 'n09_member_jet_sections'
assert not output.exists(), '正式结果存在，需先归档。'
models = ['GFDL-ESM4', 'NorESM2-LM', 'CESM2', 'CNRM-ESM2-1']
levels = [85000, 70000, 50000]


# %% 3. 读取相同窗口的12层剖面与三层空间场

with xr.open_dataset(jet_base / 'n27_endpoint10_jet_sections/jet_temperature_sections.nc') as source:
    sections = source.sel(model=models, region='North Atlantic-Europe', drop=True).load()
    section_method = dict(source.attrs)
with xr.open_dataset(jet_base / 'n23_eight_model_pressure_thermal/pressure_matched_patterns.nc') as source:
    pressure_fields = source[['ua', 'gradient']].sel(model=models, plev=levels).load()
with xr.open_dataset(base / 'n07_member_consistency/member_changes.nc') as source:
    changes = source.load()
jet = pd.read_csv(base / 'n07_member_consistency/member_jet_centroid.csv')
jet = jet[jet.region == 'North Atlantic-Europe'].copy()
assert sections.wind_profile.dims == ('model', 'period', 'plev', 'lat')
assert sections.gradient_profile.dims == ('model', 'period', 'plev', 'lat')
assert pressure_fields.ua.dims == ('model', 'period', 'plev', 'lat', 'lon')
assert pressure_fields.gradient.dims == ('model', 'period', 'plev', 'lat', 'lon')
assert list(sections.period.values) == ['early', 'later']
assert list(pressure_fields.period.values) == ['early', 'later']
assert np.array_equal(sections.plev.values, [100000, 92500, 85000, 70000, 60000, 50000, 40000, 30000, 25000, 20000, 15000, 10000])


# %% 4. 独立复核差值，分开保存风场与温度梯度变化

wind_error = float(abs(sections.wind_profile.sel(period='later') - sections.wind_profile.sel(period='early') - sections.delta_wind_profile).max())
gradient_error = float(abs(sections.gradient_profile.sel(period='later') - sections.gradient_profile.sel(period='early') - sections.delta_gradient_profile).max())
assert max(wind_error, gradient_error) < 1.e-10
pressure_changes = pressure_fields.sel(period='later', drop=True) - pressure_fields.sel(period='early', drop=True)
pressure_changes = pressure_changes.rename({'ua': 'delta_u', 'gradient': 'delta_gradient'})
u850_error = float(abs(pressure_changes.delta_u.sel(plev=85000, drop=True) - changes.u850).max())
assert u850_error < 1.e-8
# 保留各自既有支持；n01温度梯度与风的联合掩膜可能比n23原梯度更严格。
gradient850_error = float(abs(pressure_changes.delta_gradient.sel(plev=85000, drop=True) - changes.gradient850).max())
assert gradient850_error < 1.e-10
assert float(abs(sections.upper_wind - sections.lower_wind - sections.actual_shear).max()) < 1.e-10


# %% 5. 保存绘图输入和维序口径，不重新读取原始数据

output.mkdir(parents=True)
sections.to_netcdf(output / 'na_member_sections.nc')
pressure_changes.to_netcdf(output / 'pressure_change_maps.nc')
changes.to_netcdf(output / 'member_surface_changes.nc')
jet.to_csv(output / 'na_jet_centroids.csv', index=False)
report = {'models': models, 'groups': {'recovered': models[:2], 'not_recovered': models[2:]},
          'period': 'laterCDR131-140 minus early1pctCO2 1-10; continuous271-280 minus1-10; 30DJF calendar months/dayweighted',
          'sections_source': str(jet_base / 'n27_endpoint10_jet_sections/jet_temperature_sections.nc'),
          'maps_source': str(jet_base / 'n23_eight_model_pressure_thermal/pressure_matched_patterns.nc'),
          'source_section_method': section_method,
          'dimensions': {'wind_profile_and_gradient_profile': ['model', 'period', 'plev', 'lat'], 'section_changes': ['model', 'plev', 'lat'], 'pressure_maps': ['model', 'plev', 'lat', 'lon'], 'surface_maps': ['model', 'lat', 'lon']},
          'pressure_units': 'Pa in output; hPa only in plotting coordinates',
          'gradient': '-dT/dy in K per1000km; positive change means stronger poleward cooling at the same pressure level; not surface-temperature gradient',
          'section_support': 'Existing n27: each model uses common temperature-gradient/wind support across both periods, >=80% regional longitudes per level/latitude; no filling. Longitude support can differ by model; no group section mean.',
          'mechanism_limit': 'Temperature gradient constrains vertical wind shear, not absolute jet latitude; sections and same-level maps are descriptive checks, not causal attribution.',
          'errors': {'wind_section_identity': wind_error, 'gradient_section_identity': gradient_error, 'u850_vs_existing': u850_error, 'gradient850_on_overlap_vs_existing': gradient850_error},
          'selection': 'Fixed current four models; UK omitted at prior user request; no new selection, smoothing or infilling.'}
(output / 'qc.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(report['errors'], indent=2), flush=True)
