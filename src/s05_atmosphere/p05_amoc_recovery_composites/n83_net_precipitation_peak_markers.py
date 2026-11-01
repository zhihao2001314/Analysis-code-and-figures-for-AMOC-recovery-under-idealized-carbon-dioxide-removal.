# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import xarray as xr


# %% 2. 参数与路径设置

root = Path('F:/sea/Antic sea ice constraint')
base = root / 'output/s05_atmosphere/p05_amoc_recovery_composites'
source = base / 'n81_six_model_net_precipitation/group_climatology.nc'
destination = base / 'n83_net_precipitation_peak_markers'
products, periods = ['recovered', 'not_recovered'], ['early', 'later']
tropical_bounds, marker_spacing_deg = [-20., 20.], 10.
peak_rows, shift_rows = [], []
assert not destination.exists(), '正式结果存在，须先归档。'


# %% 3. 从初末期绝对净降水场提取热带峰值

with xr.open_dataset(source) as dataset:
    climate = dataset.pe.sel(product=products, period=periods).load()
assert climate.dims == ('product', 'period', 'lat', 'lon')
for product in products:
    for period in periods:
        field = climate.sel(product=product, period=period, lat=slice(*tropical_bounds))
        assert field.dims == ('lat', 'lon') and np.isfinite(field).all()
        maxima = field.values.argmax(axis=0)
        for longitude_index, longitude in enumerate(field.lon.values):
            column = field.values[:, longitude_index]
            latitude_index = maxima[longitude_index]
            peak_rows.append({'product': product, 'period': period, 'lon': float(longitude),
                              'peak_lat': float(field.lat.values[latitude_index]),
                              'peak_pe_mm_day': float(column[latitude_index]),
                              'tied_maxima': int(np.sum(column == column[latitude_index])),
                              'peak_at_search_boundary': latitude_index in [0, len(field.lat) - 1],
                              'plotted': bool(np.isclose(np.mod(longitude, marker_spacing_deg), 0.))})
peaks = pd.DataFrame(peak_rows)
for product in products:
    early = peaks.loc[peaks['product'].eq(product) & peaks.period.eq('early')].set_index('lon')
    later = peaks.loc[peaks['product'].eq(product) & peaks.period.eq('later')].set_index('lon')
    for longitude in early.index:
        shift_rows.append({'product': product, 'lon': float(longitude),
                           'initial_peak_lat': float(early.loc[longitude, 'peak_lat']),
                           'later_peak_lat': float(later.loc[longitude, 'peak_lat']),
                           'shift_deg_north': float(later.loc[longitude, 'peak_lat'] - early.loc[longitude, 'peak_lat']),
                           'plotted': bool(early.loc[longitude, 'plotted'])})


# %% 4. 保存结果、维序和标记定义

destination.mkdir(parents=True)
peaks.to_csv(destination / 'net_precipitation_peak_gridpoints.csv', index=False)
pd.DataFrame(shift_rows).to_csv(destination / 'net_precipitation_peak_shifts.csv', index=False)
metadata = {'status': 'passed', 'groups': products, 'periods': periods,
            'input_dimensions': 'pe(product,period,lat,lon), 2-degree grid; saved by n81',
            'peak_definition': 'Latitude of the maximum absolute group-mean P-E climatology at each longitude within20S-20N; not maximum P-E change.',
            'order': 'First calculate each model P-E climatology, equal-weight group mean, then locate the group-mean peak.',
            'marker_spacing_deg': marker_spacing_deg, 'bounds_degN': tropical_bounds,
            'ties': 'First (southernmost) grid latitude if exact ties; tie counts saved.',
            'initial_final': 'Experiment years1-10 and271-280; early maps to Initial and later maps to Later.',
            'direction': 'Positive peak shift is northward; negative is southward.',
            'interpretation': 'P-E maximum is a net precipitation rainband indicator, not identical to precipitation-maximum or median ITCZ index; evaporation and competing rainbands can alter it.',
            'peak_units': 'Latitude degrees north; peak value mm day-1.',
            'display': 'b/d only: initial blue open circles, later magenta crosses; f has no absolute-climatology peak markers.',
            'source_files': [{'path': str(source), 'sha256': hashlib.sha256(source.read_bytes()).hexdigest()}]}
(destination / 'metadata.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
print(peaks.groupby(['product', 'period'])[['plotted', 'peak_at_search_boundary']].sum().to_string(), flush=True)
print(destination.resolve(), flush=True)
