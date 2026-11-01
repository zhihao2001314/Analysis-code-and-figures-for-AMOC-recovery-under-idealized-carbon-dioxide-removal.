# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import numpy as np
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
from shapely.geometry import shape


# %% 2. 参数与最终图件路径

root = Path(__file__).resolve().parents[3]
source = root / 'output/s05_atmosphere/p05_amoc_recovery_composites/n86_three_by_three_hydrology'
destination = root / 'figure/s05_atmosphere/p05_amoc_recovery_composites/n87_three_by_three_hydrology'
png_path = destination / 'three_by_three_precipitation_evaporation_net_precipitation.png'
pdf_path = destination / 'three_by_three_precipitation_evaporation_net_precipitation.pdf'
svg_path = destination / 'three_by_three_precipitation_evaporation_net_precipitation.svg'
figure_paths = [png_path, pdf_path, svg_path]
assert not destination.exists(), '正式图件存在，须先归档。'
figure_size, dpi = (12.6, 8.8), 450
column_labels = [r'$\Delta P$', r'$\Delta E$', r'$\Delta(P-E)$']
row_labels = ['Recovery (n = 3)', 'Non-recovery (n = 3)', 'Non-recovery − Recovery']
records = []


# %% 3. 仅读取三对三分量及各自检验

with xr.open_dataset(source / 'hydrology_display.nc') as dataset:
    fields = dataset.load()
metadata = json.loads((source / 'metadata.json').read_text(encoding='utf-8'))
coast = [shape(item['geometry']) for item in json.loads((source / 'coastline.geojson').read_text(encoding='utf-8'))['features']]
lon_grid, lat_grid = np.meshgrid(fields.lon, fields.lat)


# %% 4. 三行三列全球地图

with plt.rc_context({'font.family': 'Arial', 'font.size': 8, 'pdf.fonttype': 42, 'svg.fonttype': 'none',
                     'hatch.color': '#777777', 'hatch.linewidth': .2}):
    figure = plt.figure(figsize=figure_size)
    grid = figure.add_gridspec(3, 3, left=.065, right=.985, bottom=.12, top=.945, wspace=.16, hspace=.26)
    for row, product in enumerate(metadata['products']):
        for column, variable in enumerate(metadata['variables']):
            axis = figure.add_subplot(grid[row, column], projection=ccrs.Robinson(central_longitude=180.))
            axis.set_global()
            selected = fields.sel(variable=variable, product=product)
            cyclic, longitude = add_cyclic_point(selected.delta.values, coord=fields.lon.values)
            image = axis.pcolormesh(longitude, fields.lat, cyclic, shading='auto', cmap='BrBG',
                norm=Normalize(-3., 3.), transform=ccrs.PlateCarree(), rasterized=True)
            failed, mask_longitude = add_cyclic_point((~selected.passed.values).astype(float), coord=fields.lon.values)
            axis.contourf(mask_longitude, fields.lat, failed, levels=[.5, 1.5], colors='none', hatches=['xx'],
                transform=ccrs.PlateCarree(), zorder=3)
            mask = selected.passed.values.copy()
            mask[np.arange(mask.shape[0]) % 3 != 0] = False
            mask[:, np.arange(mask.shape[1]) % 3 != 0] = False
            axis.scatter(lon_grid[mask], lat_grid[mask], s=.65, marker='.', color='#333333', linewidths=0,
                alpha=.7, transform=ccrs.PlateCarree(), zorder=4)
            axis.add_geometries(coast, ccrs.PlateCarree(), facecolor='none', edgecolor='#6B7073', linewidth=.35, zorder=5)
            ticks = axis.gridlines(xlocs=[60, 120, 180, -120, -60], ylocs=[-60, -30, 0, 30, 60],
                draw_labels=True, x_inline=False, y_inline=False, rotate_labels=False,
                xlabel_style={'size': 6}, ylabel_style={'size': 6}, linewidth=0, xpadding=2, ypadding=2)
            ticks.top_labels, ticks.right_labels, ticks.geo_labels = False, False, False
            ticks.left_labels, ticks.bottom_labels = column == 0, True
            axis.set_title(chr(97 + row * 3 + column), loc='left', weight='bold', fontsize=10, pad=7)
            if row == 0:
                axis.set_title(column_labels[column], fontsize=11, pad=9)
            if column == 0:
                axis.text(-.16, .5, row_labels[row], rotation=90, ha='center', va='center',
                    transform=axis.transAxes, fontsize=9)
            records.append({'variable': variable, 'product': product,
                'field_sha256': hashlib.sha256(selected.delta.values.tobytes()).hexdigest(),
                'passed_sha256': hashlib.sha256(selected.passed.values.tobytes()).hexdigest()})
    cax = figure.add_axes([.28, .055, .47, .018])
    bar = figure.colorbar(image, cax=cax, orientation='horizontal', extend='both', ticks=[-3, -1.5, 0, 1.5, 3])
    bar.set_label('Change (mm day$^{-1}$)', fontsize=9)
    destination.mkdir(parents=True)
    for path in figure_paths:
        figure.savefig(path, dpi=dpi, facecolor='white')
        print(path.resolve(), flush=True)
    plt.close(figure)


# %% 5. 保存图注、真实显示值与来源

(destination / 'caption.txt').write_text(metadata['caption'], encoding='utf-8')
manifest = {'panel_count': 9, 'figure_size_inches': figure_size, 'dpi': dpi, 'records': records,
    'caption': metadata['caption'], 'source_files': [{'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
        for p in source.iterdir() if p.is_file()]}
(destination / 'figure_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
