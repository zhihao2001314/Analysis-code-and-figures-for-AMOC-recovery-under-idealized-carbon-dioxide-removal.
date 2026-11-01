# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm
import cartopy.crs as ccrs
from shapely.geometry import shape


# %% 2. 参数与最终图件路径

root = Path(__file__).resolve().parents[2]
source = root / 'output/s03_tem/n19_arctic_partition_tas_ice'
destination = root / 'figure/s03_tem/n20_arctic_partition_tas_ice'
png_path = destination / 'arctic_partition_temperature_ice_nine_models.png'
pdf_path = destination / 'arctic_partition_temperature_ice_nine_models.pdf'
svg_path = destination / 'arctic_partition_temperature_ice_nine_models.svg'
figure_paths = [png_path, pdf_path, svg_path]
assert not destination.exists(), '正式图件存在，须先归档。'
figure_size, dpi = (10.4, 11.2), 450
region_colors = {'CA': '#526E89', 'NA': '#C47D20', 'NP': '#317A65'}
edge_colors, edge_styles = ['#7B3DBA', '#27985B', '#E000C8'], ['--', '--', '-']
region_names = {'CA': 'Central Arctic', 'NA': 'North Atlantic', 'NP': 'North Pacific'}
records = []


# %% 3. 仅读取正式计算产品

data = np.load(source / 'arctic_partition_tas_ice.npz')
metadata = json.loads((source / 'metadata.json').read_text(encoding='utf-8'))
summary = pd.read_csv(source / 'regional_tas_changes.csv').set_index(['model', 'region'])
obs = np.load(source / 'observed_sic_1988_2020_climatology_edge.npz')
coast = [shape(item['geometry']) for item in json.loads((source / 'coastline.geojson').read_text(encoding='utf-8'))['features']]
longitude = np.append(data['longitude'], data['longitude'][0] + 360.)
angle = np.linspace(0., 2. * np.pi, 361)
circle = mpath.Path(np.column_stack([.5 + .5 * np.sin(angle), .5 + .5 * np.cos(angle)]))
projection, geographic = ccrs.NorthPolarStereo(central_longitude=0.), ccrs.PlateCarree()
observation_projection = ccrs.epsg(3411)


# %% 4. 九模式温度分区与海冰边缘

with plt.rc_context({'font.family': 'Arial', 'font.size': 8, 'pdf.fonttype': 42,
                     'svg.fonttype': 'none', 'axes.linewidth': .5}):
    figure = plt.figure(figsize=figure_size)
    figure.subplots_adjust(left=.035, right=.965, top=.965, bottom=.15, wspace=.16, hspace=.20)
    for index, model in enumerate(data['models']):
        axis = figure.add_subplot(3, 3, index + 1, projection=projection)
        axis.set_extent([-180., 180., 45., 90.], crs=geographic)
        axis.set_boundary(circle, transform=axis.transAxes)
        field = np.concatenate([data['tas_change'][index], data['tas_change'][index, :, :1]], axis=1)
        image = axis.pcolormesh(longitude, data['latitude'], field, transform=geographic,
            cmap='RdBu_r', norm=TwoSlopeNorm(vmin=-8., vcenter=0., vmax=8.), rasterized=True, shading='nearest')
        axis.add_geometries(coast, geographic, facecolor='none', edgecolor='#6A7175', linewidth=.35, zorder=3)
        axis.gridlines(xlocs=[-180, -120, -60, 0, 60, 120], ylocs=[60, 70, 80],
                       linewidth=.3, color='#8C9498', alpha=.4, linestyle=':')
        segments = []
        for key, color, style in zip(['initial_sic', 'final_sic'], edge_colors[:2], edge_styles[:2]):
            ice = np.concatenate([data[key][index], data[key][index, :, :1]], axis=1)
            contour = axis.contour(longitude, data['latitude'], ice, levels=[.15], colors=[color],
                linestyles=style, linewidths=.7, transform=geographic, zorder=4)
            segments.append(sum(len(part) > 1 for part in contour.allsegs[0]))
        contour = axis.contour(obs['projected_x_m'], obs['projected_y_m'], obs['sic_smoothed_for_edge_fraction'],
            levels=[.15], colors=[edge_colors[2]], linewidths=.7, transform=observation_projection, zorder=4)
        segments.append(sum(len(part) > 1 for part in contour.allsegs[0]))
        for region, bounds in metadata['regions'].items():
            south, north, west, east = bounds
            pieces = [(west, east)] if west < east else [(west, 180.), (-180., east)]
            for left, right in pieces:
                x = np.linspace(left, right, 181)
                axis.plot(x, np.full_like(x, south), color=region_colors[region], lw=.85, transform=geographic, zorder=5)
                if region != 'CA':
                    axis.plot(x, np.full_like(x, north), color=region_colors[region], lw=.85, transform=geographic, zorder=5)
            if region != 'CA':
                y = np.linspace(south, north, 81)
                for x in [west, east]:
                    axis.plot(np.full_like(y, x), y, color=region_colors[region], lw=.85, transform=geographic, zorder=5)
        for region, xy in [('CA', (.50, .57)), ('NA', (.47, .20)), ('NP', (.57, .87))]:
            value = summary.loc[(model, region_names[region]), 'raw_delta_tas_degC']
            text = axis.text(*xy, f'{region} {value:+.2f}°C', transform=axis.transAxes,
                ha='center', va='center', fontsize=6.9, color=region_colors[region], zorder=7)
            text.set_path_effects([pe.withStroke(linewidth=2., foreground='white')])
        axis.set_title(f'({chr(97+index)}) {model}', loc='left', fontsize=9, pad=11)
        for x, y, text in [(.5, -.028, '0°'), (.5, 1.025, '180°'), (-.035, .5, '90°W'), (1.035, .5, '90°E')]:
            axis.text(x, y, text, transform=axis.transAxes, ha='center', va='center', fontsize=6.5, color='#68747B')
        records.append({'model': str(model), 'nonempty_ice_segments': segments,
            'field_sha256': hashlib.sha256(data['tas_change'][index].tobytes()).hexdigest()})
    color_axis = figure.add_axes([.24, .095, .52, .014])
    bar = figure.colorbar(image, cax=color_axis, orientation='horizontal', ticks=[-8, -4, 0, 4, 8], extend='both')
    bar.set_label('Near-surface air-temperature change (°C)', fontsize=9)
    handles = [Line2D([], [], color=color, ls=style, lw=1., label=label) for color, style, label in
        zip(edge_colors, edge_styles, ['Initial ice edge', 'Final ice edge', 'Observed ice edge (1988–2020)'])]
    figure.legend(handles=handles, loc='center', bbox_to_anchor=(.5, .032), ncol=3, frameon=False, fontsize=8)
    destination.mkdir(parents=True)
    for path in figure_paths:
        figure.savefig(path, dpi=dpi, facecolor='white')
        print(path.resolve(), flush=True)
    plt.close(figure)


# %% 5. 保存图注和来源

(destination / 'caption.txt').write_text(metadata['caption'], encoding='utf-8')
manifest = {'models': metadata['models'], 'panel_count': 9, 'figure_size_inches': figure_size, 'dpi': dpi,
    'color_limits_degC': [-8, 8], 'records': records, 'regions': metadata['regions'],
    'caption': metadata['caption'], 'source_files': [{'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
        for p in source.iterdir() if p.is_file()]}
(destination / 'figure_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
