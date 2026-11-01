# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
from matplotlib.ticker import FixedLocator
import cartopy.crs as ccrs
from shapely.geometry import shape


# %% 2. 参数与最终图件完整路径

root = Path(__file__).resolve().parents[3]
source = root / 'output/s02_sia/s02_arctic_region/p24_mld_method_comparison/n095_mld_method_difference_patterns'
figure2_source = root / 'output/s02_sia/s02_arctic_region/p23_stratification_mld_composite/n093_compact_regional_hydrography'
destination = root / 'figure/s02_sia/s02_arctic_region/p24_mld_method_comparison/n096_mld_method_difference_patterns'
quantities = ['initial', 'late', 'change']
figure_paths = {quantity: [destination / ('march_mld_' + quantity + '_derived_minus_mlotst.' + suffix)
                           for suffix in ['png', 'pdf', 'svg']] for quantity in quantities}
colorbar_labels = {'initial': 'Initial March MLD difference: derived − mlotst (m)',
                   'late': 'Final March MLD difference: derived − mlotst (m)',
                   'change': 'March MLD response difference: derived − mlotst (m)'}
figure_size, dpi = (15.2, 8.0), 400
palette_anchors = ['#145D9E', '#347FB5', '#70A6CA', '#B6D0DF', '#F5F3EF',
                   '#ECC8B4', '#DF936E', '#CE5144', '#B60026']
projection = ccrs.LambertConformal(central_longitude=-25, central_latitude=65, standard_parallels=(55, 75))
geographic = ccrs.PlateCarree()
exports = []
assert not destination.exists(), '正式图件已存在，须先归档。'


# %% 3. 仅读取已保存差值与底图

metadata = json.loads((source / 'metadata.json').read_text(encoding='utf-8'))
with np.load(source / 'mld_method_difference_patterns.npz') as saved:
    models = saved['models'].tolist()
    np.testing.assert_array_equal(saved['quantities'], quantities)
    difference = saved['difference'].copy()
    longitude, latitude = saved['longitude'].copy(), saved['latitude'].copy()
with np.load(figure2_source / 'grouped_mld_sic_fields.npz') as saved:
    # 完全沿用Figure2的色标与分级；实际方法差值不变。
    levels, ticks = saved['color_levels_m'].copy(), saved['color_ticks_m'].copy()
assert difference.shape == (7, 3, 141, 341)
land = [shape(feature['geometry']) for feature in json.loads((source / 'land_50m.geojson').read_text(encoding='utf-8'))['features']]
coast = [shape(feature['geometry']) for feature in json.loads((source / 'coastline_50m.geojson').read_text(encoding='utf-8'))['features']]


# %% 4. 七模式两方法差值图，统一色标

with plt.rc_context({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.linewidth': .6,
                     'pdf.fonttype': 42, 'svg.fonttype': 'none', 'savefig.facecolor': 'white'}):
    colormap = LinearSegmentedColormap.from_list('mld_residual', palette_anchors, N=256)
    colormap.set_bad('#EDF1F5')
    normalization = BoundaryNorm(levels, colormap.N)
    for quantity_index, quantity in enumerate(quantities):
        figure = plt.figure(figsize=figure_size, facecolor='white')
        grid = figure.add_gridspec(2, 4, left=.04, right=.99, bottom=.15, top=.95, wspace=.075, hspace=.17)
        records = []
        for model_index, model in enumerate(models):
            row, column = divmod(model_index, 4)
            axis = figure.add_subplot(grid[row, column], projection=projection)
            axis.set_extent(metadata['extent'], crs=geographic)
            axis.set_facecolor('#E2E5E7')
            axis.set_rasterization_zorder(2)
            field = difference[model_index, quantity_index]
            filled = axis.contourf(longitude, latitude, np.ma.masked_invalid(field), levels=levels,
                                   cmap=colormap, norm=normalization, extend='both', transform=geographic, antialiased=False, zorder=1)
            axis.add_geometries(land, geographic, facecolor='#E2E5E7', edgecolor='none', zorder=3)
            axis.add_geometries(coast, geographic, facecolor='none', edgecolor='#656565', linewidth=.4, zorder=4)
            lines = axis.gridlines(crs=geographic, draw_labels=True, x_inline=False, y_inline=False,
                                   linewidth=.3, color='#A5AFB5', alpha=.5, linestyle=':')
            lines.xlocator = FixedLocator([-60, -40, -20, 0, 20])
            lines.ylocator = FixedLocator([50, 60, 70, 80])
            lines.top_labels, lines.right_labels = False, False
            lines.left_labels, lines.bottom_labels = column == 0, True
            lines.rotate_labels = False
            lines.xlabel_style, lines.ylabel_style = {'size': 7}, {'size': 7}
            axis.spines['geo'].set_edgecolor('#777777')
            axis.spines['geo'].set_linewidth(.5)
            axis.set_title(chr(97 + model_index) + '  ' + model, loc='left', fontsize=10, weight='bold', pad=7)
            records.append({'panel': chr(97 + model_index), 'model': model,
                            'field_sha256': hashlib.sha256(np.asarray(field, dtype=float).tobytes()).hexdigest()})
        legend_axis = figure.add_subplot(grid[1, 3])
        legend_axis.axis('off')
        legend_axis.text(.5, .63, 'Temperature–salinity diagnostic\nminus native mlotst', ha='center', va='center',
                          fontsize=10, color='#424C53', linespacing=1.7)
        legend_axis.text(.5, .39, {'initial': 'Difference in initial MLD', 'late': 'Difference in final MLD',
                                  'change': 'Difference in late-minus-initial response'}[quantity],
                          ha='center', va='center', fontsize=9, color='#424C53', wrap=True)
        cax = figure.add_axes([.23, .065, .54, .023])
        colorbar = figure.colorbar(filled, cax=cax, orientation='horizontal', ticks=ticks, extend='both')
        colorbar.set_label(colorbar_labels[quantity], fontsize=10, labelpad=7)
        colorbar.ax.tick_params(labelsize=8)
        colorbar.ax.minorticks_off()
        figure.canvas.draw()
        if quantity_index == 0:
            destination.mkdir(parents=True)
        for path in figure_paths[quantity]:
            figure.savefig(path, dpi=dpi, facecolor='white')
            print(path.resolve(), flush=True)
        plt.close(figure)
        exports.append({'quantity': quantity, 'panels': records, 'figure_size_inches': list(figure_size),
                         'dpi': dpi, 'paths': [str(path) for path in figure_paths[quantity]],
                         'caption': metadata['captions'][quantity].replace('beyond600m', 'beyond1000m')
                                    + ' The colormap, levels and ticks match Figure2 exactly.'})


# %% 5. 保存图注与来源记录

source_files = [source / name for name in ['mld_method_difference_patterns.npz', 'metadata.json', 'land_50m.geojson', 'coastline_50m.geojson']]
source_files.append(figure2_source / 'grouped_mld_sic_fields.npz')
manifest = {'exports': exports, 'models': models, 'panel_count_per_figure': 7, 'overall_title': None,
            'color_limits_m': [float(levels[0]), float(levels[-1])], 'color_ticks_m': ticks.tolist(),
            'color_levels_m': levels.tolist(), 'palette_anchors': palette_anchors,
            'normalization': 'Same BoundaryNorm and256-colour LinearSegmentedColormap as Figure2.',
            'difference_sign': metadata['difference_sign'],
            'source_files': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in source_files]}
(destination / 'figure_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
for quantity in quantities:
    caption = metadata['captions'][quantity].replace('beyond600m', 'beyond1000m') + ' The colormap, levels and ticks match Figure2 exactly.'
    (destination / ('caption_' + quantity + '.txt')).write_text(caption + '\n', encoding='utf-8')
