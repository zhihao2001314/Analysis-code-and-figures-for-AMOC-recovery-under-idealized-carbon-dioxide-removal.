# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, to_rgb, to_hex
from matplotlib.lines import Line2D
import cartopy.crs as ccrs
from shapely.geometry import shape


# %% 2. 参数与最终图件完整路径

root = Path(__file__).resolve().parents[3]
source = root / 'output/s05_atmosphere/p07_toa_energy_balance/n18_five_model_figure3'
significance_source = root / 'output/s05_atmosphere/p07_toa_energy_balance/n22_figure3_temperature_significance'
ice_source = root / 'output/s05_atmosphere/p07_toa_energy_balance/n23_figure3_bipolar_ice_edges'
destination = root / 'figure/s05_atmosphere/p07_toa_energy_balance/n19_five_model_figure3'
png_path = destination / 'figure3_eei_temperature_groups.png'
pdf_path = destination / 'figure3_eei_temperature_groups.pdf'
svg_path = destination / 'figure3_eei_temperature_groups.svg'
figure_paths = [png_path, pdf_path, svg_path]
assert not destination.exists() or not any(destination.iterdir()), '正式图件存在，须先归档。'
figure_size, dpi = (12.6, 6.5), 450
column_lefts, panel_width = [.075, .395, .715], .255
labels = ['Non-recovery', 'Recovery', 'Non-recovery − Recovery']
component_labels = ['Net', 'Shortwave', 'Longwave']
component_colors = ['#27323A', '#C85B3A', '#2B6F96']
component_styles = ['-', '-', '-']
# 单模式仅画短波和长波：同一色相，短波深、长波浅，全部实线。
model_base_colors = {'CESM2': '#CC6677', 'CNRM-ESM2-1': '#332288', 'CAS-ESM2-0': '#B39A32',
                     'GFDL-ESM4': '#44AA99', 'NorESM2-LM': '#882255'}
model_component_colors = {}
for model_name, base_color in model_base_colors.items():
    rgb = np.asarray(to_rgb(base_color))
    model_component_colors[model_name] = [to_hex(.85 * rgb), to_hex(.65 * rgb + .35)]
model_legend_records = []
ice_colors = ['#7B3DBA', '#27985B', '#E000C8']
ice_labels = ['Initial ice edge', 'Final ice edge', 'Observed ice edge (1988–2020)']
ice_records, hatch_records = [], []
stipple_records = []
stipple_stride = 4
artist_records = []
graticule_records = []


# %% 3. 仅读取已保存的完整组序列与地图

metadata = json.loads((source / 'metadata.json').read_text(encoding='utf-8'))
groups = pd.read_csv(source / 'group_radiation.csv', float_precision='round_trip')
contrast = pd.read_csv(source / 'radiation_contrast.csv', float_precision='round_trip')
model_radiation = pd.read_csv(source / 'model_radiation.csv', float_precision='round_trip')
significance_metadata = json.loads((significance_source / 'metadata.json').read_text(encoding='utf-8'))
ice_metadata = json.loads((ice_source / 'metadata.json').read_text(encoding='utf-8'))
with np.load(significance_source / 'temperature_significance.npz') as saved:
    significant = saved['significant'].copy()
    nonsignificant = saved['nonsignificant'].copy()
    mask_longitude = np.append(saved['longitude'], saved['longitude'][0] + 360)
    cyclic_mask = np.concatenate([nonsignificant, nonsignificant[:, :, :1]], axis=2)
with np.load(ice_source / 'grouped_bipolar_sic.npz') as saved:
    ice_longitude, ice_latitude = saved['longitude_cyclic'].copy(), saved['latitude'].copy()
    ice_cyclic = saved['ice_cyclic'].copy()
observations = {}
for hemisphere in ['north', 'south']:
    with np.load(ice_source / f'observed_{hemisphere}.npz') as saved:
        observations[hemisphere] = {key: saved[key].copy() for key in saved.files}
with np.load(source / 'temperature_panels.npz') as saved:
    latitude = saved['latitude'].copy()
    longitude = saved['longitude'].copy()
    temperature = saved['temperature_panels'].copy()
coastlines = [shape(item['geometry']) for item in json.loads(
    (source / 'coastline_global_50m.geojson').read_text(encoding='utf-8'))['features']]


# %% 4. 上排两组EEI分解与不恢复减恢复

with plt.rc_context({'font.family': 'Arial', 'font.size': 9, 'axes.spines.top': False,
                     'axes.spines.right': False, 'pdf.fonttype': 42, 'svg.fonttype': 'none',
                     'hatch.color': '#707070', 'hatch.linewidth': .24,
                     'axes.linewidth': .7, 'savefig.facecolor': 'white'}):
    figure = plt.figure(figsize=figure_size)
    for column, label in enumerate(labels):
        axis = figure.add_axes([column_lefts[column], .575, panel_width, .315])
        axis.set(xlim=(1, 280), xticks=[1, 70, 140, 210, 280],
                 ylim=(-5, 7), xlabel='Year')
        if column == 0:
            axis.set_ylabel('Control-relative radiation (W m$^{-2}$)')
        axis.set_title(label, fontsize=10, pad=9)
        axis.axhline(0, color='#AEB6BD', lw=.6, zorder=0)
        axis.axvline(140.5, color='#B8BEC5', lw=.7, ls=':', zorder=0)
        axis.grid(axis='y', color='#EDF0F2', lw=.5)
        axis.tick_params(labelsize=8, length=3)
        if column < 2:
            selected = groups.loc[groups.group.eq(['non_recovery', 'recovery'][column])]
            axis.text(.97, .96, f'n = {[2, 3][column]}', transform=axis.transAxes,
                       ha='right', va='top', fontsize=8)
        else:
            selected = contrast
        if column < 2:
            group = ['non_recovery', 'recovery'][column]
            for model in metadata['groups'][group]:
                individual = model_radiation.loc[model_radiation.model.eq(model)]
                for component_index, field in enumerate(metadata['radiation_fields'][1:]):
                    for first, last in [(1, 280)]:
                        part = individual.loc[individual.simulation_year.between(first, last)]
                        line, = axis.plot(part.simulation_year, part[field + '_smooth11'],
                                           color=model_component_colors[model][component_index], ls='-',
                                           lw=.80, alpha=.80, zorder=1)
                        artist_records.append({'panel': chr(97 + column), 'kind': 'individual',
                                               'model': model, 'field': field, 'first_year': first, 'last_year': last,
                                               'linewidth': line.get_linewidth(), 'alpha': line.get_alpha(),
                                               'linestyle': line.get_linestyle(),
                                               'model_color': line.get_color(),
                                               'x_sha256': hashlib.sha256(np.asarray(line.get_xdata(), dtype=float).tobytes()).hexdigest(),
                                               'y_sha256': hashlib.sha256(np.asarray(line.get_ydata(), dtype=float).tobytes()).hexdigest()})
        for field, color, style in zip(metadata['radiation_fields'], component_colors, component_styles):
            for first, last in [(1, 280)]:
                part = selected.loc[selected.simulation_year.between(first, last)]
                line, = axis.plot(part.simulation_year, part[field + '_smooth11'],
                                   color=color, ls=style, lw=1.55, zorder=3)
                artist_records.append({'panel': chr(97 + column), 'kind': 'group_mean_or_contrast',
                                       'field': field, 'first_year': first, 'last_year': last, 'linestyle': line.get_linestyle(),
                                       'x_sha256': hashlib.sha256(np.asarray(line.get_xdata(),dtype=float).tobytes()).hexdigest(),
                                       'y_sha256': hashlib.sha256(np.asarray(line.get_ydata(),dtype=float).tobytes()).hexdigest()})
        if column < 2:
            model_handles = [Line2D([], [], color=model_component_colors[model][0], ls='-', lw=1.1)
                              for model in metadata['groups'][group]]
            model_legend = axis.legend(handles=model_handles, labels=metadata['groups'][group],
                                        loc='upper left', frameon=False,
                                        fontsize=6.5, handlelength=2.7, labelspacing=.45,
                                        borderaxespad=.65, handletextpad=.7)
            model_legend_records.append({'panel': chr(97 + column),
                                         'handle_style': 'one solid model-colour line per model; no segmented swatches',
                                         'title': model_legend.get_title().get_text(),
                                         'labels': [item.get_text() for item in model_legend.get_texts()]})
        axis.text(-.10, 1.055, chr(97 + column), transform=axis.transAxes, fontsize=12, weight='bold')
    handles = [Line2D([], [], color=color, ls=style, lw=1.55, label=label)
               for label, color, style in zip(component_labels, component_colors, component_styles)]
    figure.legend(handles=handles, loc='center', bbox_to_anchor=(.515, .955), ncol=3,
                    frameon=False, fontsize=9, handlelength=2.5, columnspacing=2.,
                    title='Group means / difference', title_fontsize=7)


# %% 5. 下排完整3对3温度残余与组差，只标经纬度不画内部网线

    for column in range(3):
        axis = figure.add_axes([column_lefts[column], .16, panel_width, .30], projection=ccrs.Robinson())
        axis.set_global()
        mesh = axis.pcolormesh(longitude, latitude, temperature[column], transform=ccrs.PlateCarree(),
                               cmap='RdBu_r', norm=Normalize(-metadata['temperature_limit_degC'], metadata['temperature_limit_degC']),
                               shading='auto', rasterized=True, zorder=1)
        axis.contourf(mask_longitude, latitude, cyclic_mask[column].astype(float),
                       levels=[.5, 1.5], colors='none', hatches=['xx'],
                       transform=ccrs.PlateCarree(), zorder=2)
        hatch_records.append({'panel': chr(100 + column), 'criterion': 'three-model sign agreement' if column < 2
                               else 'centered circular-block bootstrap, global BH q<0.05',
                               'hatched_cells': int(nonsignificant[column].sum()),
                               'mask_sha256': hashlib.sha256(nonsignificant[column].tobytes()).hexdigest()})
        # 固定每4个网格点显示一次，通过判据的原始格点打点，不改变检验掩膜。
        dot_longitude, dot_latitude = np.meshgrid(longitude[::stipple_stride], latitude[::stipple_stride])
        dot_mask = significant[column, ::stipple_stride, ::stipple_stride]
        dots = axis.scatter(dot_longitude[dot_mask], dot_latitude[dot_mask], s=1.2, marker='.',
                              color='#252525', alpha=.70, linewidths=0, edgecolors='none',
                              transform=ccrs.PlateCarree(), zorder=2.5, rasterized=True)
        stipple_records.append({'panel': chr(100 + column), 'stride': stipple_stride,
                                 'display_points': int(dot_mask.sum()), 'marker_size_pt2': 1.2,
                                 'all_passed_cells': int(significant[column].sum()),
                                 'mask_sha256': hashlib.sha256(significant[column].tobytes()).hexdigest(),
                                 'offsets_sha256': hashlib.sha256(np.asarray(dots.get_offsets(), dtype=float).tobytes()).hexdigest()})
        if column < 2:
            for period_index in range(2):
                contours = axis.contour(ice_longitude, ice_latitude, ice_cyclic[column, period_index],
                                         levels=[.15], colors=[ice_colors[period_index]], linewidths=.50,
                                         linestyles='--', transform=ccrs.PlateCarree(), zorder=5)
                ice_records.append({'panel': chr(100 + column), 'period': ['initial', 'final'][period_index],
                                    'level': .15, 'linewidth': .50, 'color': ice_colors[period_index],
                                    'segments': len(contours.allsegs[0]),
                                    'field_sha256': hashlib.sha256(ice_cyclic[column, period_index].tobytes()).hexdigest()})
            for hemisphere, observed in observations.items():
                contours = axis.contour(observed['x'], observed['y'], observed['sic'],
                                         levels=[.15], colors=[ice_colors[2]], linewidths=.50,
                                         transform=ccrs.epsg(int(observed['epsg'])), zorder=5)
                ice_records.append({'panel': chr(100 + column), 'period': 'observed', 'hemisphere': hemisphere,
                                    'level': .15, 'linewidth': .50, 'color': ice_colors[2],
                                    'segments': len(contours.allsegs[0])})
        axis.add_geometries(coastlines, ccrs.PlateCarree(), facecolor='none',
                             edgecolor='#65717A', linewidth=.35, zorder=3)
        grid = axis.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                              xlocs=[-120, -60, 0, 60, 120], ylocs=[-60, -30, 0, 30, 60],
                              x_inline=False, y_inline=False, rotate_labels=False,
                              linewidth=0, zorder=2)
        grid.xlines, grid.ylines = False, False
        grid.top_labels, grid.right_labels = False, False
        grid.bottom_labels, grid.left_labels = True, True
        grid.xlabel_style, grid.ylabel_style = {'size': 6.8}, {'size': 6.8}
        grid.xpadding, grid.ypadding = 3, 3
        axis.spines['geo'].set_linewidth(.6)
        axis.spines['geo'].set_edgecolor('#9CA5AC')
        axis.text(-.10, 1.055, chr(100 + column), transform=axis.transAxes, fontsize=12, weight='bold')
        if column < 2:
            axis.text(.98, 1.055, 'n = 3', transform=axis.transAxes,
                       fontsize=8, ha='right', va='bottom')
        graticule_records.append({'panel': chr(100 + column), 'artist': grid})
        artist_records.append({'panel': chr(100 + column),
                               'field_sha256': hashlib.sha256(temperature[column].tobytes()).hexdigest()})
    color_axis = figure.add_axes([.32, .09, .40, .021])
    limit = metadata['temperature_limit_degC']
    colorbar = figure.colorbar(mesh, cax=color_axis, orientation='horizontal', ticks=np.linspace(-limit, limit, 5))
    colorbar.set_label('Near-surface air temperature change (°C)', fontsize=9)
    colorbar.ax.tick_params(labelsize=8, length=2.5)
    colorbar.outline.set_linewidth(.6)
    ice_handles = [Line2D([], [], color=color, ls=style, lw=.65, label=label)
                   for color, style, label in zip(ice_colors, ['--', '--', '-'], ice_labels)]
    figure.legend(handles=ice_handles, loc='center', bbox_to_anchor=(.52, .017), ncol=3,
                    frameon=False, fontsize=7.2, handlelength=2.5, columnspacing=1.8)
    figure.canvas.draw()
    for record in graticule_records:
        grid = record.pop('artist')
        record.update({'internal_xlines': grid.xlines, 'internal_ylines': grid.ylines,
                       'longitude_labels': [item.get_text() for item in grid.xlabel_artists if item.get_visible()],
                       'latitude_labels': [item.get_text() for item in grid.ylabel_artists if item.get_visible()]})
        assert record['longitude_labels'] and record['latitude_labels']
    destination.mkdir(parents=True, exist_ok=True)
    for path in figure_paths:
        figure.savefig(path, dpi=dpi, facecolor='white')
        print(path.resolve(), flush=True)
    plt.close(figure)


# %% 6. 保存实际画线与完整面板记录

curve_counts = {panel: sum(record['panel'] == panel for record in artist_records if 'kind' in record)
                 for panel in ['a', 'b', 'c']}
assert curve_counts == {'a': 7, 'b': 9, 'c': 3}
manifest = {'status': metadata['status'], 'panel_count': 6, 'figure_size_inches': list(figure_size),
            'dpi': dpi, 'overall_title': None, 'groups': metadata['groups'], 'group_sizes': [2, 3],
            'temperature_groups': metadata['temperature_groups'], 'temperature_group_sizes': [3, 3],
            'graticule_records': graticule_records,
            'radiation_limits': [-5, 7], 'component_linestyles': component_styles,
            'model_component_colors': model_component_colors, 'model_legend_records': model_legend_records,
            'radiation_curve_counts': curve_counts,
            'individual_SW_LW_encoding': 'Only individual SW and LW; one hue per model, SW darker, LW lighter, all solid. No individual net curves. Three thick solid group means and difference retain black/orange/blue component colours.',
            'hatch_records': hatch_records, 'ice_records': ice_records, 'panel_titles': labels,
            'stipple_records': stipple_records,
            'significance_metadata': significance_metadata, 'ice_metadata': ice_metadata,
            'missing_radiation_panels': [], 'artist_records': artist_records,
            'CAS_alignment_limitation': metadata['CAS_alignment_limitation'],
            'caption': significance_metadata['caption'].replace(
                'Thin lines show individual models and thick lines show equal-weight group means; all components are plotted as solid lines.',
                'Thin solid lines show individual-model shortwave and longwave contributions, using a model-specific hue with darker shades for shortwave and lighter shades for longwave. Thick solid lines show equal-weight group means in (a, b) and their difference in (c), with net, shortwave, and longwave contributions in black, orange, and blue, respectively. Individual-model net radiation is omitted for clarity.')
                .replace(
                'Lines show saved 11-year centred means, computed within each forcing phase;',
                'Lines show 11-year centred means computed across the full 280-year sequence;')
                + ' Stippling marks three-model sign agreement in (d, e) and FDR-significant differences in (f); cross-hatching marks cells failing the respective criterion. Dots are displayed every fourth grid point for clarity. In (d, e), purple and green dashed lines show the initial and final 15% contours of three-model group-mean annual sea-ice concentration in both hemispheres. Magenta solid lines show NOAA/NSIDC observed annual-mean ice edges for 1988–2020.',
            'source_files': [{'path': str(source / name), 'sha256': hashlib.sha256((source / name).read_bytes()).hexdigest()}
                             for name in ['metadata.json', 'temperature_panels.npz', 'group_radiation.csv',
                                          'radiation_contrast.csv', 'coastline_global_50m.geojson', 'model_radiation.csv']]
            + [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in
               [significance_source / 'temperature_significance.npz', significance_source / 'metadata.json',
                ice_source / 'grouped_bipolar_sic.npz', ice_source / 'metadata.json',
                ice_source / 'observed_north.npz', ice_source / 'observed_south.npz']]}
(destination / 'figure_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
(destination / 'caption.txt').write_text(manifest['caption'] + '\n', encoding='utf-8')
