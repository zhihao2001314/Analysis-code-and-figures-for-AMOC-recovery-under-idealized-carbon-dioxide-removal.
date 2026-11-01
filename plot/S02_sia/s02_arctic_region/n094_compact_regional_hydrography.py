# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.ticker import FixedLocator
import cartopy.crs as ccrs
from shapely.geometry import shape
from scipy.interpolate import interp1d


# %% 2. 参数与最终图件完整路径

root = Path(__file__).resolve().parents[3]
stage = Path('s02_sia/s02_arctic_region/p23_stratification_mld_composite')
source = root / 'output' / stage / 'n093_compact_regional_hydrography'
destination = root / 'figure' / stage / 'n094_compact_regional_hydrography'
png_path = destination / 'li_nordic_nine_models_grouped_mld_ice.png'
pdf_path = destination / 'li_nordic_nine_models_grouped_mld_ice.pdf'
svg_path = destination / 'li_nordic_nine_models_grouped_mld_ice.svg'
figure_paths = [png_path, pdf_path, svg_path]
assert not destination.exists() or not any(destination.iterdir()), '正式图件存在，须先归档。'
figure_size, dpi = (11.4, 10.6), 450
column_lefts, column_width = [.085, .57], .39
profile_bottom, profile_height = .655, .295
history_bottom, history_height = .39, .162
map_lefts, map_width, map_bottom, map_height = [.055, .375, .695], .285, .095, .24
initial_ice_color, final_ice_color, observed_ice_color = '#6A3D9A', '#2E8B57', '#CC00A6'
mld_contour_label_levels_m = [200, 400, 600, 800, 1000]
projection = ccrs.LambertConformal(central_longitude=-25., central_latitude=67.5, standard_parallels=(50., 80.))
data_crs, observed_crs = ccrs.PlateCarree(), ccrs.epsg(3411)
artist_records, map_records, axis_records = [], [], []
text_checks = []
profile_axes = []


# %% 3. 仅读取保存的剖面序列及分组空间场

metadata = json.loads((source / 'metadata.json').read_text(encoding='utf-8'))
profiles = pd.read_csv(source / 'member_n2_profiles.csv', float_precision='round_trip')
observations = pd.read_csv(source / 'iap_profiles_10_1000m.csv', float_precision='round_trip')
history = pd.read_csv(source / 'member_mld_timeseries.csv', float_precision='round_trip')
membership = pd.read_csv(source / 'map_membership.csv')
with np.load(source / 'grouped_mld_sic_fields.npz') as saved:
    fields = {name: saved[name].copy() for name in saved.files}
land = [shape(item['geometry']) for item in json.loads((source / 'map_land_50m.geojson').read_text())['features']]
coastlines = [shape(item['geometry']) for item in json.loads((source / 'map_coastline_50m.geojson').read_text())['features']]
vertices = pd.read_csv(source / 'mld_rectangular_vertices.csv')


# %% 4. 两区九模式层结与MLD，不画组平均

with plt.rc_context({'font.family': 'Arial', 'font.size': 9, 'axes.labelsize': 10,
                     'axes.linewidth': .7, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
                     'pdf.fonttype': 42, 'svg.fonttype': 'none', 'savefig.facecolor': 'white'}):
    figure = plt.figure(figsize=figure_size, facecolor='white')
    for column, region in enumerate(metadata['regions']):
        profile_axis = figure.add_axes([column_lefts[column], profile_bottom, column_width, profile_height])
        history_axis = figure.add_axes([column_lefts[column], history_bottom, column_width, history_height])
        # 库提供的插值对象作正逆轴变换；不定义自编函数，不修改原始N²曲线。
        profile_axis.set_xscale('function', functions=(
            interp1d(metadata['n2_axis_transform']['raw_knots'], metadata['n2_axis_transform']['display_knots'],
                      bounds_error=False, fill_value='extrapolate', assume_sorted=True),
            interp1d(metadata['n2_axis_transform']['display_knots'], metadata['n2_axis_transform']['raw_knots'],
                      bounds_error=False, fill_value='extrapolate', assume_sorted=True)))
        profile_axes.append((chr(97 + column), profile_axis))
        obs = observations.loc[observations.region.eq(region)].sort_values('depth_m')
        profile_axis.fill_betweenx(obs.depth_m, obs.lower_1sd_n2_1e5_s2, obs.upper_1sd_n2_1e5_s2,
                                    color='#C7C7C7', alpha=.4, lw=0, zorder=1)
        for model in metadata['models']:
            profile = profiles.loc[profiles.region.eq(region) & profiles.model.eq(model)].sort_values('depth_m')
            line, = profile_axis.plot(profile.n2_1e5_s2, profile.depth_m,
                                       color=metadata['model_colors'][model], lw=1.25, ls='-', zorder=3)
            artist_records.append({'panel': chr(97 + column), 'region': region, 'quantity': 'member_n2', 'model': model,
                                   'x_sha256': hashlib.sha256(np.asarray(line.get_xdata(), dtype=float).tobytes()).hexdigest(),
                                   'y_sha256': hashlib.sha256(np.asarray(line.get_ydata(), dtype=float).tobytes()).hexdigest()})
            values = history.loc[history.region.eq(region) & history.model.eq(model)].sort_values('simulation_year')
            line, = history_axis.plot(values.simulation_year, values.march_mld_11yr_m,
                                       color=metadata['model_colors'][model], lw=1.2, ls='-', zorder=3)
            artist_records.append({'panel': chr(99 + column), 'region': region, 'quantity': 'member_mld', 'model': model,
                                   'x_sha256': hashlib.sha256(np.asarray(line.get_xdata(), dtype=float).tobytes()).hexdigest(),
                                   'y_sha256': hashlib.sha256(np.asarray(line.get_ydata(), dtype=float).tobytes()).hexdigest()})
        line, = profile_axis.plot(obs.iapv4_march_mean_n2_1e5_s2, obs.depth_m, color='#111111', lw=1.5, zorder=5)
        artist_records.append({'panel': chr(97 + column), 'region': region, 'quantity': 'IAP', 'model': 'IAP',
                               'x_sha256': hashlib.sha256(np.asarray(line.get_xdata(), dtype=float).tobytes()).hexdigest(),
                               'y_sha256': hashlib.sha256(np.asarray(line.get_ydata(), dtype=float).tobytes()).hexdigest(),
                               'depth_min_m': float(obs.depth_m.min()), 'depth_max_m': float(obs.depth_m.max())})
        profile_axis.set(xlim=metadata['n2_x_limits'], ylim=(1000, 0), xticks=metadata['n2_x_ticks'],
                         yticks=[0, 200, 400, 600, 800, 1000], xlabel=r'$N^2$ ($10^{-5}$ s$^{-2}$)', ylabel='Depth (m)')
        profile_axis.axvline(0, color='#B8BEC3', lw=.55, ls=':', zorder=0)
        profile_axis.grid(axis='y', color='#E8EAEC', lw=.5)
        title = 'Labrador–Irminger Seas' if column == 0 else 'Nordic Seas'
        profile_axis.set_title(title, loc='center', fontsize=10.5, pad=9)
        history_axis.set(xlim=(1, 280), ylim=metadata['mld_y_limits_m'], xticks=[1, 70, 140, 210, 280],
                         yticks=metadata['mld_y_ticks_m'], xlabel='Year', ylabel='March MLD (m)')
        history_axis.axvline(140.5, color='#A5ABB1', lw=.8, ls=':', zorder=0)
        history_axis.grid(axis='y', color='#E8EAEC', lw=.5)
        for letter, axis in [(chr(97 + column), profile_axis), (chr(99 + column), history_axis)]:
            axis.spines[['top', 'right']].set_visible(False)
            axis.tick_params(length=3)
            letter_artist = axis.text(-.10, 1.02, letter, transform=axis.transAxes, fontsize=12, weight='bold')
            text_checks.extend([(letter + '_letter', letter_artist), (letter + '_xlabel', axis.xaxis.label)])
            axis_records.append({'panel': letter, 'bounds': list(axis.get_position().bounds),
                                 'xlim': list(axis.get_xlim()), 'ylim': list(axis.get_ylim()), 'xticks': axis.get_xticks().tolist()})
        if column == 0:
            inset = profile_axis.inset_axes([.48, .12, .48, .43], projection=projection)
            inset.set_extent([-66, 24, 49, 83], crs=data_crs)
            inset.set_facecolor('#F8FBFC')
            inset.add_geometries(land, data_crs, facecolor='#DEDEDE', edgecolor='none', zorder=1)
            inset.add_geometries(coastlines, data_crs, facecolor='none', edgecolor='#666666', linewidth=.35, zorder=3)
            for name, color, abbreviation, label_x, label_y in [('Labrador-Irminger', '#D55E00', 'LI', -42, 57),
                                                               ('Nordic Seas', '#0072B2', 'NS', 2, 72)]:
                points = vertices.loc[vertices.region.eq(name)].sort_values('vertex_order')
                inset.add_patch(MplPolygon(points[['longitude', 'latitude']].to_numpy(), closed=True,
                                            facecolor=color, alpha=.15, edgecolor=color, lw=.8,
                                            transform=data_crs, zorder=2))
                inset.add_patch(MplPolygon(points[['longitude', 'latitude']].to_numpy(), closed=True,
                                            facecolor='none', edgecolor=color, lw=.65,
                                            transform=data_crs, zorder=4))
                inset.text(label_x, label_y, abbreviation, transform=data_crs, color=color,
                            fontsize=6.5, weight='bold', ha='center', va='center', zorder=5)
            inset.spines['geo'].set_linewidth(.5)
    handles = [Line2D([], [], color=metadata['model_colors'][model], lw=1.5, label=model) for model in metadata['models']]
    handles.append(Line2D([], [], color='#111111', lw=1.5, label='IAP mean ±1 SD'))
    model_legend = figure.legend(handles=handles, loc='center', bbox_to_anchor=(.53, .595), ncol=5,
                                  frameon=False, fontsize=7.6, handlelength=2.2,
                                  columnspacing=1.7, handletextpad=.6, labelspacing=.6)


# %% 5. 三组三模式MLD平均及组平均海冰边缘

    cmap = LinearSegmentedColormap.from_list('mld_residual',
           ['#145D9E', '#347FB5', '#70A6CA', '#B6D0DF', '#F5F3EF', '#ECC8B4', '#DF936E', '#CE5144', '#B60026'], N=256)
    cmap.set_bad('#EDF1F5')
    norm = BoundaryNorm(fields['color_levels_m'], cmap.N)
    for column, (group, name) in enumerate(zip(metadata['groups'], metadata['group_names'])):
        axis = figure.add_axes([map_lefts[column], map_bottom, map_width, map_height], projection=projection)
        axis.set_extent(metadata['map_extent'], crs=data_crs)
        axis.set_facecolor('#EDF1F5')
        grid = axis.gridlines(draw_labels=True, x_inline=False, y_inline=False, rotate_labels=False,
                              linewidth=.35, color='#999999', linestyle=':', alpha=.55, zorder=0)
        grid.xlocator = FixedLocator([-60, -40, -20, 0, 20])
        grid.ylocator = FixedLocator([50, 60, 70, 80])
        grid.top_labels, grid.right_labels, grid.left_labels = False, False, column == 0
        grid.xlabel_style, grid.ylabel_style = {'size': 7}, {'size': 7}
        filled = axis.contourf(fields['longitude'], fields['latitude'], np.ma.masked_invalid(fields[group + '_delta_mld_m']),
                                levels=fields['color_levels_m'], cmap=cmap, norm=norm, extend='both',
                                transform=data_crs, antialiased=False, corner_mask=False, zorder=1)
        contours = axis.contour(fields['longitude'], fields['latitude'], np.ma.masked_invalid(fields[group + '_initial_mld_m']),
                                levels=fields['contour_levels_m'], colors='#222222', linewidths=.6,
                                linestyles='-', transform=data_crs, corner_mask=False, zorder=3)
        contour_labels = axis.clabel(contours, levels=mld_contour_label_levels_m, fmt='%d',
                                     fontsize=6, inline=True, inline_spacing=3)
        edge_counts = {}
        for period, color in [('initial', initial_ice_color), ('later', final_ice_color)]:
            ice_contour = axis.contour(fields['ice_longitude'], fields['ice_latitude'], np.ma.masked_invalid(fields[group + '_ice_' + period]),
                                        levels=[metadata['ice_threshold']], colors=[color], linewidths=1.05,
                                        linestyles='--', transform=data_crs, zorder=7)
            edge_counts[period] = sum(len(part) > 1 for part in ice_contour.allsegs[0])
            assert edge_counts[period] > 0
        observed = axis.contour(fields['observed_x_m'], fields['observed_y_m'], np.ma.masked_invalid(fields['observed_sic_fraction']),
                                levels=[metadata['ice_threshold']], colors=[observed_ice_color], linewidths=1.05,
                                transform=observed_crs, zorder=8)
        observed.set_path_effects([path_effects.Stroke(linewidth=1.8, foreground='white'), path_effects.Normal()])
        axis.add_geometries(land, data_crs, facecolor='#DEDEDE', edgecolor='none', zorder=9)
        axis.add_geometries(coastlines, data_crs, facecolor='none', edgecolor='#555555', linewidth=.4, zorder=10)
        axis.spines['geo'].set_linewidth(.65)
        axis.spines['geo'].set_edgecolor('#777777')
        axis.set_title(name, fontsize=10, pad=7)
        axis.annotate(chr(101 + column), xy=(0, 1), xycoords='axes fraction', xytext=(-5, 5),
                       textcoords='offset points', fontsize=12, weight='bold', annotation_clip=False)
        map_records.append({'panel': chr(101 + column), 'group': group,
                            'members': membership.loc[membership.group.eq(group)].sort_values('member_order').model.tolist(),
                            'initial_sha256': hashlib.sha256(fields[group + '_initial_mld_m'].tobytes()).hexdigest(),
                            'delta_sha256': hashlib.sha256(fields[group + '_delta_mld_m'].tobytes()).hexdigest(),
                            'initial_SIC_sha256': hashlib.sha256(fields[group + '_ice_initial'].tobytes()).hexdigest(),
                            'later_SIC_sha256': hashlib.sha256(fields[group + '_ice_later'].tobytes()).hexdigest(),
                            'mld_contour_levels_m': contours.levels.tolist(),
                            'mld_contour_labels_m': sorted(set(int(label.get_text()) for label in contour_labels)),
                            'ice_contour_segments': edge_counts})
    colorbar_axis = figure.add_axes([.26, .073, .51, .013])
    colorbar = figure.colorbar(filled, cax=colorbar_axis, orientation='horizontal', ticks=fields['color_ticks_m'])
    colorbar.set_label('March MLD change (m)', fontsize=9, labelpad=3)
    colorbar.ax.tick_params(labelsize=8, length=3, pad=2)
    colorbar.outline.set_linewidth(.6)
    map_handles = [Line2D([], [], color='#222222', lw=.7, label='Initial group-mean MLD'),
                   Line2D([], [], color=initial_ice_color, ls='--', lw=1.05, label='Initial group-mean ice edge'),
                   Line2D([], [], color=final_ice_color, ls='--', lw=1.05, label='Final group-mean ice edge'),
                   Line2D([], [], color=observed_ice_color, lw=1.05, label='Observed ice edge (1988–2020)')]
    map_legend = figure.legend(handles=map_handles, loc='center', bbox_to_anchor=(.52, .021), ncol=4,
                                frameon=False, fontsize=7.2, handlelength=2.5, columnspacing=1.4)
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    tick_spacing_records = []
    for panel, profile_axis in profile_axes:
        tick_values = np.asarray(metadata['n2_x_ticks'], dtype=float)
        tick_pixels = profile_axis.transData.transform(np.column_stack([tick_values, np.zeros(len(tick_values))]))[:, 0]
        np.testing.assert_allclose(np.diff(tick_pixels), np.repeat(np.diff(tick_pixels)[0], len(tick_values) - 1), atol=1e-6, rtol=0)
        tick_spacing_records.append({'panel': panel, 'tick_values': tick_values.tolist(), 'display_x_pixels': tick_pixels.tolist(),
                                     'equal_spacing_verified': True})
    legend_box = model_legend.get_window_extent(renderer)
    overlaps = [name for name, artist in text_checks if legend_box.overlaps(artist.get_window_extent(renderer))]
    assert not overlaps, ('Shared legend overlaps axes text', overlaps)
    legend_records = {'model_legend_bounds': list(model_legend.get_window_extent(renderer).transformed(figure.transFigure.inverted()).bounds),
                       'map_legend_bounds': list(map_legend.get_window_extent(renderer).transformed(figure.transFigure.inverted()).bounds),
                       'axis_text_overlap': overlaps}
    destination.mkdir(parents=True, exist_ok=True)
    for path in figure_paths:
        figure.savefig(path, dpi=dpi, facecolor='white')
        print(path.resolve(), flush=True)
    plt.close(figure)


# %% 6. 保存实际绘图清单

manifest = {'panel_count': 7, 'figure_size_inches': list(figure_size), 'dpi': dpi, 'overall_title': None,
            'artist_records': artist_records, 'map_records': map_records, 'axes': axis_records,
            'legends': legend_records, 'no_group_mean_curves': True, 'group_sizes': [3, 3, 3],
            'n2_axis_transform': metadata['n2_axis_transform'], 'n2_tick_spacing': tick_spacing_records,
            'caption': metadata['caption'], 'model_colors': metadata['model_colors'],
            'mld_color_limits_m': [float(fields['color_levels_m'][0]), float(fields['color_levels_m'][-1])],
            'mld_color_ticks_m': fields['color_ticks_m'].tolist(),
            'source_files': [{'path': str(source / name), 'sha256': hashlib.sha256((source / name).read_bytes()).hexdigest()}
                             for name in ['metadata.json', 'member_n2_profiles.csv', 'iap_profiles_10_1000m.csv', 'member_mld_timeseries.csv',
                                          'map_membership.csv', 'grouped_mld_sic_fields.npz', 'map_land_50m.geojson', 'map_coastline_50m.geojson',
                                          'mld_rectangular_vertices.csv']],
            'alt_text': 'Seven-panel compact hydrography figure. Top: LI and Nordic Seas initial March N2, all9models and IAP reference only below10m. Middle: same9models regional March MLD with no group mean. Bottom: non-recovery, recovery and larger-bias MLD residual maps,3models each, with initial MLD contours, group-mean initial/final sea-ice edges and common observed edge. UK belongs to non-recovery and CAS to recovery.'}
(destination / 'figure_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
