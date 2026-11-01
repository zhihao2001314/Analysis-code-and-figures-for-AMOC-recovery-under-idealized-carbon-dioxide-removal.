# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
from shapely.geometry import shape


# %% 2. 参数与最终图件完整路径

root = Path(__file__).resolve().parents[3]
base = root / 'output/s05_atmosphere/p05_amoc_recovery_composites'
hydrology_source = base / 'n81_six_model_net_precipitation'
peak_source = base / 'n83_net_precipitation_peak_markers'
flow_source = base / 'n74_six_model_circulation_composites'
destination = root / 'figure/s05_atmosphere/p05_amoc_recovery_composites/n82_six_model_circulation_net_precipitation'
png_path = destination / 'annual_circulation_global_net_precipitation_six_models.png'
pdf_path = destination / 'annual_circulation_global_net_precipitation_six_models.pdf'
svg_path = destination / 'annual_circulation_global_net_precipitation_six_models.svg'
figure_paths = [png_path, pdf_path, svg_path]
assert not destination.exists(), '正式图件存在，须先归档。'
products = ['recovered', 'not_recovered', 'contrast']
row_names = ['Recovery (n = 3)', 'Non-recovery (n = 3)', 'Non-recovery − Recovery']
psi_levels, pr_levels = np.linspace(-40, 40, 258), np.linspace(-3, 3, 258)
psi_cmap = ListedColormap(np.load(base / 'check_red_white_blue_style/red_white_blue_rgb.npy'))
psi_cmap.set_bad('#E7E7E7')
pr_cmap = plt.get_cmap('BrBG', 257)
psi_norm, pr_norm = BoundaryNorm(psi_levels, 257), BoundaryNorm(pr_levels, 257)
coast_path = root / 'output/s05_atmosphere/p07_toa_energy_balance/n04_toa_spatial_patterns/coastline_global_110m.geojson'
figure_size, dpi = (10.4, 9.6), 450
central_longitude = 180
phase_colors = {'early': '#225EA8', 'later': '#9E1B78'}
artist_records, significance_records, label_records = [], [], []
marker_records, connector_records = [], []


# %% 3. 仅读取完整三对三结果与净降水掩膜

metadata = json.loads((flow_source / 'metadata.json').read_text(encoding='utf-8'))
assert metadata['status'] == 'passed' and metadata['group_sizes'] == [3, 3]
with xr.open_dataarray(flow_source / 'group_circulation_climatology.nc') as dataset:
    climate = dataset.load()
with xr.open_dataarray(flow_source / 'group_circulation_changes.nc') as dataset:
    changes = dataset.load()
with xr.open_dataset(flow_source / 'psi_significance.nc') as dataset:
    psi_statistics = dataset.load()
with xr.open_dataset(hydrology_source / 'group_changes.nc') as dataset:
    pr_changes = dataset.pe.load()
with xr.open_dataset(hydrology_source / 'pe_display_criteria.nc') as dataset:
    pr_display = dataset.load()
peaks = pd.read_csv(peak_source / 'net_precipitation_peak_gridpoints.csv')
peak_metadata = json.loads((peak_source / 'metadata.json').read_text(encoding='utf-8'))
coast = [shape(feature['geometry']) for feature in json.loads(coast_path.read_text(encoding='utf-8'))['features']]


# %% 4. 六面板环流与净降水，保留点状标记和网纹

with plt.rc_context({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.linewidth': .65,
                     'pdf.fonttype': 42, 'svg.fonttype': 'none', 'savefig.facecolor': 'white',
                     'hatch.color': '#777777', 'hatch.linewidth': .25}):
    figure = plt.figure(figsize=figure_size, facecolor='white')
    grid = figure.add_gridspec(3, 2, left=.075, right=.985, bottom=.14, top=.91,
                               width_ratios=[1., 1.22], wspace=.10, hspace=.26)
    flow_handles = [Line2D([], [], color='#414141', lw=.7, ls='-', label='Initial: positive'),
                    Line2D([], [], color='#414141', lw=.7, ls='--', label='Initial: negative')]
    figure.legend(handles=flow_handles, loc='upper center', bbox_to_anchor=(.270, .986), ncol=2,
                   frameon=False, fontsize=8, columnspacing=1.5)
    peak_handles = [Line2D([], [], color=phase_colors['early'], marker='o', mfc='white', ls='none', ms=4.5,
                           label='Initial P−E peak'),
                    Line2D([], [], color=phase_colors['later'], marker='x', ls='none', ms=4.5,
                           label='Later P−E peak')]
    figure.legend(handles=peak_handles, loc='upper center', bbox_to_anchor=(.722, .986), ncol=2,
                   frameon=False, fontsize=8, columnspacing=1.5)
    figure.text(.270, .944, r'Mass streamfunction ($10^9$ kg s$^{-1}$)', ha='center', fontsize=10, weight='semibold')
    figure.text(.722, .944, r'Net precipitation change (mm day$^{-1}$)', ha='center', fontsize=10, weight='semibold')
    for row, product in enumerate(products):
        axis = figure.add_subplot(grid[row, 0])
        field = changes.sel(product=product, plev=slice(10000, 100000))
        image = axis.contourf(field.lat, field.plev / 100, field, levels=psi_levels,
                             cmap=psi_cmap, norm=psi_norm, extend='both', antialiased=False, zorder=1)
        axis.set_rasterization_zorder(2)
        axis.set_facecolor('#E7E7E7')
        artist_records.append({'kind': 'psi_change', 'product': product,
                               'sha256': hashlib.sha256(np.asarray(field.values, dtype=float).tobytes()).hexdigest()})
        if product == 'contrast':
            diagnostic = psi_statistics.sel(product=product)
            np.testing.assert_allclose(field, diagnostic.delta, atol=1e-10, rtol=0)
            axis.contourf(field.lat, field.plev / 100, diagnostic.nonsignificant.astype(float),
                          levels=[.5, 1.5], colors='none', hatches=['xx'], zorder=3)
            latitude_grid, pressure_grid = np.meshgrid(field.lat.values, field.plev.values / 100)
            selected = diagnostic.significant.values.copy()
            selected[:, np.arange(selected.shape[1]) % 3 != 0] = False
            dots = axis.scatter(latitude_grid[selected], pressure_grid[selected], marker='.', s=3.,
                                color='#222222', linewidths=0, zorder=5)
            significance_records.append({'kind': 'psi', 'product': product, 'n_displayed_dots': int(selected.sum()),
                                         'dot_offsets_sha256': hashlib.sha256(np.asarray(dots.get_offsets(), dtype=float).tobytes()).hexdigest(),
                                         'significant_mask_sha256': hashlib.sha256(diagnostic.significant.values.tobytes()).hexdigest(),
                                         'nonsignificant_mask_sha256': hashlib.sha256(diagnostic.nonsignificant.values.tobytes()).hexdigest()})
        else:
            initial = climate.sel(product=product, period='early', plev=field.plev)
            contour_levels = [-200, -150, -100, -75, -50, -25, 25, 50, 75, 100, 150, 200]
            contours = axis.contour(initial.lat, initial.plev / 100, initial, levels=contour_levels, colors='#414141',
                                    linestyles=['--' if value < 0 else '-' for value in contour_levels], linewidths=.55)
            axis.clabel(contours, inline=True, inline_spacing=2, fontsize=6, fmt='%g')
            axis.contour(initial.lat, initial.plev / 100, initial, levels=[0], colors='#8C8C8C', linewidths=.65)
            artist_records.append({'kind': 'psi_initial', 'product': product,
                                   'sha256': hashlib.sha256(np.asarray(initial.values, dtype=float).tobytes()).hexdigest()})
        axis.axvline(0, color='#B3B3B3', lw=.5)
        axis.set_yscale('log')
        axis.yaxis.set_minor_formatter(NullFormatter())
        axis.tick_params(which='minor', left=False, right=False)
        axis.set_ylim(1000, 100)
        axis.set_xlim(-90, 90)
        axis.set_xticks([-90, -60, -30, 0, 30, 60, 90], labels=['90°S', '60°S', '30°S', '0°', '30°N', '60°N', '90°N'])
        axis.set_yticks([100, 200, 300, 500, 700, 1000], labels=['100', '200', '300', '500', '700', '1000'])
        axis.tick_params(labelsize=7.5, length=3)
        axis.set_ylabel('Pressure (hPa)', fontsize=9)
        axis.set_title(chr(97 + row * 2), loc='left', weight='bold', fontsize=11, pad=8)
        axis.set_title(row_names[row], loc='right', fontsize=9.3, pad=8)
        if row == 2:
            axis.set_xlabel('Latitude', fontsize=9)
            cax = figure.add_axes([.105, .070, .30, .017])
            colorbar = figure.colorbar(image, cax=cax, orientation='horizontal', ticks=[-40, -20, 0, 20, 40],
                                       label=r'$\Delta\Psi$ ($10^9$ kg s$^{-1}$)')
            colorbar.ax.minorticks_off()
            colorbar.solids.set_edgecolor('face')
        map_axis = figure.add_subplot(grid[row, 1], projection=ccrs.Robinson(central_longitude=central_longitude))
        map_axis.set_anchor('W')
        map_axis.set_global()
        field = pr_changes.sel(product=product)
        cyclic_values, cyclic_lon = add_cyclic_point(field.values, coord=field.lon.values)
        map_image = map_axis.pcolormesh(cyclic_lon, field.lat, cyclic_values, transform=ccrs.PlateCarree(),
                                       cmap=pr_cmap, norm=pr_norm, shading='auto', rasterized=True)
        artist_records.append({'kind': 'pe_change', 'product': product,
                               'sha256': hashlib.sha256(np.asarray(field.values, dtype=float).tobytes()).hexdigest()})
        display_criterion = pr_display.sel(product=product)
        np.testing.assert_array_equal(field, display_criterion.delta)
        criterion_passed = display_criterion.passed.values.copy()
        criterion_failed = display_criterion.not_passed.values.copy()
        cyclic_mask, mask_lon = add_cyclic_point(criterion_failed.astype(float), coord=field.lon.values)
        map_axis.contourf(mask_lon, field.lat, cyclic_mask, levels=[.5, 1.5], colors='none', hatches=['xx'],
                          transform=ccrs.PlateCarree(), zorder=3)
        longitude_grid, latitude_grid = np.meshgrid(field.lon.values, field.lat.values)
        selected = criterion_passed.copy()
        selected[np.arange(selected.shape[0]) % 3 != 0, :] = False
        selected[:, np.arange(selected.shape[1]) % 3 != 0] = False
        dots = map_axis.scatter(longitude_grid[selected], latitude_grid[selected], marker='.', s=1.2,
                                color='#222222', linewidths=0, alpha=.7, transform=ccrs.PlateCarree(), zorder=4)
        significance_records.append({'kind': 'pe', 'product': product, 'n_displayed_dots': int(selected.sum()),
                                     'criterion': str(display_criterion.criterion.item()),
                                     'dot_offsets_sha256': hashlib.sha256(np.asarray(dots.get_offsets(), dtype=float).tobytes()).hexdigest(),
                                     'significant_mask_sha256': hashlib.sha256(criterion_passed.tobytes()).hexdigest(),
                                     'nonsignificant_mask_sha256': hashlib.sha256(criterion_failed.tobytes()).hexdigest()})
        map_axis.add_geometries(coast, ccrs.PlateCarree(), facecolor='none', edgecolor='#666666', linewidth=.4)
        gridliner = map_axis.gridlines(xlocs=[60, 120, 180, -120, -60], ylocs=[-60, -30, 0, 30, 60],
                                       draw_labels=True, x_inline=False, y_inline=False, rotate_labels=False,
                                       xlabel_style={'size': 7, 'color': '#555555'}, xpadding=3,
                                       linewidth=.4, color='#A0A0A0', alpha=.4, linestyle=':')
        gridliner.top_labels, gridliner.right_labels, gridliner.left_labels = False, False, False
        gridliner.bottom_labels, gridliner.geo_labels = True, False
        for latitude_tick, label in zip([-60, -30, 0, 30, 60], ['60°S', '30°S', '0°', '30°N', '60°N']):
            boundary_x, boundary_y = map_axis.projection.transform_point(central_longitude - 180, latitude_tick, ccrs.PlateCarree())
            annotation = map_axis.annotate(label, xy=(boundary_x, boundary_y), xycoords='data',
                                            xytext=(-3, 0), textcoords='offset points', fontsize=7,
                                            color='#555555', ha='right', va='center', annotation_clip=False)
            label_records.append({'product': product, 'latitude': latitude_tick, 'label': annotation.get_text()})
        map_axis.set_title(chr(98 + row * 2), loc='left', weight='bold', fontsize=11, pad=8)
        if product != 'contrast':
            sampled = peaks.loc[peaks['product'].eq(product) & peaks.plotted]
            early = sampled.loc[sampled.period.eq('early')].set_index('lon').sort_index()
            later = sampled.loc[sampled.period.eq('later')].set_index('lon').sort_index()
            for longitude in early.index:
                if early.loc[longitude, 'peak_lat'] != later.loc[longitude, 'peak_lat']:
                    line, = map_axis.plot([longitude, longitude], [early.loc[longitude, 'peak_lat'], later.loc[longitude, 'peak_lat']],
                                           color='#555555', lw=.55, alpha=.8, marker='None',
                                           transform=ccrs.PlateCarree(), zorder=5)
                    connector_records.append({'product': product, 'lon': float(longitude),
                                               'x': list(line.get_xdata()), 'y': list(line.get_ydata())})
            for period in ['early', 'later']:
                points = sampled.loc[sampled.period.eq(period)].sort_values('lon')
                markers = map_axis.scatter(points.lon, points.peak_lat, s=18 if period == 'early' else 8,
                                           marker='o' if period == 'early' else 'x',
                                           facecolors='white' if period == 'early' else phase_colors[period],
                                           edgecolors=phase_colors[period] if period == 'early' else None,
                                           color=phase_colors[period] if period == 'later' else None,
                                           linewidths=.6, transform=ccrs.PlateCarree(), zorder=6)
                marker_records.append({'product': product, 'period': period, 'n': len(points),
                                       'marker': 'o' if period == 'early' else 'x',
                                       'offsets_sha256': hashlib.sha256(np.asarray(markers.get_offsets(), dtype=float).tobytes()).hexdigest()})
        if row == 2:
            cax = figure.add_axes([.54, .070, .37, .017])
            colorbar = figure.colorbar(map_image, cax=cax, orientation='horizontal', extend='both', ticks=[-3, -1.5, 0, 1.5, 3],
                                       label=r'$\Delta(P-E)$ (mm day$^{-1}$)')
            colorbar.ax.minorticks_off()
            colorbar.solids.set_edgecolor('face')
    figure.canvas.draw()
    destination.mkdir(parents=True)
    for path in figure_paths:
        figure.savefig(path, dpi=dpi, facecolor='white')
        print(path.resolve(), flush=True)
    plt.close(figure)


# %% 5. 保存图注、来源与真实显示记录

caption = ('Atmospheric circulation and net precipitation responses across AMOC recovery groups. '
           'Rows show recovery (GFDL-ESM4, NorESM2-LM, CAS-ESM2-0), non-recovery '
           '(CESM2, CNRM-ESM2-1, UKESM1-0-LL), and non-recovery minus recovery. '
           'Each group is an equal-weight mean of three fixed models. '
           '(a, c) Shading shows annual zonal-mean meridional mass streamfunction changes; '
           'contours show initial circulation, solid positive and dashed negative. Positive and '
           'negative streamfunction denote clockwise and counterclockwise circulation, respectively, '
           'in the displayed latitude-pressure sections; units are 10^9 kg s^-1. '
           '(e) Difference between group-mean streamfunction changes. '
           '(b, d) Annual net precipitation changes, P-E, where E is upward surface evaporation, '
           'including sublimation and transpiration. Green/brown indicate increased/decreased '
           'net precipitation in mm day^-1. Positive changes indicate a shift toward a more positive '
           'surface freshwater balance, not necessarily positive absolute P-E. '
           '(f) Non-recovery minus recovery net precipitation changes. '
           'Changes are experiment years 271-280 minus 1-10. '
           'In (b, d), dots indicate that all three individual-model P-E changes have the same sign; '
           'cross-hatching indicates no unanimous sign agreement, not a failed statistical test. '
           'In (e, f), dots indicate bootstrap-FDR significance and cross-hatching indicates failure: '
           '9999 centered two-sided circular block draws, maximum p across 2/3/4-year blocks, BH q<0.05. '
           'Dots are subsampled for clarity. Circulation fields and their masks are unchanged from '
           'the precipitation version; P-E contrast significance is recomputed using P-E samples. '
           'Blue open circles and magenta crosses in (b, d) mark the initial and later latitude of '
           'maximum group-mean absolute P-E within 20S-20N, at longitudes sampled every 10 degrees. '
           'Thin grey connectors show latitude differences at the same longitude. These markers '
           'indicate net-precipitation rainband positions, not significance; they differ from the '
           'precipitation-maximum or median ITCZ indices. Initial and later periods are years1-10 '
           'and271-280, respectively.')
input_paths = [flow_source / name for name in ['group_circulation_climatology.nc', 'group_circulation_changes.nc', 'psi_significance.nc', 'metadata.json']]
input_paths += [hydrology_source / 'group_changes.nc', hydrology_source / 'pe_display_criteria.nc',
                hydrology_source / 'metadata.json', peak_source / 'net_precipitation_peak_gridpoints.csv',
                peak_source / 'metadata.json', coast_path, base / 'check_red_white_blue_style/red_white_blue_rgb.npy']
manifest = {'panel_count': 6, 'figure_size_inches': list(figure_size), 'dpi': dpi, 'overall_title': None,
            'groups': metadata['groups'], 'group_sizes': [3, 3], 'contrast': 'non-recovery minus recovery',
            'artist_records': artist_records, 'significance_records': significance_records,
            'geographic_labels': label_records, 'peak_markers': True,
            'marker_records': marker_records, 'connector_records': connector_records,
            'peak_definition': peak_metadata['peak_definition'], 'marker_spacing_deg': 10,
            'significance': 'b/d: all3 strictly same-sign P-E change; e/f: max-L234-block-bootstrap BH q<0.05. a/c have no test overlay.',
            'palette_limits': {'psi': [-40, 40], 'pe': [-3, 3]}, 'caption': caption,
            'source_files': [{'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in input_paths]}
(destination / 'figure_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
(destination / 'caption.txt').write_text(caption + '\n', encoding='utf-8')
