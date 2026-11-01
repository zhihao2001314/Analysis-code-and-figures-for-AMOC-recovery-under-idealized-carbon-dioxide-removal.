# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory


# %% 2. 参数与最终图件路径

root = Path(__file__).resolve().parents[2]
source = root / 'output/s01_amoc/n77_seven_panel_hysteresis'
destination = root / 'figure/s01_amoc/n78_seven_panel_hysteresis'
png_path = destination / 'amoc_hysteresis_seven_panels_26p5.png'
pdf_path = destination / 'amoc_hysteresis_seven_panels_26p5.pdf'
svg_path = destination / 'amoc_hysteresis_seven_panels_26p5.svg'
figure_paths = [png_path, pdf_path, svg_path]
assert not destination.exists() or not any(destination.iterdir()), '正式图件存在，须先归档。'
figure_size, dpi = (9.6, 12.6), 450
overview_bounds = [.16, .735, .68, .235]
lower_y = [.505, .285, .065]
schematic_left, model_left, lower_width, lower_height = .095, .555, .365, .170
schematic_color, annotation_color = '#364F60', '#526977'
forcing_color, forcing_width = '#73818C', 1.15
overview_band_alpha, overview_legend_alpha = .18, .15
artist_records, arrow_records, axis_records, annotation_boxes = [], [], [], []
model_legend_artists = []


# %% 3. 仅读取已保存的绘图坐标

metadata = json.loads((source / 'metadata.json').read_text(encoding='utf-8'))
overview = pd.read_csv(source / 'overview_anomaly.csv', float_precision='round_trip')
forcing = pd.read_csv(source / 'co2_forcing.csv', float_precision='round_trip')
curves = pd.read_csv(source / 'model_paths.csv', float_precision='round_trip')
schematics = pd.read_csv(source / 'schematic_paths.csv', float_precision='round_trip')
model_arrows = pd.read_csv(source / 'model_arrows.csv', float_precision='round_trip')
schematic_arrows = pd.read_csv(source / 'schematic_arrows.csv', float_precision='round_trip')
panels = pd.read_csv(source / 'panels.csv', keep_default_na=False)
area_annotations = pd.read_csv(source / 'area_annotations.csv', float_precision='round_trip')
legend_values = pd.read_csv(source / 'initial_legend_values.csv').set_index('model')


# %% 4. 绘制适中宽度的九模式时间序列

with plt.rc_context({'font.family': 'Arial', 'font.size': 9, 'text.color': '#25313A',
                     'axes.labelcolor': '#25313A', 'xtick.color': '#56616B', 'ytick.color': '#56616B',
                     'axes.linewidth': .7, 'axes.spines.top': False, 'axes.spines.right': False,
                     'pdf.fonttype': 42, 'svg.fonttype': 'none', 'savefig.facecolor': 'white'}):
    figure = plt.figure(figsize=figure_size, facecolor='white')
    axis = figure.add_axes(overview_bounds)
    overview_band = axis.fill_between(overview.simulation_year, overview.lower, overview.upper,
                                       color='#C7CDD2', alpha=overview_band_alpha, lw=0, zorder=1)
    for model in metadata['models']:
        line, = axis.plot(overview.simulation_year, overview[model], color=metadata['model_colors'][model],
                          lw=1.15, label=legend_values.loc[model, 'label'], zorder=3)
        artist_records.append({'panel': 'a', 'model': model, 'color': line.get_color(),
                               'x_sha256': hashlib.sha256(np.asarray(line.get_xdata(), dtype=float).tobytes()).hexdigest(),
                               'y_sha256': hashlib.sha256(np.asarray(line.get_ydata(), dtype=float).tobytes()).hexdigest()})
    line, = axis.plot(overview.simulation_year, overview.mmm, color='#27323A', lw=2., label=legend_values.loc['mmm', 'label'], zorder=4)
    artist_records.append({'panel': 'a', 'model': 'mmm', 'color': line.get_color(),
                           'x_sha256': hashlib.sha256(np.asarray(line.get_xdata(), dtype=float).tobytes()).hexdigest(),
                           'y_sha256': hashlib.sha256(np.asarray(line.get_ydata(), dtype=float).tobytes()).hexdigest()})
    axis.axhline(0, color='#A5AFB7', lw=.6, zorder=2)
    axis.axvline(140, color='#BEC6CC', lw=.7, zorder=2)
    axis.set(xlim=(1, 280), xticks=[1, 70, 140, 210, 280],
             ylim=metadata['overview_y_limits_sv'], yticks=metadata['overview_y_ticks_sv'],
             xlabel='Year', ylabel='AMOC anomaly (Sv)')
    axis.grid(axis='y', color='#EDF0F2', lw=.5)
    axis.tick_params(labelsize=8, length=3)
    forcing_axis = axis.twinx()
    # 保留细虚线；淡化灰带与图例背景，使强迫线可见而不增加线宽。
    forcing_line, = axis.plot(forcing.simulation_year, forcing.co2_ppm, color=forcing_color, lw=forcing_width,
                              ls=(0, (4, 2.5)), zorder=6,
                              transform=blended_transform_factory(axis.transData, forcing_axis.transData))
    forcing_axis.set(ylim=(250, 1200), yticks=[300, 600, 900, 1200])
    forcing_axis.tick_params(labelsize=7.5, colors='#66727B', length=3)
    forcing_axis.spines['right'].set_visible(True)
    forcing_axis.spines['right'].set_color('#A5AFB7')
    forcing_axis.text(.74, .96, 'CO$_2$ forcing (ppm)', transform=forcing_axis.transAxes,
                      ha='center', va='top', fontsize=7.7, color=forcing_color)
    axis.set_zorder(forcing_axis.get_zorder() + 1)
    axis.patch.set_visible(False)
    overview_legend = axis.legend(loc='lower left', bbox_to_anchor=(.017, .025), ncol=1, frameon=True,
                                  facecolor='white', edgecolor='none', framealpha=overview_legend_alpha, fontsize=6.8,
                                  title='Model (initial AMOC, Sv)', title_fontsize=7.0,
                                  handlelength=1.6, columnspacing=1., handletextpad=.5, labelspacing=.03,
                                  borderaxespad=0, borderpad=0)
    overview_legend.set_zorder(10)
    overview_axis = axis
    forcing_report = {'x_sha256': hashlib.sha256(np.asarray(forcing_line.get_xdata(), dtype=float).tobytes()).hexdigest(),
                      'y_sha256': hashlib.sha256(np.asarray(forcing_line.get_ydata(), dtype=float).tobytes()).hexdigest(),
                      'line_zorder': forcing_line.get_zorder(), 'band_zorder': overview_band.get_zorder(),
                      'line_color': forcing_line.get_color(), 'line_width': forcing_line.get_linewidth(),
                      'white_halo': False, 'legend_background_alpha': overview_legend.get_frame().get_alpha(),
                      'band_alpha': overview_band.get_alpha(), 'rendered_in_main_axes': forcing_line.axes is axis}
    axis.text(-.08, 1.025, 'a', transform=axis.transAxes, fontsize=12, weight='bold')
    axis_records.append({'panel': 'a', 'bounds': list(axis.get_position().bounds)})


# %% 5. 三行示意图与组内模式叠加

    for record in panels.loc[panels.panel.ne('a')].itertuples(index=False):
        left = schematic_left if record.kind == 'schematic' else model_left
        axis = figure.add_axes([left, lower_y[record.row - 1], lower_width, lower_height])
        axis.set(xlim=(250, 1200), ylim=(0, 1.2), xticks=[300, 600, 900, 1200], yticks=[0, .4, .8, 1.2])
        axis.tick_params(labelsize=8, length=3, pad=4)
        axis.spines[['left', 'bottom']].set_color('#7C8993')
        axis.axhline(1., color='#DBE1E6', lw=.6, zorder=0)
        axis.text(-.13, 1.035, record.panel, transform=axis.transAxes, fontsize=12, weight='bold')
        if record.kind == 'schematic':
            branches = {}
            for phase in ['increase', 'decrease']:
                part = schematics.loc[schematics.panel.eq(record.panel) & schematics.phase.eq(phase)].sort_values('node')
                branches[phase] = part
                line, = axis.plot(part.co2_ppm, part.relative_amoc, color=schematic_color, lw=1.7, ls='-', zorder=3)
                artist_records.append({'panel': record.panel, 'phase': phase, 'color': line.get_color(),
                                       'x_sha256': hashlib.sha256(np.asarray(line.get_xdata(), dtype=float).tobytes()).hexdigest(),
                                       'y_sha256': hashlib.sha256(np.asarray(line.get_ydata(), dtype=float).tobytes()).hexdigest()})
            axis.fill_between(branches['increase'].co2_ppm, branches['increase'].relative_amoc,
                               branches['decrease'].relative_amoc, facecolor='#E3EBF0', lw=0, zorder=1)
            axis.plot(metadata['reference_co2_ppm'], 1., marker='o', ls='none', ms=4.3,
                       mfc='white', mec=schematic_color, mew=1., zorder=5)
            if record.panel != 'b':
                axis.plot(metadata['reference_co2_ppm'], branches['decrease'].relative_amoc.iloc[0],
                           marker='o', ls='none', ms=3.8, mfc=schematic_color, mec='white', mew=.5, zorder=5)
            for arrow in schematic_arrows.loc[schematic_arrows.panel.eq(record.panel)].itertuples(index=False):
                axis.annotate('', xy=(arrow.end_x, arrow.end_y), xytext=(arrow.start_x, arrow.start_y),
                              arrowprops={'arrowstyle': '-|>', 'color': schematic_color, 'lw': 1., 'mutation_scale': 9}, zorder=5)
            axis.text(.98, .96, 'Schematic', transform=axis.transAxes, ha='right', va='top', fontsize=8, color='#6B7A85')
            label = record.label.replace(' A', r' $A$')
            axis.text(.98, .035, label, transform=axis.transAxes, ha='right', va='bottom', fontsize=8.5)
            axis.set_ylabel('Relative AMOC strength', fontsize=9)
            if record.panel == 'd':
                annotation = area_annotations.loc[area_annotations.panel.eq(record.panel)].iloc[0]
                artist = axis.annotate('Hysteresis\narea ($A$)',
                                        xy=(annotation.target_co2_ppm, annotation.target_relative_amoc), xycoords='data',
                                        xytext=(annotation.text_axes_x, annotation.text_axes_y), textcoords='axes fraction',
                                        ha='center', va='center', fontsize=8.5, color=annotation_color,
                                        arrowprops={'arrowstyle': '-|>', 'color': annotation_color, 'lw': .85,
                                                    'shrinkA': 6, 'shrinkB': 3, 'mutation_scale': 9},
                                        zorder=6)
                annotation_boxes.append((record.panel, artist))
        else:
            for model in metadata['members'][record.row - 1]:
                part = curves.loc[curves.model.eq(model)].sort_values('simulation_year')
                line, = axis.plot(part.co2_ppm, part.relative_amoc, color=metadata['model_colors'][model],
                                  lw=1.65, ls='-', label=model, solid_capstyle='round', zorder=3)
                artist_records.append({'panel': record.panel, 'model': model, 'color': line.get_color(),
                                       'x_sha256': hashlib.sha256(np.asarray(line.get_xdata(), dtype=float).tobytes()).hexdigest(),
                                       'y_sha256': hashlib.sha256(np.asarray(line.get_ydata(), dtype=float).tobytes()).hexdigest()})
                for arrow in model_arrows.loc[model_arrows.model.eq(model)].itertuples(index=False):
                    axis.annotate('', xy=(arrow.end_x, arrow.end_y), xytext=(arrow.start_x, arrow.start_y),
                                  arrowprops={'arrowstyle': '-|>', 'color': metadata['model_colors'][model],
                                              'lw': .9, 'mutation_scale': 8}, zorder=5)
                    arrow_records.append({'panel': record.panel, 'model': model, 'color': metadata['model_colors'][model],
                                          'start': [arrow.start_x, arrow.start_y], 'end': [arrow.end_x, arrow.end_y]})
            legend_position = 'upper right' if record.row == 3 else 'lower left'
            legend_anchor = (.985, .97) if record.row == 3 else (.015, .035)
            model_legend = axis.legend(loc=legend_position, bbox_to_anchor=legend_anchor, ncol=2, frameon=False,
                                       fontsize=7.2, handlelength=1.7, columnspacing=1., labelspacing=.5, borderaxespad=0)
            model_legend_artists.append((record.panel, axis, model_legend))
        axis.set_xlabel('Atmospheric CO$_2$ (ppm)', fontsize=9, labelpad=5)
        axis_records.append({'panel': record.panel, 'bounds': list(axis.get_position().bounds)})
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    model_legend_records = []
    for panel, legend_axis, model_legend in model_legend_artists:
        model_legend_box = model_legend.get_window_extent(renderer)
        overlapping = []
        for line in legend_axis.lines:
            if line.get_label() in metadata['models']:
                if line.get_path().transformed(line.get_transform()).intersects_bbox(model_legend_box, filled=False):
                    overlapping.append(line.get_label())
        assert not overlapping, (panel, 'model legend overlaps data', overlapping)
        model_legend_records.append({'panel': panel, 'labels': [text.get_text() for text in model_legend.get_texts()],
                                     'overlapping_models': overlapping})
    legend_box = overview_legend.get_window_extent(renderer)
    intersecting_models = []
    for line in overview_axis.lines:
        if line.get_label() in legend_values.label.tolist():
            if line.get_path().transformed(line.get_transform()).intersects_bbox(legend_box, filled=False):
                intersecting_models.append(line.get_label())
    assert not intersecting_models, ('Initial AMOC legend overlaps scientific curves', intersecting_models)
    legend_report = {'title': overview_legend.get_title().get_text(),
                     'labels': [item.get_text() for item in overview_legend.get_texts()],
                     'bounds_figure': list(legend_box.transformed(figure.transFigure.inverted()).bounds),
                     'overlapping_models': intersecting_models}
    annotation_records = [{'panel': panel, 'bbox_figure': list(artist.get_window_extent(renderer).transformed(figure.transFigure.inverted()).bounds)}
                          for panel, artist in annotation_boxes]
    destination.mkdir(parents=True, exist_ok=True)
    for path in figure_paths:
        figure.savefig(path, dpi=dpi, facecolor='white')
        print(path.resolve(), flush=True)
    plt.close(figure)


# %% 6. 保存实际绘图记录

manifest = {'panel_count': 7, 'figure_size_inches': list(figure_size), 'dpi': dpi, 'overall_title': None,
            'overview_target_latitude_degN': 26.5, 'pathway_target_latitude_degN': 26.5,
            'group_names': metadata['group_names'], 'model_counts': [2, 3, 4],
            'model_color_policy': 'Same panel-a model color throughout the complete trajectory; solid paths, no phase color legend.',
            'artist_records': artist_records, 'model_arrow_records': arrow_records, 'axes': axis_records,
            'forcing_display': forcing_report,
            'initial_legend': legend_report, 'model_legends': model_legend_records,
            'area_annotations': annotation_records, 'caption': metadata['caption'],
            'source_files': [{'path': str(source / name), 'sha256': hashlib.sha256((source / name).read_bytes()).hexdigest()}
                             for name in ['overview_anomaly.csv', 'co2_forcing.csv', 'model_paths.csv', 'schematic_paths.csv',
                                          'model_arrows.csv', 'schematic_arrows.csv', 'panels.csv', 'area_annotations.csv', 'initial_legend_values.csv', 'metadata.json']],
            'alt_text': 'Seven panels at26.5N. A moderately wide top panel shows nine model AMOC anomalies and CO2 forcing. Three rows below pair a neutral-color schematic with2,3,4overlaid model loops. Actual model colors match panel a, with one solid color throughout both forcing phases. Only the middle schematic d labels hysteresis area with an arrow; f has no area label. Panel-a SD shading is semitransparent and the thin gray dashed CO2 line retains its original width without a halo. The SD band and legend background are translucent (alpha0.18 and0.15). Legend numbers give raw experiment-years1-10 mean AMOC strength.'}
(destination / 'figure_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
