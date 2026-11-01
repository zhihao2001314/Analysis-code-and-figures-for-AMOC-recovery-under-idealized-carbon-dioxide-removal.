# %% 1. 导入依赖

from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# %% 2. 参数与最终图件完整路径

root = Path(__file__).resolve().parents[2]
source = root / 'output/s03_tem/n17_nine_model_global_tas_series'
destination = root / 'figure/s03_tem/n18_nine_model_global_tas_series'
png_path = destination / 'nine_model_global_tas_280years.png'
pdf_path = destination / 'nine_model_global_tas_280years.pdf'
svg_path = destination / 'nine_model_global_tas_280years.svg'
figure_paths = [png_path, pdf_path, svg_path]
assert not destination.exists(), '正式图件存在，须先归档。'
figure_size, dpi = (9.4, 5.8), 450
artist_records = []


# %% 3. 仅读取已保存的年平均异常与配色

annual = pd.read_csv(source / 'global_annual_tas.csv', float_precision='round_trip')
metadata = json.loads((source / 'metadata.json').read_text(encoding='utf-8'))


# %% 4. 九模式全球气温演变图

with plt.rc_context({'font.family': 'Arial', 'font.size': 11, 'axes.spines.top': False,
                     'axes.spines.right': False, 'axes.linewidth': .8,
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'}):
    figure, axis = plt.subplots(figsize=figure_size)
    figure.subplots_adjust(left=.105, right=.975, bottom=.12, top=.97)
    for model in metadata['models']:
        selected = annual.loc[annual.model.eq(model)].sort_values('simulation_year')
        line, = axis.plot(selected.simulation_year, selected.tas_anomaly_smooth11_degC,
                          color=metadata['model_colors'][model], lw=1.65, ls='-', label=model)
        artist_records.append({'model': model, 'line_style': line.get_linestyle(),
            'color': line.get_color(), 'points': len(selected),
            'x_sha256': hashlib.sha256(np.asarray(line.get_xdata(), dtype=float).tobytes()).hexdigest(),
            'y_sha256': hashlib.sha256(np.asarray(line.get_ydata(), dtype=float).tobytes()).hexdigest()})
    axis.axhline(0., color='#A6AEB5', lw=.65, zorder=0)
    axis.axvline(140.5, color='#9BA6AF', lw=.8, ls=':', zorder=0)
    axis.set(xlim=(1, 280), xticks=[1, 70, 140, 210, 280],
             ylim=metadata['shared_y_limits_degC'], xlabel='Year',
             ylabel='Global near-surface air-temperature anomaly (°C)')
    axis.tick_params(length=4, width=.8)
    axis.grid(False)
    legend = axis.legend(loc='upper right', bbox_to_anchor=(.99, .985), ncol=1,
                         frameon=False, fontsize=9, handlelength=2.1,
                         columnspacing=1.5, labelspacing=.55)
    figure.canvas.draw()
    destination.mkdir(parents=True)
    for path in figure_paths:
        figure.savefig(path, dpi=dpi, facecolor='white')
        print(path.resolve(), flush=True)
    plt.close(figure)


# %% 5. 保存画线与来源记录

manifest = {'panel_count': 1, 'model_line_count': 9, 'figure_size_inches': list(figure_size),
    'dpi': dpi, 'overall_title': None, 'artist_records': artist_records,
    'legend_labels': [item.get_text() for item in legend.get_texts()], 'caption': metadata['caption'],
    'reference': metadata['reference'], 'smoothing': metadata['smoothing'],
    'source_files': [{'path': str(source / name), 'sha256': hashlib.sha256((source / name).read_bytes()).hexdigest()}
                     for name in ['global_annual_tas.csv', 'metadata.json']],
    'exports': [str(path.resolve()) for path in figure_paths]}
(destination / 'figure_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
(destination / 'caption.txt').write_text(metadata['caption'], encoding='utf-8')
