"""Plot frozen calibration summaries, spread across seeds is NOT a confidence band.
Usage: plot_calib.py RESULT_DIRECTORY
"""
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

root = Path(sys.argv[1])
files = sorted(root.glob('*.curve.tsv'))
fig, axes = plt.subplots(1, len(files), figsize=(5 * len(files), 4.5), squeeze=False)
for ax, file in zip(axes[0], files):
    frame = pd.read_csv(file, sep='\t')
    for arm, color in [('never', 'tab:blue'), ('always', 'tab:orange')]:
        group = frame[frame.arm == arm].sort_values('alpha')
        x = group.alpha.to_numpy() * 100
        for method, style in [('lower', '--'), ('combined', '-')]:
            ax.plot(x, group[f'{method}_median'] * 100, style, color=color, label=f'{arm}: {method}')
            ax.fill_between(x, group[f'{method}_min'].to_numpy()*100,
                            group[f'{method}_max'].to_numpy()*100, color=color, alpha=.08)
    ax.plot([.1, 5], [.1, 5], ':', color='black', label='identity')
    ax.set(title=file.stem.replace('.curve', ''), xlabel='Reported q threshold (%)', ylabel='Estimated FDP (%)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(alpha=.2)
axes[0, 0].legend(fontsize=7)
fig.suptitle('Median and min/max across five fits of the same raw, not confidence intervals')
fig.tight_layout()
fig.savefig(root / 'calibration_curves.png', dpi=160)
fig.savefig(root / 'calibration_curves.svg')
