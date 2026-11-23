#!/usr/bin/env python3
"""Plot descriptive full-library calibration, without confidence intervals."""
import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

root = Path(sys.argv[1])
frame = pd.read_csv(root / 'calibration_curve.tsv', sep='\t', dtype={'raw': str})
fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharex=True, sharey=True)
for ax, raw in zip(axes, ['01', '02', '03']):
    rows = frame[frame.raw.str.zfill(2) == raw].sort_values('alpha')
    for metric, style in [('lower', '--'), ('combined', '-')]:
        x = rows.alpha * 100
        ax.plot(x, rows[metric+'_median'] * 100, style, label=metric)
        ax.fill_between(x, rows[metric+'_min'] * 100, rows[metric+'_max'] * 100, alpha=.15)
    ax.plot([.1, 5], [.1, 5], 'k:', label='identity')
    ax.set(xscale='log', yscale='log', xlabel='Reported q threshold (%)', title='RAW '+raw)
    ax.grid(alpha=.2)
axes[0].set_ylabel('Entrapment FDP estimate (%)')
axes[0].legend()
fig.suptitle('Median and min/max across scoring seeds, not confidence intervals')
fig.tight_layout()
fig.savefig(root/'calibration_curves.png', dpi=160)
fig.savefig(root/'calibration_curves.svg')
