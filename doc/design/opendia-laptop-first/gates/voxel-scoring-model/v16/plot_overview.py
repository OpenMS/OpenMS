"""V16 overview figure from frozen threshold audit outputs only."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
bands = pd.read_csv(OUT/'threshold_bands.tsv', sep='\t')
churn = pd.read_csv(OUT/'historical_churn.tsv', sep='\t')
labels = ['(-0.001, 0.01]', '(0.01, 0.05]', '(0.05, 0.15]', '(0.15, 0.5]', '(0.5, 1.0]', 'lost_pair_competition']
short = ['Already q<=1%', 'q 1-5%', 'q 5-15%', 'q 15-50%', 'q 50-100%', 'Lost to paired decoy']
colors = ['#2a9d8f', '#457b9d', '#e9c46a', '#f4a261', '#e76f51', '#9d4edd']
names = {'D11FWT':'Base', 'D11FPWT':'Mass (V5)', 'D11FXWT':'Cross-run (V7)', 'D11FNWT':'Native (V9)',
         'D11FEWT':'Emp. null (V10)', 'D11FUWT':'Unique (V11)', 'D11FIWT':'Corrected (V12)'}
fig, axes = plt.subplots(2, 2, figsize=(13.5, 9), layout='constrained')
for column, raw in enumerate(['raw01','raw02']):
    data = bands[bands.raw == raw].set_index('band').reindex(labels)
    ax = axes[0, column]
    bars = ax.barh(np.arange(6), data.n, color=colors)
    ax.set_yticks(np.arange(6), short)
    ax.invert_yaxis()
    ax.bar_label(bars, padding=4)
    ax.set_xlim(0, 430)
    ax.set_xlabel('Historical Stage D reference IDs')
    ax.set_title(f'{raw}: Stage D is not one failure type')
    ax = axes[1, column]
    tab = churn[churn.raw == raw]
    for _, row in tab.iterrows():
        ax.scatter(row.D_new_vs_base, 100*row.empirical_fdp, s=70 if row.arm == 'D11FWT' else 42,
                   color='#264653' if row.arm == 'D11FWT' else '#e76f51')
        ax.annotate(names[row.arm], (row.D_new_vs_base, 100*row.empirical_fdp), xytext=(4,5),
                    textcoords='offset points', fontsize=8)
    ax.axhline(1.5, ls='--', color='#2a9d8f', label='Existing FDP_c gate')
    ax.set_xlabel('New historical Stage D IDs at q<=1%, versus base')
    ax.set_ylabel('Entrapment FDP_c at q<=1% (%)')
    ax.set_title(f'{raw}: prior recovery often trades against specificity')
    ax.set_xlim(-15, 320)
    ax.set_ylim(0.7, 3.6)
    ax.legend(loc='upper left', fontsize=8)
fig.suptitle('V16 diagnosis, unchanged D11FWT and historical experiments\nStage D is defined at empirical FDP_c 1%, not operational q 1%', fontsize=14)
fig.savefig(OUT/'stageD_overview.png', dpi=160)
fig.savefig(OUT/'stageD_overview.pdf')
