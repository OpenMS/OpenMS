"""Raw01-only pre-plan descriptive diagnostic for interference-aware profile features.
Labels classify diagnostics only. D11FWT scores define the score-matched comparison strata; no model is trained.
"""
import sys
import numpy as np
import pandas as pd
sys.path.insert(0, 'v4')
from dd_band import auc, ent
from interference_profiles import FEATURES

pop = pd.read_csv('raw01/population.tsv', sep='\t')
FI = np.load('raw01/interference_features.npy', mmap_mode='r')
s = np.load('raw01/scores_D11FWT.npy').astype(float)
w = pop.assign(s=s, row=np.arange(len(pop))).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id')
w = w[w.s > -1e20].copy()
st = pd.read_csv('v3seq/misses_D11FWT_raw01.tsv', sep='\t').set_index('precursor_id').stage
w['grp'] = np.where(w.precursor_id.map(st) == 'D right window, below cut', 'D',
                    np.where(w.is_decoy == 1, 'DEC', np.where(w.precursor_id.map(ent) == 1, 'ENT', 'HUM')))
for i, name in enumerate(FEATURES):
    w[name] = FI[w.row.values, i]
# Use the same high-score region and D-quantile strata as the prior E/U pre-plan looks.
cut = w[w.grp == 'DEC'].s.quantile(0.9)
top = w[w.s >= cut].copy()
print('raw01 winners above decoy-score q90:', top.grp.value_counts().to_dict())
print('feature medians by group:')
print(top.groupby('grp')[FEATURES].median().round(3).to_string())
for ref, title in (('D', 'D vs decoy; entrapment vs decoy'), ('HUM', 'all human vs decoy; entrapment vs decoy')):
    g = top[top.grp.isin([ref, 'DEC', 'ENT'])].copy()
    edges = np.unique(np.r_[-np.inf, g[g.grp == ref].s.quantile(np.linspace(0.1, 0.9, 9)).values, np.inf])
    g['stratum'] = pd.cut(g.s, edges, labels=False)
    print('==', title)
    for feature in FEATURES:
        aa, ee, weights = [], [], []
        for _, group in g.groupby('stratum', observed=True):
            x = group[group.grp == ref][feature]
            d = group[group.grp == 'DEC'][feature]
            e = group[group.grp == 'ENT'][feature]
            if min(len(x), len(d), len(e)) >= 5:
                aa.append(auc(x, d)); ee.append(auc(e, d)); weights.append(len(x))
        print(f'{feature:26s} human/D={np.average(aa, weights=weights):.3f} ENT={np.average(ee, weights=weights):.3f}')
