"""PLAN.md metrics for the raw 01 search+boundary diagnostic: counts at pooled paired q, fold weights, null-band symmetry."""
import sys, re
import numpy as np, pandas as pd
sys.path.insert(0, '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-design-raw01-20261107/cnn-topaz-design')
from evaluate import paired_pooled  # same pooled paired TDC as the tool (reproduces 717 on search)
import sqlite3

D = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-design-raw01-20261107'
LIB = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-entrapment-20261107/lib/library.pqp'
c = sqlite3.connect('file:' + LIB + '?mode=ro', uri=True)
p = pd.read_sql('select ID precursor_id, TRAML_ID, DECOY from PRECURSOR', c)
p['pair_key'] = np.where(p.DECOY == 1, p.TRAML_ID.str.replace('^DECOY_', '', regex=True), p.TRAML_ID)
pairs = p[['precursor_id', 'pair_key']]
for name, path in [('search (frozen)', D + '/stage2-search-design/search'), ('search+boundary', D + '/stage2-boundary-design/boundary')]:
    df = pd.read_csv(path + '/crossfit.tsv', sep='\t').merge(pairs, on='precursor_id', how='left')
    best, _ = paired_pooled(df, 'stage2_score')
    t = best[best.is_decoy == 0]
    print(f'{name:18s}', ' '.join(f'q<={q}: {int((t.q <= q).sum())}' for q in (0.005, 0.01, 0.02, 0.05)))
log = open(D + '/stage2-boundary-design/boundary/stdout.log').read()
for line in re.findall(r'second stage fold=.*', log):
    print(' ', line[:400])
f = pd.read_csv(D + '/stage2-boundary-design/boundary/features.tsv', sep='\t')
cr = pd.read_csv(D + '/stage2-boundary-design/boundary/crossfit.tsv', sep='\t', usecols=['sample_id', 'precursor_best'])
w = f.merge(cr, on='sample_id')
w = w[w.precursor_best == 1]
cols = [c for c in f.columns if c.startswith('boundary_')]
for lo, hi in [(-99, 0), (4, 6), (6, 99)]:
    b = w[(w.cnn_z > lo) & (w.cnn_z <= hi)]
    med = b.groupby('is_decoy')[cols].median().round(3)
    print(f'\nCNN z in ({lo},{hi}] winners: targets {int((b.is_decoy == 0).sum())}, decoys {int((b.is_decoy == 1).sum())} (median target / decoy)')
    for col in cols:
        print(f'  {col:30s} {med.loc[0, col]:8.3f} {med.loc[1, col]:8.3f}')
