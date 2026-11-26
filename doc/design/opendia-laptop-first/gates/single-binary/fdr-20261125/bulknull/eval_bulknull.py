#!/usr/bin/env python3
"""SCREEN_BULKNULL evaluator. Entrapment labels are read here only, from a key table exported read-only from the library
with the DIAG_STEPA label correction (exact human sequences relabelled human, as rtnull/eval_pdec_relabel.py).
Usage: eval_bulknull.py <keys_relabel.tsv: id traml decoy e> <set> <replay.tsv> ...   (bulknull_replay outputs, one per seed)
Per run and q column (tool_q = default plain TDC as exported, bn_q = bulk-null): at q <= 0.01, N, E, lower E/N, combined
E(1+1/r)/N. The library ratio r = entrapment / human targets after relabelling. N@FDP_c 1 % is the same for both columns because
the ranking is identical. It is computed on the exported score order, as in the gate evaluator. Pooled over seeds: combined with
a 90 % Poisson interval on E."""
import sys, numpy as np, pandas as pd
from scipy.stats import chi2
K = pd.read_csv(sys.argv[1], sep='\t', header=None, names=['id', 't', 'd', 'e'], dtype={'t': str})
tg = K[K.d == 0]; r = (tg.e == 1).sum() / (tg.e == 0).sum(); k = 1 + 1 / r
ent = tg.set_index('id').e
print(f'library r = {r:.4f}, k = {k:.4f}')
tot = {}
for fp in sys.argv[3:]:
    t = pd.read_csv(fp, sep='\t')
    t = t[t.is_decoy == 0].copy(); t['e'] = t.precursor_id.map(ent).fillna(0).values
    o = t.sort_values(['tool_q', 'score'], ascending=[True, False]).e.values
    fdp = np.cumsum(o) * k / np.arange(1, len(o) + 1); ok = np.where(fdp <= 0.01)[0]
    n1 = int(ok.max() + 1) if len(ok) else 0
    for col in ['tool_q', 'bn_q']:
        a = t[t[col] <= 0.01]; N, E = len(a), int(a.e.sum())
        g = tot.setdefault(col, [0, 0, []]); g[0] += N; g[1] += E; g[2].append(n1)
        print(f'{sys.argv[2]}\t{fp.split("/")[-1]}\t{col}\tN {N}\tE {E}\tlower {E / N:.4f}\tcombined {E * k / N:.4f}\tN@FDPc1% {n1}')
for col, (N, E, n1) in tot.items():
    lo, hi = chi2.ppf(0.05, 2 * E) / 2 if E else 0.0, chi2.ppf(0.95, 2 * E + 2) / 2
    print(f'{sys.argv[2]}\tPOOLED\t{col}\tN {N}\tE {E}\tlower {E / N:.4f}\tcombined {E * k / N:.4f} (90 % {lo * k / N:.4f}-{hi * k / N:.4f})'
          f'\tN@FDPc1% sum {sum(n1)}')
