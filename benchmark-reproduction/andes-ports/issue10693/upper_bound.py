"""Expected precursors (identified in 1-2 runs) reaching >= 3 runs if every precursor x run pair with real MS1
signal got a candidate that is then accepted like today's candidates. p_signal: share of no-candidate pairs with
real co-eluting signal (from xic_check.py, background- and sensitivity-corrected); p_accept: share of today's
candidates that are accepted. Monte Carlo over the pairs. Usage: upper_bound.py cand_pairs.tsv.gz p_signal p_accept"""
import sys, numpy as np, pandas as pd
g = pd.read_csv(sys.argv[1], sep='\t'); ps, pa = float(sys.argv[2]), float(sys.argv[3])
g = g[g.n_ident <= 2]
rng = np.random.default_rng(1)
acc = (g.status == 'accepted').values
cand = (g.status == 'candidate').values
nocand = g.status.isin(['no_feature', 'outside_rt']).values
res = []
for _ in range(50):
    gain = acc | (cand & (rng.random(len(g)) < pa)) | (nocand & (rng.random(len(g)) < ps * pa))
    s = pd.Series(gain).groupby(g.key.values).sum()
    n = g.groupby('key').n_ident.first()
    res.append(int(((n + s.reindex(n.index)) >= 3).sum()))
now = int(((g.groupby('key').n_ident.first() + pd.Series(acc).groupby(g.key.values).sum()) >= 3).sum())
print(f'precursors identified in 1-2 runs: {g.key.nunique()}; reach >= 3 today: {now}; '
      f'expected with candidates at every real signal: {np.mean(res):.0f} (+{np.mean(res) - now:.0f})')
