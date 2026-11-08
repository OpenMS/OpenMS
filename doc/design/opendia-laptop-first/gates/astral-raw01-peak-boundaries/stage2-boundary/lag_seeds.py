"""Seed stability of the MLP combiner result in lag.py (search+boundary only; same code path, seeds 1-3)."""
import runpy, sys
import numpy as np, torch
src = open('lag.py').read().split('res = {}')[0]
g = {}
exec(src, g)
for seed in (1, 2, 3):
    torch.manual_seed(seed); np.random.seed(seed)
    sc = g['crossfit'](g['SEARCH'] + g['BOUND'], 'mlp')
    tmp = g['df'][['sample_id', 'precursor_id', 'is_decoy', 'fold', 'in_rt_window', 'pair_key']].assign(score=sc)
    best, w = g['paired_pooled'](tmp, 'score')
    t = best[best.is_decoy == 0]
    print(f'seed {seed}: q<=0.005/0.01/0.02/0.05 =', [int((t.q <= q).sum()) for q in (0.005, 0.01, 0.02, 0.05)],
          '| OSW IDs among 1%', int(t[t.q <= 0.01].precursor_id.isin(g['o'].precursor_id).sum()))
