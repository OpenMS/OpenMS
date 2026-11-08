"""Does averaging MLP seeds (per fold, rank-free raw logits) stabilise the 1% count? Two disjoint 5-seed ensembles."""
import numpy as np, torch
src = open('lag.py').read().split('res = {}')[0]
g = {}
exec(src, g)
cols = g['SEARCH'] + g['BOUND']
for name, seeds in [('ens A', range(10, 15)), ('ens B', range(20, 25))]:
    acc = 0
    for s in seeds:
        torch.manual_seed(s); np.random.seed(s)
        sc = g['crossfit'](cols, 'mlp')
        # standardise per fold so seeds contribute on one scale
        for k in g['folds']:
            m = g['fold_of'] == k
            sc[m] = (sc[m] - sc[m].mean()) / sc[m].std()
        acc = acc + sc
    tmp = g['df'][['sample_id', 'precursor_id', 'is_decoy', 'fold', 'in_rt_window', 'pair_key']].assign(score=acc)
    best, w = g['paired_pooled'](tmp, 'score')
    t = best[best.is_decoy == 0]
    nt, nd = (w.is_decoy == 0).sum(), (w.is_decoy == 1).sum()
    low = w[w.z < w[w.is_decoy == 1].z.median()]
    R = ((low.is_decoy == 0).sum() / (low.is_decoy == 1).sum()) / (nt / nd)
    print(f'{name}: q<=0.005/0.01/0.02/0.05 =', [int((t.q <= q).sum()) for q in (0.005, 0.01, 0.02, 0.05)],
          '| OSW IDs among 1%', int(t[t.q <= 0.01].precursor_id.isin(g['o'].precursor_id).sum()), f'| R {R:.3f}', flush=True)
