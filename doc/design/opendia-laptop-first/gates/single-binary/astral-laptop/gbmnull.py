"""PRECHECK_GBMNULL statistic. Reads only targets, decoys and the library pairing (no entrapment labels).
Usage: gbmnull.py <library.pqp> <set dir> -> one line per seed: set, seed, region n, b, c, p, decision.
<set dir> holds <name>_net_s<k>/ and <name>_gbm_s<k>/ run directories (seed-matched; only the GBM differs)."""
import sqlite3, sys, glob, os, math, numpy as np, pandas as pd
from scipy.stats import binom

L, D = sys.argv[1:3]
LN2 = math.log(2.0)
c = sqlite3.connect(f'file:{L}?mode=ro', uri=True)
pq = pd.read_sql('select ID precursor_id, TRAML_ID t, DECOY d from PRECURSOR', c)
key = np.where(pq.d == 1, pq.t.str.replace('^DECOY_', '', regex=True), pq.t)
pq = pq.assign(key=key)
pairs = pq[pq.d == 0].set_index('key')[['precursor_id']].rename(columns={'precursor_id': 'tg'}).join(
    pq[pq.d == 1].set_index('key')[['precursor_id']].rename(columns={'precursor_id': 'dc'}), how='inner')


def scores(d):
    f = [x for x in glob.glob(f'{d}/*.precursors.tsv*') if '/precursors.tsv' not in x][0]
    t = pd.read_csv(f, sep='\t', usecols=['precursor_id', 'score'])
    t = t[t.score > -1e29].sort_values('score', ascending=False).drop_duplicates('precursor_id')
    return t.set_index('precursor_id').score


print('set\tseed\tregion_pairs\tb_decoy_to_target\tc_target_to_decoy\tp_one_sided\tnet_target_win\tgbm_target_win\tdecision')
for nd in sorted(glob.glob(f'{D}/*_net_s*')):
    name, seed = os.path.basename(nd).split('_net_s')
    gd = f'{D}/{name}_gbm_s{seed}'
    if not os.path.isdir(gd): continue
    sn, sg = scores(nd), scores(gd)
    x = pairs.assign(n_t=sn.reindex(pairs.tg).values, n_d=sn.reindex(pairs.dc).values,
                     g_t=sg.reindex(pairs.tg).values, g_d=sg.reindex(pairs.dc).values).dropna()
    x = x[(x.n_t <= LN2) & (x.n_d <= LN2)]
    wn, wg = x.n_t > x.n_d, x.g_t > x.g_d
    b, cc = int((~wn & wg).sum()), int((wn & ~wg).sum())
    p = float(binom.sf(b - 1, b + cc, 0.5)) if b + cc else 1.0
    print(f'{name}\t{seed}\t{len(x)}\t{b}\t{cc}\t{p:.3g}\t{wn.mean():.4f}\t{wg.mean():.4f}\t{"skip" if p < 0.01 else "use"}')
