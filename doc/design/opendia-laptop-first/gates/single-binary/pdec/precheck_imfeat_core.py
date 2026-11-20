"""PRECHECK_IMFEAT step 4: does the C++-feasible variant (core 1/K0 histogram, c_*) carry the band6 gain?
Same sample, labels and GBM ceilings as precheck_imfeat.py. Arms: 36 | 36+band6 | 36+core5 | 36+core5+band6, 2 seeds."""
import os, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
O = os.environ['FAINT2']; FEAT = os.environ['FEATS']
S = pd.read_csv(f'{O}/sample.tsv', sep='\t').set_index('pid')
D = pd.read_parquet(f'{O}/imfeat.parquet')
F = pd.read_csv(FEAT, sep='\t'); F['rta'] = F.rt_apex.round(3)
F = F.merge(D.drop_duplicates(['precursor_id', 'rta']), on=['precursor_id', 'rta'])
BAND = ['b_offset', 'b_frag_sd', 'b_cos_apex', 'b_support', 'b_contrast', 'b_local_sn']
CORE = ['c_offset', 'c_frag_sd', 'c_cos', 'c_support', 'c_contrast']
meta = ['candidate_row', 'precursor_id', 'is_decoy', 'fold', 'rt_apex', 'score', 'q_value', 'rta', 'p_cos_apex', 'p_local_sn']
X36 = [c for c in F.columns if c not in meta + BAND + CORE]
dec = F.is_decoy.values == 1
pos = (F.precursor_id.map(S.cls).values == 'dnn') & (np.abs(F.rt_apex.values - F.precursor_id.map(S.rt_s).values) <= 6)
fold = F.fold.values; qq = S.qq; dec_ids = S.index[S.cls == 'decoy']
print(f'joined {len(F)}, X36 {len(X36)}, right windows {pos.sum()}', flush=True)
def gbm(M, seed):
    s = np.zeros(len(F))
    for f in np.unique(fold):
        tr = (fold != f) & (pos | dec); y = pos[tr].astype(int); w = np.where(y == 1, (y == 0).sum() / max(y.sum(), 1), 1.0)
        s[fold == f] = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=31, random_state=seed).fit(M[tr], y, sample_weight=w).decision_function(M[fold == f])
    return s
def tpr(s):
    F['_r'] = s; best = F.groupby('precursor_id')._r.max(); rb = F[pos].groupby('precursor_id')._r.max(); neg = best.reindex(dec_ids).dropna().values
    return [(rb > np.quantile(neg, 1 - f)).mean() for f in (1e-3, 3e-3, 1e-2)], (rb[rb.index.map(qq) == 'Q1'] > np.quantile(neg, 1 - 3e-3)).mean()
res = {}
for name, cols in (('36', X36), ('36+band6', X36 + BAND), ('36+core5', X36 + CORE), ('36+core5+band6', X36 + CORE + BAND)):
    r = [tpr(gbm(F[cols].values.astype(np.float32), sd)) for sd in (0, 1)]; res[name] = r
    print(f'{name:16s} ' + ' | '.join(f'seed{sd} {a[0]:.3f}/{a[1]:.3f}/{a[2]:.3f} Q1 {q:.3f}' for sd, (a, q) in enumerate(r)), flush=True)
for name in res:
    if name != '36':
        d = [res[name][sd][0][1] - res['36'][sd][0][1] for sd in (0, 1)]
        print(f'  {name:16s} delta@3e-3 {d[0]:+.3f} / {d[1]:+.3f}')
print('done')
