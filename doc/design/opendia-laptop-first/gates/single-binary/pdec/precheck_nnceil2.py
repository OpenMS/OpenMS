"""PRECHECK_NNCEIL step 2: do the label-free CNN logit and the PRECHECK_FAINT trace features carry the same information?
DIA-NN-labelled GBM ceilings on the faint sample: 36 / 36+faint7 / 36+CNNlf / 36+faint7+CNNlf, two GBM seeds.
Labels for evaluation only."""
import glob, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
G = '/home/singjc/.jcode/scratch/dia2/gap7'; O = f'{G}/faint'
S = pd.read_csv(f'{O}/sample.tsv', sep='\t').set_index('pid')
Z = np.load(f'{O}/nnceil_scores.npz')
C = pd.DataFrame({'precursor_id': Z['pid'], 'rta': Z['rta'], 'cnn_lf': Z['cnn_lf']})
N = pd.read_parquet(f'{O}/new_features.parquet'); N['rta'] = N.rt_apex.round(3)
NF = ['ps_cos_w', 'ps_cos_apex', 'ps_n4', 'apex_spread', 'frag_run', 'local_sn', 'rank_corr']
ids = set(S.index); parts = []
for ch in pd.read_csv(glob.glob(f'{G}/ev/v08m8/*.features.tsv')[0], sep='\t', chunksize=1_000_000):
    parts.append(ch[ch.precursor_id.isin(ids)])
F = pd.concat(parts); F['rta'] = F.rt_apex.round(3)
F = F.merge(C.drop_duplicates(['precursor_id', 'rta']), on=['precursor_id', 'rta']).merge(N[['precursor_id', 'rta'] + NF].drop_duplicates(['precursor_id', 'rta']), on=['precursor_id', 'rta'])
meta = ['candidate_row', 'precursor_id', 'is_decoy', 'fold', 'rt_apex', 'score', 'q_value', 'rta', 'cnn_lf'] + NF
X36 = [c for c in F.columns if c not in meta]
dec = F.is_decoy.values == 1
pos = (F.precursor_id.map(S.cls).values == 'dnn') & (np.abs(F.rt_apex.values - F.precursor_id.map(S.rt_s).values) <= 6)
fold = F.fold.values; qq = S.qq; dec_ids = S.index[S.cls == 'decoy']
print(f'windows {len(F)}', flush=True)
def gbm(cols, seed):
    Xv = F[cols].values.astype(np.float32); s = np.zeros(len(F))
    for f in np.unique(fold):
        tr = (fold != f) & (pos | dec); y = pos[tr].astype(int); w = np.where(y == 1, (y == 0).sum() / max(y.sum(), 1), 1.0)
        s[fold == f] = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=31, random_state=seed).fit(Xv[tr], y, sample_weight=w).decision_function(Xv[fold == f])
    return s
def tpr(s):
    F['_r'] = s; best = F.groupby('precursor_id')._r.max(); rb = F[pos].groupby('precursor_id')._r.max(); neg = best.reindex(dec_ids).dropna().values
    return [(rb > np.quantile(neg, 1 - f)).mean() for f in (1e-3, 3e-3, 1e-2)], [(rb[rb.index.map(qq) == 'Q1'] > np.quantile(neg, 1 - 3e-3)).mean()]
for name, cols in (('36', X36), ('36+faint7', X36 + NF), ('36+CNNlf', X36 + ['cnn_lf']), ('36+faint7+CNNlf', X36 + NF + ['cnn_lf'])):
    r = [tpr(gbm(cols, sd)) for sd in (0, 1)]
    print(f'{name:16s} TPR @1e-3/3e-3/1e-2: ' + ' | '.join('/'.join(f'{v:.3f}' for v in a) + f' (Q1@3e-3 {b[0]:.3f})' for a, b in r), flush=True)
