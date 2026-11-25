#!/usr/bin/env python3
"""Mechanism check: are decoys distinguishable from false targets? Entrapment-free test plus entrapment reference.

Null targets: target precursor winners whose network score is below the decoy winners' median in their fold. Under
target-decoy exchangeability these are (almost all) false targets and their feature distribution should match the
decoys'. A cross-fitted classifier decoy vs null target therefore has AUC 0.5 if decoys mimic false targets.
Reference (entrapment, diagnosis only): entrapment vs decoy winners near the threshold."""
import sqlite3, sys, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score

lib = sys.argv[1]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID precursor_id, p.DECOY d,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID
  group by p.ID""", c)
ent = pq[pq.d == 0].set_index('precursor_id').e
k = 2.0
META = ('candidate_row', 'precursor_id', 'is_decoy', 'fold', 'rt_apex', 'score', 'q_value', 'network_score', 'network_q_value', 'e')

def n_at_fdpc(b, s, q):
    t = b[b.is_decoy == 0].sort_values([q, s], ascending=[True, False])
    fdp = np.cumsum(t.e.values) * k / np.arange(1, len(t) + 1)
    ok = np.where(fdp <= 0.01)[0]
    return int(ok.max() + 1) if len(ok) else 0

for path in sys.argv[2:]:
    f = pd.read_csv(path, sep='\t')
    f['e'] = f.precursor_id.map(ent).fillna(0).astype(int)
    feats = [x for x in f.columns if x not in META]
    print(f'\n== {path}')
    for tag, s, q in [('net', 'network_score', 'network_q_value'), ('gbm', 'score', 'q_value')]:
        b = f[f[s] > -1e29].sort_values(s, ascending=False).drop_duplicates('precursor_id')
        print(f'{tag}: N@FDP_c1% {n_at_fdpc(b, s, q)}')
    w = f[f.network_score > -1e29].sort_values('network_score', ascending=False).drop_duplicates('precursor_id').copy()
    med = w[w.is_decoy == 1].groupby('fold').network_score.median()
    w['null_t'] = (w.is_decoy == 0) & (w.network_score < w.fold.map(med))
    D, T = w[w.is_decoy == 1], w[w.null_t]
    print(f'decoy winners {len(D)}, null target winners {len(T)} (of which entrapment {T.e.mean():.3f})')
    smd_null = {x: (T[x].mean() - D[x].mean()) / (np.sqrt(0.5 * (T[x].var() + D[x].var())) or 1.0) for x in feats}
    near = w[w.network_q_value <= 0.05]
    E2, D2 = near[(near.is_decoy == 0) & (near.e == 1)], near[near.is_decoy == 1]
    smd_ent = {x: (E2[x].mean() - D2[x].mean()) / (np.sqrt(0.5 * (E2[x].var() + D2[x].var())) or 1.0) for x in feats}
    top = sorted(feats, key=lambda x: -abs(smd_null[x]))[:15]
    print('feature: null-target minus decoy (bulk) | entrapment minus decoy (near threshold), SD units')
    for x in top: print(f'  {x:28s} {smd_null[x]:+.2f} | {smd_ent[x]:+.2f}')
    a, b2 = np.array([smd_null[x] for x in feats]), np.array([smd_ent[x] for x in feats])
    print(f'correlation over {len(feats)} features: {np.corrcoef(a, b2)[0, 1]:.2f}')
    # Cross-fitted decoy-vs-null-target AUC (by library fold).
    X = pd.concat([D, T]); y = X.is_decoy.values
    p = np.zeros(len(X))
    for fo in sorted(X.fold.unique()):
        tr, te = X.fold.values != fo, X.fold.values == fo
        m = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.05, random_state=0).fit(X[feats].values[tr], y[tr])
        p[te] = m.predict_proba(X[feats].values[te])[:, 1]
    print(f'cross-fitted AUC decoy vs null target: {roc_auc_score(y, p):.3f}')
    # Per-feature AUC
    aucs = sorted(((x, roc_auc_score(y, X[x].values)) for x in feats), key=lambda v: -abs(v[1] - 0.5))[:10]
    print('single-feature AUC (decoy vs null target):', ', '.join(f'{x} {v:.3f}' for x, v in aucs))
