#!/usr/bin/env python3
"""Where does the GBM treat decoys and entrapments differently? Diagnosis only.
For precursor winners in the network tail, compare the GBM-induced rank change of decoys vs entrapments at matched network score,
then which features separate tail decoys from tail entrapments, and whether the item's m/z twin (decoy: its paired target) is a confident ID."""
import sqlite3, sys, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score
lib = sys.argv[1]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID precursor_id, p.TRAML_ID t, p.DECOY d,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID
  group by p.ID""", c)
pq['key'] = np.where(pq.d == 1, pq.t.str.replace('^DECOY_', '', regex=True), pq.t)
P = pq.set_index('precursor_id')
META = ('candidate_row', 'precursor_id', 'is_decoy', 'fold', 'rt_apex', 'score', 'q_value', 'network_score', 'network_q_value', 'e', 'key', 'cls')
for path in sys.argv[2:]:
    f = pd.read_csv(path, sep='\t')
    f['e'] = f.precursor_id.map(P.e).fillna(0).astype(int)
    f['key'] = f.precursor_id.map(P.key)
    feats = [x for x in f.columns if x not in META]
    # GBM-kept window per precursor carries score/q; network-kept window carries network_score.
    g = f[f.score > -1e29].set_index('precursor_id')
    n = f[f.network_score > -1e29].set_index('precursor_id')
    b = n[['network_score', 'is_decoy', 'e', 'key', 'fold'] + feats].join(g[['score', 'q_value']], how='inner')
    b['cls'] = np.where(b.is_decoy == 1, 'D', np.where(b.e == 1, 'E', 'T'))
    tgt_q = g[g.is_decoy == 0].groupby('key').q_value.min()
    print(f'== {path.split("/")[0]}: same window kept by both scores for {np.mean(g.loc[b.index].candidate_row.values == n.loc[b.index].candidate_row.values):.3f} of precursors')
    # Ranks among decoy winners: network rank and GBM rank in decoy-tail units (-log10 survival among decoys).
    for s in ('network_score', 'score'):
        dv = np.sort(b[b.cls == 'D'][s].values)
        b[s + '_u'] = -np.log10((len(dv) - np.searchsorted(dv, b[s].values, side='left') + 1) / (len(dv) + 1))
    b['dsh'] = b['score_u'] - b['network_score_u']
    print('  GBM shift in decoy-tail units (log10), by network tail bin; mean and n for D and E (2E/D excess = rows above that bin moving up):')
    for lo, hi in ((1, 1.5), (1.5, 2), (2, 2.5), (2.5, 3), (3, 6)):
        m = (b.network_score_u >= lo) & (b.network_score_u < hi)
        dd, ee = b[m & (b.cls == 'D')].dsh, b[m & (b.cls == 'E')].dsh
        print(f'    net u [{lo},{hi}): D n {len(dd):5d} mean {dd.mean():+.3f} | E n {len(ee):5d} mean {ee.mean():+.3f} | diff E-D {ee.mean() - dd.mean():+.3f}')
    # Twin presence: decoys whose paired target is a confident (GBM q<=0.01) ID.
    for lab, m in (('network tail u>=2', b.network_score_u >= 2), ('GBM tail u>=2', b.score_u >= 2)):
        dd = b[m & (b.cls == 'D')]
        print(f'  {lab}: decoys {len(dd)}, paired target is a q<=0.01 ID: {np.mean(dd.key.map(tgt_q).fillna(1).values <= 0.01):.3f}; entrapments {int((m & (b.cls == "E")).sum())}')
    # Tail separability decoy vs entrapment (network tail u >= 1.5: ~top 3 % of decoys), cross-fitted by fold.
    X = b[(b.network_score_u >= 1.5) & (b.cls != 'T')]
    y = (X.cls == 'D').values.astype(int)
    p = np.zeros(len(X))
    for fo in sorted(X.fold.unique()):
        tr, te = X.fold.values != fo, X.fold.values == fo
        p[te] = HistGradientBoostingClassifier(max_iter=100, learning_rate=0.05, min_samples_leaf=40, random_state=0).fit(X[feats].values[tr], y[tr]).predict_proba(X[feats].values[te])[:, 1]
    print(f'  tail (net u>=1.5) decoy {y.sum()} vs entrapment {(1 - y).sum()}: cross-fitted AUC {roc_auc_score(y, p):.3f}')
    sm = []
    for v in feats:
        a, e_ = X[X.cls == 'D'][v], X[X.cls == 'E'][v]
        sd = np.sqrt((a.var() + e_.var()) / 2) or 1
        sm.append((v, (e_.mean() - a.mean()) / sd, roc_auc_score(y, X[v].values)))
    sm.sort(key=lambda t: -abs(t[1]))
    print('  tail SMD entrapment minus decoy:', ', '.join(f'{v} {d:+.2f}' for v, d, _ in sm[:14]))
