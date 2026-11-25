#!/usr/bin/env python3
"""Do the features separate decoys from false (entrapment) targets? Cross-fitted classifier, diagnosis only.
If yes, any learned target-vs-decoy scorer can push decoys below false targets, so TDC under-counts false targets."""
import sqlite3, sys, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score
from sklearn.inspection import permutation_importance
lib = sys.argv[1]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID precursor_id, p.DECOY d,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID
  group by p.ID""", c)
ent = pq[pq.d == 0].set_index('precursor_id').e
META = ('candidate_row', 'precursor_id', 'is_decoy', 'fold', 'rt_apex', 'score', 'q_value', 'network_score', 'network_q_value', 'e')
for path in sys.argv[2:]:
    f = pd.read_csv(path, sep='\t')
    f['e'] = f.precursor_id.map(ent).fillna(0).astype(int)
    feats = [x for x in f.columns if x not in META]
    # Every window (the GBM trains on rows), and best windows (what TDC competes).
    for lab, x in (('all windows', f), ('best windows', f[f.network_score > -1e29].sort_values('network_score', ascending=False).drop_duplicates('precursor_id'))):
        X = x[(x.is_decoy == 1) | (x.e == 1)]
        y = X.is_decoy.values
        p = np.zeros(len(X))
        for fo in sorted(X.fold.unique()):
            tr, te = X.fold.values != fo, X.fold.values == fo
            m = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.05, max_leaf_nodes=31, random_state=0).fit(X[feats].values[tr], y[tr])
            p[te] = m.predict_proba(X[feats].values[te])[:, 1]
        print(f'{path.split("/")[0]} {lab}: n decoy {int(y.sum())} entrapment {int((1 - y).sum())}; cross-fitted AUC decoy vs entrapment {roc_auc_score(y, p):.3f}')
        if lab == 'best windows':
            single = sorted(((v, roc_auc_score(y, X[v].values)) for v in feats), key=lambda t: -abs(t[1] - 0.5))[:12]
            print('   single-feature AUC:', ', '.join(f'{a} {b:.3f}' for a, b in single))
            # drop-one: which feature group removal brings cross-fitted AUC to 0.5
            for drop in [[v] for v, _ in single[:6]]:
                keep = [v for v in feats if v not in drop]
                p2 = np.zeros(len(X))
                for fo in sorted(X.fold.unique()):
                    tr, te = X.fold.values != fo, X.fold.values == fo
                    m = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.05, random_state=0).fit(X[keep].values[tr], y[tr])
                    p2[te] = m.predict_proba(X[keep].values[te])[:, 1]
                print(f'   without {drop}: AUC {roc_auc_score(y, p2):.3f}')
