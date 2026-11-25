#!/usr/bin/env python3
"""Tail exchangeability per score: entrapment winners vs decoy winners (both false). Diagnosis only.
Under exchangeability a decoy-quantile threshold at 1 - a passes fraction a of entrapment winners too."""
import sqlite3, sys, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
lib = sys.argv[1]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID precursor_id, p.DECOY d,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID
  group by p.ID""", c)
ent = pq[pq.d == 0].set_index('precursor_id').e
for path in sys.argv[2:]:
    f = pd.read_csv(path, sep='\t', usecols=['precursor_id', 'is_decoy', 'score', 'network_score'])
    f['e'] = f.precursor_id.map(ent).fillna(0).astype(int)
    print(f'== {path}')
    for s in ('network_score', 'score'):
        b = f[f[s] > -1e29].sort_values(s, ascending=False).drop_duplicates('precursor_id')
        D, E = b[b.is_decoy == 1][s].values, b[(b.is_decoy == 0) & (b.e == 1)][s].values
        y = np.r_[np.zeros(len(D)), np.ones(len(E))]
        auc = roc_auc_score(y, np.r_[D, E])
        out = []
        for a in (0.01, 0.003, 0.001):
            t = np.quantile(D, 1 - a)
            out.append(f'{a:.3f}: entrapment {np.mean(E > t):.4f} ({np.mean(E > t) / a:.2f}x)')
        print(f'  {s}: AUC entrapment vs decoy {auc:.3f}; ' + '; '.join(out))
