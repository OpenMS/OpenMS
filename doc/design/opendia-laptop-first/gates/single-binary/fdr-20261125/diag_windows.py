#!/usr/bin/env python3
"""Window-count asymmetry: the GBM keeps max-score window per precursor. If false targets have more candidate windows
than decoys, max-of-more inflates them. Diagnosis only."""
import sqlite3, sys, numpy as np, pandas as pd
lib = sys.argv[1]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID precursor_id, p.DECOY d,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID
  group by p.ID""", c)
E = pq.set_index('precursor_id').e
for path in sys.argv[2:]:
    f = pd.read_csv(path, sep='\t', usecols=['precursor_id', 'is_decoy', 'score', 'network_score'])
    f['cls'] = np.where(f.is_decoy == 1, 'D', np.where(f.precursor_id.map(E) == 1, 'E', 'T'))
    w = f.groupby('precursor_id').agg(nw=('score', 'size'), cls=('cls', 'first'))
    print(f'== {path.split("/")[0]}')
    for k in 'DET':
        x = w[w.cls == k].nw
        print(f'  {k}: precursors {len(x)}, windows mean {x.mean():.3f}, median {x.median():.0f}, P(1) {np.mean(x == 1):.3f}, P(>=5) {np.mean(x >= 5):.3f}')
    for s in ('network_score', 'score'):
        b = f[f[s] > -1e29].set_index('precursor_id')
        b = b.join(w.nw)
        d = np.sort(b[b.cls == 'D'][s].values)[::-1]
        for nd in (20, 80, 320):
            t = d[nd - 1]
            top = b[(b[s] >= t) & (b.cls != 'T')]
            print(f'  {s:14s} top D{nd}: D windows mean {top[top.cls == "D"].nw.mean():.2f}, E n {int((top.cls == "E").sum())} windows mean {top[top.cls == "E"].nw.mean():.2f}')
