#!/usr/bin/env python3
"""Extreme-tail counts: decoy winners above t vs 2 x entrapment winners above t, at t set by decoy count.
Also the same under paired competition (pair winners only). Diagnosis only."""
import sqlite3, sys, numpy as np, pandas as pd
lib = sys.argv[1]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID precursor_id, p.TRAML_ID t, p.DECOY d,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID
  group by p.ID""", c)
pq['key'] = np.where(pq.d == 1, pq.t.str.replace('^DECOY_', '', regex=True), pq.t)
ent = pq[pq.d == 0].set_index('precursor_id').e
key = pq.set_index('precursor_id').key
tot = {s: {'D': 0, 'E': 0} for s in ()}
for path in sys.argv[2:]:
    f = pd.read_csv(path, sep='\t', usecols=['precursor_id', 'is_decoy', 'score', 'network_score'])
    f['e'] = f.precursor_id.map(ent).fillna(0).astype(int)
    f['key'] = f.precursor_id.map(key)
    print(f'== {path}')
    for s in ('network_score', 'score'):
        b = f[f[s] > -1e29].sort_values(s, ascending=False).drop_duplicates('precursor_id')
        pw = b.drop_duplicates('key')  # pair winner (best member of each pair)
        for lab, x in (('precursor winners', b), ('pair winners', pw)):
            D = np.sort(x[x.is_decoy == 1][s].values)[::-1]
            E = x[(x.is_decoy == 0) & (x.e == 1)][s].values
            row = []
            for nd in (5, 10, 20, 40, 80, 160):
                if nd > len(D): break
                t = D[nd - 1]
                row.append(f'D{nd}:2E={2 * int((E >= t).sum())}')
            print(f'  {s:14s} {lab:18s} ' + '  '.join(row))
