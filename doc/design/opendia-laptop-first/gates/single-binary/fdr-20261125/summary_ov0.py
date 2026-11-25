#!/usr/bin/env python3
"""Own q<=0.01 entrapment counts split by fragment overlap with the isobaric human parent (0, 1, 2+), with the
class share of all entrapments, so each class gives its own combined bound E_c (1 + 1/r_c) / N, r_c = r * share_c.
Post hoc, descriptive only."""
import sqlite3, sys, numpy as np, pandas as pd
lib = sys.argv[1]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID precursor_id, min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
pq2 = pd.read_sql("select ID precursor_id, DECOY d from PRECURSOR", c)
E = pq.set_index('precursor_id').e
ENT = E[(E == 1).values & (pq2.set_index('precursor_id').d.reindex(E.index).values == 0)].index  # entrapment targets only
ov = pd.read_csv('overlap.tsv', sep='\t', index_col=0).ov
eo = ov.reindex(ENT).fillna(0).clip(lower=0, upper=2).astype(int)
share = eo.value_counts(normalize=True).reindex([0, 1, 2]).fillna(0)
print('entrapment share by overlap class 0/1/2+:', ' '.join(f'{v:.4f}' for v in share.values))
acc = {('net', k): 0 for k in range(3)} | {('gbm', k): 0 for k in range(3)}
N = {'net': 0, 'gbm': 0}
for path in sys.argv[2:]:
    f = pd.read_csv(path, sep='\t', usecols=['precursor_id', 'is_decoy', 'score', 'q_value', 'network_score', 'network_q_value'])
    f['e'] = f.precursor_id.map(E).fillna(0).astype(int)
    f['oc'] = f.precursor_id.map(ov).fillna(0).clip(lower=0, upper=2).astype(int)
    for s, sc, q in (('net', 'network_score', 'network_q_value'), ('gbm', 'score', 'q_value')):
        a = f[(f[sc] > -1e29) & (f.is_decoy == 0) & (f[q] <= 0.01)]
        N[s] += len(a)
        for k in range(3): acc[(s, k)] += int(((a.e == 1) & (a.oc == k)).sum())
print(f'runs {len(sys.argv) - 2}')
for s in ('net', 'gbm'):
    parts = []
    for k in range(3):
        e = acc[(s, k)]; rc = share.loc[k]
        parts.append(f'ov{k}{"+" if k == 2 else ""}: E {e} bound {e * (1 + 1 / rc) / N[s]:.4f}')
    tot = sum(acc[(s, k)] for k in range(3))
    print(f'  {s}: N {N[s]} E {tot} all-class bound {2 * tot / N[s]:.4f} | ' + ' | '.join(parts))
