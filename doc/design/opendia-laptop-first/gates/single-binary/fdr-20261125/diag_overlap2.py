#!/usr/bin/env python3
"""At each score's own acceptance threshold (q<=0.01) and at fixed decoy depths: counts of accepted entrapments by fragment
overlap with their isobaric human parent, and decoys above the threshold, split into decoys of entrapments (DE, with
overlap) and decoys of human targets (DH). Pooled over runs. Uses overlap.tsv from diag_overlap.py. Diagnosis only."""
import sqlite3, sys, numpy as np, pandas as pd
from collections import defaultdict
lib = sys.argv[1]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID precursor_id, min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
E = pq.set_index('precursor_id').e  # for decoys: 1 if decoy of an entrapment
ov = pd.read_csv('overlap.tsv', sep='\t', index_col=0).ov
tot = defaultdict(lambda: defaultdict(int))
for path in sys.argv[2:]:
    f = pd.read_csv(path, sep='\t', usecols=['precursor_id', 'is_decoy', 'score', 'q_value', 'network_score', 'network_q_value'])
    f['e'] = f.precursor_id.map(E).fillna(0).astype(int)
    f['ov'] = f.precursor_id.map(ov).fillna(-1).astype(int)
    for s, q in (('network_score', 'network_q_value'), ('score', 'q_value')):
        b = f[f[s] > -1e29]
        dv = np.sort(b[b.is_decoy == 1][s].values)[::-1]
        for lab, t in (('q<=0.01', b[(b.is_decoy == 0) & (b[q] <= 0.01)][s].min()), ('D50', dv[49]), ('D200', dv[199])):
            a = b[b[s] >= t]
            T = tot[(s, lab)]
            T['N'] += int((a.is_decoy == 0).sum())
            for cl, m in (('E ov0-1', (a.is_decoy == 0) & (a.e == 1) & (a.ov <= 1)), ('E ov2+', (a.is_decoy == 0) & (a.e == 1) & (a.ov >= 2)),
                          ('DE ov0-1', (a.is_decoy == 1) & (a.e == 1) & (a.ov <= 1)), ('DE ov2+', (a.is_decoy == 1) & (a.e == 1) & (a.ov >= 2)),
                          ('DH', (a.is_decoy == 1) & (a.e == 0))):
                T[cl] += int(m.sum())
print(f'runs {len(sys.argv) - 2}')
for k, T in tot.items():
    print(f'{k[0]:14s} {k[1]:8s} ' + '  '.join(f'{a} {v}' for a, v in T.items()))
