#!/usr/bin/env python3
"""q sweep: at own q thresholds, accepted targets N, decoys above the threshold D, entrapment-estimated false 2E (r=1).
Under valid TDC, 2E ~ D in expectation. Conditional test: E | E+D ~ Bin(E+D, 1/3). Pooled over runs. Diagnosis only."""
import sqlite3, sys, numpy as np, pandas as pd
from scipy.stats import binom
lib = sys.argv[1]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID precursor_id, min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
E = pq.set_index('precursor_id').e
Q = (0.005, 0.01, 0.02, 0.05, 0.1, 0.2)
tot = {}
for path in sys.argv[2:]:
    f = pd.read_csv(path, sep='\t', usecols=['precursor_id', 'is_decoy', 'score', 'q_value', 'network_score', 'network_q_value'])
    f['e'] = f.precursor_id.map(E).fillna(0).astype(int)
    for s, q in (('network_score', 'network_q_value'), ('score', 'q_value')):
        b = f[f[s] > -1e29]
        for qq in Q:
            acc = b[(b.is_decoy == 0) & (b[q] <= qq)]
            t = acc[s].min()
            d = int(((b.is_decoy == 1) & (b[s] >= t)).sum())
            r = tot.setdefault((s, qq), np.zeros(3, int)); r += [len(acc), d, int(acc.e.sum())]
print(f'runs {len(sys.argv) - 2}')
for s in ('network_score', 'score'):
    print(s)
    for qq in Q:
        n, d, e = tot[(s, qq)]
        p = binom.sf(e - 1, e + d, 1 / 3)
        print(f'   q<={qq:<5}  N {n:7d}  D {d:6d}  2E {2 * e:6d}  2E/D {2 * e / max(d, 1):.2f}  FDP_c {2 * e / n:.4f}  P(E>=obs) {p:.3g}')
