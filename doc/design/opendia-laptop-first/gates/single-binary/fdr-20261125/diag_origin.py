#!/usr/bin/env python3
"""Where do the GBM's accepted false targets come from? For the GBM acceptance set (own q <= 0.01) and the
network acceptance set, count decoys above the threshold (D, the TDC estimate of false targets) and 2E (the entrapment
estimate), split by the precursor's network tail units (log10 decoy survival). Pooled over the given runs. Diagnosis only."""
import sqlite3, sys, numpy as np, pandas as pd
lib = sys.argv[1]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID precursor_id, p.DECOY d,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
E = pq.set_index('precursor_id').e
bins = [-9, 1, 2, 2.5, 3, 3.5, 99]
tot = {}
for path in sys.argv[2:]:
    f = pd.read_csv(path, sep='\t', usecols=['precursor_id', 'is_decoy', 'score', 'q_value', 'network_score', 'network_q_value'])
    f['e'] = f.precursor_id.map(E).fillna(0).astype(int)
    n = f[f.network_score > -1e29].set_index('precursor_id')
    g = f[f.score > -1e29].set_index('precursor_id')
    b = n[['network_score', 'network_q_value', 'is_decoy', 'e']].join(g[['score', 'q_value']], how='inner')
    dv = np.sort(b[b.is_decoy == 1].network_score.values)
    b['u'] = -np.log10((len(dv) - np.searchsorted(dv, b.network_score.values, side='left') + 1) / (len(dv) + 1))
    b['ub'] = pd.cut(b.u, bins)
    for s, q in (('network_score', 'network_q_value'), ('score', 'q_value')):
        tq = b[(b.is_decoy == 0) & (b[q] <= 0.01)][s].min()  # acceptance threshold in that score
        above = b[b[s] >= tq]
        for ub, x in above.groupby('ub', observed=False):
            t = tot.setdefault((s, ub), [0, 0, 0])
            t[0] += int(((x.is_decoy == 0)).sum()); t[1] += int((x.is_decoy == 1).sum()); t[2] += int(((x.is_decoy == 0) & (x.e == 1)).sum())
print(f'runs: {len(sys.argv) - 2}')
for s in ('network_score', 'score'):
    print(f'{s} acceptance set, by network tail units u (log10 decoy survival under the network):')
    T = D = E2 = 0
    for ub in pd.cut(pd.Series([0.0]), bins).cat.categories:
        t, d, e = tot.get((s, ub), [0, 0, 0])
        T += t; D += d; E2 += 2 * e
        print(f'   u {str(ub):12s} targets {t:6d}  decoys above {d:4d}  2E {2 * e:4d}')
    print(f'   total targets {T}, decoys {D}, 2E {E2}: TDC FDP est {D / T:.4f}, entrapment est {E2 / T:.4f}')
