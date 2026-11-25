#!/usr/bin/env python3
"""Step A control for recur.py: do counted decoys (score >= the lowest accepted target in the run, i.e. the decoys TDC counts
at q <= 0.01) recur across the two acquisitions as accepted entrapments do? Same rule: >= 3 of 5 seeds per acquisition.
Systematic recurrence of a false hit means it matches something present in both samples. If decoys stand in for false
targets, the recurrent fraction should be similar for entrapments and decoys. Diagnosis only.
Usage: recur_decoy.py <lib.pqp> <D glob> <H glob>"""
import sqlite3, sys, glob, collections, numpy as np, pandas as pd
lib, gD, gH = sys.argv[1:4]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
P = pd.read_sql("""select p.ID id, p.DECOY d, pe.UNMODIFIED_SEQUENCE seq,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c).set_index('id')
def runs(g):
    out = []
    for f in sorted(glob.glob(g)):
        t = pd.read_csv(f, sep='\t', usecols=['precursor_id', 'is_decoy', 'score', 'q_value'])
        t = t[t.score > -1e29].sort_values('score', ascending=False).drop_duplicates('precursor_id')
        a = t[(t.is_decoy == 0) & (t.q_value <= 0.01)]; thr = a.score.min()
        e = P.seq.reindex(a.precursor_id[a.precursor_id.map(P.e) == 1]).tolist()
        d = P.seq.reindex(t.precursor_id[(t.is_decoy == 1) & (t.score >= thr)]).tolist()
        out.append((e, d))
    return out
RD, RH = runs(gD), runs(gH)
def rec(R, k):
    cnt = collections.Counter(s for x in R for s in set(x[k])); return {s for s, v in cnt.items() if v >= 3}
for k, lab in ((0, 'accepted entrapments'), (1, 'counted decoys')):
    sD, sH = rec(RD, k), rec(RH, k)
    for name, R, other in (('D', RD, sH), ('H', RH, sD)):
        n = sum(len(x[k]) for x in R); o = sum(sum(s in other for s in x[k]) for x in R)
        print(f'{lab:21s} {name}: hits {n}, whose peptide recurs (>= 3 of 5 seeds) in the other acquisition {o} ({o / max(n, 1):.3f}); '
              f'recurrent sets D {len(sD)} H {len(sH)} overlap {len(sD & sH)}', flush=True)
