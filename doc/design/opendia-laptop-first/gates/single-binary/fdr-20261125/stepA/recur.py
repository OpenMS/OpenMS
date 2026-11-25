#!/usr/bin/env python3
"""Step A: cross-acquisition recurrence of accepted entrapments (HeLa D and H are different acquisitions). Diagnosis only.
Per acquisition: entrapment peptides accepted (own q <= QC) in >= 3 of 5 seeds. Overlap D x H vs the chance expectation
|D| |H| / (#entrapment peptides). Then, for H, accepted entrapments split into D-recurrent vs not, with the bound if the
D-recurrent ones were present peptides (label-informed, diagnostic only). The recurrent set is written for the FASTA check.
Usage: recur.py <lib.pqp> <QC> <D glob> <H glob> <out.tsv>"""
import sqlite3, sys, glob, collections, numpy as np, pandas as pd
lib, QC, gD, gH, out = sys.argv[1], float(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
P = pd.read_sql("""select p.ID id, p.DECOY d, pe.UNMODIFIED_SEQUENCE seq,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c).set_index('id')
T = P[P.d == 0]; r = (T.e == 1).sum() / (T.e == 0).sum(); kf = 1 + 1 / r
nent = T[T.e == 1].seq.nunique()
def runs(g):
    res = []
    for f in sorted(glob.glob(g)):
        t = pd.read_csv(f, sep='\t', usecols=['precursor_id', 'is_decoy', 'score', 'q_value'])
        t = t[t.score > -1e29].sort_values('score', ascending=False).drop_duplicates('precursor_id')
        a = t[(t.is_decoy == 0) & (t.q_value <= QC)]
        res.append((len(a), P.seq.reindex(a.precursor_id[a.precursor_id.map(P.e) == 1]).tolist()))
    return res
RD, RH = runs(gD), runs(gH)
def rec(R):
    cnt = collections.Counter(s for _, l in R for s in set(l)); return {s for s, v in cnt.items() if v >= 3}, cnt
sD, cD = rec(RD); sH, cH = rec(RH)
exp = len(sD) * len(sH) / nent
print(f'q <= {QC}: entrapment peptides in library {nent}; recurrent (>= 3 of 5 seeds) D {len(sD)}, H {len(sH)}; '
      f'overlap {len(sD & sH)} vs chance {exp:.2f}', flush=True)
for lab, R, other in (('H', RH, sD), ('D', RD, sH)):
    N = sum(n for n, _ in R); E = sum(len(l) for _, l in R); Eo = sum(sum(s in other for s in l) for _, l in R)
    print(f'{lab}: pooled N {N} E {E} lower {E / N:.4f} combined {E * kf / N:.4f}; E whose peptide recurs in the other acquisition '
          f'{Eo} ({Eo / max(E, 1):.2f}); without them lower {(E - Eo) / N:.4f} combined {(E - Eo) * kf / N:.4f}', flush=True)
pd.DataFrame({'seq': sorted(sD & sH)}).to_csv(out, sep='\t', index=False)
