#!/usr/bin/env python3
"""Near-isomer exposure by precursor class, apples to apples: share of precursors that have an isobaric HUMAN target
(same composition, charge, mod count; not itself) sharing >= 1 / >= 2 library fragments (10 ppm).
Classes: human targets, entrapment targets, decoys of human targets, decoys of entrapments. Diagnosis only."""
import sqlite3, numpy as np, pandas as pd
from collections import defaultdict
lib = '/home/singjc/.jcode/scratch/perf/lib/search.pqp'
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID precursor_id, p.DECOY d, p.CHARGE z, pe.UNMODIFIED_SEQUENCE seq, pe.MODIFIED_SEQUENCE mseq,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
pq['comp'] = pq.seq.map(lambda s: ''.join(sorted(s))) + '|' + pq.z.astype(str) + '|' + pq.mseq.str.count('UniMod').astype(str)
frag = defaultdict(list)
for pid, mz in c.execute('select tpm.PRECURSOR_ID, tr.PRODUCT_MZ from TRANSITION tr join TRANSITION_PRECURSOR_MAPPING tpm on tpm.TRANSITION_ID=tr.ID'):
    frag[pid].append(mz)
frag = {k: np.array(v) for k, v in frag.items()}
groups = pq[(pq.d == 0) & (pq.e == 0)].groupby('comp').precursor_id.apply(list)
def best(pid, cp):
    b_ = -1
    for t in groups.get(cp, []):
        if t == pid: continue
        a, b = frag[pid], frag[t]
        b_ = max(b_, int((np.abs(a[:, None] - b[None, :]) <= b[None, :] * 1e-5).any(axis=1).sum()))
    return b_
rng = np.random.default_rng(0)
for lab, m in (('human targets', (pq.d == 0) & (pq.e == 0)), ('entrapment targets', (pq.d == 0) & (pq.e == 1)),
               ('decoys of human', (pq.d == 1) & (pq.e == 0)), ('decoys of entrapment', (pq.d == 1) & (pq.e == 1))):
    x = pq[m]
    x = x.iloc[rng.choice(len(x), size=min(40000, len(x)), replace=False)]
    v = np.array([best(p, cp) for p, cp in zip(x.precursor_id, x.comp)])
    print(f'{lab:22s} n {len(v)}: has isobaric human {np.mean(v >= 0):.4f}  shares >=1 {np.mean(v >= 1):.4f}  >=2 {np.mean(v >= 2):.4f}  >=3 {np.mean(v >= 3):.4f}', flush=True)
