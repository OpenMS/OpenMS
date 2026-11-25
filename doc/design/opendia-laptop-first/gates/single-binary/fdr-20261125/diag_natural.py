#!/usr/bin/env python3
"""How common are near-isomer competitors in a natural (human-only) library? For each human target: does another human
target with the same composition, charge and mod count exist, and how many of its 6 library fragments match (10 ppm)?
Compare with entrapment vs its human parent (7.4 % share 2+ under 'overlap.tsv' classes 1/2+ as computed there)."""
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
hum = pq[(pq.d == 0) & (pq.e == 0)]
groups = hum.groupby('comp').precursor_id.apply(list)
best = []
for pid, cp, seq in zip(hum.precursor_id, hum.comp, hum.seq):
    others = [t for t in groups[cp] if t != pid]
    if not others: best.append(-1); continue
    a = np.array(frag[pid]); b_ = 0
    for t in others:
        b = np.array(frag[t])
        b_ = max(b_, int((np.abs(a[:, None] - b[None, :]) <= b[None, :] * 1e-5).any(axis=1).sum()))
    best.append(b_)
best = np.array(best)
print(f'human targets {len(best)}: with an isobaric human (same composition) {np.mean(best >= 0):.4f}; sharing >=1 fragment {np.mean(best >= 1):.4f}; >=2 {np.mean(best >= 2):.4f}')
ov = pd.read_csv('overlap.tsv', sep='\t', index_col=0).ov
eo = ov.reindex(pq[(pq.d == 0) & (pq.e == 1)].precursor_id)
print(f'entrapments vs their parent: >=1 {np.mean(eo >= 1):.4f}; >=2 {np.mean(eo >= 2):.4f}')
