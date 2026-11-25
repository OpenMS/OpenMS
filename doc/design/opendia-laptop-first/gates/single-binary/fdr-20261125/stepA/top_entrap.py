#!/usr/bin/env python3
"""Step A: the ultra-confident entrapments (own q <= 0.001). Are they the same peptides across seeds and across the two HeLa
runs (D and H, different acquisitions)? Recurrence across independent acquisitions suggests real signal (a present
peptide or a near-identical human/contaminant species), not random matches. Diagnosis only.
Usage: top_entrap.py <lib.pqp> <qcut> <arm>=<glob> ..."""
import sqlite3, sys, glob, collections, numpy as np, pandas as pd
lib, QC = sys.argv[1], float(sys.argv[2])
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
P = pd.read_sql("""select p.ID id, p.DECOY d, p.CHARGE z, pe.UNMODIFIED_SEQUENCE seq, pe.MODIFIED_SEQUENCE mseq,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c).set_index('id')
hits = {}
for spec in sys.argv[3:]:
    arm, g = spec.split('=', 1)
    per = []
    for f in sorted(glob.glob(g)):
        t = pd.read_csv(f, sep='\t', usecols=['precursor_id', 'is_decoy', 'score', 'q_value'])
        t = t[t.score > -1e29].sort_values('score', ascending=False).drop_duplicates('precursor_id')
        a = t[(t.is_decoy == 0) & (t.q_value <= QC)]
        e = a[a.precursor_id.map(P.e) == 1]
        dd = t[(t.is_decoy == 1) & (t.q_value <= QC)] if 'q_value' in t else None
        per.append(set(P.seq.reindex(e.precursor_id)))
        print(f'{arm} {f.split("/")[-2]}: targets q<={QC} {len(a)}, entrapment precursors {len(e)}, peptides {len(per[-1])}', flush=True)
    cnt = collections.Counter(s for x in per for s in x)
    hits[arm] = cnt
    print(f'  {arm}: distinct entrapment peptides {len(cnt)}; in >= 3 of {len(per)} seeds {sum(v >= 3 for v in cnt.values())}; in all {sum(v == len(per) for v in cnt.values())}', flush=True)
arms = list(hits)
for i in range(len(arms)):
    for j in range(i + 1, len(arms)):
        a, b = hits[arms[i]], hits[arms[j]]
        sa = {k for k, v in a.items() if v >= 3}; sb = {k for k, v in b.items() if v >= 3}
        print(f'overlap {arms[i]} x {arms[j]} (peptides in >= 3 seeds each): {len(sa & sb)} of {len(sa)} / {len(sb)}', flush=True)
allc = collections.Counter()
for k in hits.values(): allc.update(k)
print('most recurrent:', ', '.join(f'{s}:{v}' for s, v in allc.most_common(40)))
