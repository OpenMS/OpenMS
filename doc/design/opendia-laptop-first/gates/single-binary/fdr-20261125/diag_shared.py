#!/usr/bin/env python3
"""For accepted entrapments (each score's q<=0.01) and for tail decoys of entrapments (DE above the same threshold, or the
top DE): does the isobaric human parent T have ANY candidate window with apex within 2 s of the item's kept apex, and
what is T's best GBM/network raw rank there. Split by fragment overlap class. Diagnosis only."""
import sqlite3, sys, numpy as np, pandas as pd
lib = sys.argv[1]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID precursor_id, p.TRAML_ID t, p.DECOY d, p.CHARGE z, pe.UNMODIFIED_SEQUENCE seq, pe.MODIFIED_SEQUENCE mseq,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
pq['key'] = np.where(pq.d == 1, pq.t.str.replace('^DECOY_', '', regex=True), pq.t)
pq['comp'] = pq.seq.map(lambda s: ''.join(sorted(s))) + '|' + pq.z.astype(str) + '|' + pq.mseq.str.count('UniMod').astype(str)
ecomp = pq[(pq.d == 0) & (pq.e == 1)].set_index('key').comp
pq['pcomp'] = np.where((pq.e == 1) | (pq.d == 1), pq.key.map(ecomp), None)
P = pq.set_index('precursor_id')
hum = set(pq[(pq.d == 0) & (pq.e == 0)].precursor_id)
hum_comp = pq[(pq.d == 0) & (pq.e == 0)].groupby('comp').precursor_id.apply(list)
ov = pd.read_csv('overlap.tsv', sep='\t', index_col=0).ov
for path in sys.argv[2:]:
    f = pd.read_csv(path, sep='\t', usecols=['precursor_id', 'is_decoy', 'rt_apex', 'score', 'q_value', 'network_score', 'network_q_value'])
    allw = f[f.precursor_id.isin(hum)].groupby('precursor_id').rt_apex.apply(np.array)
    print(f'== {path.split("/")[0]}')
    for s, q in (('network_score', 'network_q_value'), ('score', 'q_value')):
        k = f[f[s] > -1e29].set_index('precursor_id').join(P[['e', 'pcomp']])
        k['ov'] = k.index.map(ov)
        thr = k[(k.is_decoy == 0) & (k[q] <= 0.01)][s].min()
        dv = np.sort(k[k.is_decoy == 1][s].values)[::-1]
        def shared(rows):
            out = []
            for cp, rt in zip(rows.pcomp, rows.rt_apex):
                hit = False
                for t in hum_comp.get(cp, []):
                    w = allw.get(t)
                    if w is not None and np.min(np.abs(w - rt)) <= 2: hit = True
                out.append(hit)
            return np.array(out, bool)
        for lab, rows in (('accepted E', k[(k.e == 1) & (k.is_decoy == 0) & (k[s] >= thr)]),
                          ('DE above D200', k[(k.e == 1) & (k.is_decoy == 1) & (k[s] >= dv[199])]),
                          ('E above D200', k[(k.e == 1) & (k.is_decoy == 0) & (k[s] >= dv[199])])):
            sh = shared(rows)
            hi = rows.ov.values >= 2
            print(f'  {s:14s} {lab:14s} n {len(rows):4d}: parent has a window at the apex {sh.sum():3d}; ov2+ {hi.sum():3d} (of those parent window {int((sh & hi).sum())})')
