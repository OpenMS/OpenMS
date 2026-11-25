#!/usr/bin/env python3
"""Fragment-overlap test. Entrapment E is a K/R- and C-term-preserving shuffle of a human parent T (same precursor m/z).
Its decoy DE is the pseudo-reverse of E (also isobaric). If E shares more library fragments with T than DE does, E can
take T's real signal where DE cannot, and decoys under-represent that kind of false target.
For every (E, DE) pair: ov_E = #E fragments within 10 ppm of a T fragment (T = isobaric human parent, same charge and mod
count), ov_DE likewise. Then tail and acceptance counts by overlap class under the network and the GBM. Diagnosis only."""
import sqlite3, sys, numpy as np, pandas as pd
from collections import defaultdict
lib = sys.argv[1]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID precursor_id, p.TRAML_ID t, p.DECOY d, p.CHARGE z, pe.UNMODIFIED_SEQUENCE seq, pe.MODIFIED_SEQUENCE mseq,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
pq['key'] = np.where(pq.d == 1, pq.t.str.replace('^DECOY_', '', regex=True), pq.t)
pq['comp'] = pq.seq.map(lambda s: ''.join(sorted(s))) + '|' + pq.z.astype(str) + '|' + pq.mseq.str.count('UniMod').astype(str)
frag = defaultdict(list)
for pid, mz in c.execute('select tpm.PRECURSOR_ID, tr.PRODUCT_MZ from TRANSITION tr join TRANSITION_PRECURSOR_MAPPING tpm on tpm.TRANSITION_ID=tr.ID'):
    frag[pid].append(mz)
hum = pq[(pq.d == 0) & (pq.e == 0)]
hum_by_comp = hum.groupby('comp').precursor_id.apply(list)
def overlap(pid, comp):
    ts = hum_by_comp.get(comp)
    if ts is None: return -1
    a = np.array(frag[pid])
    best = 0
    for t in ts:
        b = np.array(frag[t])
        best = max(best, int((np.abs(a[:, None] - b[None, :]) <= b[None, :] * 1e-5).any(axis=1).sum()))
    return best
E = pq[(pq.d == 0) & (pq.e == 1)].set_index('key')
D = pq[pq.d == 1].set_index('key')
pairs = E[['precursor_id', 'comp']].join(D[['precursor_id']].rename(columns={'precursor_id': 'did'}), how='inner')
pairs['ovE'] = [overlap(p, cp) for p, cp in zip(pairs.precursor_id, pairs.comp)]
pairs['ovD'] = [overlap(p, cp) for p, cp in zip(pairs.did, pairs.comp)]
print('pairs with an isobaric human parent:', int((pairs.ovE >= 0).sum()), 'of', len(pairs))
x = pairs[pairs.ovE >= 0]
print('overlap distribution (of 6 library fragments), E:', np.bincount(x.ovE, minlength=7).tolist(), ' DE:', np.bincount(x.ovD, minlength=7).tolist())
ov = pd.concat([pairs.set_index('precursor_id').ovE.rename('ov'), pairs.set_index('did').ovD.rename('ov')])
ov.to_csv('overlap.tsv', sep='\t')
cls = lambda v: 'none' if v < 0 else ('0-1' if v <= 1 else ('2' if v == 2 else '3+'))
tot = defaultdict(lambda: np.zeros(2, int))
for path in sys.argv[2:]:
    f = pd.read_csv(path, sep='\t', usecols=['precursor_id', 'is_decoy', 'score', 'q_value', 'network_score', 'network_q_value'])
    f = f[f.precursor_id.isin(ov.index)]
    for s, q in (('network_score', 'network_q_value'), ('score', 'q_value')):
        full = pd.read_csv(path, sep='\t', usecols=['is_decoy', s, q]) if False else None
        k = f[f[s] > -1e29].copy()
        k['oc'] = k.precursor_id.map(ov).map(cls)
        dv = np.sort(k[k.is_decoy == 1][s].values)[::-1]
        for lab, t in (('top40', dv[39]), ('top160', dv[159]), ('top640', dv[639])):
            a = k[k[s] >= t]
            for oc, y in a.groupby('oc'):
                tot[(s, lab, oc)] += [int((y.is_decoy == 0).sum()), int((y.is_decoy == 1).sum())]
print('E vs DE above the k-th decoy of the (E, DE) set, pooled, by fragment overlap with the isobaric human parent:')
for s in ('network_score', 'score'):
    for lab in ('top40', 'top160', 'top640'):
        print(f'  {s:14s} {lab:7s} ' + '  '.join(f'{oc}: E {tot[(s, lab, oc)][0]} DE {tot[(s, lab, oc)][1]}' for oc in ('none', '0-1', '2', '3+')))
