#!/usr/bin/env python3
"""Isobaric co-elution test. For every entrapment E (and its decoy DE as control) find the human parent T with the same
composition, charge and mod count (identical precursor m/z for E, mass-shifted for DE). Co-eluting = T's network-kept
window apex within 3 s of the item's kept apex (T at any score). Tail counts above the k-th decoy under each score,
split by co-elution, pooled over runs. If the GBM excess is isobaric signal sharing, it sits in the co-eluting E group."""
import sqlite3, sys, numpy as np, pandas as pd
lib = sys.argv[1]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID precursor_id, p.TRAML_ID t, p.DECOY d, p.CHARGE z, pe.UNMODIFIED_SEQUENCE seq, pe.MODIFIED_SEQUENCE mseq,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
pq['key'] = np.where(pq.d == 1, pq.t.str.replace('^DECOY_', '', regex=True), pq.t)
pq['comp'] = pq.seq.map(lambda s: ''.join(sorted(s))) + '|' + pq.z.astype(str) + '|' + pq.mseq.str.count('UniMod').astype(str)
P = pq.set_index('precursor_id')
ecomp = pq[(pq.d == 0) & (pq.e == 1)].set_index('key').comp  # entrapment key -> composition (DE shares the key)
ks = (10, 20, 40, 80, 160)
tot = {}
for path in sys.argv[2:]:
    f = pd.read_csv(path, sep='\t', usecols=['precursor_id', 'is_decoy', 'rt_apex', 'score', 'q_value', 'network_score', 'network_q_value'])
    nk = f[f.network_score > -1e29].set_index('precursor_id').join(P[['e', 'comp', 'd']])
    hum = nk[(nk.d == 0) & (nk.e == 0)].reset_index().groupby('comp').rt_apex.apply(np.array)
    for s, q in (('network_score', 'network_q_value'), ('score', 'q_value')):
        k = f[f[s] > -1e29].set_index('precursor_id').join(P[['e', 'key', 'd']])
        x = k[(k.d == 1) | (k.e == 1)].copy()
        x['comp'] = x.key.map(ecomp)
        x = x[x.comp.notna()]
        x['co'] = [cp in hum.index and np.min(np.abs(hum[cp] - rt)) <= 3 for cp, rt in zip(x.comp, x.rt_apex)]
        dv = np.sort(k[k.d == 1][s].values)[::-1]
        thr = b = k[(k.d == 0) & (k[q] <= 0.01)][s].min()
        for lab, t in [(f'D{n}', dv[n - 1]) for n in ks] + [('accepted', thr)]:
            a = x[x[s] >= t]
            r = tot.setdefault((s, lab), np.zeros(4, int))
            r += [int(((a.d == 0) & a.co).sum()), int(((a.d == 1) & a.co).sum()), int(((a.d == 0) & ~a.co).sum()), int(((a.d == 1) & ~a.co).sum())]
        r = tot.setdefault((s, 'all'), np.zeros(4, int))
        r += [int(((x.d == 0) & x.co).sum()), int(((x.d == 1) & x.co).sum()), int(((x.d == 0) & ~x.co).sum()), int(((x.d == 1) & ~x.co).sum())]
print(f'runs {len(sys.argv) - 2}; E = entrapments, DE = their own decoys; co = isobaric human parent co-elutes within 3 s')
for s in ('network_score', 'score'):
    print(s)
    for lab in ['all'] + [f'D{n}' for n in ks] + ['accepted']:
        a, b_, c_, d_ = tot[(s, lab)]
        print(f'   {lab:9s} co: E {a:6d} DE {b_:6d}   not co: E {c_:6d} DE {d_:6d}')
