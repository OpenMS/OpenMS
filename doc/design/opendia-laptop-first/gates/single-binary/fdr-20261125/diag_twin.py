#!/usr/bin/env python3
"""Twin test. Each entrapment E is a K/R- and C-term-preserving shuffle of a human peptide T (same composition, so same
precursor m/z at the same charge). E competes against its own decoy DE (mass-shifted pseudo-reverse). If T is present
and abundant, E can borrow T's MS1/fragment signal while DE cannot, so E and DE stop being exchangeable.
Per score: E-vs-DE win share and tail counts, split by whether the isobaric human twin is an ID (network q <= 0.01,
fixed for both scores) and by |RT E - RT T|. Diagnosis only."""
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
for path in sys.argv[2:]:
    f = pd.read_csv(path, sep='\t', usecols=['precursor_id', 'is_decoy', 'rt_apex', 'score', 'q_value', 'network_score', 'network_q_value'])
    print(f'== {path.split("/")[0]}')
    nk = f[f.network_score > -1e29].set_index('precursor_id').join(P[['e', 'comp']])
    hid = nk[(nk.is_decoy == 0) & (nk.e == 0) & (nk.network_q_value <= 0.01)]
    twin_rt = hid.reset_index().groupby('comp').rt_apex.apply(np.array)
    for s, q in (('network_score', 'network_q_value'), ('score', 'q_value')):
        k = f[f[s] > -1e29].set_index('precursor_id').join(P[['e', 'comp', 'key', 'd']])
        E = k[(k.d == 0) & (k.e == 1)][['key', 'comp', s, 'rt_apex', q]]
        DE = k[k.d == 1].set_index('key')[[s]].rename(columns={s: 'sd'})
        x = E.join(DE, on='key', how='inner')
        x['twin'] = x.comp.isin(twin_rt.index)
        x['dt'] = [np.min(np.abs(twin_rt[cp] - rt)) if tw else np.nan for cp, rt, tw in zip(x.comp, x.rt_apex, x.twin)]
        dv = np.sort(k[k.d == 1][s].values)[::-1]
        t80, t320 = dv[79], dv[319]
        print(f'  {s}:')
        for lab, m in (('twin is an ID', x.twin), ('  twin |dt|<=10 s', x.twin & (x.dt <= 10)), ('  twin 10<|dt|<=60 s', x.twin & (x.dt > 10) & (x.dt <= 60)),
                       ('  twin |dt|>60 s', x.twin & (x.dt > 60)), ('no twin ID', ~x.twin)):
            y = x[m]
            ew = np.mean(y[s] > y.sd)
            e80, d80 = int((y[s] >= t80).sum()), int((y.sd >= t80).sum())
            e320, d320 = int((y[s] >= t320).sum()), int((y.sd >= t320).sum())
            print(f'    {lab:22s} pairs {len(y):6d}  E beats DE {ew:.3f}  above D80: E {e80:3d} DE {d80:3d}  above D320: E {e320:3d} DE {d320:3d}  E q<=0.01 {int((y[q] <= 0.01).sum())}')
