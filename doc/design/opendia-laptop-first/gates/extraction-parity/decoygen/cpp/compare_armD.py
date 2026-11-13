"""Acceptance check: C++ OpenSwathDecoyGenerator -predict_decoys (cppD.pqp) against the arm D library of PLAN_DECOYGEN
(DIALibGen-predicted pseudo-reverse + K/R swap decoys). Decoys are paired through their target (DECOY_<target TRAML_ID>
mapped to the target's modified sequence and charge). Reports: shared decoys, identical decoy sequence and precursor m/z,
RT difference, and per-decoy fragment-set overlap and intensity cosine."""
import sqlite3, numpy as np, pandas as pd


def load(path):
    c = sqlite3.connect(f'file:{path}?mode=ro', uri=True)
    pr = pd.read_sql("""select p.ID id, p.TRAML_ID traml, p.CHARGE z, p.PRECURSOR_MZ mz, p.LIBRARY_RT rt, p.DECOY d,
      pe.MODIFIED_SEQUENCE ms from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
      join PEPTIDE pe on pe.ID=m.PEPTIDE_ID""", c)
    tr = pd.read_sql("""select m.PRECURSOR_ID id, t.TYPE ion, t.ORDINAL o, t.CHARGE fz, t.LIBRARY_INTENSITY li
      from TRANSITION t join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID where t.DECOY=1""", c)
    tgt = pr[pr.d == 0].set_index('traml')
    dec = pr[pr.d == 1].copy()
    t = tgt.loc[dec.traml.str.replace('^DECOY_', '', regex=True)]
    dec['key'] = t.ms.values + '/' + t.z.astype(str).values
    return dec.set_index('key'), {i: g for i, g in tr.groupby('id')}


import sys
a, at = load(sys.argv[1] if len(sys.argv) > 1 else "cppD.pqp")
b, bt = load('../armD/home/lib/library.pqp')
com = a.index.intersection(b.index)
print(f'decoys: C++ {len(a)}, arm D {len(b)}, shared targets {len(com)}')
print('identical decoy sequence %.4f, precursor m/z within 1e-3 %.4f' % (
    (a.loc[com, 'ms'] == b.loc[com, 'ms']).mean(), (abs(a.loc[com, 'mz'] - b.loc[com, 'mz']) < 1e-3).mean()))
d = abs(a.loc[com, 'rt'] - b.loc[com, 'rt'])
print('normalised RT |difference|: median %.2e, p99 %.2e, max %.2e' % (np.median(d), np.quantile(d, .99), d.max()))
rng = np.random.default_rng(0)
same, cos = [], []
for k in com[rng.choice(len(com), min(20000, len(com)), replace=False)]:
    x, y = at[a.loc[k, 'id']], bt[b.loc[k, 'id']]
    kx = dict(zip(zip(x.ion, x.o, x.fz), x.li)); ky = dict(zip(zip(y.ion, y.o, y.fz), y.li))
    same.append(len(set(kx) & set(ky)))
    keys = sorted(set(kx) | set(ky))
    u = np.array([kx.get(q, 0.0) for q in keys]); v = np.array([ky.get(q, 0.0) for q in keys])
    cos.append(u @ v / np.linalg.norm(u) / np.linalg.norm(v))
same, cos = np.array(same), np.array(cos)
print('shared fragments (of 6), share of decoys:', dict(enumerate(np.round(np.bincount(same, minlength=7) / len(same), 4))))
print('intensity cosine over the union: median %.4f, q01 %.4f' % (np.median(cos), np.quantile(cos, .01)))
