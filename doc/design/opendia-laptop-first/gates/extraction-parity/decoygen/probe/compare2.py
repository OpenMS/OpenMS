"""Probe 2 evaluation: share of library targets whose 6 transitions equal DIALibGen's top 6 (per setting), and intensity cosine."""
import sqlite3, glob, pandas as pd, numpy as np, pyarrow.parquet as pq
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
c = sqlite3.connect(f'file:{B}/pairlib-20261113/lib/library.pqp?mode=ro', uri=True)
prec = pd.read_sql("""select p.ID id,p.CHARGE z,pe.MODIFIED_SEQUENCE ms from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m
  on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID where p.DECOY=0""", c)
tr = pd.read_sql("""select m.PRECURSOR_ID id,t.TYPE ion,t.ORDINAL ord,t.CHARGE fz,t.LIBRARY_INTENSITY li from TRANSITION t
  join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID where t.DECOY=0""", c)
prec['key'] = prec.ms + '/' + prec.z.astype(str)
norm = lambda s: s.replace('(Carbamidomethyl)', '(UniMod:4)').replace('(Oxidation)', '(UniMod:35)')
for f in sorted(glob.glob(f'{B}/decoygen-20261113/probe/q_*.parquet')):
    p = pq.read_table(f).to_pandas()
    p['key'] = p['Modified.Sequence'].map(norm) + '/' + p['Precursor.Charge'].astype(str); p = p.set_index('key')
    s = prec[prec.key.isin(p.index)]; g = tr[tr.id.isin(s.id)].groupby('id')
    eq = []; cos = []
    for r in s.itertuples():
        q = p.loc[r.key]
        d = {({2: 'b', 5: 'y'}[int(a)], int(o), int(z)): float(i) for a, o, z, i in
             zip(q['Fragment.Type'], q['Fragment.Series.Number'], q['Fragment.Charge'], q['Relative.Intensity'])}
        t = g.get_group(r.id); k = list(zip(t.ion, t.ord, t.fz))
        eq.append(len(set(d) & set(k)))
        lib = t.li.values.astype(float); pr = np.array([d.get(x, 0.0) for x in k])
        if lib.any() and pr.any(): cos.append(lib @ pr / np.linalg.norm(lib) / np.linalg.norm(pr))
    eq = np.array(eq)
    print(f.split('/')[-1], 'n', len(eq), 'all6 %.3f' % (eq == 6).mean(), '>=5 %.3f' % (eq >= 5).mean(), 'cos median %.4f' % np.median(cos))
