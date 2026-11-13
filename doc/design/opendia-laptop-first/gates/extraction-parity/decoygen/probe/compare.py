"""Compare ProteoBench library target intensities with DIALibGen predictions under several settings (probe)."""
import sqlite3, sys, glob, pandas as pd, numpy as np, pyarrow.parquet as pq
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
c = sqlite3.connect(f'file:{B}/pairlib-20261113/lib/library.pqp?mode=ro', uri=True)
prec = pd.read_sql("""select p.ID id,p.CHARGE z,pe.MODIFIED_SEQUENCE ms from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m
  on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID where p.DECOY=0""", c)
tr = pd.read_sql("""select m.PRECURSOR_ID id,t.TYPE ion,t.ORDINAL ord,t.CHARGE fz,t.LIBRARY_INTENSITY li from TRANSITION t
  join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID where t.DECOY=0""", c)
prec['key'] = prec.ms + '/' + prec.z.astype(str)
norm = lambda s: s.replace('(Carbamidomethyl)', '(UniMod:4)').replace('(Oxidation)', '(UniMod:35)')
for f in sorted(glob.glob(f'{B}/decoygen-20261113/probe/p_*.parquet')) + [f'{B}/pairlib-20261113/lib/pred12.parquet']:
    p = pq.read_table(f).to_pandas()
    p['key'] = p['Modified.Sequence'].map(norm) + '/' + p['Precursor.Charge'].astype(str)
    p = p.drop_duplicates('key').set_index('key')
    s = prec[prec.key.isin(p.index)]
    if len(s) > 3000: s = s.sample(3000, random_state=0)
    g = tr[tr.id.isin(s.id)].groupby('id')
    cos = []; top6 = []; dif = []
    for r in s.itertuples():
        q = p.loc[r.key]
        keys = [({2: 'b', 5: 'y'}[int(a)], int(o), int(z)) for a, o, z in zip(q['Fragment.Type'], q['Fragment.Series.Number'], q['Fragment.Charge'])]
        d = dict(zip(keys, q['Relative.Intensity'].astype(float)))
        t = g.get_group(r.id)
        lib = t.li.values.astype(float); pr = np.array([d.get(k, 0.0) for k in zip(t.ion, t.ord, t.fz)])
        cos.append(lib @ pr / np.linalg.norm(lib) / max(np.linalg.norm(pr), 1e-12))
        top6.append(len(set(sorted(d, key=d.get, reverse=True)[:6]) & set(zip(t.ion, t.ord, t.fz))))
        m = pr > 0
        dif.append(np.median(abs(lib[m] - pr[m])) if m.any() else np.nan)
    cos = np.array(cos); top6 = np.array(top6)
    print(f.split('/')[-1], 'n', len(s), 'cos median %.4f q10 %.4f' % (np.nanmedian(cos), np.nanquantile(cos, .1)),
          'top6==6 %.3f' % (top6 == 6).mean(), 'median |dI| %.4f' % np.nanmedian(dif))
