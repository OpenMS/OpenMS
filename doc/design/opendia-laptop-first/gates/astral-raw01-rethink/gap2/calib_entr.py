"""Why do entrapment precursors have fewer empty fragment traces and more candidates than human targets and decoys?
Design diagnostic. Compares precursor m/z, charge, length, library RT and fragment m/z/type across groups, and the
zero-trace fraction conditional on them. If entrapment is not exchangeable with false human targets, the entrapment
FDP estimate is biased (likely upward), and the 'anti-conservative' reading in CALIB.md must be qualified."""
import os, sqlite3, numpy as np, pandas as pd
RUN = os.environ.get('OPENDIA_RUN', '01')
S = os.environ['JCODE_SCRATCH_DIR'] + '/xic' + RUN
meta = pd.read_parquet(S + '/meta.parquet'); lib = np.load(S + '/lib.npy'); frag = np.load(S + '/frag.npy', mmap_mode='r')
om = pd.read_parquet(S + '/oswlike_meta.parquet')
c = sqlite3.connect('/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-entrapment-20261107/lib/library.pqp')
p = pd.read_sql("""select p.ID pid, p.PRECURSOR_MZ mz, p.CHARGE z, p.LIBRARY_RT lrt, pe.UNMODIFIED_SEQUENCE seq, pe.MODIFIED_SEQUENCE mseq
 from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID""", c).drop_duplicates('pid').set_index('pid')
tr = pd.read_sql("""select m.PRECURSOR_ID pid, t.PRODUCT_MZ fmz, t.TYPE ty, t.ORDINAL ord, t.CHARGE fz, t.LIBRARY_INTENSITY li from TRANSITION t
 join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID""", c)
ent = om.groupby('precursor_id').ent.first().reindex(meta.precursor_id).fillna(0).astype(bool).values
dec = meta.decoy.values.astype(bool)
g = pd.Series(np.where(dec, 'decoy', np.where(ent, 'entrapment', 'human')), index=meta.precursor_id.values)
zero = np.zeros(len(meta))
for b in range(0, len(meta), 50000):
    f = np.asarray(frag[b:b + 50000]); l = lib[b:b + 50000]
    zero[b:b + 50000] = ((f.max(2) <= 0) & (l > 0)).sum(1) / np.maximum((l > 0).sum(1), 1)
d = p.reindex(meta.precursor_id.values).copy(); d['grp'] = g.values; d['zero'] = zero; d['len'] = d.seq.str.len()
d['oxM'] = d.mseq.str.contains('UniMod:35'); d['C'] = d.seq.str.contains('C')
print(f'raw {RUN}')
print(d.groupby('grp')[['mz', 'z', 'len', 'lrt', 'oxM', 'C']].mean().round(3).to_string())
print('charge share:'); print(pd.crosstab(d.grp, d.z, normalize='index').round(3).to_string())
t = tr.merge(g.rename('grp'), left_on='pid', right_index=True)
print('fragments:'); print(t.groupby('grp').agg(fmz=('fmz', 'mean'), y_share=('ty', lambda s: (s == 'y').mean()), ord=('ord', 'mean'), fz2=('fz', lambda s: (s > 1).mean())).round(3).to_string())
d['mzbin'] = pd.cut(d.mz, [400, 500, 600, 700, 800, 1000]); d['lbin'] = pd.cut(d.len, [6, 10, 14, 18, 24, 36])
print('zero-trace fraction by precursor m/z bin:'); print(d.pivot_table(index='mzbin', columns='grp', values='zero', observed=True).round(3).to_string())
print('zero-trace fraction by length bin:'); print(d.pivot_table(index='lbin', columns='grp', values='zero', observed=True).round(3).to_string())
print('zero-trace fraction by charge:'); print(d.pivot_table(index='z', columns='grp', values='zero', observed=True).round(3).to_string())
