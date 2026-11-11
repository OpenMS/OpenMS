"""Compare a DIALibGen library (compact parquet) with the original 6-fragment PQP on the same peptides.
Usage: lib_compare.py <new.parquet> [<fasta subset>]  (descriptive)."""
import sys, re, sqlite3, numpy as np, pandas as pd, pyarrow.parquet as pq
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
df = pq.read_table(sys.argv[1]).to_pandas()
def nd(s): return re.sub(r'\((Carbamidomethyl|UniMod:4)\)', '', re.sub(r'\((Oxidation|UniMod:35)\)', '(ox)', s))
df['k'] = df['Modified.Sequence'].map(nd) + '/' + df['Precursor.Charge'].astype(str)
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
keys = tuple(set(df['Modified.Sequence'].str.replace(r'\(.*?\)', '', regex=True)))
q = """select p.CHARGE z, pe.MODIFIED_SEQUENCE ms, pe.UNMODIFIED_SEQUENCE us, t.LIBRARY_INTENSITY li, t.PRODUCT_MZ fmz
 from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
 join TRANSITION_PRECURSOR_MAPPING tm on tm.PRECURSOR_ID=p.ID join TRANSITION t on t.ID=tm.TRANSITION_ID where p.DECOY=0"""
o = pd.read_sql(q, c); o = o[o.us.isin(set(keys))]; o['k'] = o.ms.map(nd) + '/' + o.z.astype(str); o = o[o.k.isin(set(df.k))]
nf = df.set_index('k'); hit = tot = 0; pos = []; err = []; top6 = 0; ntop = 0
for k, g in o.groupby('k'):
    r = nf.loc[k]; new = {round(mz, 3): (i, it) for i, (mz, it) in enumerate(zip(r['Product.Mz'], r['Relative.Intensity']))}
    ranks = []
    for _, x in g.iterrows():
        tot += 1; kk = round(x.fmz, 3)
        if kk in new: hit += 1; pos.append(new[kk][0]); err.append(abs(new[kk][1] - x.li)); ranks.append(new[kk][0])
    ntop += 1; top6 += int(len(ranks) == 6 and max(ranks) < 6)
print(f'{sys.argv[1]}: original fragments found {hit}/{tot} ({hit / tot:.3f}); within new top 6 {np.mean(np.array(pos) < 6):.3f}; '
      f'precursors with identical top-6 set {top6}/{ntop}; median |intensity diff| {np.median(err):.4f}')
