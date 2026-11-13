"""PLAN_DECOYGEN acceptance arm E: run root for the C++ OpenSwathDecoyGenerator -predict_decoys library (cpp/cppD.pqp).
lib/library.pqp is cppD.pqp; lib/pred12.parquet is the PAIRLIB target prediction; lib/decoy12.tsv gives each C++ decoy its
own 12-fragment prediction from lib/decoys_p12.parquet (the pred12 settings), keyed by its decoy modified sequence and
charge. Decoys whose sequence is not in decoys_p12 keep their 6 library transitions in build_pqp12.py (count printed)."""
import os, shutil, sqlite3, pandas as pd, pyarrow.parquet as pq
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
D = B + '/decoygen-20261113'
O = D + '/armE/home/lib'
os.makedirs(O, exist_ok=True)
shutil.copyfile(D + '/cpp/cppD.pqp', O + '/library.pqp')
shutil.copyfile(B + '/pairlib-20261113/lib/pred12.parquet', O + '/pred12.parquet')
c = sqlite3.connect(f'file:{O}/library.pqp?mode=ro', uri=True)
dec = pd.read_sql("""select p.TRAML_ID traml, p.CHARGE z, pe.MODIFIED_SEQUENCE ms from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m
  on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID where p.DECOY=1""", c)
norm = lambda s: s.replace('(Carbamidomethyl)', '(UniMod:4)').replace('(Oxidation)', '(UniMod:35)')
p = pq.read_table(D + '/lib/decoys_p12.parquet').to_pandas()
p['key'] = p['Modified.Sequence'].map(norm) + '/' + p['Precursor.Charge'].astype(str)
p = p.drop_duplicates('key').set_index('key')
ION = {2: 'b', 5: 'y'}
rows = []; miss = 0
for r in dec.itertuples(index=False):
    k = f'{r.ms}/{r.z}'
    if k not in p.index:
        miss += 1; continue
    q = p.loc[k]
    for i, (mz, li, ft, fz, o) in enumerate(zip(q['Product.Mz'], q['Relative.Intensity'], q['Fragment.Type'], q['Fragment.Charge'],
                                                q['Fragment.Series.Number'])):
        rows.append((r.traml, r.ms, i, ION[int(ft)], int(o), int(fz), float(mz), float(li)))
pd.DataFrame(rows, columns=['traml', 'decoy_ms', 'rank', 'ion', 'ord', 'fz', 'mz', 'li']).to_csv(O + '/decoy12.tsv', sep='\t', index=False)
msg = f'decoys {len(dec)}, with a 12-fragment prediction {len(dec) - miss}, without {miss}'
print(msg); open(O + '/build_lib.txt', 'w').write(msg + '\n')
