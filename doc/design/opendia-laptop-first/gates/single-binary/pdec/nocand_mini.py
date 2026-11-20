"""Build a mini PQP of sampled picker misses + found controls (targets only) for an OpenDIA search re-run with trace dump.
Diagnostic. Writes gap7/nocand/mini.pqp, ids.txt (precursor ids), sample.tsv (category, DIA-NN RT/IM)."""
import os, sqlite3, numpy as np, pandas as pd, runpy, sys
G = '/home/singjc/.jcode/scratch/dia2/gap7'; O = f'{G}/nocand'; os.makedirs(O, exist_ok=True)
src = open(f'{G}/precheck_nocand.py').read().split("print(f'DIA-NN IDs")[0]
g = {}; exec(compile(src, 'nocand_head', 'exec'), g)
d, tol, im_tol = g['d'], g['tol'], g['im_tol']
rng = np.random.default_rng(0)
inwin = (d.rt_err.abs() <= tol) & (d.im_err.abs() <= im_tol)
miss = d[~d.right & inwin & (d.n < 5)]; cap = d[~d.right & inwin & (d.n >= 5)]; found = d[d.right]
S = pd.concat([miss.sample(600, random_state=0).assign(cat='picker'), cap.sample(150, random_state=0).assign(cat='cap'),
               found.sample(300, random_state=0).assign(cat='found')])
S[['cat', 'rt_s', 'pred', 'IM', 'im_pred', 'n', 'qq', 'Precursor.Quantity']].to_csv(f'{O}/sample.tsv', sep='\t')
ids = S.index.astype(int).tolist()
L = g['L']
if os.path.exists(f'{O}/mini.pqp'): os.remove(f'{O}/mini.pqp')
c = sqlite3.connect(f'{O}/mini.pqp')
c.execute(f"attach database '{L}' as s")
for (n, sql) in c.execute("select name, sql from s.sqlite_master where type='table'").fetchall(): c.execute(sql)
c.execute('create temp table sel(id integer primary key)'); c.executemany('insert into sel values (?)', [(i,) for i in ids])
c.execute('insert into PRECURSOR select * from s.PRECURSOR where ID in (select id from sel)')
c.execute('insert into TRANSITION_PRECURSOR_MAPPING select * from s.TRANSITION_PRECURSOR_MAPPING where PRECURSOR_ID in (select id from sel)')
c.execute('insert into TRANSITION select * from s.TRANSITION where ID in (select TRANSITION_ID from TRANSITION_PRECURSOR_MAPPING)')
c.execute('insert into PRECURSOR_PEPTIDE_MAPPING select * from s.PRECURSOR_PEPTIDE_MAPPING where PRECURSOR_ID in (select id from sel)')
c.execute('insert into PEPTIDE select * from s.PEPTIDE where ID in (select PEPTIDE_ID from PRECURSOR_PEPTIDE_MAPPING)')
c.execute('insert into PEPTIDE_PROTEIN_MAPPING select * from s.PEPTIDE_PROTEIN_MAPPING where PEPTIDE_ID in (select ID from PEPTIDE)')
c.execute('insert into PROTEIN select * from s.PROTEIN where ID in (select PROTEIN_ID from PEPTIDE_PROTEIN_MAPPING)')
c.execute('insert into PEPTIDE_GENE_MAPPING select * from s.PEPTIDE_GENE_MAPPING where PEPTIDE_ID in (select ID from PEPTIDE)')
for t in ('VERSION', 'GENE'): c.execute(f'insert into {t} select * from s.{t}')
c.commit()
print({t: c.execute(f'select count(*) from {t}').fetchone()[0] for t in ('PRECURSOR', 'TRANSITION', 'PEPTIDE')})
open(f'{O}/ids.txt', 'w').write('\n'.join(str(i) for i in ids) + '\n')
print(S.cat.value_counts().to_dict())
