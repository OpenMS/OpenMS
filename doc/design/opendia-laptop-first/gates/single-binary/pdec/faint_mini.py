"""PRECHECK_FAINT step 1: sample for the faint-peak evidence test (r~1.56 D s0, valley 0.8 / max 8 arm ev/v08m8).
6,000 DIA-NN non-entrapment q<=0.01 targets (labels for evaluation only) + 30,000 random decoy precursors.
Writes gap7/faint/mini.pqp, ids.txt, sample.tsv."""
import os, glob, sqlite3, numpy as np, pandas as pd
G = '/home/singjc/.jcode/scratch/dia2/gap7'; O = f'{G}/faint'; os.makedirs(O, exist_ok=True)
L = '/home/singjc/.jcode/scratch/dia2/ctest_cache/c/hela_entrap_codon.predicted_decoys.pqp'
DN = '/home/singjc/.jcode/scratch/dia2/entlib12/D/report.parquet'
c = sqlite3.connect(f'file:{L}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID pid, p.DECOY d, p.CHARGE z, pe.MODIFIED_SEQUENCE ms, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
 from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
 join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
c.close()
lab = pq[(pq.d == 0) & (pq.e == 0)].assign(k=lambda x: x.ms + '/' + x.z.astype(str)).drop_duplicates('k').set_index('k')
d = pd.read_parquet(DN, columns=['Modified.Sequence', 'Precursor.Charge', 'Q.Value', 'RT', 'Precursor.Quantity'])
d = d[d['Q.Value'] <= 0.01]; k = d['Modified.Sequence'] + '/' + d['Precursor.Charge'].astype(str)
d = d[k.isin(lab.index).values]; d['pid'] = lab.pid.reindex(k[k.isin(lab.index)]).values; d = d.drop_duplicates('pid').set_index('pid')
d['qq'] = pd.qcut(np.log10(d['Precursor.Quantity'].clip(lower=1)), 4, labels=['Q1', 'Q2', 'Q3', 'Q4'])
F = pd.read_csv(glob.glob(f'{G}/ev/v08m8/*.precursors.tsv')[0], sep='\t', usecols=['precursor_id', 'is_decoy'])
searched = set(F.precursor_id)
tg = d[d.index.isin(searched)].sample(6000, random_state=0)
dec = pq[(pq.d == 1) & pq.pid.isin(searched)].sample(30000, random_state=0).pid
S = pd.concat([pd.DataFrame({'pid': tg.index, 'cls': 'dnn', 'rt_s': tg.RT.values * 60, 'qq': tg.qq.astype(str).values}),
               pd.DataFrame({'pid': dec.values, 'cls': 'decoy', 'rt_s': np.nan, 'qq': ''})])
S.to_csv(f'{O}/sample.tsv', sep='\t', index=False)
ids = S.pid.astype(int).tolist()
if os.path.exists(f'{O}/mini.pqp'): os.remove(f'{O}/mini.pqp')
c = sqlite3.connect(f'{O}/mini.pqp'); c.execute(f"attach database '{L}' as s")
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
print({t: c.execute(f'select count(*) from {t}').fetchone()[0] for t in ('PRECURSOR', 'TRANSITION')})
open(f'{O}/ids.txt', 'w').write('\n'.join(map(str, ids)) + '\n')
print(S.cls.value_counts().to_dict(), tg.qq.value_counts().to_dict())
