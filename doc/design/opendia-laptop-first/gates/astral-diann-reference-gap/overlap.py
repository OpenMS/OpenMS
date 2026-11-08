"""How many ProteoBench DIA-NN 1% IDs exist in our PQP library, and what our raw 02 pipeline did with them."""
import sqlite3, pandas as pd, numpy as np, re
B='/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
LIB=B+'/astral-entrapment-20261107/lib/library.pqp'
D='/home/singjc/Documents/datasets/proteobench-astral-entrapment/proteobench_diann_v231_res_20260717_123332/'
con=sqlite3.connect(LIB)
lib=pd.read_sql('select pr.ID precursor_id, pr.DECOY decoy, pr.CHARGE charge, pr.LIBRARY_RT lib_rt, pr.PRECURSOR_MZ mz, p.UNMODIFIED_SEQUENCE seq, p.MODIFIED_SEQUENCE modseq from PRECURSOR pr join PRECURSOR_PEPTIDE_MAPPING pp on pp.PRECURSOR_ID=pr.ID join PEPTIDE p on p.ID=pp.PEPTIDE_ID',con)
nt=pd.read_sql('select PRECURSOR_ID precursor_id, count(*) n from TRANSITION_PRECURSOR_MAPPING group by PRECURSOR_ID',con)
con.close()
lib=lib.merge(nt,on='precursor_id')
t=lib[lib.decoy==0]
print('library targets',len(t),'unique seq',t.seq.nunique(),'charges',t.charge.value_counts().sort_index().to_dict())
print('transitions per target precursor', t.n.describe()[['min','50%','max']].to_dict())
print('modseq examples', t.modseq.head(3).tolist())
d=pd.read_csv(D+'result_performance.csv')
d=d[d['Q-Value']<=0.01]
t=t.assign(key=t.seq+'/'+t.charge.astype(str))
d['key']=d.Peptide+'/'+d.Charge.astype(str)
keys=set(t.key)
for rf,g in d.groupby('Raw file'):
    inl=g.key.isin(keys)
    print(rf,'DIA-NN 1%',len(g),'in our library',int(inl.sum()), 'entrapment in lib', int((inl&(g['Target or Entrapment']=='entrapment')).sum()))
d.to_parquet('diann_1pct.parquet'); t.to_parquet('lib_targets.parquet')
