import sqlite3, pandas as pd, numpy as np
H='/home/singjc/Documents/github/testing_auto_irt/diaPASEF_evosep_published'
lib=sqlite3.connect(f'file:{H}/200623_hela_evosep_fractions_library_ptypic_decoy.pqp?mode=ro',uri=True)
pre=pd.read_sql("select p.ID pid, p.CHARGE z, pe.MODIFIED_SEQUENCE ms from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID",lib)
key=dict(zip(pre.pid,pre.ms+'/'+pre.z.astype(str)))
dn=pd.read_parquet('diann_hela/diann_report.parquet',columns=['Modified.Sequence','Precursor.Charge','Q.Value'])
D=set((dn[dn['Q.Value']<=0.01]['Modified.Sequence']+'/'+dn[dn['Q.Value']<=0.01]['Precursor.Charge'].astype(str)))
for a in ['cur','new']:
    for s in [0]:
        t=pd.read_csv(f'pgeval/D_{a}_s{s}/precursors.tsv',sep='\t'); O=set(t[t.q_value<=0.01].precursor_id.map(key))
        print(a,s,'N',len(O),'in DIA-NN',len(O&D),'not',len(O-D))
c=sqlite3.connect('file:/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/openswath-hela-published/svm.osw?mode=ro',uri=True)
o=pd.read_sql("select f.PRECURSOR_ID pid, s.QVALUE q from SCORE_MS2 s join FEATURE f on f.ID=s.FEATURE_ID join PRECURSOR p on p.ID=f.PRECURSOR_ID where s.RANK=1 and p.DECOY=0",c)
O=set(o[o.q<=0.01].pid.map(key)); print('OpenSWATH SVM N',len(O),'agree',round(len(O&D)/len(O),4),'not',len(O-D))
