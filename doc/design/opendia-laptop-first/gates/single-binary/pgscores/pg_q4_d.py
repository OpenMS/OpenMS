"""q4 precheck on D: among IDs NEW adds over CUR (seed 0), compare added-in-DIA-NN vs added-not-in-DIA-NN, and decoy winners, per pg score."""
import sqlite3, glob, pandas as pd, numpy as np
from sklearn.metrics import roc_auc_score
H='/home/singjc/Documents/github/testing_auto_irt/diaPASEF_evosep_published'
lib=sqlite3.connect(f'file:{H}/200623_hela_evosep_fractions_library_ptypic_decoy.pqp?mode=ro',uri=True)
pre=pd.read_sql("select p.ID pid, p.CHARGE z, pe.MODIFIED_SEQUENCE ms from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID",lib)
key=dict(zip(pre.pid,pre.ms+'/'+pre.z.astype(str)))
dn=pd.read_parquet('diann_hela/diann_report.parquet',columns=['Modified.Sequence','Precursor.Charge','Q.Value'])
dn=dn[dn['Q.Value']<=0.01]; DN=set(dn['Modified.Sequence']+'/'+dn['Precursor.Charge'].astype(str))
cur=pd.read_csv('rtc/D_imgrid_s0/precursors.tsv',sep='\t'); C=set(cur[cur.q_value<=0.01].precursor_id)
F=pd.read_csv(glob.glob('pg/D_new_s0/*.features.tsv')[0],sep='\t')
w=F.sort_values('score',ascending=False).drop_duplicates('precursor_id')
w['key']=w.precursor_id.map(key)
acc=w[(w.q_value<=0.01)&(w.is_decoy==0)]
add=acc[~acc.precursor_id.isin(C)]
A_in=add[add.key.isin(DN)]; A_out=add[~add.key.isin(DN)]
old=acc[acc.precursor_id.isin(C)&acc.key.isin(DN)]
dec=w[w.is_decoy==1]; dtop=dec.nlargest(len(add),'score')
print(f'added {len(add)}: in DIA-NN {len(A_in)}, not {len(A_out)}; top decoys {len(dtop)}')
pg=[c for c in F.columns if c.startswith('pg_')]; base=[c for c in F.columns if c not in pg and c not in ('candidate_row','precursor_id','is_decoy','fold','rt_apex','score','q_value')]
print(f'{"feature":24s} {"AUC in>out":>10s} {"AUC out>topdec":>14s} {"med in":>8s} {"med out":>8s} {"med topdec":>10s} {"med old":>8s}')
for c in pg+base:
    a1=roc_auc_score(np.r_[np.ones(len(A_in)),np.zeros(len(A_out))],np.r_[A_in[c],A_out[c]])
    a2=roc_auc_score(np.r_[np.ones(len(A_out)),np.zeros(len(dtop))],np.r_[A_out[c],dtop[c]])
    print(f'{c:24s} {a1:10.3f} {a2:14.3f} {A_in[c].median():8.3f} {A_out[c].median():8.3f} {dtop[c].median():10.3f} {old[c].median():8.3f}')
