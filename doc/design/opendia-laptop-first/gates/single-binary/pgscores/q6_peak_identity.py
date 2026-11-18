"""q6 (diagnostic, D seed 0 NEW features): does each feature separate identity, or only 'a real peak is here'?
Groups (precursor winners):
  T    true: CUR-accepted targets confirmed by DIA-NN 1 % at the same RT (+-3 s)
  Dpk  peak present, wrong identity: decoy winners whose apex is within 3 s of their paired target's DIA-NN 1 % RT
  Doff no real peak: decoy winners away (> 30 s) from their paired target's DIA-NN RT, or target not in DIA-NN
  X    NEW's 694 additions not in DIA-NN 1 % (mostly false by PLAUS_D)
identity AUC = AUC(T > Dpk); peak AUC = AUC(Dpk > Doff) (feature rewards a real peak regardless of identity)."""
import sqlite3, glob, pandas as pd, numpy as np
from sklearn.metrics import roc_auc_score
H='/home/singjc/Documents/github/testing_auto_irt/diaPASEF_evosep_published'
lib=sqlite3.connect(f'file:{H}/200623_hela_evosep_fractions_library_ptypic_decoy.pqp?mode=ro',uri=True)
pre=pd.read_sql("select p.ID pid, p.TRAML_ID t, p.CHARGE z, p.DECOY decoy, pe.MODIFIED_SEQUENCE ms from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID",lib)
pre['key']=pre.ms+'/'+pre.z.astype(str); tid=pre.set_index('t').pid
pre['pair']=np.where(pre.decoy==1,pre.t.str.replace('^DECOY_','',regex=True).map(tid),pre.pid)
P=pre.set_index('pid')
dn=pd.read_parquet('diann_hela/diann_report.parquet',columns=['Modified.Sequence','Precursor.Charge','Q.Value','RT'])
dn=dn[dn['Q.Value']<=0.01]; dn['key']=dn['Modified.Sequence']+'/'+dn['Precursor.Charge'].astype(str); drt=dn.drop_duplicates('key').set_index('key').RT*60
F=pd.read_csv(glob.glob('pg/D_new_s0/*.features.tsv')[0],sep='\t')
w=F.sort_values('score',ascending=False).drop_duplicates('precursor_id').copy()
w['tkey']=w.precursor_id.map(P.pair).map(P.key); w['dnrt']=w.tkey.map(drt); dt=np.abs(w.rt_apex-w.dnrt)
cur=pd.read_csv('rtc/D_imgrid_s0/precursors.tsv',sep='\t'); C=set(cur[cur.q_value<=0.01].precursor_id)
T=w[(w.is_decoy==0)&w.precursor_id.isin(C)&(dt<=3)]
Dpk=w[(w.is_decoy==1)&(dt<=3)]; Doff=w[(w.is_decoy==1)&((dt>30)|w.dnrt.isna())]
X=w[(w.is_decoy==0)&(w.q_value<=0.01)&~w.precursor_id.isin(C)&w.dnrt.isna()]
D=w[w.is_decoy==1]
print(f'T {len(T)}  Dpk {len(Dpk)}  Doff {len(Doff)}  X {len(X)} | share of decoy winners on their target peak: {len(Dpk)/len(D):.3f}')
# Where does each group sit relative to the cut? X and decoys in the score band of X
lo=X.score.min(); band=D[D.score>=lo]
print(f'decoys above the lowest X score: {len(band)}, of which on target peak: {(np.abs(band.rt_apex-band.dnrt)<=3).mean():.2f}')
meta=('candidate_row','precursor_id','is_decoy','fold','rt_apex','score','q_value','tkey','dnrt')
cols=[c for c in F.columns if c not in meta]
def auc(a,b,c):
    return roc_auc_score(np.r_[np.ones(len(a)),np.zeros(len(b))],np.r_[a[c].fillna(0),b[c].fillna(0)])
rows=[(c,auc(T,Dpk,c),auc(Dpk,Doff,c),auc(T,Doff,c),auc(X,Dpk,c)) for c in cols]
df=pd.DataFrame(rows,columns=['feature','identity T>Dpk','peak Dpk>Doff','overall T>Doff','X>Dpk'])
fold=lambda s:(s-0.5).abs()
df['identity']=fold(df['identity T>Dpk']); df['peak']=fold(df['peak Dpk>Doff'])
df['kind']=np.where(df.identity>=df.peak,'identity','PEAK-ONLY')
df['pg']=df.feature.str.startswith('pg_')
print(df.sort_values(['pg','identity'],ascending=[False,False]).round(3).to_string(index=False))
# Score itself
for name,g in (('T',T),('Dpk',Dpk),('Doff',Doff),('X',X)): print(f'score median {name}: {g.score.median():.3f}')
