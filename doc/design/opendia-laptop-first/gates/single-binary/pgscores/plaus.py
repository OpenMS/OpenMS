"""Are the IDs NEW adds on HeLa D (seed 0), and that DIA-NN does not report at 1 %, plausible? Evidence per ID:
DIA-NN at any q (and at which RT), OpenSWATH SVM q<=0.01 at the same RT, DIA-NN on the other HeLa run (H), same peptide
at another charge, and protein support. Compared with the added IDs DIA-NN does confirm, and with top decoys as a null."""
import sqlite3, glob, pandas as pd, numpy as np
H='/home/singjc/Documents/github/testing_auto_irt/diaPASEF_evosep_published'
lib=sqlite3.connect(f'file:{H}/200623_hela_evosep_fractions_library_ptypic_decoy.pqp?mode=ro',uri=True)
pre=pd.read_sql("select p.ID pid, p.CHARGE z, p.DECOY decoy, pe.MODIFIED_SEQUENCE ms, pe.UNMODIFIED_SEQUENCE seq from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID",lib)
pre['key']=pre.ms+'/'+pre.z.astype(str); P=pre.set_index('pid')
def dn(path):
    d=pd.read_parquet(path,columns=['Modified.Sequence','Precursor.Charge','Stripped.Sequence','Q.Value','RT','Protein.Group'])
    d['key']=d['Modified.Sequence']+'/'+d['Precursor.Charge'].astype(str); return d.sort_values('Q.Value').drop_duplicates('key').set_index('key')
DD=dn('diann_hela/diann_report.parquet'); DH=dn('diann_ho/diann_report.parquet')
D1=DD[DD['Q.Value']<=0.01]; pep1=set(D1['Stripped.Sequence']); prot1=set(D1['Protein.Group'])
c=sqlite3.connect('file:/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/openswath-hela-published/svm.osw?mode=ro',uri=True)
osw=pd.read_sql("select f.PRECURSOR_ID pid, f.EXP_RT rt, s.QVALUE q from SCORE_MS2 s join FEATURE f on f.ID=s.FEATURE_ID where s.RANK=1",c).set_index('pid')
cur=pd.read_csv('rtc/D_imgrid_s0/precursors.tsv',sep='\t'); C=set(cur[cur.q_value<=0.01].precursor_id)
w=pd.read_csv(glob.glob('pg/D_new_s0/*.precursors.tsv')[0],sep='\t').sort_values('score',ascending=False).drop_duplicates('precursor_id')
w['key']=w.precursor_id.map(P.key)
acc=w[(w.q_value<=0.01)&(w.is_decoy==0)]; add=acc[~acc.precursor_id.isin(C)]
groups={'added, in DIA-NN 1 %':add[add.key.isin(D1.index)],'added, NOT in DIA-NN 1 %':add[~add.key.isin(D1.index)],
        'old (CUR) accepted':acc[acc.precursor_id.isin(C)]}
# null: best decoys just below the cut (as many as the not-in-DIA-NN group), mapped to their target's key for lookups
dec=w[w.is_decoy==1].copy(); n0=len(groups['added, NOT in DIA-NN 1 %'])
rows=[]
def ev(g,name):
    k=g.key; d=DD.reindex(k)
    has=d['Q.Value'].notna().values; rt_ok=has&(np.abs(d.RT.values*60-g.rt_apex.values)<=10)
    o=osw.reindex(g.precursor_id); o_ok=(o.q.values<=0.01)&(np.abs(o.rt.values-g.rt_apex.values)<=10)
    h=DH.reindex(k)['Q.Value'].values<=0.01
    seq=g.precursor_id.map(P.seq); other_z=seq.isin(pep1).values & ~k.isin(D1.index).values
    rows.append(dict(group=name,n=len(g),
        dn_any_q=has.mean(), dn_same_rt=rt_ok.mean(), dn_q05_same_rt=(rt_ok&(d['Q.Value'].values<=0.05)).mean(),
        dn_med_q=np.nanmedian(d['Q.Value'].values),
        osw_same_rt=o_ok.mean(), diann_H_1pct=h.mean(), pep_other_charge=other_z.mean()))
for n,g in groups.items(): ev(g,n)
print(pd.DataFrame(rows).round(3).to_string(index=False))
nd=groups['added, NOT in DIA-NN 1 %']; d=DD.reindex(nd.key)
print('\nNOT-in-DIA-NN group, DIA-NN q distribution (where reported):'); print(pd.cut(d['Q.Value'].dropna(),[0,0.01,0.02,0.05,0.1,0.2,0.5]).value_counts().sort_index().to_string())
rtd=(nd.rt_apex.values-d.RT.values*60); rtd=rtd[~np.isnan(rtd)]
print('|our apex - DIA-NN RT| (s): median %.1f, share <=10 s %.2f, share >60 s %.2f'%(np.median(np.abs(rtd)),(np.abs(rtd)<=10).mean(),(np.abs(rtd)>60).mean()))
nd[['precursor_id','key','rt_apex','score','q_value']].assign(dn_q=d['Q.Value'].values,dn_rt_s=d.RT.values*60).to_csv('pg/added_not_diann_D_s0.tsv',sep='\t',index=False)

# Null groups: our target winners that NEW rejects (run q 0.05-0.5), and decoy winners mapped to their target's key at the
# decoy's own apex (a decoy has no real signal, so any DIA-NN 'confirmation' at that RT is chance).
w['dq']=w.q_value
rej=w[(w.is_decoy==0)&(w.q_value>0.05)&(w.q_value<=0.5)]
rows.clear()
ev(rej.sample(min(len(rej),2000),random_state=0),'null: our rejected targets (q .05-.5)')
dd=w[w.is_decoy==1].nlargest(694,'score').copy()
pk=pre[pre.decoy==0].set_index('seq')
tgt_of={r.pid:r.key for r in pre.itertuples()}
# decoy -> paired target key via the library pair (DECOY_ prefix on TRAML id)
tr=pd.read_sql("select ID pid, TRAML_ID t from PRECURSOR",lib).set_index('t').pid
dpid=pd.read_sql("select ID pid, TRAML_ID t from PRECURSOR where DECOY=1",lib)
dpid['tp']=dpid.t.str.replace('^DECOY_','',regex=True).map(tr)
m=dict(zip(dpid.pid,dpid.tp))
dd['precursor_id_t']=dd.precursor_id.map(m); dd=dd.dropna(subset=['precursor_id_t'])
dd['key']=dd.precursor_id_t.map(P.key); dd['precursor_id']=dd.precursor_id_t.astype(int)
ev(dd,'null: top decoys at their own apex')
print(pd.DataFrame(rows).round(3).to_string(index=False))
