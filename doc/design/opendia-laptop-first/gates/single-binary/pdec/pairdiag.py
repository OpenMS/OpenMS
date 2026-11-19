"""Diagnostic (PLAN_PDEC D, design data): where does P(entrapment beats own predicted decoy | best q<=0.05) = 0.54 come from?
Split the entrapment pairs by charge, length, decoy fragment count < target's, internal K/R, RT shift decoy-target."""
import sqlite3, glob, numpy as np, pandas as pd
c=sqlite3.connect('file:hela_pdec5.pqp?mode=ro',uri=True)
pq=pd.read_sql("""select p.ID precursor_id, p.TRAML_ID t, p.DECOY d, p.CHARGE z, p.LIBRARY_RT rt, pe.UNMODIFIED_SEQUENCE s,
  min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""",c)
nfc=pd.read_sql("select PRECURSOR_ID precursor_id, count(*) nf from TRANSITION_PRECURSOR_MAPPING group by 1",c)
pq=pq.merge(nfc,on='precursor_id',how='left')
pq['key']=np.where(pq.d==1,pq.t.str.replace('^DECOY_','',regex=True),pq.t)
T=pq[pq.d==0].set_index('key'); D=pq[pq.d==1].set_index('key')
P=T.join(D[['precursor_id','nf','rt']],rsuffix='_d',how='inner'); P=P[P.e==1]
P['len']=P.s.str.len(); P['ikr']=P.s.str[:-1].str.contains('[KR]'); P['short']=P.nf_d<P.nf; P['drt']=(P.rt_d-P.rt).abs()
rows=[]
for f in sorted(glob.glob('runs/D_*_s*/*.gz')):
    t=pd.read_csv(f,sep='\t'); t=t[t.score>-1e29].sort_values('score',ascending=False).drop_duplicates('precursor_id').set_index('precursor_id')
    x=P.assign(st=t.score.reindex(P.precursor_id).values, sd=t.score.reindex(P.precursor_id_d).values,
               qt=t.q_value.reindex(P.precursor_id).values, qd=t.q_value.reindex(P.precursor_id_d).values).dropna(subset=['st','sd'])
    x=x[np.minimum(x.qt,x.qd)<=0.05]; x['win']=x.st>x.sd; x['arm']=f.split('/')[1].split('_')[1]; rows.append(x)
x=pd.concat(rows); print('pairs',len(x),'P(win)',round(x.win.mean(),3))
def split(col,lab=None):
    g=x.groupby(col).win.agg(['mean','size']); g['se']=np.sqrt(0.25/g['size']); print(f'-- {lab or col}'); print(g.round(3).to_string())
split('z'); split(pd.cut(x.len,[0,9,12,16,44]),'length'); split('short','decoy has fewer fragments'); split('ikr','internal K/R')
split(pd.qcut(x.drt,4),'|RT decoy - RT target| quartile'); split(x.nf<6,'target <6 fragments')
x['rel']=np.sign(x.nf_d-x.nf).map({-1:'decoy fewer',0:'equal',1:'decoy more'})
split('rel','decoy vs target fragment count')
