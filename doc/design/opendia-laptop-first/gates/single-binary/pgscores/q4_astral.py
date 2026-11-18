"""q4 precheck (diagnostic, entrapment labels for evaluation only). Per feature, on precursor winners of NEW seed 0:
AUC(entrapment > decoy) in the band around the cut. Decoys model false targets on a feature only if this is ~0.5.
Also AUC(human target > entrapment) = discrimination of true from false targets."""
import sqlite3, pandas as pd, numpy as np, sys
from sklearn.metrics import roc_auc_score
L='/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114/armE_raw01/library.pqp'
c=sqlite3.connect(f'file:{L}?mode=ro',uri=True)
ent=pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""",c).set_index('precursor_id').e
for R in ('01','02'):
    w=pd.read_csv(f'win{R}.tsv.gz',sep='\t')
    w=w[w.score>-1e29].copy(); w['e']=w.precursor_id.map(ent).fillna(0).astype(int)
    w=w.sort_values('score',ascending=False).reset_index(drop=True)
    d=np.cumsum(w.is_decoy.values); t=np.cumsum(1-w.is_decoy.values); q=np.minimum.accumulate((d/np.maximum(t,1))[::-1])[::-1]
    w['q']=q
    band=w[(w.q>=0.002)&(w.q<=0.05)]
    E=band[(band.is_decoy==0)&(band.e==1)]; D=band[band.is_decoy==1]; T=band[(band.is_decoy==0)&(band.e==0)]
    print(f'raw{R}: band q 0.002-0.05: human {len(T)}, entrapment {len(E)}, decoy {len(D)}')
    cols=[x for x in w.columns if x not in ('candidate_row','precursor_id','is_decoy','fold','rt_apex','score','q_value','e','q')]
    rows=[]
    for x in cols:
        ed=roc_auc_score(np.r_[np.ones(len(E)),np.zeros(len(D))],np.r_[E[x],D[x]])
        te=roc_auc_score(np.r_[np.ones(len(T)),np.zeros(len(E))],np.r_[T[x],E[x]])
        rows.append((x,ed,te))
    df=pd.DataFrame(rows,columns=['feature','AUC_E>D','AUC_T>E']); df['bias']=(df['AUC_E>D']-0.5).abs()
    print(df.sort_values('bias',ascending=False).round(3).to_string(index=False))
