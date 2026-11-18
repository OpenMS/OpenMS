"""q4: AUC(entrapment > decoy) of the final score among precursor winners in the band q 0.002-0.05, CUR vs NEW, seeds 0-4."""
import sqlite3, pandas as pd, numpy as np
from sklearn.metrics import roc_auc_score
L='/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114/armE_raw01/library.pqp'
c=sqlite3.connect(f'file:{L}?mode=ro',uri=True)
ent=pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""",c).set_index('precursor_id').e
def auc(path):
    t=pd.read_csv(path,sep='\t'); t=t[t.score>-1e29]
    w=t.sort_values('score',ascending=False).drop_duplicates('precursor_id').reset_index(drop=True)
    d=np.cumsum(w.is_decoy.values); tt=np.cumsum(1-w.is_decoy.values); w['q']=np.minimum.accumulate((d/np.maximum(tt,1))[::-1])[::-1]
    b=w[(w.q>=0.002)&(w.q<=0.05)]; e=b.precursor_id.map(ent).fillna(0).values
    E=b[(b.is_decoy==0)&(e==1)].score; D=b[b.is_decoy==1].score
    return roc_auc_score(np.r_[np.ones(len(E)),np.zeros(len(D))],np.r_[E,D]), len(E), len(D)
for R in ('01','02','03'):
    for arm,root in (('cur','ms1' if R!='03' else 'pg'),('new','pg')):
        r=[auc(f'{root}/raw{R}_s{s}_{arm}/p.tsv.gz') for s in range(5)]
        print(f'raw{R} {arm}: AUC E>D', ' '.join(f'{a:.3f}' for a,_,_ in r), f'| median {np.median([a for a,_,_ in r]):.3f} | E/D n ~{r[0][1]}/{r[0][2]}')
