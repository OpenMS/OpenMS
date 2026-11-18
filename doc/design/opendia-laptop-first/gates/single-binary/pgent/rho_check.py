"""Diagnostic for q9 (labels for evaluation only). Per run/arm/seed 0: among precursor winners, by run-q band,
rho = estimated false targets from entrapments (E * (1 + 1/r)) / decoy winners in the same band. rho ~ 1 means the
decoys model false targets; rho > 1 means the decoys under-count them. Also AUC(entrapment > decoy) on the score near the cut."""
import sqlite3, glob, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
c=sqlite3.connect('file:hela_entrap_dec.pqp?mode=ro',uri=True)
ent=pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
 from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by 1""",c).set_index('precursor_id').e
r=1.1225; k=1+1/r
for run in ('D','H'):
  for arm in ('cur','new'):
    f=glob.glob(f'runs/runs/{run}_{arm}_s0/*.gz')[0]
    t=pd.read_csv(f,sep='\t'); t=t[t.score>-1e29].sort_values('score',ascending=False).drop_duplicates('precursor_id')
    t['e']=t.precursor_id.map(ent).fillna(0).values*(t.is_decoy==0)
    s=t.score.values; dec=t.is_decoy.values; e=t.e.values
    # band edges by score quantile of the cut region: use the tool's q of targets
    q=t.q_value.values
    out=[]
    for lo,hi in ((0,0.005),(0.005,0.01),(0.01,0.02),(0.02,0.05),(0.05,0.2)):
        cut_hi=s[(dec==0)&(q<=lo)].min() if lo>0 else np.inf
        cut_lo=s[(dec==0)&(q<=hi)].min()
        b=(s<cut_hi)&(s>=cut_lo)
        E=e[b].sum(); D=dec[b].sum()
        out.append(f'q {lo}-{hi}: E {int(E)} D {int(D)} rho {E*k/max(D,1):.2f}')
    band=(q<=0.05)|(dec==1)&(s>=s[(dec==0)&(q<=0.05)].min())
    E_=s[band&(e==1)]; D_=s[band&(dec==1)]
    auc=roc_auc_score(np.r_[np.ones(len(E_)),np.zeros(len(D_))],np.r_[E_,D_]) if len(E_) and len(D_) else float('nan')
    print(f'{run} {arm}: '+' | '.join(out)+f' | AUC(E>D, q<=0.05 band) {auc:.3f}')
