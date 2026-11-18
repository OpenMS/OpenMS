"""Among entrapment-vs-own-decoy pairs whose better member has run q <= 0.05 (D s0 CUR, PGENT2 lib, -features_out),
which feature difference (entrapment minus its decoy, each at its own winning window) predicts which one wins?
A feature with mean diff > 0 and high AUC for 'entrapment wins' is a source of the asymmetry."""
import sqlite3, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
c=sqlite3.connect('file:hela_ath2.pqp?mode=ro',uri=True)
P=pd.read_sql("""select p.ID precursor_id, p.TRAML_ID t, p.DECOY d, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
 from PRECURSOR p left join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID left join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID
 left join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""",c)
P['key']=np.where(P.d==1,P.t.str.replace('^DECOY_','',regex=True),P.t)
tg=P[(P.d==0)&(P.e==1)].set_index('key')[['precursor_id']]; dc=P[P.d==1].set_index('key')[['precursor_id']].rename(columns={'precursor_id':'dec'})
pairs=tg.join(dc,how='inner')
w=pd.read_csv('win_D_cur_s0.tsv.gz',sep='\t'); w=w[w.score>-1e29].set_index('precursor_id')
meta=('candidate_row','is_decoy','fold','rt_apex','score','q_value')
cols=[x for x in w.columns if x not in meta]
pairs=pairs[pairs.precursor_id.isin(w.index)&pairs.dec.isin(w.index)]
E=w.loc[pairs.precursor_id]; D=w.loc[pairs.dec]
best_q=np.minimum(E.q_value.values,D.q_value.values); m=best_q<=0.05
win=(E.score.values>D.score.values)[m]
print(f'pairs with best q<=0.05: {m.sum()}, entrapment wins {win.mean():.3f}')
rows=[]
for x in cols:
    diff=(E[x].values-D[x].values)[m]
    a=roc_auc_score(win,diff)
    rows.append((x,np.nanmean(diff),a))
df=pd.DataFrame(rows,columns=['feature','mean diff E-D','AUC diff -> E wins']); df['dev']=(df.iloc[:,2]-0.5).abs()
print(df.sort_values('dev',ascending=False).head(15).round(3).to_string(index=False))
# all pairs (no q restriction) for contrast
diffs={x:(E[x].values-D[x].values) for x in ['f12_libcos','abs_library_rt_error','boundary_library_correlation','i_libcos_resid']}
print('all pairs mean diff:',{k:round(float(np.nanmean(v)),4) for k,v in diffs.items()})
