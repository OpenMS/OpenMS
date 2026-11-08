"""Raw 01 design ablation (T/D only): add library/search feature groups to the frozen stage-2 LDA one at a time."""
import pandas as pd, numpy as np
from lda import stage2, evaluate, BASE
df=pd.read_parquet('raw01_rows.parquet')
df['abs_drt']=df.apex_rt_error_seconds.abs().clip(upper=120)
df['drt2']=(df.apex_rt_error_seconds/20.7).clip(-6,6)**2
df['log_apex']=np.log1p(df.apex_intensity)
df['log_area']=np.log1p(df.area)
df['support']=df.transition_support
df['log_scans']=np.log1p(df.matched_scans)
df['log_width']=np.log1p(df.peak_width_seconds)
for z in (1,2,3,4): df[f'chg{z}']=(df.charge==z).astype(float)
groups={'base':[],'+delta_rt':['abs_drt','drt2'],'+intensity':['log_apex','log_area'],'+support/scans/width':['support','log_scans','log_width'],'+charge':['chg1','chg2','chg3','chg4']}
res={}
for name,extra in groups.items():
    s,W=stage2(df,BASE+extra); r=evaluate(df,s); res[name]=r
    print(f'{name:22s}', r, 'weights fold0', dict(zip(BASE+extra,np.round(W[0],2))) if extra else '', flush=True)
allx=sum(groups.values(),[])
s,W=stage2(df,BASE+allx); print(f'{"all groups":22s}', evaluate(df,s)); print(dict(zip(BASE+allx,np.round(W.mean(0),2))))
