"""Raw 01 (T/D only): how much does the CNN add on top of hand features? Feature 0 seeds the first ranking."""
import pandas as pd, numpy as np
from lda import stage2, evaluate, BASE
df=pd.read_parquet('raw01_rows.parquet')
df['support']=df.transition_support; df['log_scans']=np.log1p(df.matched_scans); df['log_width']=np.log1p(df.peak_width_seconds)
df['abs_drt']=df.apex_rt_error_seconds.abs().clip(upper=120); df['drt2']=(df.apex_rt_error_seconds/20.7).clip(-6,6)**2
S=['support','log_scans','log_width','abs_drt','drt2']
H=[f for f in BASE if f!='cnn_z']
for name,feats in {'cnn_z only':['cnn_z'],'cnn_z + search':['cnn_z']+S,'hand only (coelution seeds)':['fragment_coelution']+[f for f in H if f!='fragment_coelution']+S,'all':BASE+S}.items():
    s,W=stage2(df,feats); print(f'{name:30s}', evaluate(df,s), flush=True)
