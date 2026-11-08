"""Raw 01 (T/D only): exact native candidate feature sets (abs dRT unclipped, no squared term)."""
import pandas as pd, numpy as np, warnings
warnings.filterwarnings('ignore')
from lda import stage2, evaluate, BASE
df=pd.read_parquet('raw01_rows.parquet')
df['support']=df.transition_support.astype(float); df['log_scans']=np.log1p(df.matched_scans); df['log_width']=np.log1p(df.peak_width_seconds)
df['abs_drt']=df.apex_rt_error_seconds.abs()
print('abs_drt max', df.abs_drt.max())
S=['support','log_scans','log_width','abs_drt']
for name,feats in {'search set: cnn_z + 4':['cnn_z']+S,'extraction+search: 8 + 4':BASE+S}.items():
    s,W=stage2(df,feats); print(f'{name:28s}', evaluate(df,s), np.round(W,2).tolist(), flush=True)
