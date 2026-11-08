"""Raw 01 (T/D only): starting from cnn_z + search shape features, add each extraction feature singly."""
import pandas as pd, numpy as np, warnings
warnings.filterwarnings('ignore')
from lda import stage2, evaluate, BASE
df=pd.read_parquet('raw01_rows.parquet')
df['support']=df.transition_support; df['log_scans']=np.log1p(df.matched_scans); df['log_width']=np.log1p(df.peak_width_seconds)
df['abs_drt']=df.apex_rt_error_seconds.abs().clip(upper=120); df['drt2']=(df.apex_rt_error_seconds/20.7).clip(-6,6)**2
core=['cnn_z','support','log_scans','log_width']
print(f'{"core":28s}', evaluate(df,stage2(df,core)[0]), flush=True)
for f in ['abs_drt','fragment_fraction','fragment_coelution','library_correlation','ms2_abs_ppm','ms1_abs_ppm','ms1_coelution','ms1_isotope']:
    print(f'{"core + "+f:28s}', evaluate(df,stage2(df,core+[f])[0]), flush=True)
