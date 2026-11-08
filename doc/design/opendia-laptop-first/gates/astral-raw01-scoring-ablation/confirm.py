"""Diagnostic only: share of accepted raw 01 targets that DIA-NN (full library) also reports at 1% with RT within 10 s.
Not used to fit or select; it flags features that exploit target/decoy asymmetry."""
import pandas as pd, numpy as np
from lda import stage2, winners, qvals, BASE
G='/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-diann-gap-20261108/'
ref=pd.read_parquet(G+'ref_exact.parquet'); ref=ref[ref.Run.str.endswith('_01')][['precursor_id','RT']]
df=pd.read_parquet('raw01_rows.parquet')
df['support']=df.transition_support; df['log_scans']=np.log1p(df.matched_scans); df['log_width']=np.log1p(df.peak_width_seconds)
df['abs_drt']=df.apex_rt_error_seconds.abs().clip(upper=120); df['drt2']=(df.apex_rt_error_seconds/20.7).clip(-6,6)**2
def accepted(score):
    d=df.assign(s2=score); win=d.loc[winners(d,score)]
    def st(g):
        x=g[g.is_decoy==1].s2; m=x.median(); s=1.4826*(x-m).abs().median(); return (g.s2-m)/s
    win['z']=win.groupby('fold',group_keys=False).apply(st); q=qvals(win,'z'); return q[(q.is_decoy==0)]
for name,feats in {'base':BASE,'+support/scans/width':BASE+['support','log_scans','log_width'],'+delta_rt':BASE+['abs_drt','drt2'],'+both':BASE+['support','log_scans','log_width','abs_drt','drt2']}.items():
    s,_=stage2(df,feats); t=accepted(s)
    for qq in (0.01,0.05):
        a=t[t.q<=qq].merge(ref,on='precursor_id',how='left')
        conf=((a.rt_apex-a.RT*60).abs()<=10)
        print(f'{name:22s} q<={qq}: accepted {len(a)}, DIA-NN-confirmed at same RT {int(conf.sum())} ({conf.mean():.3f}), in DIA-NN at all {int(a.RT.notna().sum())}', flush=True)
