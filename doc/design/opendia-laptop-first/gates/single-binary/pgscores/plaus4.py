"""Decoy null for the PLAUS_D indicators (+-3 s): decoy winners scored like the extra IDs (same score band as NEW's
non-DIA-NN additions), mapped to their paired target's key at the decoy's own apex."""
import sqlite3, glob, pandas as pd, numpy as np
exec(open('plaus3.py').read().split("rows={")[0])
tr=pd.read_sql("select ID pid, TRAML_ID t from PRECURSOR",lib).set_index('t').pid
dp=pd.read_sql("select ID pid, TRAML_ID t from PRECURSOR where DECOY=1",lib); dp['tp']=dp.t.str.replace('^DECOY_','',regex=True).map(tr)
m=dict(zip(dp.pid,dp.tp))
nd=add[~add.key.isin(D1.index)]; lo,hi=nd.score.min(),nd.score.max()
dec=w[(w.is_decoy==1)].copy()
for name,g in (('decoys in the score band of the 694',dec[(dec.score>=lo)&(dec.score<=hi)]),('top 694 decoys',dec.nlargest(694,'score')),('all decoy winners (random 3000)',dec.sample(3000,random_state=0))):
    g=g.copy(); g['precursor_id']=g.precursor_id.map(m); g=g.dropna(subset=['precursor_id']); g['precursor_id']=g.precursor_id.astype(int); g['key']=g.precursor_id.map(P.key)
    tgt=w[(w.is_decoy==0)].set_index('precursor_id').rt_apex
    dt=np.abs(g.rt_apex.values-tgt.reindex(g.precursor_id).values)
    print(f'{name:38s}', {k:round(v,3) for k,v in sup(g).items()}, f'| decoy apex within 3 s of its target winner apex: {np.nanmean(dt<=3):.2f}')
print('NEW added, not in DIA-NN              ', {k:round(v,3) for k,v in sup(nd).items()})
