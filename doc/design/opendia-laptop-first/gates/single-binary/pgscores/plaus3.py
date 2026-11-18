"""Compare support of NEW's extra non-DIA-NN IDs with CUR's own non-DIA-NN IDs (same proxy, same run, seed 0)."""
import sqlite3, glob, pandas as pd, numpy as np
exec(open('plaus.py').read().split("groups={")[0])
ws=pd.read_csv(glob.glob('rtc/D_imgrid_s0/*.precursors.tsv')[0],sep='\t').sort_values('score',ascending=False).drop_duplicates('precursor_id')
ws['key']=ws.precursor_id.map(P.key); cacc=ws[(ws.q_value<=0.01)&(ws.is_decoy==0)]
N=set(acc.precursor_id)
def sup(g):
    d=DD.reindex(g.key); o=osw.reindex(g.precursor_id)
    a=(np.abs(d.RT.values*60-g.rt_apex.values)<=3)&(d['Q.Value'].values<=0.05)
    b=(np.abs(o.rt.values-g.rt_apex.values)<=3)&(o.q.values<=0.01)
    h=DH.reindex(g.key)['Q.Value'].values<=0.01
    return dict(n=len(g), dn05_same_rt=a.mean(), osw1_same_rt=b.mean(), either=(a|b).mean(), diann_other_run=h.mean())
rows={'CUR accepted, not in DIA-NN':sup(cacc[~cacc.key.isin(D1.index)]),
      ' of which NEW keeps':sup(cacc[~cacc.key.isin(D1.index)&cacc.precursor_id.isin(N)]),
      ' of which NEW drops':sup(cacc[~cacc.key.isin(D1.index)&~cacc.precursor_id.isin(N)]),
      'NEW added, not in DIA-NN':sup(add[~add.key.isin(D1.index)]),
      'NEW added, in DIA-NN':sup(add[add.key.isin(D1.index)]),
      'CUR accepted, in DIA-NN':sup(cacc[cacc.key.isin(D1.index)])}
print(pd.DataFrame(rows).T.round(3).to_string())
