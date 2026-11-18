"""Expected true count of NEW's additions not in DIA-NN 1 %, from DIA-NN's own posterior (1 - PEP) at the same peak (+-3 s),
compared with CUR's own non-DIA-NN IDs and with TDC (decoys in the same score band). D seed 0."""
import sqlite3, glob, pandas as pd, numpy as np
exec(open('plaus.py').read().split("groups={")[0])
d=pd.read_parquet('diann_hela/diann_report.parquet',columns=['Modified.Sequence','Precursor.Charge','Q.Value','PEP','RT'])
d['key']=d['Modified.Sequence']+'/'+d['Precursor.Charge'].astype(str); d=d.sort_values('Q.Value').drop_duplicates('key').set_index('key')
ws=pd.read_csv(glob.glob('rtc/D_imgrid_s0/*.precursors.tsv')[0],sep='\t').sort_values('score',ascending=False).drop_duplicates('precursor_id')
ws['key']=ws.precursor_id.map(P.key); cacc=ws[(ws.q_value<=0.01)&(ws.is_decoy==0)]
def post(g,name):
    x=d.reindex(g.key); same=(np.abs(x.RT.values*60-g.rt_apex.values)<=3)
    p=np.where(same,1-x.PEP.values,np.nan)
    seen=np.isfinite(p)
    print(f'{name:30s} n {len(g):5d} | DIA-NN same peak {seen.mean():.2f} | sum(1-PEP) {np.nansum(p):7.1f} '
          f'| expected true if unseen are all false: {np.nansum(p)/len(g):.2f} | mean (1-PEP) where seen {np.nanmean(p):.2f}')
post(add[~add.key.isin(D1.index)],'NEW added, not in DIA-NN 1 %')
post(cacc[~cacc.key.isin(D1.index)],'CUR accepted, not in DIA-NN 1 %')
post(add[add.key.isin(D1.index)],'NEW added, in DIA-NN 1 %')
# TDC view: in the score range of NEW's non-DIA-NN additions, targets vs decoys
X=add[~add.key.isin(D1.index)]
for lo,hi in ((X.score.quantile(q),X.score.quantile(q2)) for q,q2 in ((0,.25),(.25,.5),(.5,.75),(.75,1))):
    b=w[(w.score>=lo)&(w.score<=hi)]
    nx=((X.score>=lo)&(X.score<=hi)).sum()
    print(f'score {lo:7.2f}-{hi:7.2f}: targets {int((b.is_decoy==0).sum()):5d}, decoys {int(b.is_decoy.sum()):4d}, of targets in DIA-NN 1 % {b[b.is_decoy==0].key.isin(D1.index).mean():.2f}, X {nx}')
