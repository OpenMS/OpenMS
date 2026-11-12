"""Descriptive: our search-run geometry vs DIA-NN's peak, raw01. Labels only define groups."""
import numpy as np, pandas as pd
B='/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; V=B+'/astral-voxel-pilot-20261110/v2'
P=pd.read_csv(V+'/raw01/population.tsv',sep='\t'); H=np.load(V+'/raw01/heuristics.npy')
W=pd.read_csv(V+'/v4/stageD_v16/score/reference_winners_raw01.tsv.gz',sep='\t')
W['width']=np.expm1(H[W.row.values,2]); W['support']=H[W.row.values,0]; W['scans']=np.expm1(H[W.row.values,1])
dn=pd.read_csv(B+'/diann-ourlib-20261112/diann_ourlib_raw01.tsv.gz',sep='\t'); dn=dn[dn['Q.Value']<=0.01].set_index('precursor_id')
W=W.join(dn[['RT','RT.Start','RT.Stop','FWHM']],on='precursor_id')
W['at_rt']=(W.rt_apex-W.RT*60).abs()<=10
g=np.full(len(W),'',object)
g[(W.RT.notna()&W.at_rt&(W.full_pct<0.999)).values]='DN'; g[(W.RT.notna()&W.at_rt&(W.full_pct>=0.999)).values]='BOTH'
g[(W.group=='DEC').values]='DEC'; g[(W.group=='ENT').values]='ENT'; W['g']=g; W=W[W.g!=''].copy()
W['diann_fwhm_s']=W.FWHM*60; W['diann_span_s']=(W['RT.Stop']-W['RT.Start'])*60
print('medians'); print(W.groupby('g')[['width','scans','support','diann_fwhm_s','diann_span_s']].median().round(2).to_string())
print('q10'); print(W.groupby('g')[['width','scans']].quantile(0.1).round(1).to_string())
print('q90'); print(W.groupby('g')[['width','scans']].quantile(0.9).round(1).to_string())
x=W[W.g.isin(['DN','BOTH'])]
print('Spearman our width vs DIA-NN span (DN+BOTH):', round(x[['width','diann_span_s']].corr('spearman').iloc[0,1],3))
for gg in ['DN','BOTH']:
    y=W[W.g==gg]; print(gg,'our width / DIA-NN span median', round((y.width/y.diann_span_s).median(),2), '| frac our width < 0.5x span', round((y.width<0.5*y.diann_span_s).mean(),3))
print('scans per second of width, median:', W.assign(r=W.scans/W.width.clip(lower=0.1)).groupby('g').r.median().round(3).to_dict())
W.to_csv('width_raw01.tsv.gz',sep='\t',index=False)
