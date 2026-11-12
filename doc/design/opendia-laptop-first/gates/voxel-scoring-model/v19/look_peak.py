"""Raw01 look at peak.npy (descriptive; labels only define groups). Raw02 not read."""
import numpy as np, pandas as pd
V='/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2'
W=pd.read_csv('width_raw01.tsv.gz',sep='\t'); Pk=np.load(V+'/raw01/peak.npy')
W['nsup']=Pk[W.row.values,0]; W['npts']=np.expm1(Pk[W.row.values,1]); W['nfwhm']=np.expm1(Pk[W.row.values,2])
x=W[W.g.isin(['DN','BOTH'])]
print('Spearman native FWHM vs DIA-NN FWHM (DN+BOTH):', round(x[['nfwhm','diann_fwhm_s']].corr('spearman').iloc[0,1],3), '| old width vs DIA-NN FWHM', round(x[['width','diann_fwhm_s']].corr('spearman').iloc[0,1],3))
print('medians'); print(W.groupby('g')[['nsup','npts','nfwhm','support','scans','width','diann_fwhm_s']].median().round(2).to_string())
def auc(a,b):
    r=pd.Series(np.r_[a,b]).rank().values; return (r[:len(a)].sum()-len(a)*(len(a)+1)/2)/(len(a)*len(b))
print('AUC vs DEC (all winners, ties half):')
for f in ['nsup','npts','nfwhm','support','scans','width']:
    d=W[W.g=='DEC'][f].values
    print(f.ljust(8), {g: round(auc(W[W.g==g][f].values,d),3) for g in ['DN','BOTH','ENT']})
