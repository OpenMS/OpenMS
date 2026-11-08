"""Raw 02: features of DIA-NN-confirmed precursors whose stage-2 winner is at the right RT but not accepted, vs decoy winners."""
import pandas as pd, numpy as np
B='/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
t=pd.read_parquet('lib_targets.parquet')
r=pd.read_parquet('diann_runs_1pct.parquet'); r=r[r.Run.str.endswith('_02')]
print('unique keys DIA-NN raw02 in library', r.key.isin(set(t.key)).sum(), '| library rows per key max', t.key.value_counts().max(), 'keys with >1 precursor', int((t.key.value_counts()>1).sum()))
ref=t.merge(r[['key','RT','entrap']],on='key'); ref['rt_s']=ref.RT*60
f=pd.read_csv(B+'/astral-raw02-stage2-check-20261108/crossfit/features.tsv',sep='\t')
c=pd.read_csv(B+'/astral-raw02-stage2-check-20261108/crossfit/crossfit.tsv',sep='\t',usecols=['sample_id','precursor_best','q_value'])
w=pd.read_csv(B+'/astral-entrapment-raw02-20261108/all_windows.tsv',sep='\t',usecols=['sample_id','rt_apex'])
f=f.merge(c,on='sample_id').merge(w,on='sample_id'); f['q']=pd.to_numeric(f.q_value,errors='coerce')
win=f[f.precursor_best==1]
m=win.merge(ref[['precursor_id','rt_s']],on='precursor_id'); m['drt']=m.rt_apex-m.rt_s
right=m[m.drt.abs()<=10]
cols=['cnn_z','fragment_fraction','fragment_coelution','library_correlation','ms2_abs_ppm','ms1_abs_ppm','ms1_coelution','ms1_isotope','stage2_score']
groups={'ref winner right RT, accepted':right[right.q<=0.01],'ref winner right RT, q>0.01':right[right.q>0.01],'ref winner wrong RT':m[m.drt.abs()>10],'decoy winners':win[win.is_decoy==1],'decoy winners top 1%':win[win.is_decoy==1].nlargest(2700,'stage2_score')}
for k,g in groups.items(): print(f'{k:32s} n={len(g):6d}', g[cols].median().round(2).to_dict())
print('right-RT unaccepted q quantiles', right[right.q>0.01].q.quantile([.25,.5,.75]).round(3).to_dict())
print('wrong-RT |drt| quantiles', m[m.drt.abs()>10].drt.abs().quantile([.25,.5,.75]).round(0).to_dict())
