"""Raw 02 (spent, descriptive): funnel against DIA-NN 1% precursors that exactly match a library precursor (seq, charge, m/z incl. mods)."""
import pandas as pd, numpy as np, json
B='/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
x=pd.read_parquet('ref_exact.parquet'); ref=x[x.Run.str.endswith('_02')].copy(); ref['rt_s']=ref.RT*60
print('reference', len(ref), 'entrapment', int(ref.entrap.sum()))
w=pd.read_csv(B+'/astral-entrapment-raw02-20261108/all_windows.tsv',sep='\t',usecols=['sample_id','rt_apex'])
c=pd.read_csv(B+'/astral-raw02-stage2-check-20261108/crossfit/crossfit.tsv',sep='\t').merge(w,on='sample_id')
c['q']=pd.to_numeric(c.q_value,errors='coerce')
m=c[c.is_decoy==0].merge(ref[['precursor_id','rt_s']],on='precursor_id'); m['right']=(m.rt_apex-m.rt_s).abs()<=10
g=m.groupby('precursor_id')
print('with a candidate window', g.ngroups, '| with a candidate within 10 s', int(g.right.any().sum()))
win=m[m.precursor_best==1]
print('winner at the right RT', int(win.right.sum()), '| accepted q<=0.01', int((win.q<=0.01).sum()), 'right RT and accepted', int((win.right&(win.q<=0.01)).sum()))
for q in (0.01,0.05,0.1): print(' accepted at q<=',q, int((win.q<=q).sum()))
# RT: intrinsic scatter of library RT vs DIA-NN RT (robust, local linear in 20 bins)
ref=ref.sort_values('lib_rt'); bins=pd.qcut(ref.lib_rt,20,labels=False)
res=ref.rt_s-ref.groupby(bins).rt_s.transform('median')
print('DIA-NN RT vs library RT: residual 1.4826*MAD (s)', round(1.4826*(res-res.median()).abs().median(),1), 'q90 |res|', round(res.abs().quantile(.9),1))
cal=json.load(open(B+'/astral-entrapment-raw02-20261108/cal/cal.json'))
print('our calibration rt block', {k:v for k,v in cal['rt'].items() if not isinstance(v,(list,dict))})
