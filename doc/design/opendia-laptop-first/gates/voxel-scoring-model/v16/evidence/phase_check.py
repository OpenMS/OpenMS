"""Post-primary clarification, not tuned: isolate bin-pair origin from endpoint cropping.
The primary origin1 diagnostic discarded outer bins. This check retains every input
point and gives both origins 17 output bins, with zero-padding [0,2] vs [1,1].
No parameters selected from outcomes. Physical 0.9375-second original bins fixed.
"""
import numpy as np,pandas as pd
from diagnose import OUT,V2,paircorr,stratified

def main():
 summaries=[];aucs=[]
 for raw in ['raw01','raw02']:
  t=pd.read_csv(OUT/f'window_measurements_{raw}.tsv',sep='\t',usecols=['row','s','fold','D','E','HUMAN','ENT','DEC'])
  x=np.load(V2/raw/'xic.npy',mmap_mode='r');rr=t.row.to_numpy();parts=[]
  for k in range(0,len(rr),10000):
   z=np.asarray(x[rr[k:k+10000]],float)
   a=np.pad(z,((0,0),(0,0),(0,2))).reshape(len(z),6,17,2).sum(3)
   b=np.pad(z,((0,0),(0,0),(1,1))).reshape(len(z),6,17,2).sum(3)
   assert np.allclose(a.sum((1,2)),z.sum((1,2))) and np.allclose(b.sum((1,2)),z.sum((1,2)))
   parts.append(abs(paircorr(a)-paircorr(b)))
  t['mass_preserving_phase_sensitivity']=np.concatenate(parts)
  for grp in ['D','E','HUMAN','ENT','DEC']:
   g=t[t[grp]].mass_preserving_phase_sensitivity
   summaries.append(dict(raw=raw,group=grp,n=len(g),finite=int(g.notna().sum()),median=g.median(),q90=g.quantile(.9),share_gt_0_1=(g>.1).mean()))
  a,_=stratified(t,['mass_preserving_phase_sensitivity'],raw);aucs.append(a)
 pd.DataFrame(summaries).to_csv(OUT/'mass_preserving_phase_summary.tsv',sep='\t',index=False)
 pd.concat(aucs).to_csv(OUT/'mass_preserving_phase_auc.tsv',sep='\t',index=False)
 print('PASS: all point intensities preserved for both origins, both raws')
 print(pd.DataFrame(summaries).round(4).to_string(index=False))
 print(pd.concat(aucs)[['raw','left','right','auc_high','matched_left']].round(4).to_string(index=False))
if __name__=='__main__':main()
