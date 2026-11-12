from pathlib import Path
import pandas as pd
O=Path(__file__).resolve().parent
out=[]
for raw in ['raw01','raw02']:
    W=pd.read_csv(O/raw/'coverage_winners.tsv.gz',sep='\t')
    for group,mask in [('D',W.stage.eq('D right window, below cut')),('E',W.stage.eq('E accepted'))]:
        for category,G in W[mask].groupby('exact_all_positive_category'):
            out.append(dict(raw=raw,group=group,category=category,n=len(G),full_linear_far=int((G.final_full_linear_rt_distance>10).sum()),deep=int((G.full_pct<.9).sum()),deep_and_far=int(((G.full_pct<.9)&(G.final_full_linear_rt_distance>10)).sum()),support_outside_median=G.features_outside_training_positive_1_99.median()))
pd.DataFrame(out).to_csv(O/'coverage_selector_overlap.tsv',sep='\t',index=False)
print(pd.DataFrame(out).to_string(index=False))
