"""Execute the frozen label-free rules. No evaluation imports or biological annotations."""
from pathlib import Path
import numpy as np
import pandas as pd
from paircal_helpers import pair_winners, tdc_q, self_test

ROOT = Path('/home/singjc/.jcode/scratch/paircal-20261126')
INPUT = Path('/home/singjc/.jcode/scratch/blend-20261125/eval/ref/entlib-gbt-20261122/runs')
self_test()
meta = pd.read_pickle(ROOT/'library_label_free.pkl')[['precursor_id','key','decoy']]
first = True
summaries = []
for dataset in ['D','H']:
    for seed in range(5):
        nf = list((INPUT/f'{dataset}_cur_s{seed}').glob('*.precursors.tsv*'))
        gf = list((INPUT/f'{dataset}_gbt_s{seed}').glob('*.precursors.tsv*'))
        assert len(nf)==len(gf)==1, (nf,gf)
        n = pair_winners(nf[0],meta).set_index('key')
        g = pair_winners(gf[0],meta).reset_index(drop=True)
        # Quantile pool is all NET-kept decoy precursor winners, not only pair winners.
        nd = pd.read_csv(nf[0],sep='\t',usecols=['precursor_id','score','is_decoy'])
        nd = nd[np.isfinite(nd.score)&(nd.score>-1e29)].sort_values('score',ascending=False).drop_duplicates('precursor_id')
        cut = np.quantile(nd.loc[nd.is_decoy==1,'score'],.99)
        net_id = g.key.map(n.precursor_id)
        net_max = g.key.map(n.score)
        ag = np.where(net_id.isna(),'missing',np.where(net_id.to_numpy()==g.precursor_id.to_numpy(),'agree','disagree'))
        hi = np.where(net_max.isna(),'missing',np.where(net_max>cut,'high','low'))
        for rule,groups in [('baseline',np.repeat('all',len(g))),('A_agreement',ag),('B_net99',hi)]:
            q = np.ones(len(g))
            for group in sorted(set(groups)):
                ix = np.flatnonzero(groups==group)
                q[ix] = tdc_q(g.score.to_numpy()[ix],g.decoy.to_numpy()[ix])
                summaries.append(dict(dataset=dataset,seed=seed,rule=rule,group=group,pairs=len(ix),targets=int((g.decoy.to_numpy()[ix]==0).sum()),decoys=int(g.decoy.to_numpy()[ix].sum()),net_decoy_q99=cut))
            output = g[['precursor_id','score','decoy']].copy()
            output['q_value'] = q
            output['dataset'],output['seed'],output['rule'],output['group'] = dataset,seed,rule,groups
            output.to_csv(ROOT/'selected.tsv',sep='\t',index=False,mode='w' if first else 'a',header=first)
            first = False
            print(dataset,seed,rule,'pairs',len(g),'targets_q01',int(((g.decoy==0)&(q<=.01)).sum()),'strata',dict(zip(*np.unique(groups,return_counts=True))),flush=True)
        print(f'JCODE_PROGRESS {{"dataset":"{dataset}","seed":{seed}}}',flush=True)
pd.DataFrame(summaries).to_csv(ROOT/'strata.tsv',sep='\t',index=False)
