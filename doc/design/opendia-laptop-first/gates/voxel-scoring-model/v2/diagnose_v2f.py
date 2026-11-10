"""Part 6: does the round-0 training set allow an intensity shortcut? AUC of single window descriptors, anchors vs
round-0 negatives (paired + random decoys), and V0 top-decile share of bright windows among all decoy windows."""
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
O = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2'
for raw in ('raw01', 'raw02'):
    pop = pd.read_csv(f'{O}/{raw}/population.tsv', sep='\t'); tw = pd.read_csv(f'{O}/{raw}/train_windows.tsv', sep='\t')
    tw = tw[tw.group.isin(['anchor', 'anchor_decoy', 'decoy'])]
    row = pd.Series(np.arange(len(pop)), index=pop.sample_id).loc[tw.sample_id].values
    xic = np.load(f'{O}/{raw}/xic.npy', mmap_mode='r'); o = np.argsort(row); x = np.empty((len(row), 6, 32), np.float32); x[o] = xic[row[o]]
    a = x.sum(2); tot = np.log1p(a.sum(1)); nfr = (a > 0.05 * a.sum(1, keepdims=True)).sum(1)
    y = tw.label.values
    print(raw, f'anchors {y.sum()} neg {(y == 0).sum()} | AUC log MS2 sum {roc_auc_score(y, tot):.3f} | AUC n frags>5% {roc_auc_score(y, nfr):.3f}',
          f'| median MS2 sum anchors {np.expm1(np.median(tot[y == 1])):.3g} neg {np.expm1(np.median(tot[y == 0])):.3g}')
    s = np.load(f'{O}/{raw}/scores_V0.npy'); d = np.where(pop.is_decoy.values == 1)[0]; d = np.random.default_rng(1).choice(d, 200000, replace=False)
    od = np.argsort(d); xd = np.empty((len(d),), np.float32); xd[od] = np.log1p(np.asarray(xic[d[od]]).sum((1, 2)))
    hi = xd >= np.quantile(xd, 0.99)
    print(raw, f'V0 decoy windows: logit median bright top 1 % {np.median(s[d][hi]):.2f} vs rest {np.median(s[d][~hi]):.2f}; share of bright among decoy logit top 0.1 %',
          f'{hi[s[d] >= np.quantile(s[d], 0.999)].mean():.2f}')
