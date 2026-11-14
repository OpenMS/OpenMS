"""Label-free tell check for featI (and featA/G): in null windows (|library RT error| > 3 sigma, no peptide expected),
AUC(target > decoy) per column over windows and over per-precursor null winners. A column that separates targets from
decoys where no peptide is present encodes a library/sequence tell. raw01 arm E arrays."""
import json, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
D = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114/armE_raw01'
pop = pd.read_csv(f'{D}/population.tsv', sep='\t'); H = np.load(f'{D}/heuristics.npy')
sig = json.load(open(f'{D}/cal.json'))['rt']['sigma_s']; null = H[:, 3] > 3 * sig
y = pop.is_decoy.values == 0
s = np.load(f'{D}/scores_D11FDWT.npy')
d = pop[null].assign(s=s[null], r=np.where(null)[0]).sort_values('s', ascending=False).drop_duplicates('precursor_id')
print(f'null windows {null.sum()}, null winners {len(d)}')
for feat in ('featA', 'featG', 'featI'):
    F = np.load(f'{D}/{feat}.npy')
    for j in range(F.shape[1]):
        a1 = roc_auc_score(y[null], F[null, j]); a2 = roc_auc_score(d.is_decoy.values == 0, F[d.r.values, j])
        print(f'  {feat} col {j}: null-window AUC(T>D) {a1:.4f}  null-winner {a2:.4f}')
