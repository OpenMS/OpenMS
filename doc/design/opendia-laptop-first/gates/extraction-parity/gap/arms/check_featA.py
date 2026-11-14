"""Pre-plan look (raw01 only, evaluation labels): featA columns, rejected DIA-NN IDs vs score-matched decoys, entrapment
in the upper band vs its matched decoys (as evidence_band.py)."""
import numpy as np, os, sys
FEAT = sys.argv[1] if len(sys.argv) > 1 else 'featA'
src = open('/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/gap-20261114/evidence.py').read().split('XIC = np.load')[0]
exec(src.replace("E = w[(w.is_decoy == 0) & w.precursor_id.map(ent).eq(1) & (w.s < cut)]",
 "R0 = wr.loc[sorted(near_is_winner - accepted)]\nE = w[(w.is_decoy == 0) & w.precursor_id.map(ent).eq(1) & (w.s < cut) & (w.s >= R0.s.median())]"))
from sklearn.metrics import roc_auc_score
A = np.load(f'{D}/{FEAT}.npy')
def auc(a, b, j):
    va, vb = A[a, j], A[b, j]; return roc_auc_score(np.r_[np.ones(len(va)), np.zeros(len(vb))], np.r_[va, vb])
print(f'{FEAT}: R {len(R)}, E band {len(E)}')
for j in range(A.shape[1]):
    print(f'  col {j}: AUC R>dec {auc(grp["R"], grp["R_dec"], j):.3f}  E>dec {auc(grp["E"], grp["E_dec"], j):.3f}  '
          f'all-target-winner>decoy-winner {auc(w[w.is_decoy==0].row.values, w[w.is_decoy==1].row.values, j):.3f}')
