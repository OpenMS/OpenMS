"""HeLa (diagnostic). Uses the features table (every scored candidate window).
1. Candidates per precursor and the share at the max_peaks cap.
2. Wrong peak: is DIA-NN's apex among our candidates (then the W rule / scorer chose another one) or never picked?
3. Right peak below cut: per feature AUC of R vs score-matched decoys, and the same for accepted vs decoys, to see which
   evidence is weak. Also library RT error and run-time span of the RT window."""
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
B = '/home/singjc/.jcode/scratch/e2e/hela'
F = pd.read_csv(f'{B}/feat/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A1_1_2737.features.tsv', sep='\t')
f = pd.read_csv(f'{B}/funnel_hela.tsv', sep='\t')
n = F.groupby('precursor_id').size()
print('candidates per precursor:', n.describe().round(2).to_dict(), '| at cap (5):', round((n >= 5).mean(), 3))
tol = 3.43
wp = f[f.grp == 'wrong_peak'][['pid', 'RT']]
m = F.merge(wp, left_on='precursor_id', right_on='pid')
m['hit'] = (m.rt_apex - m.RT * 60).abs() <= tol
h = m.groupby('precursor_id').hit.any()
print(f'wrong peak {len(wp)}: DIA-NN apex among our candidates {int(h.sum())} ({h.mean():.2f}); never a candidate {int((~h).sum())}')
hit = m[m.hit].sort_values('score', ascending=False).drop_duplicates('precursor_id').set_index('precursor_id')
win = F.sort_values('score', ascending=False).drop_duplicates('precursor_id').set_index('precursor_id')
print(f'  when present: right-candidate score median {hit.score.median():.2f} vs chosen {win.loc[hit.index].score.median():.2f}; right one would pass the cut: {(hit.score >= 5.539).mean():.3f}')
nh = wp[wp.pid.isin(h[~h].index)]
print('  never a candidate: candidates per such precursor', n.loc[nh.pid].value_counts().sort_index().to_dict())
feats = [c for c in F.columns if c not in ('candidate_row', 'precursor_id', 'is_decoy', 'fold', 'rt_apex', 'score', 'q_value')]
rb = f[f.grp == 'right_below_cut']; ac = f[f.grp == 'accepted']
def right_rows(grp):
    mm = F.merge(grp[['pid', 'RT']], left_on='precursor_id', right_on='pid')
    mm = mm[(mm.rt_apex - mm.RT * 60).abs() <= tol]
    return mm.sort_values('score', ascending=False).drop_duplicates('precursor_id')
R, A = right_rows(rb), right_rows(ac)
dec = F[F.is_decoy == 1].sort_values('score', ascending=False).drop_duplicates('precursor_id').sort_values('score')
ds = dec.score.values
ix = np.clip(np.searchsorted(ds, R.score.values), 0, len(ds) - 1)
Dm = dec.iloc[np.unique(ix)]
print(f'R {len(R)} A {len(A)} score-matched decoys {len(Dm)}')
rows = []
for c in feats:
    y = np.r_[np.ones(len(R)), np.zeros(len(Dm))]
    a1 = roc_auc_score(y, np.r_[R[c].values, Dm[c].values])
    y2 = np.r_[np.ones(len(A)), np.zeros(len(dec))]
    a2 = roc_auc_score(y2, np.r_[A[c].values, dec[c].values])
    rows.append((c, round(a1, 3), round(a2, 3), round(R[c].median(), 3), round(A[c].median(), 3), round(Dm[c].median(), 3)))
print(pd.DataFrame(rows, columns=['feature', 'AUC R vs matched D', 'AUC A vs all D', 'R med', 'A med', 'D med']).to_string(index=False))
