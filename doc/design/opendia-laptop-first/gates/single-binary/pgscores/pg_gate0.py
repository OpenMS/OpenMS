"""PLAN_PGSCORES gate 0: labelled GBM ceiling on the common set (ablate.py method), our 24 vs 24 + 10 pg scores (diagnostic)."""
import pandas as pd, numpy as np, glob
from sklearn.ensemble import HistGradientBoostingClassifier
fs = pd.read_csv('img/funnel_imgrid_svm.tsv', sep='\t'); fd = pd.read_csv('img/funnel_imgrid_diann.tsv', sep='\t'); pos = set(fd.pid) & set(fs.pid)
F = pd.read_csv(glob.glob('pg/D_new_s0/*.features.tsv')[0], sep='\t')
meta = ('candidate_row', 'precursor_id', 'is_decoy', 'fold', 'rt_apex', 'score', 'q_value')
pgc = [x for x in F.columns if x.startswith('pg_')]; oc = [x for x in F.columns if x not in meta and x not in pgc]
assert len(oc) == 24 and len(pgc) == 10, (len(oc), len(pgc))
print(F[pgc].describe().T[['mean', 'std', 'min', 'max']].round(3).to_string())
# Winners per precursor ranked by the OpenDIA 24-feature score of the tip run, as in ablate.py (same common set rule:
# winners of precursors present in OpenSWATH's rank-1 table).
import sqlite3
R = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/openswath-hela-published'
c = sqlite3.connect(f'file:{R}/svm.osw?mode=ro', uri=True)
osw = set(pd.read_sql("select distinct f.PRECURSOR_ID pid from SCORE_MS2 s join FEATURE f on f.ID=s.FEATURE_ID where s.RANK=1", c).pid)
w = F.sort_values('score', ascending=False).drop_duplicates('precursor_id')
M = w[w.precursor_id.isin(osw)].reset_index(drop=True)
yy = np.where(M.is_decoy == 1, 0, np.where(M.precursor_id.isin(pos), 1, -1)); dec = M.is_decoy.values
def n_at(s):
    o = np.argsort(-s); d = np.cumsum(dec[o]); t = np.cumsum(1 - dec[o]); q = d / np.maximum(t, 1); qm = np.minimum.accumulate(q[::-1])[::-1]
    return int(((qm <= 0.01) & (dec[o] == 0)).sum())
def ceil(cc):
    s = np.zeros(len(M))
    for k in range(3):
        tr = (M.fold != k).values & (yy >= 0); te = (M.fold == k).values
        s[te] = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, random_state=0).fit(M.loc[tr, cc].fillna(0).values, yy[tr]).decision_function(M.loc[te, cc].fillna(0).values)
    return n_at(s)
b = ceil(oc); n = ceil(oc + pgc)
print(f'common set {len(M)} | 24 alone {b} | 24 + pg {n} | gain {n - b:+d} | gate 0 (>= +300): {"PASS" if n - b >= 300 else "FAIL"}')
from sklearn.metrics import roc_auc_score
m = yy >= 0
for x in pgc: print(f'  AUC {x:24s} {roc_auc_score(yy[m], M.loc[m, x].fillna(0)):.3f}')
