"""Can an in-run, teacher-free RT refit sharpen RT? (raw 01 design diagnostic.)
Train: our own accepted IDs (tool stage-2 search+boundary, q<=0.01, 746) -> observed apex RT. Model: residual
(observed - calibrated library RT) ~ ridge on sequence composition (20 AA counts, N/C-terminal AA one-hot, length,
charge, cubic in calibrated RT). Cross-fitted by the tool's folds. Evaluate on DIA-NN-confirmed precursors NOT in
the training set (DIA-NN RT as truth): robust sigma before vs after refit."""
import sqlite3, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
D = B + '/astral-design-raw01-20261107'
cr = pd.read_csv(D + '/stage2-boundary-design/boundary/crossfit.tsv', sep='\t', usecols=['sample_id', 'precursor_id', 'is_decoy', 'fold', 'precursor_best', 'q_value'])
s = pd.read_csv(D + '/pipeline/search.features.tsv', sep='\t', usecols=['candidate_row', 'rt_apex', 'apex_rt_error_seconds']).rename(columns={'candidate_row': 'sample_id'})
w = cr[cr.precursor_best == 1].merge(s, on='sample_id')
w['q'] = pd.to_numeric(w.q_value, errors='coerce'); w['cal_rt'] = w.rt_apex - w.apex_rt_error_seconds
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
pep = pd.read_sql('select pr.ID precursor_id, pr.CHARGE charge, p.UNMODIFIED_SEQUENCE seq from PRECURSOR pr join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=pr.ID join PEPTIDE p on p.ID=m.PEPTIDE_ID', c)
cal = cr[['precursor_id', 'fold']].drop_duplicates('precursor_id').merge(s.merge(cr[['sample_id', 'precursor_id']], on='sample_id').assign(cal_rt=lambda d: d.rt_apex - d.apex_rt_error_seconds).drop_duplicates('precursor_id')[['precursor_id', 'cal_rt']], on='precursor_id')
AA = 'ACDEFGHIKLMNPQRSTVWY'
def feats(t):
    t = t.merge(pep, on='precursor_id')
    cnt = np.array([[q.count(a) for a in AA] for q in t.seq], float)
    nt = np.array([[q[0] == a for a in AA] for q in t.seq], float); ct = np.array([[q[-1] == a for a in AA] for q in t.seq], float)
    n2 = np.array([[q[1] == a for a in AA] for q in t.seq], float); c2 = np.array([[q[-2] == a for a in AA] for q in t.seq], float)
    L = t.seq.str.len().values[:, None].astype(float)
    r = (t.cal_rt.values[:, None] - 700) / 400
    X = np.hstack([cnt / L, nt, ct, n2, c2, np.log(L), t.charge.values[:, None] == 2, r, r ** 2, r ** 3, np.ones_like(L)])
    return t, X
tr = w[(w.is_decoy == 0) & (w.q <= 0.01)][['precursor_id', 'fold', 'rt_apex', 'cal_rt']]
ref = pd.read_parquet(B + '/astral-diann-gap-20261108/ref_exact.parquet'); ref = ref[ref.Run.str.endswith('_01')]
ref = ref[['precursor_id']].assign(diann_rt=ref.RT.values * 60).merge(cal, on='precursor_id')
ev = ref[~ref.precursor_id.isin(tr.precursor_id)]
tr, Xtr = feats(tr); ev, Xev = feats(ev)
ytr = tr.rt_apex.values - tr.cal_rt.values
pred = np.zeros(len(ev))
for k in tr.fold.unique():
    m = tr.fold.values != k
    lam = 3.0
    A = Xtr[m]; wgt = np.linalg.solve(A.T @ A + lam * np.eye(A.shape[1]), A.T @ ytr[m])
    e = ev.fold.values == k
    pred[e] = Xev[e] @ wgt
rs = lambda r: round(1.4826 * np.median(np.abs(r - np.median(r))), 1)
before = ev.diann_rt.values - ev.cal_rt.values; after = before - pred
print(f'train IDs {len(tr)} | eval confirmed not in train {len(ev)}')
print('robust sigma before', rs(before), 'after', rs(after), '| q90 |res| before', round(np.quantile(np.abs(before - np.median(before)), .9), 1), 'after', round(np.quantile(np.abs(after - np.median(after)), .9), 1))
