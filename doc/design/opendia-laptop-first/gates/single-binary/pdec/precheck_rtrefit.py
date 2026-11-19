"""PRECHECK_RTREFIT (diagnostic): refit the library-RT -> run-RT map from the run's own first-pass IDs (tool q <= 0.01
target winners, label-free, as DIA-NN does) and recompute |RT error| for every window. Then: AUC of the old and new
RT error, and label-free GBM (same as precheck_scorer B, 2 rounds) with the new error replacing the old.
DIA-NN and entrapment labels for evaluation only."""
import sqlite3, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score
from sklearn.isotonic import IsotonicRegression
G = '/home/singjc/.jcode/scratch/dia2/gap7'
L = '/home/singjc/.jcode/scratch/dia2/entlib10/cut/hela_r5.predicted_decoys.pqp'
DN = '/home/singjc/.jcode/scratch/dia2/entlib12/D/report.parquet'
F = pd.read_parquet(f'{G}/D_feat.parquet')
meta = ['candidate_row', 'precursor_id', 'is_decoy', 'fold', 'rt_apex', 'score', 'q_value']
X = [c for c in F.columns if c not in meta]
c = sqlite3.connect(f'file:{L}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID pid, p.TRAML_ID t, p.DECOY d, p.CHARGE z, p.LIBRARY_RT lrt, p.LIBRARY_DRIFT_TIME lim, pe.MODIFIED_SEQUENCE ms,
  min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
pq['key'] = np.where(pq.d == 1, pq.t.str.replace('^DECOY_', '', regex=True), pq.t)
tg = pq[pq.d == 0]; r = (tg.e == 1).sum() / (tg.e == 0).sum(); k = 1 + 1 / r
P = pq.set_index('pid')
pair_id = pd.Series(pd.factorize(pq.key)[0], index=pq.pid); ent = P.e
d = pd.read_parquet(DN, columns=['Modified.Sequence', 'Precursor.Charge', 'Q.Value', 'RT', 'Precursor.Quantity'])
d = d[d['Q.Value'] <= 0.01]
tgs = tg.assign(sk=tg.ms + '/' + tg.z.astype(str)).drop_duplicates('sk').set_index('sk')
d['pid'] = tgs.pid.reindex(d['Modified.Sequence'] + '/' + d['Precursor.Charge'].astype(str)).values
d = d.dropna(subset=['pid']); d['pid'] = d.pid.astype(int); d = d[ent.reindex(d.pid).values == 0].drop_duplicates('pid').set_index('pid')
pid = F.precursor_id.values; dec = F.is_decoy.values.astype(bool); fold = F.fold.values; pair = pair_id.reindex(pid).values
lrt = P.lrt.reindex(pid).values; rt = F.rt_apex.values

def best_rows(score):
    o = np.lexsort((-score, pid)); first = np.r_[True, pid[o][1:] != pid[o][:-1]]
    rows = o[first]
    b = pd.DataFrame({'row': rows, 'pid': pid[rows], 's': score[rows], 'dec': dec[rows], 'pair': pair[rows]})
    b = b.sort_values(['pair', 's', 'dec'], ascending=[True, False, False]).drop_duplicates('pair')
    b = b.sort_values(['s', 'dec'], ascending=[False, False])
    nd = np.cumsum(b.dec.values); nt = np.cumsum(~b.dec.values)
    b['q'] = np.minimum.accumulate((nd / np.maximum(nt, 1))[::-1])[::-1]
    return b

def evaluate(score, name):
    b = best_rows(score); t = b[~b.dec.values].copy(); t['e'] = ent.reindex(t.pid).values
    acc = t[t.q <= 0.01]; N, E = len(acc), int(acc.e.sum())
    fdp = np.cumsum(t.e.values) * k / np.arange(1, len(t) + 1); ok = np.where(fdp <= 0.01)[0]
    print(f'{name:40s} q<=0.01 N {N:6d} lower {E / max(N, 1):.4f} combined {E * k / max(N, 1):.4f} | '
          f'N@FDP_c 1 % {int(ok.max() + 1) if len(ok) else 0:6d} | in DIA-NN {acc.pid.isin(d.index).mean():.3f}', flush=True)

kept = F.score.values > -1e29
s0 = np.where(kept, F.score.values, -1e30).astype(np.float64)
b0 = best_rows(s0); ids = b0[(~b0.dec.values) & (b0.q <= 0.01)]
x, y = lrt[ids.row.values], rt[ids.row.values]
print(f'first-pass IDs used for the refit: {len(ids)}')
# Isotonic fit on first-pass IDs, then a robust trim (residual > 3 MAD) and refit.
iso = IsotonicRegression(out_of_bounds='clip').fit(x, y)
res = y - iso.predict(x); mad = np.median(np.abs(res - np.median(res))) * 1.4826
keep = np.abs(res) <= 3 * mad
iso = IsotonicRegression(out_of_bounds='clip').fit(x[keep], y[keep])
res2 = y[keep] - iso.predict(x[keep])
print(f'refit residual on kept IDs: |res| p50 {np.median(np.abs(res2)):.2f} s, p95 {np.quantile(np.abs(res2), 0.95):.2f} s (kept {keep.sum()})')
new_err = np.abs(rt - iso.predict(lrt)).astype(np.float32)
old_err = F.abs_library_rt_error.values

dn_rt = (d.RT * 60).reindex(pid).values
posw = (~dec) & (np.abs(rt - dn_rt) <= 6)
Pw = pd.DataFrame({'pid': pid[posw], 'old': old_err[posw], 'new': new_err[posw], 's': s0[posw]}).sort_values('s', ascending=False).drop_duplicates('pid')
Pw['qq'] = pd.qcut(np.log10(d['Precursor.Quantity'].reindex(Pw.pid).clip(lower=1).values), 4, labels=['Q1', 'Q2', 'Q3', 'Q4'])
dsel = np.where(dec & kept)[0]; rng = np.random.default_rng(0); dsel = rng.choice(dsel, 400000, replace=False)
for name, v in (('old', old_err), ('new', new_err)):
    print(f'{name} RT error: DIA-NN-confirmed windows p50 {np.median(Pw[name]):.2f} p95 {np.quantile(Pw[name], 0.95):.2f} s; '
          f'decoy windows p50 {np.median(v[dsel]):.2f} s; AUC by quartile ' +
          ' '.join(f'{q} {roc_auc_score(np.r_[np.ones((Pw.qq == q).sum()), np.zeros(len(dsel))], -np.r_[Pw.loc[Pw.qq == q, name].values, v[dsel]]):.3f}' for q in ['Q1', 'Q2', 'Q3', 'Q4']), flush=True)

def gbm(Xm, pos_mask, neg_mask, seed=0, max_neg=1_500_000):
    s = np.zeros(len(F)); rng = np.random.default_rng(seed)
    for f in range(3):
        tp = np.where(pos_mask & (fold != f))[0]; tn = np.where(neg_mask & (fold != f))[0]
        if len(tn) > max_neg: tn = rng.choice(tn, max_neg, replace=False)
        idx = np.r_[tp, tn]; yy = np.r_[np.ones(len(tp)), np.zeros(len(tn))]
        w = np.r_[np.full(len(tp), len(tn) / max(len(tp), 1)), np.ones(len(tn))]
        m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=31, random_state=seed).fit(Xm[idx], yy, sample_weight=w)
        te = np.where(fold == f)[0]; s[te] = m.decision_function(Xm[te])
    return s

def winners(score):
    b = best_rows(score); m = np.zeros(len(F), bool); m[b.row.values[(b.q <= 0.01) & ~b.dec.values]] = True; return m

evaluate(s0, 'A tool score')
for label, cols, extra in (('old RT error', X, None), ('refit RT error', [x for x in X if x != 'abs_library_rt_error'], new_err)):
    Xm = F[cols].values.astype(np.float32)
    if extra is not None: Xm = np.c_[Xm, extra]
    s = s0
    for it in range(2):
        s = gbm(Xm, winners(s), dec, seed=it)
    evaluate(s, f'B label-free GBM, {label}')
    if extra is not None: np.save(f'{G}/score_B_rtrefit.npy', s)
