"""PRECHECK_SCORER (diagnostic, r~7 D seed 0 CUR, PLAN_ODP code path with -features_out).
Where is the ID gap: scorer capacity or feature information?
  A  tool score as written (per-precursor best window, paired TDC as recomputed here)
  B  label-free: cross-fitted GBM on the same 24 features, Percolator-style self-training from the tool score
     (positives = target winners at q <= 0.01 of the previous round, negatives = decoy windows), folds as the tool's
  C  ceiling (uses DIA-NN labels, NOT a valid ID method): positives = windows within 6 s of DIA-NN RT of DIA-NN
     non-entrapment q <= 0.01 precursors, negatives = decoy windows
Metric: N at entrapment FDP_c 1 % (combined bound, r from the library) and the lower/combined bound at own q <= 0.01.
Entrapment labels are read for evaluation only. Usage: precheck_scorer.py"""
import os, sqlite3, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
G = '/home/singjc/.jcode/scratch/dia2/gap7'
L = '/home/singjc/.jcode/scratch/dia2/entlib10/cut/hela_r5.predicted_decoys.pqp'
DN = '/home/singjc/.jcode/scratch/dia2/entlib12/D/report.parquet'
if not os.path.exists(f'{G}/D_feat.parquet'):
    F = pd.read_csv(f'{G}/D_feat.tsv.gz', sep='\t', engine='c', dtype={'candidate_row': np.int32, 'precursor_id': np.int32,
                    'is_decoy': np.int8, 'fold': np.int8, 'rt_apex': np.float64, 'q_value': np.float64}, low_memory=False)
    for c in F.columns:
        if F[c].dtype == np.float64 and c not in ('rt_apex', 'q_value'): F[c] = F[c].astype(np.float32)
    F.to_parquet(f'{G}/D_feat.parquet')
F = pd.read_parquet(f'{G}/D_feat.parquet')
meta = ['candidate_row', 'precursor_id', 'is_decoy', 'fold', 'rt_apex', 'score', 'q_value']
X = [c for c in F.columns if c not in meta]
print(f'windows {len(F)}, precursors {F.precursor_id.nunique()}, features {len(X)}')

c = sqlite3.connect(f'file:{L}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID pid, p.TRAML_ID t, p.DECOY d, p.CHARGE z, pe.MODIFIED_SEQUENCE ms,
  min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
pq['key'] = np.where(pq.d == 1, pq.t.str.replace('^DECOY_', '', regex=True), pq.t)
tg = pq[pq.d == 0]; r = (tg.e == 1).sum() / (tg.e == 0).sum(); k = 1 + 1 / r
pair_id = pd.Series(pd.factorize(pq.key)[0], index=pq.pid)
ent = pq.set_index('pid').e
print(f'r = {r:.4f}')

# DIA-NN labels for the ceiling only.
d = pd.read_parquet(DN, columns=['Modified.Sequence', 'Precursor.Charge', 'Q.Value', 'RT'])
d = d[d['Q.Value'] <= 0.01]
tg = tg.assign(sk=tg.ms + '/' + tg.z.astype(str)).drop_duplicates('sk').set_index('sk')
d['pid'] = tg.pid.reindex(d['Modified.Sequence'] + '/' + d['Precursor.Charge'].astype(str)).values
d = d.dropna(subset=['pid']); d['pid'] = d.pid.astype(int); d = d[ent.reindex(d.pid).values == 0].drop_duplicates('pid')
dn_rt = d.set_index('pid').RT * 60
print(f'DIA-NN non-entrapment q<=0.01 mapped: {len(dn_rt)}')

pid = F.precursor_id.values; dec = F.is_decoy.values.astype(bool); fold = F.fold.values
pair = pair_id.reindex(pid).values

def evaluate(score, name):
    """Per-precursor best window, paired TDC (higher of target/decoy precursor wins, ties to decoy), q, bounds."""
    o = np.lexsort((-score, pid)); first = np.r_[True, pid[o][1:] != pid[o][:-1]]
    best = pd.DataFrame({'pid': pid[o][first], 's': score[o][first], 'dec': dec[o][first], 'pair': pair[o][first]})
    best = best.sort_values(['pair', 's', 'dec'], ascending=[True, False, False]).drop_duplicates('pair')
    best = best.sort_values(['s', 'dec'], ascending=[False, False])
    nd = np.cumsum(best.dec.values); nt = np.cumsum(~best.dec.values)
    q = np.minimum.accumulate((nd / np.maximum(nt, 1))[::-1])[::-1]
    t = best[~best.dec.values].assign(q=q[~best.dec.values]); t['e'] = ent.reindex(t.pid).values
    acc = t[t.q <= 0.01]; N, E = len(acc), int(acc.e.sum())
    fdp = np.cumsum(t.e.values) * k / np.arange(1, len(t) + 1); ok = np.where(fdp <= 0.01)[0]
    n1 = int(ok.max() + 1) if len(ok) else 0
    dn = acc.pid.isin(dn_rt.index).mean()
    print(f'{name:44s} q<=0.01 N {N:6d} lower {E / max(N, 1):.4f} combined {E * k / max(N, 1):.4f} | N@FDP_c 1 % {n1:6d} | in DIA-NN {dn:.3f}')
    return t

kept = F.score.values > -1e29
s_tool = np.where(kept, F.score.values, -1e30).astype(np.float64)
evaluate(s_tool, 'A tool score (kept windows)')

Xv = F[X].values.astype(np.float32)
def gbm_fit_predict(pos_mask, neg_mask, rounds=300, seed=0, max_neg=1_500_000):
    s = np.zeros(len(F))
    rng = np.random.default_rng(seed)
    for f in range(3):
        tr_pos = np.where(pos_mask & (fold != f))[0]; tr_neg = np.where(neg_mask & (fold != f))[0]
        if len(tr_neg) > max_neg: tr_neg = rng.choice(tr_neg, max_neg, replace=False)
        idx = np.r_[tr_pos, tr_neg]; y = np.r_[np.ones(len(tr_pos)), np.zeros(len(tr_neg))]
        w = np.r_[np.full(len(tr_pos), len(tr_neg) / max(len(tr_pos), 1)), np.ones(len(tr_neg))]
        m = HistGradientBoostingClassifier(max_iter=rounds, learning_rate=0.05, max_leaf_nodes=31, random_state=seed)
        m.fit(Xv[idx], y, sample_weight=w)
        te = np.where(fold == f)[0]; s[te] = m.decision_function(Xv[te])
    return s

def winners_at(score, qcut):
    t = evaluate(score, f'  (positives from previous round, q<={qcut})') if False else None
    o = np.lexsort((-score, pid)); first = np.r_[True, pid[o][1:] != pid[o][:-1]]
    rows = o[first]
    best = pd.DataFrame({'row': rows, 's': score[rows], 'dec': dec[rows], 'pair': pair[rows]})
    best = best.sort_values(['pair', 's', 'dec'], ascending=[True, False, False]).drop_duplicates('pair')
    best = best.sort_values(['s', 'dec'], ascending=[False, False])
    nd = np.cumsum(best.dec.values); nt = np.cumsum(~best.dec.values)
    q = np.minimum.accumulate((nd / np.maximum(nt, 1))[::-1])[::-1]
    m = np.zeros(len(F), bool); m[best.row.values[(q <= qcut) & ~best.dec.values]] = True
    return m

s = s_tool
for it, qc in enumerate((0.01, 0.01, 0.01)):
    pos = winners_at(s, qc)
    s = gbm_fit_predict(pos, dec, seed=it)
    evaluate(s, f'B label-free GBM round {it + 1} ({pos.sum()} positives)')
np.save(f'{G}/score_B.npy', s)

rt = F.rt_apex.values; drt = np.abs(rt - dn_rt.reindex(pid).values)
pos_c = (~dec) & (drt <= 6)
s_c = gbm_fit_predict(pos_c, dec, seed=9)
evaluate(s_c, f'C DIA-NN-labelled ceiling ({pos_c.sum()} positive windows)')
np.save(f'{G}/score_C.npy', s_c)
