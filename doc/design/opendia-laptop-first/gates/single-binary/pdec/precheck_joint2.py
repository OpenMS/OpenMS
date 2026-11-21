"""PRECHECK_JOINT2: precheck_joint.py plus arms Bt (hard decoys) and Bts (hard decoys + siblings); optional 6th arg 'nolabel' stops after them.
PRECHECK_JOINT (diagnostic, 21 Nov 2026): do sibling (charge-state) features or interference explain the missed IDs?
Usage: precheck_joint.py <feat.parquet> <lib.pqp> <diann report.parquet> <tag> [r]"""
import sys, time, sqlite3, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
PQ, L, DN, TAG = sys.argv[1:5]; R = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0; T0 = time.time()
F = pd.read_parquet(PQ).reset_index(drop=True)
meta = ['candidate_row', 'precursor_id', 'is_decoy', 'fold', 'rt_apex', 'score', 'q_value']
X0 = [c for c in F.columns if c not in meta]
c = sqlite3.connect(f'file:{L}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID pid, p.TRAML_ID t, p.DECOY d, p.CHARGE z, pe.ID pep, pe.MODIFIED_SEQUENCE ms,
  min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
pq['key'] = np.where(pq.d == 1, pq.t.str.replace('^DECOY_', '', regex=True), pq.t)
tg = pq[pq.d == 0]; k = 1 + 1 / R if R > 0 else np.nan
pair_id = pd.Series(pd.factorize(pq.key)[0], index=pq.pid); ent = pq.set_index('pid').e; pep_of = pq.set_index('pid').pep
d = pd.read_parquet(DN, columns=['Modified.Sequence', 'Precursor.Charge', 'Q.Value', 'RT'])
d = d[d['Q.Value'] <= 0.01]
tgs = tg.assign(sk=tg.ms + '/' + tg.z.astype(str)).drop_duplicates('sk').set_index('sk')
d['pid'] = tgs.pid.reindex(d['Modified.Sequence'] + '/' + d['Precursor.Charge'].astype(str)).values
d = d.dropna(subset=['pid']); d['pid'] = d.pid.astype(np.int64)
d = d[ent.reindex(d.pid).values == 0].sort_values('Q.Value').drop_duplicates('pid').set_index('pid')
pid = F.precursor_id.values; dec = F.is_decoy.values.astype(bool); fold = F.fold.values; pair = pair_id.reindex(pid).values
kept = F.score.values > -1e29
s_tool = np.where(kept, F.score.values, -1e30).astype(np.float64)

# ---- sibling features (same PEPTIDE row = same modified sequence, other charge)
F['pep'] = pep_of.reindex(pid).values
npep = pq.groupby('pep').pid.nunique(); F['s_n'] = (npep.reindex(F.pep).values - 1).astype(np.float32)
smin = float(np.min(s_tool[kept])) if kept.any() else 0.0
W = pd.DataFrame({'i': np.arange(len(F)), 'pep': F.pep.values, 'pid': pid, 'rt': F.rt_apex.values, 's': np.where(kept, s_tool, smin)})
W = W[F.s_n.values > 0]
M = W.merge(W, on='pep', suffixes=('', '_o'))
M = M[(M.pid != M.pid_o) & ((M.rt - M.rt_o).abs() <= 3)]
M['dd'] = (M.rt - M.rt_o).abs()
g = M.sort_values('s_o', ascending=False).drop_duplicates('i').set_index('i')
co = M.groupby('i').pid_o.nunique()
F['s_co'] = co.reindex(np.arange(len(F))).fillna(0).values.astype(np.float32)
F['s_best'] = g.s_o.reindex(np.arange(len(F))).fillna(smin).values.astype(np.float32)
F['s_best_d'] = g.dd.reindex(np.arange(len(F))).fillna(3.0).values.astype(np.float32)
del M, W
XS = X0 + ['s_n', 's_co', 's_best', 's_best_d']
print(f'{TAG}: windows {len(F)}, base features {len(X0)}, with siblings {len(XS)}; windows with a co-eluting sibling '
      f'target {(F.s_co.values[~dec] > 0).mean():.3f} decoy {(F.s_co.values[dec] > 0).mean():.3f}  [{time.time() - T0:.0f}s]', flush=True)

def best_rows(score):
    o = np.lexsort((-score, pid)); first = np.r_[True, pid[o][1:] != pid[o][:-1]]; rows = o[first]
    b = pd.DataFrame({'row': rows, 'pid': pid[rows], 's': score[rows], 'dec': dec[rows], 'pair': pair[rows]})
    b = b.sort_values(['pair', 's', 'dec'], ascending=[True, False, False]).drop_duplicates('pair')
    b = b.sort_values(['s', 'dec'], ascending=[False, False])
    nd = np.cumsum(b.dec.values); nt = np.cumsum(~b.dec.values)
    b['q'] = np.minimum.accumulate((nd / np.maximum(nt, 1))[::-1])[::-1]
    return b
RES = {}
def evaluate(score, name):
    b = best_rows(score); t = b[~b.dec.values].copy(); t['e'] = ent.reindex(t.pid).values
    acc = t[t.q <= 0.01]; N, E = len(acc), int(acc.e.sum()); hela = acc[acc.e == 0]
    s = f'{name:34s} q<=0.01 N {N:6d} HeLa {len(hela):6d} (DIA-NN covered {d.index.isin(hela.pid).mean():.3f})'
    n1 = np.nan
    if R > 0:
        fdp = np.cumsum(t.e.values) * k / np.arange(1, len(t) + 1); ok = np.where(fdp <= 0.01)[0]; n1 = int(ok.max() + 1) if len(ok) else 0
        s += f' lower {E / max(N, 1):.4f} combined {E * k / max(N, 1):.4f} N@FDP_c1% {n1}'
    print(s + f'  [{time.time() - T0:.0f}s]', flush=True)
    RES[name] = dict(hela=len(hela), n1=n1, lower=E / max(N, 1)); return b
def gbm(Xv, pos_mask, neg_mask, seed=0, max_neg=1_500_000):
    s = np.zeros(len(F)); rng = np.random.default_rng(seed)
    for f in range(3):
        tp = np.where(pos_mask & (fold != f))[0]; tn = np.where(neg_mask & (fold != f))[0]
        if len(tn) > max_neg: tn = rng.choice(tn, max_neg, replace=False)
        idx = np.r_[tp, tn]; y = np.r_[np.ones(len(tp)), np.zeros(len(tn))]
        w = np.r_[np.full(len(tp), len(tn) / max(len(tp), 1)), np.ones(len(tn))]
        m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=31, random_state=seed).fit(Xv[idx], y, sample_weight=w)
        te = np.where(fold == f)[0]; s[te] = m.decision_function(Xv[te])
    return s
def winners(score):
    b = best_rows(score); m = np.zeros(len(F), bool); m[b.row.values[(b.q <= 0.01) & ~b.dec.values]] = True; return m
bA = evaluate(s_tool, 'A tool score')
for name, cols in (('B label-free GBM', X0), ('Bs label-free GBM + siblings', XS)):
    Xv = F[cols].values.astype(np.float32); s = s_tool
    for it in range(3): s = gbm(Xv, winners(s), dec, seed=it)
    evaluate(s, name)
def hard_decoys(score, frac=0.02):
    b = best_rows(score); bd = b[b.dec.values]; top = set(bd.pid.values[:max(1, int(frac * len(bd)))])
    return dec & np.isin(pid, list(top))
for name, cols in (('Bt hard decoys', X0), ('Bts hard decoys + siblings', XS)):
    Xv = F[cols].values.astype(np.float32); s = s_tool
    for it in range(3): s = gbm(Xv, winners(s), hard_decoys(s), seed=20 + it)
    evaluate(s, name)
print('JOINT2 ' + TAG + ' ' + ' '.join(f"{n.split()[0]}={v['n1']}/{v['lower']:.4f}/{v['hela']}" for n, v in RES.items()), flush=True)
if len(sys.argv) > 6 and sys.argv[6] == 'nolabel': sys.exit(0)
drt = np.abs(F.rt_apex.values - (d.RT * 60).reindex(pid).values); right = (~dec) & (drt <= 6)
for name, cols in (('C DIA-NN ceiling', X0), ('Cs DIA-NN ceiling + siblings', XS)):
    evaluate(gbm(F[cols].values.astype(np.float32), right, dec, seed=9), name)

# ---- interference and sibling profile by class (tool score)
tb = bA[~bA.dec.values].set_index('pid'); db = bA[bA.dec.values]
acc_t = tb[tb.q <= 0.01]
missed = d.index[~d.index.isin(acc_t.index)]
best_row = tb.row.reindex(missed)
A_low = best_row.dropna().astype(int); A_low = A_low[right[A_low.values]]
thr = acc_t.s.min()
near_dec = db[(db.s > thr - 2) & (db.s <= thr + 2)].row.values      # decoys around the 1 % cut
acc_rows = acc_t[acc_t.index.map(ent).fillna(1) == 0].row.values
acc_set_rt = pd.Series(F.rt_apex.values[acc_t.row.values], index=acc_t.index)  # accepted precursors (targets) -> rt
acc_dec = bA[bA.dec.values & (bA.q <= 0.01)]  # decoys never 'accepted' as IDs; use decoy winners above the target cut
acc_dec_rt = pd.Series(F.rt_apex.values[db[db.s >= thr].row.values], index=db[db.s >= thr].pid.values)
pep_of_pid = pq.set_index('pid').pep
def prof(rows, label, decoy_side=False):
    rows = np.asarray(rows)
    ie = F.i_explained.values[rows] if 'i_explained' in F else np.full(len(rows), np.nan)
    pool = acc_dec_rt if decoy_side else acc_set_rt
    P = pd.DataFrame({'pid_o': pool.index.values, 'rt_o': pool.values}); P['pep'] = pep_of_pid.reindex(P.pid_o).values
    Q = pd.DataFrame({'r': rows, 'pid': pid[rows], 'rt': F.rt_apex.values[rows], 'pep': F.pep.values[rows]})
    J = Q.merge(P, on='pep'); J = J[(J.pid != J.pid_o) & ((J.rt - J.rt_o).abs() <= 3)]
    hit = J.r.nunique()
    print(f'  {label:44s} n {len(rows):6d}  i_explained median {np.nanmedian(ie):.3f}  share>=0.5 {np.nanmean(ie >= 0.5):.3f}  '
          f'i_nclaim median {np.nanmedian(F.i_nclaim.values[rows]):.2f}  co-eluting accepted sibling {hit / max(len(rows), 1):.3f}  '
          f'has sibling in library {(F.s_n.values[rows] > 0).mean():.3f}', flush=True)
    return np.nanmean(ie >= 0.5)
print(f'\nprofiles (tool score; cut {thr:.3f}):')
a = prof(A_low.values, 'A_low (right window best, not accepted)')
b = prof(acc_rows, 'accepted HeLa target best windows')
prof(near_dec, 'decoy best windows within +-2 of the cut', decoy_side=True)
nb, ns = RES['B label-free GBM'], RES['Bs label-free GBM + siblings']; cb, cs = RES['C DIA-NN ceiling'], RES['Cs DIA-NN ceiling + siblings']
print(f'\nSUMMARY {TAG}: HeLa Bs/B {ns["hela"] / nb["hela"] - 1:+.2%}, Cs/C {cs["hela"] / cb["hela"] - 1:+.2%}'
      + (f'; N@FDP_c 1 % Bs/B {ns["n1"] / nb["n1"] - 1:+.2%}, lower {nb["lower"]:.4f} -> {ns["lower"]:.4f}' if R > 0 else '')
      + f'; A_low i_explained>=0.5 {a:.3f} vs accepted {b:.3f} (ratio {a / max(b, 1e-9):.2f})')
