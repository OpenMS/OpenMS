"""PRECHECK_FAINT step 2: do new per-peak trace features carry information the 36 features lack?
Sample: 6,000 DIA-NN targets + 30,000 decoys (faint_mini.py), r~1.56 D s0, search valley 0.8 / max 8 (= ev/v08m8 arm).
New features per candidate peak, from the search's own trace matrix (6 fragments x native scans) over [rt_start, rt_end]:
  ps_cos_w      consensus-weighted mean over scans of cosine(fragment vector, sqrt library intensities)
  ps_cos_apex   the same cosine on apex +-1 scans summed
  ps_n4         share of peak scans with >= 4 fragments nonzero
  apex_spread   std (scans) of per-fragment apex positions, fragments with signal
  frag_run      mean over fragments of the longest nonzero run inside the peak / peak scans
  local_sn      log apex consensus / median consensus of the trace outside the peak
  rank_corr     Spearman of apex+-1 fragment sums vs library intensities
Test: cross-fitted (feature 'fold') GBM on DIA-NN labels (positives: windows within 6 s of DIA-NN RT; negatives: decoy
windows), 36 features vs 36 + new. Precursor-level max score; TPR of DIA-NN targets at decoy FPR 1e-3 / 1e-2, by
quartile. A ceiling comparison: gain means new information, not a method. Labels for evaluation only."""
import glob, sqlite3, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from scipy.stats import spearmanr
G = '/home/singjc/.jcode/scratch/dia2/gap7'; O = f'{G}/faint'
S = pd.read_csv(f'{O}/sample.tsv', sep='\t').set_index('pid')
P = pd.read_csv(f'{O}/peaks.tsv', sep='\t')
c = sqlite3.connect(f'{O}/mini.pqp')
lib = pd.read_sql('select m.PRECURSOR_ID pid, t.PRODUCT_MZ mz, t.LIBRARY_INTENSITY li from TRANSITION t join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID', c)
libg = {p: g for p, g in lib.groupby('pid')}
T = pd.read_csv(f'{O}/traces.tsv', sep='\t', dtype={'precursor_id': np.int32, 'scan': np.int32, 'slot': np.int8})
rows = []
for pid, t in T.groupby('precursor_id', sort=False):
    m = t.pivot_table(index='rt', columns='slot', values='intensity', aggfunc='max').fillna(0.0)
    fm = t.drop_duplicates('slot').set_index('slot').fragment_mz.reindex(m.columns).values
    lg = libg[pid]; li = np.array([lg.li.values[np.argmin(np.abs(lg.mz.values - x))] for x in fm]); w = np.sqrt(np.maximum(li, 0))
    w = w / (np.linalg.norm(w) + 1e-12)
    rt = m.index.values; X = m.values; cons = X.sum(1)
    for _, pk in P[P.precursor_id == pid].iterrows():
        a, b = np.searchsorted(rt, pk.rt_start - 1e-6), np.searchsorted(rt, pk.rt_end + 1e-6)
        if b - a < 1: continue
        Y = X[a:b]; cy = cons[a:b]; ia = int(np.argmin(np.abs(rt[a:b] - pk.rt_apex)))
        nrm = np.linalg.norm(Y, axis=1); cs = np.where(nrm > 0, (Y @ w) / np.maximum(nrm, 1e-12), 0)
        ps_cos_w = float((cs * cy).sum() / max(cy.sum(), 1e-12))
        v = Y[max(ia - 1, 0):ia + 2].sum(0); ps_cos_apex = float(v @ w / max(np.linalg.norm(v), 1e-12))
        ps_n4 = float(((Y > 0).sum(1) >= 4).mean())
        has = Y.max(0) > 0; apex_spread = float(np.std(np.argmax(Y[:, has], 0))) if has.sum() > 1 else 10.0
        runs = []
        for k in range(Y.shape[1]):
            z = np.r_[0, (Y[:, k] > 0).astype(int), 0]; d = np.diff(z); st, en = np.where(d == 1)[0], np.where(d == -1)[0]
            runs.append((en - st).max() if len(st) else 0)
        frag_run = float(np.mean(runs) / (b - a))
        out = np.r_[cons[:a], cons[b:]]; bg = np.median(out) if len(out) else 0.0
        local_sn = float(np.log1p(cons[a + ia]) - np.log1p(bg))
        rc = spearmanr(v, li).correlation if (v > 0).sum() >= 3 else 0.0
        rows.append((pid, round(pk.rt_apex, 3), ps_cos_w, ps_cos_apex, ps_n4, apex_spread, frag_run, local_sn, 0.0 if np.isnan(rc) else rc))
NF = ['ps_cos_w', 'ps_cos_apex', 'ps_n4', 'apex_spread', 'frag_run', 'local_sn', 'rank_corr']
N = pd.DataFrame(rows, columns=['precursor_id', 'rta'] + NF)
print(f'new-feature rows {len(N)} for {N.precursor_id.nunique()} precursors', flush=True)
# 36 workflow features of the same arm, sample precursors only.
ids = set(S.index); parts = []
for ch in pd.read_csv(glob.glob(f'{G}/ev/v08m8/*.features.tsv')[0], sep='\t', chunksize=1_000_000):
    parts.append(ch[ch.precursor_id.isin(ids)])
F = pd.concat(parts); F['rta'] = F.rt_apex.round(3)
F = F.merge(N, on=['precursor_id', 'rta'], how='inner')
print(f'joined windows {len(F)} ({F.precursor_id.nunique()} precursors)', flush=True)
meta = ['candidate_row', 'precursor_id', 'is_decoy', 'fold', 'rt_apex', 'score', 'q_value', 'rta']
X36 = [c for c in F.columns if c not in meta + NF]
dec = F.is_decoy.values == 1; isd = F.precursor_id.map(S.cls).values == 'dnn'
F['right'] = isd & (np.abs(F.rt_apex - F.precursor_id.map(S.rt_s)) <= 6)
pos = F.right.values; fold = F.fold.values
def fit(cols, seed=0):
    Xv = F[cols].values.astype(np.float32); s = np.zeros(len(F))
    for f in np.unique(fold):
        tr = (fold != f) & (pos | dec); y = pos[tr].astype(int); wts = np.where(y == 1, (y == 0).sum() / max(y.sum(), 1), 1.0)
        m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=31, random_state=seed).fit(Xv[tr], y, sample_weight=wts)
        s[fold == f] = m.decision_function(Xv[fold == f])
    return s
qq = S.qq
def report(name, s):
    F['_s'] = s; best = F.groupby('precursor_id')._s.max()
    # right-window score for DIA-NN targets (does the right peak beat the decoy tail?)
    rb = F[F.right].groupby('precursor_id')._s.max()
    neg = best[best.index.map(S.cls) == 'decoy'].values
    out = []
    for fpr in (1e-3, 1e-2):
        thr = np.quantile(neg, 1 - fpr)
        out.append(f'FPR {fpr:g}: ' + ' '.join(f'{k} {(rb[rb.index.map(qq) == k] > thr).mean():.3f}' for k in ['Q1', 'Q2', 'Q3', 'Q4']) + f' all {(rb > thr).mean():.3f}')
    print(f'{name:28s} ' + ' | '.join(out), flush=True)

rng = np.random.default_rng(0)
def tpr(s, fpr, idx_dec, idx_pos):
    F['_s'] = s; best = F.groupby('precursor_id')._s.max(); rb = F[F.right].groupby('precursor_id')._s.max()
    return best, rb
res = {}
for seed in (0, 1):
    res[('36', seed)] = fit(X36, seed); res[('43', seed)] = fit(X36 + NF, seed)
dec_ids = np.array([p for p in S.index if S.cls[p] == 'decoy']); pos_ids = None
out = {}
for key, s in res.items():
    F['_s'] = s; best = F.groupby('precursor_id')._s.max(); rb = F[F.right].groupby('precursor_id')._s.max()
    out[key] = (best.reindex(dec_ids).values, rb)
pos_ids = out[('36', 0)][1].index.values
for fpr in (1e-3, 3e-3, 1e-2):
    for seed in (0, 1):
        diffs = []
        for b in range(300):
            di = rng.integers(0, len(dec_ids), len(dec_ids)); pi = rng.integers(0, len(pos_ids), len(pos_ids))
            v = []
            for k in ('36', '43'):
                neg, rb = out[(k, seed)]; thr = np.quantile(neg[di], 1 - fpr); v.append((rb.reindex(pos_ids).values[pi] > thr).mean())
            diffs.append(v[1] - v[0])
        neg36, rb36 = out[('36', seed)]; neg43, rb43 = out[('43', seed)]
        p36 = (rb36 > np.quantile(neg36, 1 - fpr)).mean(); p43 = (rb43 > np.quantile(neg43, 1 - fpr)).mean()
        print(f'FPR {fpr:g} seed {seed}: TPR 36 {p36:.3f} -> 43 {p43:.3f} (diff {p43 - p36:+.3f}, bootstrap 90 % {np.quantile(diffs, 0.05):+.3f} to {np.quantile(diffs, 0.95):+.3f})', flush=True)
