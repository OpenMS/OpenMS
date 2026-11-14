"""PLAN_GAP arm G: run-specific fragment intensity correction learned from the calibration anchors (label-free).
The predicted 12-fragment intensities (lib12) are corrected per fragment by a model fitted on the run's anchors
(cal/anchors.tsv: confident target precursors of the calibration, used = 1), then the 12-fragment evidence is
recomputed with the corrected intensities.
Fit (anchors only):
  - Each anchor's window is its population window with apex nearest the anchor run RT, within the calibrated peak
    extent (cal.json peak_boundary.half_widths x median population half width FWHM / 2); else skipped.
  - Central area a_f (bins 12:20) of its 12 fragments; scale = median over fragments with a_f > 0 of a_f / lib_f.
  - Target y_f = log((a_f + N) / (scale lib_f + N)), with N the deficit_features noise floor (median nonzero central
    area over 200,000 windows, numpy seed 0). Every library fragment counts, also when a_f = 0.
  - Covariates (from the 12-fragment PQP annotation): ion type x charge; fragment m/z decile (deciles of the anchor
    fragments); relative position ordinal / (peptide length); log10 predicted intensity decile.
  - Model: ridge regression on the one-hot covariates, penalty chosen by 5-fold CV over the anchors (grouped by anchor).
    The CV gain (mean squared error of y with and without the model) is printed: the arm is only useful when > 0.
Apply (every population window, targets and decoys alike):
  lib_corr_f = lib_f exp(yhat_f), renormalised to max 1.
Columns of featG.npy [n, 4] (central bins as features12.py):
  g_libcos     cosine of sqrt(a) and sqrt(lib_corr) over the 12
  g_support6   fragments among the top 6 by lib_corr with a_f > 0 and profile r > 0.5 to the sum of the others
  g_coel6      mean of those correlations over the top 6 by lib_corr
  g_dmax       max over fragments of log10(scale lib_corr_f + N) - log10(a_f + N) (deficit_features, corrected library)
Usage: features_G.py <raw01|raw02> <root> <library12_pop.pqp> <search.features.tsv>."""
import json, sqlite3, sys, warnings, numpy as np, pandas as pd
RAW, O, LIB12, FEAT = sys.argv[1:5]
C0, C1, BS, K = 12, 20, 100000, 12
warnings.filterwarnings('ignore', category=RuntimeWarning)
pop = pd.read_csv(f'{O}/{RAW}/population.tsv', sep='\t')
X = np.load(f'{O}/{RAW}/xic12.npy', mmap_mode='r'); L = np.load(f'{O}/{RAW}/lib12.npy', mmap_mode='r'); n = len(X)
cal = json.load(open(f'{O}/{RAW}/cal/cal.json'))
samp = np.sort(np.random.default_rng(0).choice(n, 200000, replace=False))
xa = np.asarray(X[samp], np.float64)[:, :, C0:C1].sum(2); la = np.asarray(L[samp]) > 0
N = float(np.median(xa[(xa > 0) & la]))

# Slot annotations (build12.py order: per precursor by transition ID).
c = sqlite3.connect(f'file:{LIB12}?mode=ro', uri=True)
t = pd.read_sql("""select m.PRECURSOR_ID pid, t.ID tid, t.ANNOTATION ann, t.PRODUCT_MZ mz from TRANSITION t
  join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID order by m.PRECURSOR_ID, t.ID""", c)
seq = pd.read_sql("""select p.ID pid, pe.UNMODIFIED_SEQUENCE seq from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
  join PEPTIDE pe on pe.ID=m.PEPTIDE_ID""", c).drop_duplicates('pid').set_index('pid').seq.str.len()
t['slot'] = t.groupby('pid').cumcount()
t = t[t.slot < K]
a = t.ann.str.extract(r'^([a-z])(\d+)\^(\d+)')
t['ion'] = a[0] + a[2]; t['ordrel'] = a[1].astype(float) / t.pid.map(seq).values
pids = pop.precursor_id.values
uniq, inv = np.unique(pids, return_inverse=True)
def slot_matrix(col, fill):
    M = np.full((len(uniq), K), fill, dtype=object if col == 'ion' else np.float64)
    i = np.searchsorted(uniq, t.pid.values); ok = (i < len(uniq)) & (uniq[np.minimum(i, len(uniq) - 1)] == t.pid.values)
    M[i[ok], t.slot.values[ok]] = t[col].values[ok]
    return M
ION, MZ, ORD = slot_matrix('ion', 'none'), slot_matrix('mz', np.nan), slot_matrix('ordrel', np.nan)

# Anchor windows.
anc = pd.read_csv(f'{O}/{RAW}/cal/anchors.tsv', sep='\t'); anc = anc[anc.used == 1]
# Anchor-window match tolerance from the data: the calibrated peak extent (peak_boundary.half_widths, in half widths)
# times the median half width (FWHM / 2) of the population windows.
fw = pd.read_csv(FEAT, sep='\t', usecols=['candidate_row', 'fwhm_seconds']).set_index('candidate_row').fwhm_seconds
tol = float(cal['peak_boundary']['half_widths']) * float(np.median(fw.reindex(pop.sample_id).values)) / 2
pa = pop.reset_index().merge(anc[['precursor_id', 'run_rt']], on='precursor_id')
pa['d'] = (pa.rt_apex - pa.run_rt).abs()
pa = pa.sort_values('d').drop_duplicates('precursor_id'); pa = pa[pa.d <= tol]
rows = np.sort(pa['index'].values)
xa = np.asarray(X[rows], np.float64)[:, :, C0:C1].sum(2); li = np.asarray(L[rows], np.float64)
pres = li > 0
scale = np.nanmedian(np.where(pres & (xa > 0), xa / np.maximum(li, 1e-12), np.nan), 1)
y = np.log((xa + N) / (scale[:, None] * li + N))
ion, mz, ordr = ION[inv[rows]], MZ[inv[rows]], ORD[inv[rows]]
m = pres & np.isfinite(y) & np.isfinite(mz) & np.isfinite(scale)[:, None]
grp_anchor = np.repeat(np.arange(len(rows)), K).reshape(len(rows), K)[m]
yv = y[m]
ions = sorted(set(ion[m]))
mz_edges = np.quantile(mz[m], np.linspace(0, 1, 11)[1:-1]); ord_edges = np.quantile(ordr[m], np.linspace(0, 1, 11)[1:-1])
li_edges = np.quantile(np.log10(li[m]), np.linspace(0, 1, 11)[1:-1])
def design(ion_, mz_, ord_, li_):
    cols = [np.asarray([i == v for i in ion_], float) for v in ions]
    for e, v in ((mz_edges, mz_), (ord_edges, ord_), (li_edges, np.log10(np.maximum(li_, 1e-12)))):
        b = np.searchsorted(e, v); cols += [(b == k).astype(float) for k in range(len(e) + 1)]
    return np.stack(cols, 1)
Dm = design(ion[m], mz[m], ordr[m], li[m])
def ridge(A, b, lam):
    mu = b.mean(); return np.linalg.solve(A.T @ A + lam * np.eye(A.shape[1]), A.T @ (b - mu)), mu
folds = grp_anchor % 5
best = None
for lam in 10.0 ** np.arange(-2, 5):
    err = 0.0; base = 0.0
    for k in range(5):
        tr, te = folds != k, folds == k
        w, mu = ridge(Dm[tr], yv[tr], lam)
        err += ((yv[te] - mu - Dm[te] @ w) ** 2).sum(); base += ((yv[te] - yv[tr].mean()) ** 2).sum()
    if best is None or err < best[0]: best = (err, lam, base)
err, lam, base = best
w, mu = ridge(Dm, yv, lam)
print(f'{RAW}: anchor match tolerance {tol:.2f} s; anchors used {len(rows)} (of {len(anc)}), fragments {len(yv)}, N {N:.0f}; ridge lambda {lam:g}; '
      f'CV MSE {err / len(yv):.4f} vs no-model {base / len(yv):.4f} (gain {1 - err / base:.3f})')
coef = dict(zip([f'ion {v}' for v in ions] + [f'mz{k}' for k in range(10)] + [f'ord{k}' for k in range(10)] + [f'li{k}' for k in range(10)], w.round(3)))
print('  coefficients', coef)

out = np.zeros((n, 4), np.float32)
def corr(a_, b_):
    a_ = a_ - a_.mean(-1, keepdims=True); b_ = b_ - b_.mean(-1, keepdims=True)
    return (a_ * b_).sum(-1) / np.sqrt((a_ * a_).sum(-1) * (b_ * b_).sum(-1) + 1e-12)
for b in range(0, n, BS):
    r = np.arange(b, min(b + BS, n)); x = np.asarray(X[b:b + BS], np.float64); li = np.asarray(L[b:b + BS], np.float64)
    pres = li > 0; ii = inv[r]
    yh = np.zeros_like(li)
    if pres.any():
        Dp = design(ION[ii][pres], MZ[ii][pres], ORD[ii][pres], li[pres]); yh[pres] = mu + Dp @ w
    lc = np.where(pres, li * np.exp(yh), 0.0); lc = lc / np.maximum(lc.max(1, keepdims=True), 1e-12)
    area = x[:, :, C0:C1].sum(2)
    cos = (np.sqrt(area) * np.sqrt(lc)).sum(1) / np.sqrt(area.sum(1) * lc.sum(1) + 1e-12)
    top6 = np.argsort(-lc, 1)[:, :6]; rr = np.arange(len(x))[:, None]
    tot = x.sum(1, keepdims=True); rc = corr(x, tot - x)
    rc6 = rc[rr, top6]; sig6 = (x[:, :, C0:C1].max(2) > 0)[rr, top6] & pres[rr, top6]
    sup = (sig6 & (rc6 > 0.5)).sum(1); coel = np.where(pres[rr, top6], rc6, 0).mean(1)
    sg = pres & (area > 0)
    sc = np.where(sg.any(1), np.nanmedian(np.where(sg, area / np.maximum(lc, 1e-12), np.nan), 1), 0.0)
    d = np.where(pres, np.log10(sc[:, None] * lc + N) - np.log10(area + N), -np.inf).max(1)
    out[b:b + BS] = np.c_[cos, sup, coel, np.where(np.isfinite(d), d, 0)]
np.save(f'{O}/{RAW}/featG.npy', out)
print(RAW, 'featG', out.shape, 'finite', bool(np.isfinite(out).all()), 'p5/p50/p95', np.percentile(out, [5, 50, 95], axis=0).round(3).T.tolist())
