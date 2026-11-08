"""Gap anatomy v3 (raw 01, OSW XICs): stable, learner-fixed information measure per data treatment.

v2 showed the semi-supervised LDA bootstrap is unstable at 1% on hand features, so counts there do not isolate the data
treatment. Here every regime is scored by the SAME supervised cross-fitted oracle (positives: DIA-NN-confirmed right
peaks; negatives: all decoy peak groups; pair-hash folds), an upper bound on what any learner can extract from those
features. Then winner per precursor, paired TDC, count at q<=0.01/0.05 and FDP_c.
Feature sets: 'simple' (v2's 11) and 'simple+xcorr' (adds OpenSWATH-style lagged cross-correlation co-elution/shape,
library-weighted shape, library dot product and Manhattan distance, MS1-vs-fragment xcorr shape) -- to test whether
the extra information OSW gets comes from these specific, cheap, well-known transforms of the same traces.
"""
import numpy as np, pandas as pd, numba as nb, torch, torch.nn as nn
import sys; sys.path.insert(0, '.')
from regimes_lib import *  # noqa: F401,F403 (data, numba kernels, lda/winners/select/auc)
torch.set_num_threads(8)
NX = 6
XN = ['xcorr_coelution', 'xcorr_shape', 'xcorr_shape_w', 'lib_dotprod', 'lib_manhattan', 'ms1_xcorr_shape']

@nb.njit(cache=True)
def xcorr_max(a, b, maxlag):
    n = len(a); ma = a.mean(); mb = b.mean(); sa = np.sqrt(((a - ma) ** 2).sum()); sb = np.sqrt(((b - mb) ** 2).sum())
    if sa <= 0 or sb <= 0: return 0.0, maxlag
    best = -2.0; bl = 0
    for lag in range(-maxlag, maxlag + 1):
        s = 0.0
        for i in range(n):
            j = i + lag
            if 0 <= j < n: s += (a[i] - ma) * (b[j] - mb)
        v = s / (sa * sb)
        if v > best: best = v; bl = lag
    return best, abs(bl)

@nb.njit(parallel=True, cache=True)
def xfeat(frag, ms1, lib, npts, rt, PID, apex, lft, rgt, mode):
    n_ = len(apex); out = np.full((n_, NX), np.nan, np.float32)
    for g in nb.prange(n_):
        i = PID[g]; n = npts[i]; l = lft[g]; r = rgt[g]
        if r - l < 2: continue
        X = frag[i, :, :n]; nf = X.shape[0]
        seg = np.empty((nf, r - l + 1), np.float32)
        for f in range(nf):
            for t in range(l, r + 1):
                if mode >= 2:
                    lo = max(0, t - 18); hi = min(n, t + 19)
                    seg[f, t - l] = max(0.0, X[f, t] - np.quantile(X[f, lo:hi], 0.2))
                else:
                    seg[f, t - l] = X[f, t]
        L = np.sqrt(np.maximum(lib[i], 0.0)); Ls = L.sum()
        coel = 0.0; shape = 0.0; shw = 0.0; ws = 0.0; cnt = 0
        for a in range(nf):
            for b in range(a + 1, nf):
                v, lg = xcorr_max(seg[a], seg[b], 3)
                coel += lg; shape += v; w = L[a] * L[b]; shw += w * v; ws += w; cnt += 1
        out[g, 0] = coel / cnt; out[g, 1] = shape / cnt; out[g, 2] = shw / ws if ws > 0 else 0.0
        area = np.zeros(nf)
        for f in range(nf): area[f] = np.sqrt(max(seg[f].sum(), 0.0))
        an = np.sqrt((area ** 2).sum()); ln = np.sqrt((L ** 2).sum())
        out[g, 3] = (area * L).sum() / (an * ln) if an > 0 and ln > 0 else 0.0
        as_ = area.sum()
        out[g, 4] = np.abs(area / as_ - L / Ls).sum() if as_ > 0 and Ls > 0 else 2.0
        s2 = np.zeros(r - l + 1, np.float32)
        for f in range(nf): s2 += seg[f]
        v, lg = xcorr_max(ms1[i, 0, l:r + 1].astype(np.float32), s2, 3)
        out[g, 5] = v
    return out

def oracle(X, right, dec, fold, kind, seed=0):
    Z = (X - np.nanmean(X, 0)) / (np.nanstd(X, 0) + 1e-9); Z = np.nan_to_num(Z)
    out = np.zeros(len(Z)); rng = np.random.default_rng(seed); torch.manual_seed(seed)
    for k in range(3):
        pos = np.where(right & (fold != k))[0]; neg = np.where(dec & (fold != k))[0]
        if kind == 'lda':
            out[fold == k] = lda(Z, pos, neg)[fold == k]
            continue
        neg = rng.choice(neg, min(len(neg), 400000), replace=False)
        Xt = torch.tensor(np.r_[Z[pos], Z[neg]], dtype=torch.float32); yt = torch.tensor(np.r_[np.ones(len(pos)), np.zeros(len(neg))], dtype=torch.float32)
        net = nn.Sequential(nn.Linear(Z.shape[1], 32), nn.GELU(), nn.Linear(32, 16), nn.GELU(), nn.Linear(16, 1))
        opt = torch.optim.AdamW(net.parameters(), lr=3e-3, weight_decay=1e-3); lf = nn.BCEWithLogitsLoss()
        for _ in range(6):
            perm = torch.randperm(len(Xt))
            for b in range(0, len(Xt), 4096):
                ix = perm[b:b + 4096]; loss = lf(net(Xt[ix]).squeeze(1), yt[ix]); opt.zero_grad(); loss.backward(); opt.step()
        with torch.no_grad():
            out[fold == k] = net(torch.tensor(Z[fold == k], dtype=torch.float32)).squeeze(1).numpy()
    return out

def tdc(score, prec, dec, pair, right, ent):
    allw = winners(np.arange(len(score)), prec, score)
    w = pd.DataFrame({'pair': pair[allw], 's': score[allw], 'dec': dec[allw], 'right': right[allw], 'ent': ent[allw]})
    w = w.sort_values(['s', 'dec'], ascending=[False, False]).drop_duplicates('pair').sort_values('s', ascending=False)
    d = w.dec.values; Tn = np.cumsum(~d); Dn = np.cumsum(d)
    q = np.minimum.accumulate(np.minimum(1, (Dn + 1) / np.maximum(Tn, 1))[::-1])[::-1]
    res = []
    for qq in (0.01, 0.05):
        a = w[(~d) & (q <= qq)]
        res.append(f'{len(a)} (conf-right {int(a.right.sum())}, FDP_c {a.ent.sum() * (1 + 1 / R_ENT) / max(1, len(a)):.3f})')
    return ' | '.join(res)

R_ENT = 148983 / 151239
cal_a = meta.cal_rt.fillna(-1e9).values.astype(np.float32)

@nb.njit(parallel=True, cache=True)
def bounds(frag, npts, rt, cal, okm, win, mode, KM):
    # re-derive apex/left/right per candidate (same logic as one())
    n_ = frag.shape[0]; A = np.full((n_, KM), -1, np.int32); Lb = np.full((n_, KM), -1, np.int32); Rb = np.full((n_, KM), -1, np.int32)
    for i in nb.prange(n_):
        if not okm[i]: continue
        n = npts[i]; X = frag[i, :, :n]; nf = X.shape[0]
        tr = np.zeros(n, np.float32)
        for f in range(nf):
            for t in range(n):
                if mode >= 2:
                    lo = max(0, t - 18); hi = min(n, t + 19)
                    tr[t] += max(0.0, X[f, t] - np.quantile(X[f, lo:hi], 0.2))
                else:
                    tr[t] += X[f, t]
        sm = np.zeros(n, np.float32)
        for t in range(n):
            a = tr[t - 1] if t > 0 else tr[t]; b = tr[t + 1] if t < n - 1 else tr[t]
            sm[t] = 0.25 * a + 0.5 * tr[t] + 0.25 * b
        cand = np.full(KM, -1); ch = np.full(KM, -1.0)
        for t in range(1, n - 1):
            if abs(rt[i, t] - cal[i]) <= win and sm[t] > 0 and sm[t] >= sm[t - 1] and sm[t] > sm[t + 1]:
                j = np.argmin(ch)
                if sm[t] > ch[j]: ch[j] = sm[t]; cand[j] = t
        k = 0
        for c_ in range(KM):
            a = cand[c_]
            if a < 0: continue
            if mode == 0:
                l = a; r = a
                while l > 0 and rt[i, a] - rt[i, l - 1] <= 22.0: l -= 1
                while r < n - 1 and rt[i, r + 1] - rt[i, a] <= 22.0: r += 1
            else:
                l = a
                while l > 0 and sm[l - 1] <= sm[l] and sm[l - 1] > 0.05 * sm[a]: l -= 1
                r = a
                while r < n - 1 and sm[r + 1] <= sm[r] and sm[r + 1] > 0.05 * sm[a]: r += 1
            if r - l < 2:
                l = max(0, a - 1); r = min(n - 1, a + 1)
            A[i, k] = a; Lb[i, k] = l; Rb[i, k] = r; k += 1
    return A, Lb, Rb

for mode, name in [(0, 'wide 45s'), (1, 'picker'), (2, 'baseline')]:
    res, ap = run(frag, ms1, lib, npts, rt, cal_a, ok, WIN, mode)
    A, Lb, Rb = bounds(frag, npts, rt, cal_a, ok, WIN, mode, KMAX)
    have = np.isfinite(res[:, :, 0]); ii, kk = np.where(have)
    assert (A[ii, kk] >= 0).all()
    PID = ii.astype(np.int64)
    Xs = res[ii, kk].astype(np.float64)
    Xx = xfeat(frag, ms1, lib, npts, rt, PID, A[ii, kk], Lb[ii, kk], Rb[ii, kk], mode).astype(np.float64)
    m = meta.iloc[ii].reset_index(drop=True); apex = ap[ii, kk]
    prec = m.precursor_id.values; dec = m.decoy.values.astype(bool); fold = m.fold.values; pair = m.pair.values; ent = m.ent.fillna(0).values
    right = (~dec) & m.dn_rt.notna().values & (np.abs(apex - m.dn_rt.values) <= 10)
    print(f'\n== {name}: groups {len(Xs)}', flush=True)
    rw = right; dw = dec
    print('   AUC xcorr feats (confirmed-right vs all decoy groups):', ' '.join(f'{n_}={auc(Xx[rw, j], Xx[dw, j]):.2f}' for j, n_ in enumerate(XN)))
    for fs, X in [('simple', Xs), ('simple+xcorr', np.c_[Xs, Xx])]:
        for kind in ('lda', 'mlp'):
            print(f'   oracle {kind} {fs:13s}: q<=0.01 / 0.05 = {tdc(oracle(X, right, dec, fold, kind), prec, dec, pair, right, ent)}', flush=True)
