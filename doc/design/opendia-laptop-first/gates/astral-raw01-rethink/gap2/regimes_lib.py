"""Gap anatomy v2 on raw 01, from the OSW-extracted XICs (same library, same raw; ~34 ppm, 1.6 s cycle, +-97 s).

Same learner for every regime (PyProphet-style semi-supervised LDA, cross-fitted on pair-hash folds, target/decoy
labels only), so differences come from the data treatment. Each regime changes one thing:
  wide     : up to 3 peaks per precursor (local maxima of the smoothed summed fragment trace within +-WIN s of the
             calibrated library RT); boundaries fixed at apex +-22 s (mimics our 45 s search groups)
  picker   : same peaks, boundaries = walk to valley / 5% of apex
  baseline : picker after per-fragment baseline subtraction (rolling 20th percentile, +-18 scans), before picking
  robust   : baseline + interference-robust features (drop the worst fragment for co-elution and library correlation)
Reported: IDs at q<=0.01/0.05, DIA-NN-confirmed with the right apex among them, entrapment FDP_c at q<=0.01, and AUC of
each feature for confirmed-right winners vs decoy winners. Also the per-precursor 'right peak among candidates' rate.
"""
import os, numpy as np, pandas as pd, sqlite3, numba as nb
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; D = B + '/astral-design-raw01-20261107'
S = os.environ['JCODE_SCRATCH_DIR'] + '/xic01'
frag = np.load(S + '/frag.npy'); lib = np.load(S + '/lib.npy'); rt = np.load(S + '/rt.npy'); ms1 = np.load(S + '/ms1.npy')
npts = np.load(S + '/npts.npy'); meta = pd.read_parquet(S + '/meta.parquet')
N, F, T = frag.shape
s = pd.read_csv(D + '/pipeline/search.features.tsv', sep='\t', usecols=['precursor_id', 'rt_apex', 'apex_rt_error_seconds'])
cal = (s.rt_apex - s.apex_rt_error_seconds).groupby(s.precursor_id).first()
meta['cal_rt'] = cal.reindex(meta.precursor_id).values
ref = pd.read_parquet(B + '/astral-diann-gap-20261108/ref_exact.parquet'); ref = ref[ref.Run.str.endswith('_01')]
meta = meta.merge(ref[['precursor_id']].assign(dn_rt=ref.RT.values * 60), on='precursor_id', how='left')
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
pk = pd.read_sql('select ID precursor_id, TRAML_ID from PRECURSOR', c)
pk['pair'] = pk.TRAML_ID.str.replace('^DECOY_', '', regex=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) ent
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c)
meta = meta.merge(pk[['precursor_id', 'pair']], on='precursor_id', how='left').merge(ent, on='precursor_id', how='left')
meta['fold'] = pd.util.hash_pandas_object(meta.pair, index=False).values % 3
WIN = 60.0; KMAX = 3; NF = 11
ok = (npts > 20) & np.isfinite(meta.cal_rt.values)
NAMES = ['coelution', 'shape_to_sum', 'lib_corr', 'apex_agree', 'n_frag_signal', 'ms1_coelution', 'ms1_isotope',
         'log_snr', 'log_width', 'log_area', 'abs_rt_err']

@nb.njit(cache=True)
def corr(a, b):
    ma = a.mean(); mb = b.mean(); va = ((a - ma) ** 2).sum(); vb = ((b - mb) ** 2).sum()
    if va <= 0 or vb <= 0: return 0.0
    return ((a - ma) * (b - mb)).sum() / np.sqrt(va * vb)

@nb.njit(cache=True)
def one(X, M, L, n, rtv, centre, win, mode, out, apx):
    nf = X.shape[0]
    Xc = np.empty((nf, n), np.float32)
    for f in range(nf):
        if mode >= 2:
            for i in range(n):
                lo = max(0, i - 18); hi = min(n, i + 19)
                Xc[f, i] = max(0.0, X[f, i] - np.quantile(X[f, lo:hi], 0.2))
        else:
            for i in range(n): Xc[f, i] = X[f, i]
    tr = np.zeros(n, np.float32)
    for i in range(n): tr[i] = Xc[:, i].sum()
    sm = np.zeros(n, np.float32)
    for i in range(n):
        a = tr[i - 1] if i > 0 else tr[i]; b = tr[i + 1] if i < n - 1 else tr[i]
        sm[i] = 0.25 * a + 0.5 * tr[i] + 0.25 * b
    # local maxima within window, top KMAX by height
    cand = np.full(KMAX, -1); ch = np.full(KMAX, -1.0)
    for i in range(1, n - 1):
        if abs(rtv[i] - centre) <= win and sm[i] > 0 and sm[i] >= sm[i - 1] and sm[i] > sm[i + 1]:
            j = np.argmin(ch)
            if sm[i] > ch[j]: ch[j] = sm[i]; cand[j] = i
    k = 0
    for c_ in range(KMAX):
        a = cand[c_]
        if a < 0: continue
        if mode == 0:
            l = a; r = a
            while l > 0 and rtv[a] - rtv[l - 1] <= 22.0: l -= 1
            while r < n - 1 and rtv[r + 1] - rtv[a] <= 22.0: r += 1
        else:
            l = a
            while l > 0 and sm[l - 1] <= sm[l] and sm[l - 1] > 0.05 * sm[a]: l -= 1
            r = a
            while r < n - 1 and sm[r + 1] <= sm[r] and sm[r + 1] > 0.05 * sm[a]: r += 1
        if r - l < 2:
            l = max(0, a - 1); r = min(n - 1, a + 1)
        seg = Xc[:, l:r + 1]; w = r - l + 1
        s2 = np.zeros(w, np.float32)
        for f in range(nf): s2 += seg[f]
        cm = np.ones((nf, nf))
        for i in range(nf):
            for j in range(i + 1, nf):
                v = corr(seg[i], seg[j]); cm[i, j] = v; cm[j, i] = v
        worst = -1
        if mode == 3:
            wv = 1e9
            for i in range(nf):
                v = cm[i].sum()
                if v < wv: wv = v; worst = i
        tot = 0.0; cnt = 0
        for i in range(nf):
            for j in range(i + 1, nf):
                if i != worst and j != worst: tot += cm[i, j]; cnt += 1
        o = out[k]
        o[0] = tot / max(cnt, 1)
        st = 0.0
        for f in range(nf): st += corr(seg[f], s2)
        o[1] = st / nf
        area = np.zeros(nf); ll = np.sqrt(np.maximum(L, 0.0))
        for f in range(nf): area[f] = np.sqrt(max(seg[f].sum(), 0.0))
        if mode == 3:
            best = -2.0
            for d in range(nf):
                msk = np.ones(nf, np.bool_); msk[d] = False
                best = max(best, corr(area[msk], ll[msk]))
            o[2] = best
        else:
            o[2] = corr(area, ll)
        ag = 0; ns = 0
        for f in range(nf):
            am = np.argmax(seg[f]) + l
            if am >= a - 1 and am <= a + 1: ag += 1
            if seg[f].max() > 0: ns += 1
        o[3] = ag / nf; o[4] = ns
        o[5] = corr(M[0, l:r + 1], s2); o[6] = corr(M[0, l:r + 1], M[1, l:r + 1])
        bgsum = 0.0; bgn = 0
        tmp = np.empty(n, np.float32)
        for i in range(n):
            if (i < l or i > r) and abs(rtv[i] - rtv[a]) <= 60: tmp[bgn] = tr[i]; bgn += 1
        noise = np.median(tmp[:bgn]) if bgn > 3 else 0.0
        o[7] = np.log1p(tr[a] / (noise + 1e-3 * tr[a] + 1e-9))
        o[8] = np.log(max(rtv[r] - rtv[l], 0.5)); o[9] = np.log1p(s2.sum()); o[10] = abs(rtv[a] - centre)
        apx[k] = rtv[a]
        k += 1
    return k

@nb.njit(parallel=True, cache=True)
def run(frag, ms1, lib, npts, rt, cal, okm, win, mode):
    n_ = frag.shape[0]
    res = np.full((n_, KMAX, NF), np.nan, np.float32); ap = np.full((n_, KMAX), np.nan, np.float32)
    for i in nb.prange(n_):
        if not okm[i]: continue
        n = npts[i]
        one(np.ascontiguousarray(frag[i, :, :n]), np.ascontiguousarray(ms1[i, :, :n]), lib[i], n, rt[i, :n], cal[i], win, mode, res[i], ap[i])
    return res, ap

def lda(Z, pos, neg):
    a, b = Z[pos], Z[neg]
    Sw = (np.cov(a.T) * (len(a) - 1) + np.cov(b.T) * (len(b) - 1)) / (len(a) + len(b) - 2) + 1e-3 * np.eye(Z.shape[1])
    return Z @ np.linalg.solve(Sw, a.mean(0) - b.mean(0))

def winners(idx, prec, score):
    o = idx[np.lexsort((idx, -score[idx]))]
    _, first = np.unique(prec[o], return_index=True)
    return o[first]

def select(idx, prec, dec, score, fdr):
    w = winners(idx, prec, score); w = w[np.argsort(-score[w], kind='stable')]
    d = dec[w]; Tn = np.cumsum(~d); Dn = np.cumsum(d)
    q = np.minimum.accumulate(np.minimum(1, (Dn + 1) / np.maximum(Tn, 1))[::-1])[::-1]
    return w[(~d) & (q <= fdr)], w[d]

def auc(a, b):
    r = pd.Series(np.r_[a, b]).rank().values
    return (r[:len(a)].sum() - len(a) * (len(a) + 1) / 2) / (len(a) * len(b))

