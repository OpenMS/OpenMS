"""OSW-like scoring on the shared XICs, teacher-free (raw 01, design diagnostic).

Question: on the SAME extracted traces OSW scored, does our own implementation of an OSW-style candidate generator,
subscore set and PyProphet-style semi-supervised LDA reproduce OSW+Percolator's ~1151 at entrapment FDP_c 1%? If yes,
the information-to-score path is understood and we can ablate which subscore groups carry it. Own code, not OSW's.
Labels: target/decoy only (DIA-NN used only for the 'confirmed-right' diagnostic column; entrapment only for FDP_c).

Candidates: per precursor, smoothed (5-point Gaussian) summed fragment trace over the whole XIC (+-97 s), local maxima,
top-5 by height, boundaries = walk to valley or 5% of apex; overlapping candidates are dropped (keep higher).
Subscores (OSW names in brackets, own implementation):
  xcorr_coel (var_xcorr_coelution), xcorr_shape (var_xcorr_shape), xcorr_coel_w / xcorr_shape_w (library-weighted),
  lib_corr (var_library_corr), lib_dotprod, lib_manhattan, lib_rmsd, lib_sangle, intensity (area / total XIC signal),
  log_sn (mean log apex S/N vs median trace), n_sig, log_width, norm_rt (|apex - calibrated library RT|),
  ms1_xcorr_coel, ms1_xcorr_shape, ms1_iso_corr, ms1_log_sn, frag_min_share, frag_corr_min.
"""
import os, time, sqlite3, numpy as np, pandas as pd, numba as nb
RUN = os.environ.get('OPENDIA_RUN', '01')
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; D = B + '/astral-design-raw01-20261107'
SEARCH = D + '/pipeline/search.features.tsv' if RUN == '01' else B + '/astral-entrapment-raw02-20261108/search.features.tsv'
S = os.environ['JCODE_SCRATCH_DIR'] + '/xic' + RUN
frag = np.load(S + '/frag.npy'); lib = np.load(S + '/lib.npy'); rt = np.load(S + '/rt.npy'); ms1 = np.load(S + '/ms1.npy')
npts = np.load(S + '/npts.npy'); meta = pd.read_parquet(S + '/meta.parquet')
N = len(meta)
s = pd.read_csv(SEARCH, sep='\t', usecols=['precursor_id', 'rt_apex', 'apex_rt_error_seconds'])
cal = (s.rt_apex - s.apex_rt_error_seconds).groupby(s.precursor_id).first()
meta['cal_rt'] = cal.reindex(meta.precursor_id).values
ref = pd.read_parquet(B + '/astral-diann-gap-20261108/ref_exact.parquet'); ref = ref[ref.Run.str.endswith('_' + RUN)]
meta = meta.merge(ref[['precursor_id']].assign(dn_rt=ref.RT.values * 60), on='precursor_id', how='left')
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
pk = pd.read_sql('select ID precursor_id, TRAML_ID from PRECURSOR', c); pk['pair'] = pk.TRAML_ID.str.replace('^DECOY_', '', regex=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) ent
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c)
meta = meta.merge(pk[['precursor_id', 'pair']], on='precursor_id', how='left').merge(ent, on='precursor_id', how='left')
meta['fold'] = pd.util.hash_pandas_object(meta.pair, index=False).values % 3
ok = (npts > 20) & np.isfinite(meta.cal_rt.values)
K = 5
NAMES = ['xcorr_coel', 'xcorr_shape', 'xcorr_coel_w', 'xcorr_shape_w', 'lib_corr', 'lib_dotprod', 'lib_manhattan',
         'lib_rmsd', 'lib_sangle', 'intensity', 'log_sn', 'n_sig', 'log_width', 'norm_rt', 'ms1_xcorr_coel',
         'ms1_xcorr_shape', 'ms1_iso_corr', 'ms1_log_sn', 'frag_min_share', 'frag_corr_min']
NF = len(NAMES)

@nb.njit(cache=True)
def xc(a, b, maxlag):
    n = len(a); ma = a.mean(); mb = b.mean(); sa = np.sqrt(((a - ma) ** 2).sum()); sb = np.sqrt(((b - mb) ** 2).sum())
    if sa <= 0 or sb <= 0: return 0.0, float(maxlag)
    best = -2.0; bl = 0
    for lag in range(-maxlag, maxlag + 1):
        t = 0.0
        for i in range(n):
            j = i + lag
            if 0 <= j < n: t += (a[i] - ma) * (b[j] - mb)
        v = t / (sa * sb)
        if v > best: best = v; bl = lag
    return best, float(abs(bl))

@nb.njit(cache=True)
def pearson(a, b):
    ma = a.mean(); mb = b.mean(); va = ((a - ma) ** 2).sum(); vb = ((b - mb) ** 2).sum()
    if va <= 0 or vb <= 0: return 0.0
    return ((a - ma) * (b - mb)).sum() / np.sqrt(va * vb)

@nb.njit(cache=True)
def score_one(X, M, L, rtv, centre, out, apx, given):
    """given[0] < 0: generate candidates from local maxima; otherwise use the apex scan indices in given (>= 0)."""
    nf, n = X.shape
    G = np.array([0.0625, 0.25, 0.375, 0.25, 0.0625])
    tr = np.zeros(n)
    for f in range(nf): tr += X[f]
    sm = np.zeros(n)
    for i in range(n):
        v = 0.0
        for d in range(-2, 3):
            j = min(max(i + d, 0), n - 1); v += G[d + 2] * tr[j]
        sm[i] = v
    # candidate maxima
    cand = np.full(K, -1); ch = np.full(K, -1.0)
    for i in range(1, n - 1):
        if sm[i] > 0 and sm[i] >= sm[i - 1] and sm[i] > sm[i + 1]:
            j = np.argmin(ch)
            if sm[i] > ch[j]: ch[j] = sm[i]; cand[j] = i
    order = np.argsort(-ch)
    if given[0] >= 0:
        cand = given.copy(); order = np.arange(K)
    used = np.zeros(n, np.bool_)
    noise = np.zeros(nf)
    for f in range(nf): noise[f] = max(np.median(X[f]), 1e-3 * max(X[f].max(), 1.0))
    mnoise = max(np.median(M[0]), 1e-3 * max(M[0].max(), 1.0))
    total = tr.sum() + 1e-9
    Ls = np.sqrt(np.maximum(L, 0.0)); Lw = Ls / (Ls.sum() + 1e-12)
    k = 0
    for oi in range(K):
        a = cand[order[oi]]
        if a < 0 or (given[0] < 0 and used[a]): continue
        l = a
        while l > 0 and sm[l - 1] <= sm[l] and sm[l - 1] > 0.05 * sm[a]: l -= 1
        r = a
        while r < n - 1 and sm[r + 1] <= sm[r] and sm[r + 1] > 0.05 * sm[a]: r += 1
        if r - l < 2: l = max(0, a - 1); r = min(n - 1, a + 1)
        for t in range(l, r + 1): used[t] = True
        seg = X[:, l:r + 1]; w = r - l + 1
        o = out[k]
        sc = 0.0; ss = 0.0; scw = 0.0; ssw = 0.0; ww = 0.0; cnt = 0; lagsum = 0.0
        cmin = 2.0
        for i in range(nf):
            for j in range(i + 1, nf):
                v, lg = xc(seg[i], seg[j], 3)
                sc += lg; ss += v; wt = Lw[i] * Lw[j]; scw += wt * lg; ssw += wt * v; ww += wt; cnt += 1
                if v < cmin: cmin = v
        o[0] = sc / cnt; o[1] = ss / cnt
        o[2] = scw / ww if ww > 0 else o[0]; o[3] = ssw / ww if ww > 0 else o[1]
        area = np.zeros(nf)
        for f in range(nf): area[f] = max(seg[f].sum(), 0.0)
        sa = np.sqrt(area); at = area.sum() + 1e-12
        o[4] = pearson(sa, Ls)
        na = np.sqrt((sa ** 2).sum()); nl = np.sqrt((Ls ** 2).sum())
        dp = (sa * Ls).sum() / (na * nl) if na > 0 and nl > 0 else 0.0
        o[5] = dp
        an = sa / (sa.sum() + 1e-12)
        o[6] = np.abs(an - Lw).sum(); o[7] = np.sqrt(((an - Lw) ** 2).mean())
        o[8] = 1.0 - 2.0 * np.arccos(min(1.0, max(-1.0, dp))) / np.pi
        o[9] = at / total
        lsn = 0.0
        for f in range(nf): lsn += np.log1p(X[f, a] / noise[f])
        o[10] = lsn / nf
        ns = 0
        for f in range(nf):
            if seg[f].max() > 3 * noise[f]: ns += 1
        o[11] = ns
        o[12] = np.log(max(rtv[r] - rtv[l], 0.5)); o[13] = abs(rtv[a] - centre)
        s2 = np.zeros(w)
        for f in range(nf): s2 += seg[f]
        m0 = M[0, l:r + 1].astype(np.float64)
        v, lg = xc(m0, s2, 3); o[14] = lg; o[15] = v
        o[16] = pearson(m0, M[1, l:r + 1].astype(np.float64))
        o[17] = np.log1p(M[0, a] / mnoise)
        o[18] = area.min() / at
        o[19] = cmin
        apx[k] = rtv[a]
        k += 1
    return k

@nb.njit(parallel=True, cache=True)
def run(frag, ms1, lib, npts, rt, cal, okm, given):
    """given: [N, K] apex scan indices (-1 = none); given[i, 0] < 0 means generate candidates from maxima."""
    n_ = frag.shape[0]
    res = np.full((n_, K, NF), np.nan, np.float32); ap = np.full((n_, K), np.nan, np.float32)
    for i in nb.prange(n_):
        if not okm[i]: continue
        n = npts[i]
        score_one(frag[i, :, :n].astype(np.float64), ms1[i, :, :n], lib[i].astype(np.float64), rt[i, :n].astype(np.float64), cal[i], res[i], ap[i], given[i])
    return res, ap

if __name__ == '__main__':
    t0 = time.time()
    res, ap = run(frag, ms1, lib, npts, rt, meta.cal_rt.fillna(-1e9).values, ok, np.full((N, K), -1, np.int64))
    ii, kk = np.where(np.isfinite(res[:, :, 0]))
    X = res[ii, kk].astype(np.float64); apex = ap[ii, kk]
    m = meta.iloc[ii].reset_index(drop=True)
    prec = m.precursor_id.values; dec = m.decoy.values.astype(bool); fold = m.fold.values; pair = m.pair.values
    entr = m.ent.fillna(0).values.astype(bool)
    right = (~dec) & m.dn_rt.notna().values & (np.abs(apex - m.dn_rt.values) <= 10)
    cov = pd.Series(right).groupby(prec).any().sum()
    print(f'candidates {len(X)} in {time.time() - t0:.0f} s | confirmed precursors with the right peak among candidates {cov} of {int(meta.dn_rt.notna().sum())}', flush=True)
    np.save(S + '/oswlike_X.npy', X); m.assign(apex=apex, right=right).to_parquet(S + '/oswlike_meta.parquet')
