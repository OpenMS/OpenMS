"""PRECHECK_MF: fragment-level, mobility-aware matched filter for faint right peaks (diagnostic, 20 Nov 2026).
Sample: faint2 (6,000 DIA-NN targets + 30,000 decoys, r~1.56 D), raw points from maps.xipm (+-15 ppm, +-0.08 1/K0).
Baseline features: the CURRENT default (10 ppm trace, -im_tolerance_scale 0.6), faintcpp D_cur_s0 features.tsv.

Per candidate peak, from the MS2 points of its 6 fragments:
  m*        consensus mobility apex: ppm-weighted 1/K0 histogram (0.004 bins, +-0.06) over |rt - apex| <= 4 s
  weights   w_ppm = N(ppm - offset; sigma_signal), w_im = N(im - m*; 0.008), points with |im - m*| <= 0.024
  traces    Xw[scan, fragment] = sum intensity * w_ppm * w_im
  mf_z_med, mf_z_trim, mf_n_z2, mf_n_z3   per-fragment robust apex z (apex +-1 scan vs 3-scan sums outside the peak,
                                          median / MAD with a sqrt floor); trim = mean without the largest fragment
  mf_loo_med, mf_loo_n                    leave-one-out shape: corr(fragment, sum of others) over apex +-3 scans
  mf_cos, mf_dom                          apex cosine with sqrt library; log max observed/library share
  mf_ppm                                  intensity-weighted |ppm - offset| at apex, |im - m*| <= 0.016 (not ppm-weighted)
  mf_llr                                  sum_k libw_k * log1p(max(z_k, 0))
  mf_gain                                 log weighted apex consensus - log unweighted pooled apex consensus

Arms (GBM, cross-fit on the tool's folds, 2 seeds):
  ceiling     DIA-NN right windows (|rt - DIA-NN| <= 6 s) vs decoys, as PRECHECK_IMFEAT
  tool        label-free proxy: positives = target windows the tool accepts (q <= 0.01), negatives = decoys
Metric: TPR of DIA-NN right windows at decoy FPR 1e-3 / 3e-3 / 1e-2, all and the 'missed' subset (tool best q > 0.01).

Frozen gate (written before any result): proceed to C++ only if 36+mf beats 36 by
  >= +0.030 TPR at decoy FPR 3e-3 in the ceiling arm on both seeds, AND
  >= +0.020 TPR at decoy FPR 3e-3 in the tool arm on the missed subset on both seeds.
core5 had +0.022 ceiling and recovered no IDs, so the bar sits above it and needs a label-free signal.

Amendment 1 (20 Nov, job 39074453 cancelled at 1:15 before any feature or result existed): add the tensor-feasible
variant mt_* = the same mf_* statistics computed from what the workflow's trace tensor already holds: 10 ppm window,
32 RT bins over +-rt_window (2 x 35.2 s / 32 = 2.2 s), 8 IM bins over +-0.6 x 0.0524 around the calibrated mobility, no
ppm weighting, IM weight = gaussian on bin centre around the tensor's own IM apex bin. The same gate applies to
36+mt. Implementation choice: mt if it passes (no extraction change), else mf if it passes (needs a new raw grid),
else neither."""
import os, sys, sqlite3, zlib, json, numpy as np, pandas as pd, pyarrow as pa, pyarrow.compute, pyarrow.parquet as pq
from multiprocessing import Pool

O = os.environ['FAINT2']; CUR = os.environ['CUR_FEATS']; CURP = os.environ['CUR_PREC']; CAL = os.environ['CAL']
NPROC = int(os.environ.get('NPROC', '16')); OUT = os.environ.get('OUT', f'{O}/mf.parquet')
cal = json.load(open(CAL))
IM_A, IM_B = cal['im']['intercept_vssc'], cal['im']['slope']
PPM_OFF, PPM_SIG = cal['ms2_mz']['offset_ppm'], cal['ms2_mz']['signal_sigmas_ppm'][0]
IM_SIG, IM_CUT, CORE_S = 0.008, 0.024, 4.0
RTW = cal['rt'].get('tolerance_half_width_s', None)
TEN_PPM, TEN_RT_BINS, TEN_IM_BINS = 10.0, 32, 8
TEN_IM_HW = 0.6 * cal['im']['tolerance_half_width_vssc']
MF = ['mf_z_med', 'mf_z_trim', 'mf_n_z2', 'mf_n_z3', 'mf_loo_med', 'mf_loo_n', 'mf_cos', 'mf_dom', 'mf_ppm', 'mf_llr', 'mf_gain']
MT = ['mt' + c[2:] for c in MF if c not in ('mf_ppm', 'mf_gain')]
print(f'ppm offset {PPM_OFF:.3f} sigma {PPM_SIG:.3f}, im {IM_A:.5f}+{IM_B:.5f}x', flush=True)

def dec(b, code):
    import pyopenms as oms
    if not b: return np.zeros(0)
    out = []; c = oms.NumpressConfig()
    c.np_compression = {5: oms.MSNumpressCoder.NumpressCompression.LINEAR, 6: oms.MSNumpressCoder.NumpressCompression.SLOF}[code]
    oms.MSNumpressCoder().decodeNPRaw(zlib.decompress(b), out, c)
    return np.asarray(out)

def cosw(v, w):
    n = np.linalg.norm(v); return float(v @ w / n) if n > 0 else 0.0

def robust_z(x, ia, a, b):
    """Apex (+-1 scan) sum against 3-scan sums outside [a, b): (apex - median) / (1.4826 MAD + sqrt(median) + 1)."""
    s = np.convolve(x, np.ones(3), 'same'); apex = s[ia]
    bg = np.r_[s[:max(a - 1, 0)], s[b + 1:]]
    if len(bg) < 3: bg = np.r_[s[:max(ia - 3, 0)], s[ia + 4:]]
    if len(bg) == 0: return 0.0
    med = float(np.median(bg)); mad = float(np.median(np.abs(bg - med)))
    return (apex - med) / (1.4826 * mad + np.sqrt(max(med, 0.0)) + 1.0)

def mf_stats(X, ia, a, b, lw, lshare, prefix, dom_win=1, loo_win=3):
    """Shared matched-filter statistics on a [scans, fragments] trace X at apex index ia, peak [a, b)."""
    nf = X.shape[1]; r = {}
    z = np.array([robust_z(X[:, k], ia, a, b) for k in range(nf)])
    r[prefix + '_z_med'] = float(np.median(z)); r[prefix + '_z_trim'] = float(np.sort(z)[:-1].mean())
    r[prefix + '_n_z2'] = int((z >= 2).sum()); r[prefix + '_n_z3'] = int((z >= 3).sum())
    lo, hi = max(ia - loo_win, 0), min(ia + loo_win + 1, X.shape[0]); W = X[lo:hi]; loo = []
    for k in range(nf):
        o = W.sum(1) - W[:, k]
        loo.append(float(np.corrcoef(W[:, k], o)[0, 1]) if W[:, k].std() > 0 and o.std() > 0 else 0.0)
    loo = np.array(loo); r[prefix + '_loo_med'] = float(np.median(loo)); r[prefix + '_loo_n'] = int((loo > 0.5).sum())
    v = X[max(ia - dom_win, 0):ia + dom_win + 1].sum(0); r[prefix + '_cos'] = cosw(v, lw)
    r[prefix + '_dom'] = float(np.log((v / v.sum()).max() / lshare[np.argmax(v)])) if v.sum() > 0 and lshare[np.argmax(v)] > 0 else 0.0
    r[prefix + '_llr'] = float((lw * np.log1p(np.maximum(z, 0))).sum())
    return r, v

G = {}  # filled before the pool forks: T, groups, PG, lib, prec

def work(pids):
    T, groups, PG, lib, prec = G['T'], G['groups'], G['PG'], G['lib'], G['prec']
    rows = []
    for pid in pids:
        idx = groups.get(pid)
        if idx is None or pid not in PG: continue
        sub = T.take(idx).to_pydict(); frags = []
        for k in range(len(idx)):
            rt = dec(sub['RT_DATA'][k], sub['RT_COMPRESSION'][k]); im = dec(sub['MOBILITY_DATA'][k], sub['MOBILITY_COMPRESSION'][k])
            it = dec(sub['INTENSITY_DATA'][k], sub['INTENSITY_COMPRESSION'][k]); mz = dec(sub['MZ_DATA'][k], sub['MZ_COMPRESSION'][k])
            tmz = sub['TARGET_MZ'][k]; ppm = (mz - tmz) / tmz * 1e6 if len(mz) == len(it) else np.zeros(len(it))
            frags.append((rt, im, it, ppm, lib.get(sub['TRANSITION_ID'][k], 0.0)))
        nf = len(frags)
        if nf < 2: continue
        lw = np.sqrt(np.maximum([f[4] for f in frags], 0)); lw = lw / (np.linalg.norm(lw) + 1e-12)
        lshare = lw ** 2 / max((lw ** 2).sum(), 1e-12)
        mcal = IM_A + IM_B * prec.get(pid, np.nan)
        if not np.isfinite(mcal): continue
        # Fragments of one frame carry slightly different RT values (float jitter), so snap to a frame grid: a new scan
        # starts where consecutive sorted RTs differ by more than 0.3 s (frames here are ~0.95 s apart).
        allrt = np.sort(np.concatenate([f[0] for f in frags]))
        if len(allrt) < 5: continue
        starts = np.r_[0, np.where(np.diff(allrt) > 0.3)[0] + 1]
        scans = np.array([allrt[s:e].mean() for s, e in zip(starts, np.r_[starts[1:], len(allrt)])])
        if len(scans) < 5: continue
        mids = 0.5 * (scans[:-1] + scans[1:])
        si = [np.searchsorted(mids, f[0]) for f in frags]
        wppm = [np.exp(-0.5 * ((f[3] - PPM_OFF) / PPM_SIG) ** 2) for f in frags]
        for _, pk in PG[pid].iterrows():
            ia = int(np.argmin(np.abs(scans - pk.rt_apex)))
            a, b = np.searchsorted(scans, pk.rt_start - 1e-6), np.searchsorted(scans, pk.rt_end + 1e-6)
            r = dict(precursor_id=pid, rta=pk.rta)
            # consensus mobility apex, ppm weighted, core RT window
            edges = np.arange(mcal - 0.06, mcal + 0.0601, 0.004); hst = np.zeros(len(edges) - 1)
            for k, (rt, im, it, ppm, _) in enumerate(frags):
                m = np.abs(rt - pk.rt_apex) <= CORE_S
                if m.any(): hst += np.histogram(im[m], edges, weights=(it * wppm[k])[m])[0]
            if hst.sum() > 0:
                hs = np.convolve(hst, [0.25, 0.5, 0.25], 'same'); ms = 0.5 * (edges[:-1] + edges[1:])[np.argmax(hs)]
            else:
                ms = mcal
            Xw = np.zeros((len(scans), nf)); Xu = np.zeros((len(scans), nf))
            ppm_num = ppm_den = 0.0
            for k, (rt, im, it, ppm, _) in enumerate(frags):
                d = im - ms; inband = np.abs(d) <= IM_CUT
                wim = np.exp(-0.5 * (d / IM_SIG) ** 2)
                np.add.at(Xw[:, k], si[k][inband], (it * wppm[k] * wim)[inband])
                pooled = np.abs(im - mcal) <= 0.0314
                np.add.at(Xu[:, k], si[k][pooled], it[pooled])
                ap = (np.abs(si[k] - ia) <= 1) & (np.abs(d) <= 0.016)
                ppm_num += float((it[ap] * np.abs(ppm[ap] - PPM_OFF)).sum()); ppm_den += float(it[ap].sum())
            z_r, v = mf_stats(Xw, ia, a, b, lw, lshare, 'mf'); r.update(z_r)
            r['mf_ppm'] = ppm_num / ppm_den if ppm_den > 0 else 15.0
            cu = Xu[max(ia - 1, 0):ia + 2].sum()
            r['mf_gain'] = float(np.log1p(v.sum()) - np.log1p(cu))
            # Tensor-feasible variant: the workflow's own 32 RT x 8 IM binning around the apex at 10 ppm.
            tw = 2.0 * RTW / TEN_RT_BINS; Xt = np.zeros((TEN_RT_BINS, TEN_IM_BINS, nf))
            for k, (rt, im, it, ppm, _) in enumerate(frags):
                rb = np.floor((rt - pk.rt_apex + RTW) / tw).astype(int)
                ib = np.floor((im - mcal + TEN_IM_HW) / (2 * TEN_IM_HW) * TEN_IM_BINS).astype(int)
                m = (np.abs(ppm) <= TEN_PPM) & (rb >= 0) & (rb < TEN_RT_BINS) & (ib >= 0) & (ib < TEN_IM_BINS)
                np.add.at(Xt[:, :, k], (rb[m], ib[m]), it[m])
            ca = TEN_RT_BINS // 2; core = Xt[max(ca - 2, 0):ca + 3].sum((0, 2))
            jb = int(np.argmax(np.convolve(core, [0.25, 0.5, 0.25], 'same'))) if core.sum() > 0 else TEN_IM_BINS // 2
            wib = np.exp(-0.5 * ((np.arange(TEN_IM_BINS) - jb) / (IM_SIG / (2 * TEN_IM_HW / TEN_IM_BINS))) ** 2)
            XT = (Xt * wib[None, :, None]).sum(1)
            a_t = int(np.clip(np.floor((pk.rt_start - pk.rt_apex + RTW) / tw), 0, TEN_RT_BINS - 1))
            b_t = int(np.clip(np.floor((pk.rt_end - pk.rt_apex + RTW) / tw) + 1, 1, TEN_RT_BINS))
            t_r, _ = mf_stats(XT, ca, a_t, b_t, lw, lshare, 'mt', dom_win=0, loo_win=2); r.update(t_r)
            rows.append(r)
    return rows

if __name__ == '__main__' and os.environ.get('STEP', 'extract') == 'extract':
    P = pd.read_csv(f'{O}/peaks.tsv', sep='\t'); P['rta'] = P.rt_apex.round(3)
    PG = {p: g for p, g in P.groupby('precursor_id')}
    c = sqlite3.connect(f'{O}/mini.pqp')
    lib = pd.read_sql('select ID tid, LIBRARY_INTENSITY li from TRANSITION', c).set_index('tid').li.to_dict()
    prec = pd.read_sql('select ID pid, LIBRARY_DRIFT_TIME lim from PRECURSOR', c).set_index('pid').lim.to_dict()
    T = pq.read_table(f'{O}/maps.xipm', columns=['MS_LEVEL', 'PRECURSOR_ID', 'TRANSITION_ID', 'TARGET_MZ', 'MZ_DATA', 'RT_DATA',
                                                 'MOBILITY_DATA', 'INTENSITY_DATA', 'MZ_COMPRESSION', 'RT_COMPRESSION',
                                                 'MOBILITY_COMPRESSION', 'INTENSITY_COMPRESSION'])
    T = T.filter(pa.compute.equal(T.column('MS_LEVEL'), 2))
    groups = pd.Series(np.arange(T.num_rows)).groupby(T.column('PRECURSOR_ID').to_numpy()).apply(list).to_dict()
    print(f'MS2 rows {T.num_rows}, precursors {len(groups)}, peaks {len(P)}', flush=True)
    pids = list(groups.keys()); chunks = [pids[i::64] for i in range(64)]
    G.update(T=T, groups=groups, PG=PG, lib=lib, prec=prec)
    with Pool(NPROC) as pool:  # fork: workers inherit G
        res = pool.map(work, chunks)
    D = pd.DataFrame([r for rr in res for r in rr]); D.to_parquet(OUT)
    print(f'windows {len(D)}', flush=True); print(D[MF + MT].describe().T.round(3).to_string(), flush=True)

if __name__ == '__main__' and os.environ.get('STEP', 'extract') == 'eval':
    from sklearn.ensemble import HistGradientBoostingClassifier
    S = pd.read_csv(f'{O}/sample.tsv', sep='\t').set_index('pid')
    keep = set(S.index)
    parts = [ch[ch.precursor_id.isin(keep)] for ch in pd.read_csv(CUR, sep='\t', chunksize=500000)]
    F = pd.concat(parts); F['rta'] = F.rt_apex.round(3)
    D = pd.read_parquet(OUT).drop_duplicates(['precursor_id', 'rta'])
    C = pd.read_parquet(f'{O}/imfeat.parquet')[['precursor_id', 'rta', 'c_offset', 'c_frag_sd', 'c_cos', 'c_support', 'c_contrast']]
    F = F.merge(D, on=['precursor_id', 'rta']).merge(C.drop_duplicates(['precursor_id', 'rta']), on=['precursor_id', 'rta'], how='left')
    CORE = ['c_offset', 'c_frag_sd', 'c_cos', 'c_support', 'c_contrast']
    F[CORE] = F[CORE].fillna(0.0)
    meta = ['candidate_row', 'precursor_id', 'is_decoy', 'fold', 'rt_apex', 'score', 'q_value', 'rta']
    X36 = [c for c in F.columns if c not in meta + MF + MT + CORE]
    Pq = pd.read_csv(CURP, sep='\t', usecols=['precursor_id', 'is_decoy', 'q_value'])
    best_q = Pq[Pq.is_decoy == 0].groupby('precursor_id').q_value.min()
    dec_ = F.is_decoy.values == 1
    pos = (F.precursor_id.map(S.cls).values == 'dnn') & (np.abs(F.rt_apex.values - F.precursor_id.map(S.rt_s).values) <= 6)
    missed_ids = set(S.index[(S.cls == 'dnn')]) - set(best_q.index[best_q <= 0.01])
    toolpos = (~dec_) & (F.q_value.values <= 0.01)
    fold = F.fold.values; dec_ids = S.index[S.cls == 'decoy']
    print(f'joined {len(F)}, X36 {len(X36)}, right windows {pos.sum()}, missed DIA-NN ids {len(missed_ids)}, tool positives {toolpos.sum()}', flush=True)
    def gbm(M, seed, y_pos):
        s = np.zeros(len(F))
        for f in np.unique(fold):
            tr = (fold != f) & (y_pos | dec_); y = y_pos[tr].astype(int); w = np.where(y == 1, (y == 0).sum() / max(y.sum(), 1), 1.0)
            s[fold == f] = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=31, random_state=seed
                                                          ).fit(M[tr], y, sample_weight=w).decision_function(M[fold == f])
        return s
    def tpr(s):
        F['_r'] = s; best = F.groupby('precursor_id')._r.max(); rb = F[pos].groupby('precursor_id')._r.max()
        neg = best.reindex(dec_ids).dropna().values; thr = [np.quantile(neg, 1 - f) for f in (1e-3, 3e-3, 1e-2)]
        rm = rb[rb.index.isin(missed_ids)]
        return [(rb > t).mean() for t in thr], [(rm > t).mean() for t in thr]
    arms = (('36', X36), ('36+core5', X36 + CORE), ('36+mf', X36 + MF), ('36+mt', X36 + MT), ('36+mf+core5', X36 + MF + CORE))
    res = {}
    for label, ypos in (('ceiling', pos), ('tool', toolpos)):
        for name, cols in arms:
            M = F[cols].values.astype(np.float32); res[(label, name)] = r = [tpr(gbm(M, sd, ypos)) for sd in (0, 1)]
            print(f'{label:7s} {name:12s} ' + ' | '.join(f's{sd} all {a[0]:.3f}/{a[1]:.3f}/{a[2]:.3f} missed {m[0]:.3f}/{m[1]:.3f}/{m[2]:.3f}'
                                                     for sd, (a, m) in enumerate(r)), flush=True)
    print('deltas at decoy FPR 3e-3 against 36 (seed 0 / seed 1):')
    for label in ('ceiling', 'tool'):
        for name, _ in arms[1:]:
            da = [res[(label, name)][sd][0][1] - res[(label, '36')][sd][0][1] for sd in (0, 1)]
            dm = [res[(label, name)][sd][1][1] - res[(label, '36')][sd][1][1] for sd in (0, 1)]
            print(f'  {label:7s} {name:12s} all {da[0]:+.3f} / {da[1]:+.3f}   missed {dm[0]:+.3f} / {dm[1]:+.3f}')
    for arm in ('36+mf', '36+mt'):
        c1 = all(res[('ceiling', arm)][sd][0][1] - res[('ceiling', '36')][sd][0][1] >= 0.030 for sd in (0, 1))
        c2 = all(res[('tool', arm)][sd][1][1] - res[('tool', '36')][sd][1][1] >= 0.020 for sd in (0, 1))
        print(f'GATE {arm}: ceiling>=+0.030: {c1}; tool missed>=+0.020: {c2}; PROCEED: {c1 and c2}')
    # single-feature view: decoy vs DIA-NN right-window AUC, and on missed right windows
    from sklearn.metrics import roc_auc_score
    rmask = pos & F.precursor_id.isin(missed_ids).values
    for c in MF + MT:
        y = np.r_[np.ones(rmask.sum()), np.zeros(dec_.sum())]; x = np.r_[F.loc[rmask, c].values, F.loc[dec_, c].values]
        print(f'  AUC missed-right vs decoy {c:11s} {roc_auc_score(y, x):.3f}')
