"""PRECHECK_MF2: which raw-point ingredient carries the extra matched-filter gain, and can the extraction pass keep it?
(diagnostic, 20 Nov 2026). Same faint2 sample as PRECHECK_MF (6,000 DIA-NN targets, 30,000 decoys, r~1.56 D).

Baseline B45: the NEW default (36 + the 9 C++ mt_ features), from mtcheck-20261120/runs/D_mt_s0 (byte-identical to the
new default, OpenMS 01c5891ce0). 'Missed' = DIA-NN sample IDs the new default does not accept at q <= 0.01.

Variants, each the same 9 matched-filter statistics (precheck_mf.mf_stats) on a different trace:
  diagnostic (raw points, native ~0.95 s frame grid):
    mf    ppm weight x Gaussian (sd 0.008) around the window's own consensus mobility apex (search +-0.06)
    p     ppm weight only, pooled +-0.0314 around the calibrated mobility
    m     own-mobility Gaussian only, hard |ppm| <= 10
  implementable in the extraction pass (the workflow's 32 RT bins of 0.94 s over +-15 s, hard |ppm| <= 10):
    tp    8 IM bins over +-0.03 (today's tensor) with ppm-weighted accumulation, Gaussian on bins around own apex
    t16   16 IM bins over +-0.06 (0.0075 bins), no ppm weight, Gaussian around own apex      (2x tensor memory)
    t16p  t16 with ppm-weighted accumulation                                                 (2x tensor memory)
ppm weight = N(ppm - run offset; run signal sigma) from the calibration.

Frozen gate (written before any result): an implementable variant V in {tp, t16, t16p} proceeds to C++ only if B45+V
beats B45 on both seeds by
  >= +0.030 TPR at decoy FPR 3e-3 in the label-free tool arm on the missed subset, AND
  >= +0.015 TPR at decoy FPR 3e-3 in the ceiling arm over all DIA-NN right windows.
If several pass, take the cheapest (tp, then t16, then t16p) unless a dearer one beats it by >= +0.015 tool-missed on
both seeds. mf, p and m only explain where the gain comes from; they cannot pass the gate."""
import os, sys, sqlite3, json, importlib.util, numpy as np, pandas as pd, pyarrow as pa, pyarrow.compute, pyarrow.parquet as pq
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('mf', os.path.join(HERE, 'precheck_mf.py')); mf = importlib.util.module_from_spec(spec)
sys.modules['mf'] = mf; spec.loader.exec_module(mf)
O = os.environ['FAINT2']; OUT = os.environ.get('OUT2', f'{O}/mf2.parquet')
RTW, PPM_HARD = 15.0, 10.0
STATS = ['z_med', 'z_trim', 'n_z2', 'n_z3', 'loo_med', 'loo_n', 'cos', 'dom', 'llr']
VARS = ['mf', 'p', 'm', 'tp', 't16', 't16p']
COLS = {v: [f'{v}_{s}' for s in STATS] for v in VARS}
G = {}

def own_apex(frags, apex_rt, mcal, weights):
    edges = np.arange(mcal - 0.06, mcal + 0.0601, 0.004); h = np.zeros(len(edges) - 1)
    for k, (rt, im, it, ppm, _) in enumerate(frags):
        m = (np.abs(rt - apex_rt) <= mf.CORE_S) & (weights[k] > 0)
        if m.any(): h += np.histogram(im[m], edges, weights=(it * weights[k])[m])[0]
    if h.sum() <= 0: return mcal
    return 0.5 * (edges[:-1] + edges[1:])[np.argmax(np.convolve(h, [0.25, 0.5, 0.25], 'same'))]

def frame_trace(frags, si, nscan, wts):
    X = np.zeros((nscan, len(frags)))
    for k, (rt, im, it, ppm, _) in enumerate(frags):
        m = wts[k] > 0
        np.add.at(X[:, k], si[k][m], (it * wts[k])[m])
    return X

def tensor_trace(frags, apex_rt, mcal, nb, hw, wppm, ppm_weighted):
    tw = 2.0 * RTW / 32; nf = len(frags); Xt = np.zeros((32, nb, nf))
    for k, (rt, im, it, ppm, _) in enumerate(frags):
        rb = np.floor((rt - apex_rt + RTW) / tw).astype(int); ib = np.floor((im - mcal + hw) / (2 * hw) * nb).astype(int)
        m = (np.abs(ppm) <= PPM_HARD) & (rb >= 0) & (rb < 32) & (ib >= 0) & (ib < nb)
        val = it * wppm[k] if ppm_weighted else it
        np.add.at(Xt[:, :, k], (rb[m], ib[m]), val[m])
    core = Xt[14:19].sum((0, 2))
    jb = int(np.argmax(np.convolve(core, [0.25, 0.5, 0.25], 'same'))) if core.sum() > 0 else nb // 2
    sd = mf.IM_SIG / (2 * hw / nb)
    return (Xt * np.exp(-0.5 * ((np.arange(nb) - jb) / sd) ** 2)[None, :, None]).sum(1)

def work(pids):
    T, groups, PG, lib, prec = G['T'], G['groups'], G['PG'], G['lib'], G['prec']
    rows = []
    for pid in pids:
        idx = groups.get(pid)
        if idx is None or pid not in PG: continue
        sub = T.take(idx).to_pydict(); frags = []
        for k in range(len(idx)):
            rt = mf.dec(sub['RT_DATA'][k], sub['RT_COMPRESSION'][k]); im = mf.dec(sub['MOBILITY_DATA'][k], sub['MOBILITY_COMPRESSION'][k])
            it = mf.dec(sub['INTENSITY_DATA'][k], sub['INTENSITY_COMPRESSION'][k]); mz = mf.dec(sub['MZ_DATA'][k], sub['MZ_COMPRESSION'][k])
            tmz = sub['TARGET_MZ'][k]; ppm = (mz - tmz) / tmz * 1e6 if len(mz) == len(it) else np.zeros(len(it))
            frags.append((rt, im, it, ppm, lib.get(sub['TRANSITION_ID'][k], 0.0)))
        nf = len(frags)
        if nf < 2: continue
        lw = np.sqrt(np.maximum([f[4] for f in frags], 0)); lw = lw / (np.linalg.norm(lw) + 1e-12); lshare = lw ** 2 / max((lw ** 2).sum(), 1e-12)
        mcal = mf.IM_A + mf.IM_B * prec.get(pid, np.nan)
        if not np.isfinite(mcal): continue
        allrt = np.sort(np.concatenate([f[0] for f in frags]))
        if len(allrt) < 5: continue
        starts = np.r_[0, np.where(np.diff(allrt) > 0.3)[0] + 1]
        scans = np.array([allrt[s:e].mean() for s, e in zip(starts, np.r_[starts[1:], len(allrt)])])
        if len(scans) < 5: continue
        mids = 0.5 * (scans[:-1] + scans[1:]); si = [np.searchsorted(mids, f[0]) for f in frags]
        wppm = [np.exp(-0.5 * ((f[3] - mf.PPM_OFF) / mf.PPM_SIG) ** 2) for f in frags]
        hard = [(np.abs(f[3]) <= PPM_HARD).astype(float) for f in frags]
        for _, pk in PG[pid].iterrows():
            ia = int(np.argmin(np.abs(scans - pk.rt_apex)))
            a, b = np.searchsorted(scans, pk.rt_start - 1e-6), np.searchsorted(scans, pk.rt_end + 1e-6)
            a_t = int(np.clip(np.floor((pk.rt_start - pk.rt_apex + RTW) / (2 * RTW / 32)), 0, 31))
            b_t = int(np.clip(np.floor((pk.rt_end - pk.rt_apex + RTW) / (2 * RTW / 32)) + 1, 1, 32))
            r = dict(precursor_id=pid, rta=pk.rta)
            ms_w = own_apex(frags, pk.rt_apex, mcal, wppm); ms_u = own_apex(frags, pk.rt_apex, mcal, hard)
            g = lambda im, c: np.where(np.abs(im - c) <= mf.IM_CUT, np.exp(-0.5 * ((im - c) / mf.IM_SIG) ** 2), 0.0)
            frame = {
                'mf': [wppm[k] * g(f[1], ms_w) for k, f in enumerate(frags)],
                'p': [wppm[k] * (np.abs(f[1] - mcal) <= 0.0314) for k, f in enumerate(frags)],
                'm': [hard[k] * g(f[1], ms_u) for k, f in enumerate(frags)],
            }
            for v, wts in frame.items():
                st, _ = mf.mf_stats(frame_trace(frags, si, len(scans), wts), ia, a, b, lw, lshare, v); r.update(st)
            for v, (nb, hw, pw) in {'tp': (8, 0.03, True), 't16': (16, 0.06, False), 't16p': (16, 0.06, True)}.items():
                st, _ = mf.mf_stats(tensor_trace(frags, pk.rt_apex, mcal, nb, hw, wppm, pw), 16, a_t, b_t, lw, lshare, v, dom_win=0, loo_win=2)
                r.update(st)
            rows.append(r)
    return rows

if __name__ == '__main__' and os.environ.get('STEP2', 'extract') == 'extract':
    P = pd.read_csv(f'{O}/peaks.tsv', sep='\t'); P['rta'] = P.rt_apex.round(3); PG = {p: g for p, g in P.groupby('precursor_id')}
    c = sqlite3.connect(f'{O}/mini.pqp')
    lib = pd.read_sql('select ID tid, LIBRARY_INTENSITY li from TRANSITION', c).set_index('tid').li.to_dict()
    prec = pd.read_sql('select ID pid, LIBRARY_DRIFT_TIME lim from PRECURSOR', c).set_index('pid').lim.to_dict()
    T = pq.read_table(f'{O}/maps.xipm', columns=['MS_LEVEL', 'PRECURSOR_ID', 'TRANSITION_ID', 'TARGET_MZ', 'MZ_DATA', 'RT_DATA', 'MOBILITY_DATA',
                                                 'INTENSITY_DATA', 'MZ_COMPRESSION', 'RT_COMPRESSION', 'MOBILITY_COMPRESSION', 'INTENSITY_COMPRESSION'])
    T = T.filter(pa.compute.equal(T.column('MS_LEVEL'), 2))
    groups = pd.Series(np.arange(T.num_rows)).groupby(T.column('PRECURSOR_ID').to_numpy()).apply(list).to_dict()
    G.update(T=T, groups=groups, PG=PG, lib=lib, prec=prec)
    pids = list(groups.keys()); chunks = [pids[i::64] for i in range(64)]
    with Pool(int(os.environ.get('NPROC', '16'))) as pool:
        res = pool.map(work, chunks)
    D = pd.DataFrame([r for rr in res for r in rr]); D.to_parquet(OUT); print(f'windows {len(D)}', flush=True)

if __name__ == '__main__' and os.environ.get('STEP2', 'extract') == 'eval':
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.metrics import roc_auc_score
    S = pd.read_csv(f'{O}/sample.tsv', sep='\t').set_index('pid'); keep = set(S.index)
    F = pd.concat([ch[ch.precursor_id.isin(keep)] for ch in pd.read_csv(os.environ['NEW_FEATS'], sep='\t', chunksize=500000)]); F['rta'] = F.rt_apex.round(3)
    F = F.merge(pd.read_parquet(OUT).drop_duplicates(['precursor_id', 'rta']), on=['precursor_id', 'rta'])
    meta = ['candidate_row', 'precursor_id', 'is_decoy', 'fold', 'rt_apex', 'score', 'q_value', 'rta']
    allv = [c for v in VARS for c in COLS[v]]
    B = [c for c in F.columns if c not in meta + allv]
    Pq = pd.read_csv(os.environ['NEW_PREC'], sep='\t', usecols=['precursor_id', 'is_decoy', 'q_value'])
    bq = Pq[Pq.is_decoy == 0].groupby('precursor_id').q_value.min()
    dec = F.is_decoy.values == 1
    pos = (F.precursor_id.map(S.cls).values == 'dnn') & (np.abs(F.rt_apex.values - F.precursor_id.map(S.rt_s).values) <= 6)
    missed = set(S.index[S.cls == 'dnn']) - set(bq.index[bq <= 0.01]); toolpos = (~dec) & (F.q_value.values <= 0.01)
    fold = F.fold.values; dec_ids = S.index[S.cls == 'decoy']
    print(f'joined {len(F)}, baseline features {len(B)} (mt_ in baseline: {sum(c.startswith("mt_") for c in B)}), right windows {pos.sum()}, '
          f'missed DIA-NN ids {len(missed)}, tool positives {toolpos.sum()}', flush=True)
    def gbm(M, sd, yp):
        s = np.zeros(len(F))
        for f in np.unique(fold):
            tr = (fold != f) & (yp | dec); y = yp[tr].astype(int); w = np.where(y == 1, (y == 0).sum() / max(y.sum(), 1), 1.0)
            s[fold == f] = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=31, random_state=sd).fit(M[tr], y, sample_weight=w).decision_function(M[fold == f])
        return s
    def tpr(s):
        F['_r'] = s; best = F.groupby('precursor_id')._r.max(); rb = F[pos].groupby('precursor_id')._r.max()
        t = np.quantile(best.reindex(dec_ids).dropna().values, 1 - 3e-3); rm = rb[rb.index.isin(missed)]
        return (rb > t).mean(), (rm > t).mean()
    res = {}
    for lab, yp in (('ceiling', pos), ('tool', toolpos)):
        res[(lab, 'B45')] = base = [tpr(gbm(F[B].values.astype(np.float32), sd, yp)) for sd in (0, 1)]
        print(f'{lab:7s} B45          all {base[0][0]:.3f} / {base[1][0]:.3f}   missed {base[0][1]:.3f} / {base[1][1]:.3f}', flush=True)
        for v in VARS:
            r = [tpr(gbm(F[B + COLS[v]].values.astype(np.float32), sd, yp)) for sd in (0, 1)]; res[(lab, v)] = r
            print(f'{lab:7s} B45+{v:5s}    all {r[0][0] - base[0][0]:+.3f} / {r[1][0] - base[1][0]:+.3f}   '
                  f'missed {r[0][1] - base[0][1]:+.3f} / {r[1][1] - base[1][1]:+.3f}', flush=True)
    d = lambda lab, v, i: [res[(lab, v)][s][i] - res[(lab, 'B45')][s][i] for s in (0, 1)]
    passing = []
    for v in ('tp', 't16', 't16p'):
        ok = all(x >= 0.030 for x in d('tool', v, 1)) and all(x >= 0.015 for x in d('ceiling', v, 0))
        print(f'GATE {v}: tool missed {d("tool", v, 1)}, ceiling all {d("ceiling", v, 0)} -> {"PASS" if ok else "fail"}')
        if ok: passing.append(v)
    choice = None
    for v in passing:
        if choice is None or all(x >= y + 0.015 for x, y in zip(d('tool', v, 1), d('tool', choice, 1))): choice = v
    print(f'CHOICE: {choice}')
    rm = pos & F.precursor_id.isin(missed).values
    y = np.r_[np.ones(rm.sum()), np.zeros(dec.sum())]
    for v in VARS:
        print(f'  AUC missed-right vs decoy {v:5s} ' + '  '.join(f'{s} {roc_auc_score(y, np.r_[F.loc[rm, f"{v}_{s}"].values, F.loc[dec, f"{v}_{s}"].values]):.3f}'
                                                         for s in ('z_med', 'n_z2', 'llr', 'loo_med', 'cos')))
    for s in ('z_med', 'n_z2', 'llr', 'loo_med', 'cos'):
        c = f'mt_{s}'
        print(f'  AUC missed-right vs decoy C++ mt {s:7s} {roc_auc_score(y, np.r_[F.loc[rm, c].values, F.loc[dec, c].values]):.3f}')
