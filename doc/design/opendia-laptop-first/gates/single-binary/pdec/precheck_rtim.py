"""PRECHECK_RTIM: Sage-style per-run linear correction of library RT and 1/K0 (diagnostic, 20 Nov 2026).

Idea (as Sage): fit, per run, a ridge regression on a featurised sequence using that run's own confident IDs, then predict
every precursor. Here it is a correction on top of the run's existing calibration, so the target is the residual.

Training set (label-free selection): the run's own OpenDIA target IDs at q <= 0.01 (new default, seed 0, r~7 library),
cross-fitted in 5 folds by stripped sequence, so every precursor is predicted by a model that never saw it.
  RT target   OpenDIA apex RT - current map(library RT)                      (OpenDIA measurement, no DIA-NN data)
  IM target   DIA-NN reported 1/K0 - calibrated library 1/K0 (im_a + im_b x lib)  (CAVEAT: OpenDIA does not write an observed
              1/K0 per ID, so the measurement is taken from DIA-NN for IDs OpenDIA itself accepted; DIA-NN chooses no IDs)
Features: residue counts (modified residues as their own tokens), N- and C-terminal residue one-hots, length, log length,
charge one-hot, precursor m/z; variant 'hyb' adds the mapped library value itself (lets the model rescale it).
Evaluation: DIA-NN non-entrapment IDs of the same run, all and the ones OpenDIA misses (not in its q <= 0.01 set).
Reported: |residual| p50 / p95 for current map, corrected map, and the best monotone map fitted on DIA-NN's own IDs
(the PRECHECK floor). Also the RT window needed to cover 95 % / 99 % of the DIA-NN IDs.

Frozen gate (written before any result): proceed to C++ only if, on H, the corrected RT p95 on missed DIA-NN IDs is at
least 15 % below the current map's, and on D it is not worse by more than 5 %. IM proceeds separately on the same
rule applied to |1/K0 residual| p95."""
import json, glob, re, zlib, sqlite3, numpy as np, pandas as pd
from sklearn.linear_model import RidgeCV
from sklearn.preprocessing import StandardScaler
C = '/nfs/research/juan/singj/opendia'
LIB = f'{C}/entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp'
c = sqlite3.connect(f'file:{LIB}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID pid, p.DECOY d, p.CHARGE z, p.PRECURSOR_MZ mz, p.LIBRARY_RT lrt, p.LIBRARY_DRIFT_TIME lim,
  pe.MODIFIED_SEQUENCE ms, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c).set_index('pid')
tg = pq[pq.d == 0].copy(); sk = (tg.ms + '/' + tg.z.astype(str)).reset_index().drop_duplicates(0).set_index(0).pid
TOK = re.compile(r'[A-Z](?:\(UniMod:\d+\))?|\(UniMod:\d+\)')
def tokens(s): return TOK.findall(s.lstrip('.'))
vocab = sorted({t for s in tg.ms.unique() for t in tokens(s)})
AA = sorted({t[0] for t in vocab if t[0].isalpha()})
print(f'library targets {len(tg)}, token vocab {len(vocab)}: {vocab}')

def featurise(df):
    seqs = df.ms.values; n = len(seqs); V = {t: i for i, t in enumerate(vocab)}; A = {a: i for i, a in enumerate(AA)}
    X = np.zeros((n, len(vocab) + 2 * len(AA) + 2 + 4 + 1), np.float32)
    for r, s in enumerate(seqs):
        t = tokens(s)
        for x in t: X[r, V[x]] += 1
        res = [x for x in t if x[0].isalpha()]
        if res:
            X[r, len(vocab) + A[res[0][0]]] = 1; X[r, len(vocab) + len(AA) + A[res[-1][0]]] = 1
        o = len(vocab) + 2 * len(AA); X[r, o] = len(res); X[r, o + 1] = np.log(max(len(res), 1))
    o = len(vocab) + 2 * len(AA) + 2
    z = np.clip(df.z.values.astype(int), 1, 4); X[np.arange(n), o + z - 1] = 1
    X[:, -1] = df.mz.values
    return X

def apply_map(cal, lrt):
    r = cal['rt']; kx = np.asarray(r['knots_library_rt']); ky = np.asarray(r['knots_run_rt_s'])
    y = np.interp(lrt, kx, ky); lo, hi = lrt < kx[0], lrt > kx[-1]
    y[lo] = ky[0] + (lrt[lo] - kx[0]) * r['low_slope_s']; y[hi] = ky[-1] + (lrt[hi] - kx[-1]) * r['high_slope_s']
    return y

def fold_of(ms): return np.array([zlib.crc32(re.sub(r'\(UniMod:\d+\)', '', s).encode()) % 5 for s in ms])

def crossfit(X, y, train_mask, folds, extra=None):
    """Ridge per fold on train rows of the other folds; returns predictions for every row."""
    Xf = X if extra is None else np.hstack([X, extra[:, None]])
    pred = np.zeros(len(y)); alphas = []
    for f in range(5):
        tr = train_mask & (folds != f)
        sc = StandardScaler().fit(Xf[tr]); m = RidgeCV(alphas=np.logspace(-2, 4, 13)).fit(sc.transform(Xf[tr]), y[tr])
        pred[folds == f] = m.predict(sc.transform(Xf[folds == f])); alphas.append(m.alpha_)
    return pred, alphas

def q(a, p): return float(np.quantile(np.abs(a), p))

out = {}
for run in ('D', 'H'):
    rd = f'{C}/entlib-mt-20261120/runs/{run}_mt_s0'
    cal = json.load(open(glob.glob(f'{rd}/*.calibration.json')[0]))
    P = pd.read_csv([f for f in glob.glob(f'{rd}/*.precursors.tsv') if not f.endswith('/precursors.tsv')][0], sep='\t',
                    usecols=['precursor_id', 'is_decoy', 'rt_apex', 'score', 'q_value'])
    P = P[P.score > -1e29].sort_values('score', ascending=False).drop_duplicates('precursor_id').set_index('precursor_id')
    ours = P[(P.is_decoy == 0) & (P.q_value <= 0.01)]
    d = pd.read_parquet(f'{C}/entlib-diann5-20261119/{run}/report.parquet', columns=['Modified.Sequence', 'Precursor.Charge', 'Q.Value', 'Decoy', 'RT', 'IM'])
    d = d[(d['Q.Value'] <= 0.01) & (d.Decoy == 0)].copy()
    d['pid'] = sk.reindex(d['Modified.Sequence'] + '/' + d['Precursor.Charge'].astype(str)).values
    d = d.dropna(subset=['pid']); d['pid'] = d.pid.astype(np.int64)
    d = d[tg.e.reindex(d.pid).values == 0].drop_duplicates('pid').set_index('pid')
    # one table over every precursor that is either an OpenDIA ID or a DIA-NN ID
    ids = sorted(set(ours.index) | set(d.index)); T = tg.loc[ids].copy()
    T['ours'] = T.index.isin(ours.index); T['dnn'] = T.index.isin(d.index)
    T['rt_obs_ours'] = ours.rt_apex.reindex(T.index).values; T['rt_dnn'] = d.RT.reindex(T.index).values * 60
    T['im_dnn'] = d.IM.reindex(T.index).values
    T['rt_map'] = apply_map(cal, T.lrt.values)
    T['im_map'] = cal['im']['intercept_vssc'] + cal['im']['slope'] * T.lim.values
    X = featurise(T); folds = fold_of(T.ms.values); tr = T.ours.values
    print(f'\n===== {run}: OpenDIA q<=0.01 IDs {tr.sum()} (training pool), DIA-NN IDs {T.dnn.sum()}, overlap {(T.ours & T.dnn).sum()}, '
          f'missed {(T.dnn & ~T.ours).sum()}; features {X.shape[1]}')
    # RT: target from OpenDIA's own apex RT
    y_rt = np.where(tr, T.rt_obs_ours.values - T.rt_map.values, 0.0)
    corr_comp, a1 = crossfit(X, y_rt, tr, folds)
    corr_hyb, a2 = crossfit(X, y_rt, tr, folds, extra=T.rt_map.values)
    ev = T.dnn.values; mi = ev & ~tr
    res = {'map': T.rt_dnn.values - T.rt_map.values}
    res['map+comp'] = res['map'] - corr_comp; res['map+hyb'] = res['map'] - corr_hyb
    from sklearn.isotonic import IsotonicRegression
    o = np.argsort(T.lrt.values[ev]); iso = IsotonicRegression(out_of_bounds='clip').fit(T.lrt.values[ev][o], T.rt_dnn.values[ev][o])
    res['floor (monotone on DIA-NN)'] = np.full(len(T), np.nan); res['floor (monotone on DIA-NN)'][ev] = T.rt_dnn.values[ev] - iso.predict(T.lrt.values[ev])
    span = np.quantile(T.rt_dnn.values[ev], .99) - np.quantile(T.rt_dnn.values[ev], .01)
    print(f'  RT (s), DIA-NN IDs: all p50/p95 | missed p50/p95 | window for 95 % / 99 % of all   (span {span:.0f} s; ridge alpha {np.median(a1):.0f}/{np.median(a2):.0f})')
    for k, v in res.items():
        print(f'    {k:28s} {q(v[ev], .5):6.1f} / {q(v[ev], .95):6.1f} | {q(v[mi], .5):6.1f} / {q(v[mi], .95):6.1f} | {q(v[ev], .95):6.1f} / {q(v[ev], .99):6.1f}')
    out[(run, 'rt')] = {k: q(v[mi], .95) for k, v in res.items()}
    # sanity: OpenDIA apex vs DIA-NN RT on the overlap (the RT target is measured as DIA-NN does)
    ov = T.ours.values & T.dnn.values
    print(f'    OpenDIA apex vs DIA-NN RT on overlap: |diff| p50 {q(T.rt_obs_ours.values[ov] - T.rt_dnn.values[ov], .5):.2f} p95 {q(T.rt_obs_ours.values[ov] - T.rt_dnn.values[ov], .95):.2f} s')
    # IM: target = DIA-NN 1/K0 for OpenDIA IDs (see caveat)
    trim = tr & np.isfinite(T.im_dnn.values) & np.isfinite(T.im_map.values)
    y_im = np.where(trim, T.im_dnn.values - T.im_map.values, 0.0)
    ic_comp, _ = crossfit(X, y_im, trim, folds); ic_hyb, _ = crossfit(X, y_im, trim, folds, extra=np.nan_to_num(T.im_map.values))
    ri = {'map': T.im_dnn.values - T.im_map.values}; ri['map+comp'] = ri['map'] - ic_comp; ri['map+hyb'] = ri['map'] - ic_hyb
    evi = ev & np.isfinite(ri['map']); mii = evi & ~tr
    print(f'  1/K0, DIA-NN IDs: all p50/p95 | missed p50/p95   (training pool with a DIA-NN 1/K0: {trim.sum()}; current tolerance '
          f'{0.6 * cal["im"]["tolerance_half_width_vssc"]:.4f})')
    for k, v in ri.items():
        print(f'    {k:28s} {q(v[evi], .5):.4f} / {q(v[evi], .95):.4f} | {q(v[mii], .5):.4f} / {q(v[mii], .95):.4f}')
    out[(run, 'im')] = {k: q(v[mii], .95) for k, v in ri.items()}

print('\n===== GATE (missed DIA-NN IDs, p95)')
for kind in ('rt', 'im'):
    best = min(('map+comp', 'map+hyb'), key=lambda k: out[('H', kind)][k])
    h = out[('H', kind)][best] / out[('H', kind)]['map'] - 1; dd = out[('D', kind)][best] / out[('D', kind)]['map'] - 1
    ok = h <= -0.15 and dd <= 0.05
    print(f'  {kind}: variant {best}; H {h * 100:+.1f} % (need <= -15 %), D {dd * 100:+.1f} % (need <= +5 %) -> {"PASS" if ok else "fail"}')
