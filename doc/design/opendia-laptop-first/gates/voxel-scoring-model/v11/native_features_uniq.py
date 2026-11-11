"""PLAN_V11 features U: native-sampling fragment evidence (empirical per-m/z null, as native_features_emp.py) on
UNIQUE fragments only: a fragment of precursor p is dropped for window w (apex t) when any other TARGET library precursor
q (q != p and q != p's target-decoy pair partner) whose precursor m/z lies in the same isolation window and whose
predicted run RT is within the run's RT tolerance of t has a library fragment within 2 x tolerance of it.
Same rule for target and decoy p (label-free, model-free; entrapment labels are never read).
Neighbour fragments: the union of the full 6-fragment library (all 300,223 targets) and the 12-fragment population
library (its targets). Predicted run RT: the run's calibration knots (cal.json rt, np.interp, constant ends as the
search that built this population). RT tolerance: cal.json rt.tolerance_half_width_s (84 s on raw 01).
Columns (natU.npy [n, 6]): the 5 columns of native_features_emp.py computed on the unique fragments, then n_unique.
Usage: native_features_uniq.py <raw01|raw02>."""
import json, sqlite3, sys, time, numpy as np, pandas as pd, pyopenms as oms
RAW = sys.argv[1]
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; V = B + '/astral-voxel-pilot-20261110/v2'
CAL = {'raw01': B + '/astral-design-raw01-20261107/window300k/cal.json', 'raw02': B + '/astral-entrapment-raw02-20261108/cal/cal.json'}[RAW]
MZML = B + f'/dd-plots-20261111/{RAW}.mzML'
cj = json.load(open(CAL)); OFF, TOL = float(cj['ms2_mz']['offset_ppm']), float(cj['ms2_mz']['tolerance_half_width_ppm'])
KX, KY, RTTOL = np.array(cj['rt']['knots_library_rt']), np.array(cj['rt']['knots_run_rt_s']), float(cj['rt']['tolerance_half_width_s'])
K, NOFF, PMIN, PMAX, LCAP, EXCL, NULLW = 12, 7, 1e-3, 0.999, 30.0, 8.0, 60.0
c12 = sqlite3.connect('file:' + B + '/astral-lib12-20261111/library12_pop.pqp?mode=ro', uri=True)
c6 = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
pre = pd.read_sql('select ID precursor_id, TRAML_ID traml, PRECURSOR_MZ pmz, LIBRARY_RT lrt, DECOY decoy from PRECURSOR', c6).set_index('precursor_id')
tid = pre[pre.decoy == 0].reset_index().set_index('traml').precursor_id
partner = pd.Series(np.where(pre.decoy == 1, pre.traml.str.replace('^DECOY_', '', regex=True).map(tid), -1), index=pre.index).fillna(-1).astype(np.int64)
tr = pd.read_sql("""select m.PRECURSOR_ID precursor_id, t.PRODUCT_MZ mz, t.LIBRARY_INTENSITY li from TRANSITION_PRECURSOR_MAPPING m
  join TRANSITION t on t.ID=m.TRANSITION_ID""", c12).sort_values(['precursor_id', 'li', 'mz'], ascending=[True, False, True])
tr['k'] = tr.groupby('precursor_id').cumcount(); tr = tr[tr.k < K]
# Neighbour fragments (targets only): full 6-fragment library + 12-fragment population library, uncalibrated library m/z.
nb6 = pd.read_sql("""select m.PRECURSOR_ID q, t.PRODUCT_MZ mz from TRANSITION_PRECURSOR_MAPPING m join TRANSITION t on t.ID=m.TRANSITION_ID
  join PRECURSOR p on p.ID=m.PRECURSOR_ID where p.DECOY=0""", c6)
nb12 = pd.read_sql("""select m.PRECURSOR_ID q, t.PRODUCT_MZ mz from TRANSITION_PRECURSOR_MAPPING m join TRANSITION t on t.ID=m.TRANSITION_ID
  join PRECURSOR p on p.ID=m.PRECURSOR_ID where p.DECOY=0""", c12)
nb = pd.concat([nb6, nb12]).drop_duplicates()
nb['pmz'] = pre.pmz.reindex(nb.q).values; nb['rt'] = np.interp(pre.lrt.reindex(nb.q).values, KX, KY)
pop = pd.read_csv(f'{V}/{RAW}/population.tsv', sep='\t'); n = len(pop)
pids = np.unique(pop.precursor_id.values); pix = pd.Series(np.arange(len(pids)), index=pids)
FM = np.full((len(pids), K), np.nan); FL = np.zeros((len(pids), K)); FR = np.full((len(pids), K), np.nan)  # calibrated / raw library m/z
t_ = tr[tr.precursor_id.isin(pix.index)]; ii = pix[t_.precursor_id].values
FM[ii, t_.k.values] = t_.mz.values * (1 + OFF * 1e-6); FL[ii, t_.k.values] = t_.li.values; FR[ii, t_.k.values] = t_.mz.values
wp = pix[pop.precursor_id].values; apex = pop.rt_apex.values; wmz = pre.pmz.reindex(pop.precursor_id).values
wpid = pop.precursor_id.values.astype(np.int64); wpart = partner.reindex(pop.precursor_id).values

od = oms.OnDiscMSExperiment(); od.openFile(MZML, False); meta = od.getMetaData(); rows = []
for i in range(meta.getNrSpectra()):
    sp = meta.getSpectrum(i)
    if sp.getMSLevel() != 2 or not sp.getPrecursors(): continue
    p = sp.getPrecursors()[0]; rows.append((i, sp.getRT(), round(p.getMZ() - p.getIsolationWindowLowerOffset(), 3), round(p.getMZ() + p.getIsolationWindowUpperOffset(), 3)))
scans = pd.DataFrame(rows, columns=['idx', 'rt', 'lo', 'hi']).sort_values('rt')
wins = scans[['lo', 'hi']].drop_duplicates().sort_values('lo').reset_index(drop=True)
wi = np.searchsorted(wins.lo.values, wmz, side='right') - 1
ok = (wi >= 0) & (wmz <= wins.hi.values[np.clip(wi, 0, None)])
print(RAW, 'MS2 scans', len(scans), 'isolation windows', len(wins), 'tol ppm', TOL, 'rt tol s', RTTOL, 'neighbour fragments', len(nb), flush=True)

L = np.zeros((n, NOFF), np.float32); M = np.zeros((n, NOFF, K), np.uint8); I = np.zeros((n, NOFF, K), np.float32)
P = np.zeros((n, K), np.float32); U = np.zeros((n, K), bool)

def ptail(p, k):
    m = p.shape[0]; dp = np.zeros((m, K + 1)); dp[:, 0] = 1
    for f in range(K):
        q = np.nan_to_num(p[:, f]); nd = dp * (1 - q)[:, None]; nd[:, 1:] += dp[:, :-1] * q[:, None]; dp = nd
    cs = np.cumsum(dp[:, ::-1], 1)[:, ::-1]
    return np.minimum(-np.log10(np.maximum(cs[np.arange(m), k], 1e-300)), LCAP)

def expand(sorted_vals, lo, hi):
    """Indices into sorted_vals within [lo, hi] per query: (query index, position)."""
    a = np.searchsorted(sorted_vals, lo); b = np.searchsorted(sorted_vals, hi, 'right'); ln = b - a
    q = np.repeat(np.arange(len(ln)), ln); pos = np.arange(ln.sum()) - np.repeat(np.cumsum(ln) - ln, ln) + np.repeat(a, ln)
    return q, pos

t0 = time.time(); n_shared = 0; n_valid = 0
for g in range(len(wins)):
    sc = scans[(scans.lo == wins.lo[g]) & (scans.hi == wins.hi[g])]; srt = sc.rt.values; sidx = sc.idx.values
    w = np.where(ok & (wi == g))[0]
    if len(w) == 0 or len(sc) == 0: continue
    # Sharing mask: neighbour targets of this isolation window, fragments sorted by raw library m/z.
    ng = nb[(nb.pmz >= wins.lo[g]) & (nb.pmz <= wins.hi[g])].sort_values('mz'); nmz, nrt, nq = ng.mz.values, ng.rt.values, ng.q.values.astype(np.int64)
    Fr = FR[wp[w]]; valid = ~np.isnan(Fr); Frz = np.where(valid, Fr, 0.0); stol = 2 * Frz * TOL * 1e-6
    q, pos = expand(nmz, (Frz - stol).ravel(), (Frz + stol).ravel()); ww_ = q // K
    hit = (np.abs(nrt[pos] - apex[w][ww_]) <= RTTOL) & (nq[pos] != wpid[w][ww_]) & (nq[pos] != wpart[w][ww_])
    shared = np.zeros(Fr.size, bool); shared[q[hit]] = True; shared = shared.reshape(Fr.shape) & valid
    uniq = valid & ~shared; U[w] = uniq; n_shared += shared.sum(); n_valid += valid.sum()
    # Empirical null (as native_features_emp.py) for the unique fragments.
    spk = [np.asarray(od.getSpectrum(int(i)).get_peaks()[0], float) for i in sidx]
    allmz = np.concatenate(spk); allj = np.repeat(np.arange(len(spk)), [len(x) for x in spk]); o = np.argsort(allmz, kind='stable'); allmz, allj = allmz[o], allj[o]
    F = FM[wp[w]]; Fz = np.where(uniq, F, 0.0); tolmz = Fz * TOL * 1e-6
    a0 = np.searchsorted(srt, apex[w] - NULLW); a1 = np.searchsorted(srt, apex[w] - EXCL); b0 = np.searchsorted(srt, apex[w] + EXCL, 'right'); b1 = np.searchsorted(srt, apex[w] + NULLW, 'right')
    nnull = (a1 - a0) + (b1 - b0)
    q, pos = expand(allmz, (Fz - tolmz).ravel(), (Fz + tolmz).ravel()); jj = allj[pos]; ww_ = q // K
    inn = ((jj >= a0[ww_]) & (jj < a1[ww_])) | ((jj >= b0[ww_]) & (jj < b1[ww_]))
    key = np.unique(q[inn] * 100000 + jj[inn]) // 100000
    cnt = np.bincount(key, minlength=F.size).reshape(F.shape)
    pf = np.clip((cnt + 0.5) / (nnull[:, None] + 1.0), PMIN, PMAX); pf = np.where(uniq, pf, np.nan); P[w] = np.nan_to_num(pf)
    j0 = np.clip(np.searchsorted(srt, apex[w]), 1, len(srt) - 1); j0 -= (apex[w] - srt[j0 - 1]) < (srt[j0] - apex[w])
    pw = np.repeat(np.arange(len(w)), NOFF); po = np.tile(np.arange(NOFF), len(w)); pj = np.repeat(j0, NOFF) + po - 3
    keep = (pj >= 0) & (pj < len(srt)); pw, po, pj = pw[keep], po[keep], pj[keep]
    order = np.argsort(pj, kind='stable'); pw, po, pj = pw[order], po[order], pj[order]
    bounds = np.r_[0, np.flatnonzero(np.diff(pj)) + 1, len(pj)]
    for a, b in zip(bounds[:-1], bounds[1:]):
        mz, it = od.getSpectrum(int(sidx[pj[a]])).get_peaks(); mz = np.asarray(mz, float); it = np.asarray(it, float)
        lw, oo = pw[a:b], po[a:b]; Fq = Fz[lw]; vq = uniq[lw]; tq = tolmz[lw]
        if len(mz) == 0: mat = np.zeros_like(Fq, bool); inten = np.zeros_like(Fq)
        else:
            j = np.clip(np.searchsorted(mz, Fq), 1, len(mz) - 1); jn = np.where(np.abs(mz[j - 1] - Fq) < np.abs(mz[j] - Fq), j - 1, j)
            mat = (np.abs(mz[jn] - Fq) <= tq) & vq; inten = np.where(mat, it[jn], 0.0)
        L[w[lw], oo] = ptail(pf[lw], mat.sum(1)); M[w[lw], oo] = mat; I[w[lw], oo] = inten
    if g % 30 == 0: print(RAW, f'isolation window {g}/{len(wins)} ({time.time() - t0:.0f} s), shared fraction so far {n_shared / max(n_valid, 1):.3f}', flush=True)

out = np.zeros((n, 6), np.float32)
out[:, 0] = L[:, 2:5].max(1); out[:, 1] = L.sum(1); out[:, 2] = (L >= 2).sum(1)
cm = M[:, 2:5].sum(1) >= 2; qq = np.where(U, 3 * P ** 2 - 2 * P ** 3, np.nan)
for a in range(0, n, 200000):
    out[a:a + 200000, 3] = ptail(np.clip(qq[a:a + 200000], 1e-9, 0.999), cm[a:a + 200000].sum(1))
ia = np.sqrt(I[:, 2:5].sum(1)); il = np.sqrt(np.where(U, FL[wp], 0.0))
out[:, 4] = (ia * il).sum(1) / np.maximum(np.linalg.norm(ia, axis=1) * np.linalg.norm(il, axis=1), 1e-12)
out[:, 5] = U.sum(1)
np.save(f'{V}/{RAW}/natU.npy', out)
dec = pop.is_decoy.values == 1
print(RAW, 'done', f'{time.time() - t0:.0f} s', 'shared fraction', round(n_shared / max(n_valid, 1), 4),
      'n_unique median target/decoy', float(np.median(out[~dec, 5])), float(np.median(out[dec, 5])),
      'mean n_unique target/decoy', round(float(out[~dec, 5].mean()), 3), round(float(out[dec, 5].mean()), 3),
      'medians target', np.median(out[~dec], 0).round(3).tolist(), 'decoy', np.median(out[dec], 0).round(3).tolist(), flush=True)
