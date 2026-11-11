"""PLAN_V10 features E: native-sampling fragment evidence with an EMPIRICAL per-m/z local null (label-free, model-free,
same rule for targets and decoys). Usage: native_features_emp.py <raw01|raw02>. Writes <raw>/natE.npy [n, 5].
Same scans, fragments, tolerance and 5 columns as native_features.py (PLAN_V9); only the chance rate p_f changes:
  p_f = (#scans of the same isolation window with a peak within tolerance of f, |RT - apex| in (EXCL, NULLW] s) + 0.5)
        / (#such scans + 1),  clipped to [1e-3, 0.999];  EXCL = 8 s, NULLW = 60 s (about 63 null scans).
Reason (NAT_NULL_CHECK_raw01.txt): the uniform-density rate predicted 0.08 matched fragments per away-from-apex decoy
scan against 1.38 observed; specific fragment m/z recur across scans (shared low-mass ions, co-isolated peptides), so the
null has to be measured per m/z, locally in RT, in the same isolation window."""
import json, sqlite3, sys, time, numpy as np, pandas as pd, pyopenms as oms
RAW = sys.argv[1]
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; V = B + '/astral-voxel-pilot-20261110/v2'
CAL = {'raw01': B + '/astral-design-raw01-20261107/window300k/cal.json', 'raw02': B + '/astral-entrapment-raw02-20261108/cal/cal.json'}[RAW]
MZML = B + f'/dd-plots-20261111/{RAW}.mzML'
cal = json.load(open(CAL))['ms2_mz']; OFF, TOL = float(cal['offset_ppm']), float(cal['tolerance_half_width_ppm'])
K, NOFF, PMIN, PMAX, LCAP, EXCL, NULLW = 12, 7, 1e-3, 0.999, 30.0, 8.0, 60.0
c = sqlite3.connect('file:' + B + '/astral-lib12-20261111/library12_pop.pqp?mode=ro', uri=True)
pmz = pd.read_sql('select ID precursor_id, PRECURSOR_MZ pmz from PRECURSOR', c).set_index('precursor_id').pmz
tr = pd.read_sql("""select m.PRECURSOR_ID precursor_id, t.PRODUCT_MZ mz, t.LIBRARY_INTENSITY li from TRANSITION_PRECURSOR_MAPPING m
  join TRANSITION t on t.ID=m.TRANSITION_ID""", c).sort_values(['precursor_id', 'li', 'mz'], ascending=[True, False, True])
tr['k'] = tr.groupby('precursor_id').cumcount(); tr = tr[tr.k < K]
pop = pd.read_csv(f'{V}/{RAW}/population.tsv', sep='\t'); n = len(pop)
pids = np.unique(pop.precursor_id.values); pix = pd.Series(np.arange(len(pids)), index=pids)
FM = np.full((len(pids), K), np.nan); FL = np.zeros((len(pids), K))
t_ = tr[tr.precursor_id.isin(pix.index)]; ii = pix[t_.precursor_id].values
FM[ii, t_.k.values] = t_.mz.values * (1 + OFF * 1e-6); FL[ii, t_.k.values] = t_.li.values
wp = pix[pop.precursor_id].values; apex = pop.rt_apex.values; wmz = pmz[pop.precursor_id].values

od = oms.OnDiscMSExperiment(); od.openFile(MZML, False); meta = od.getMetaData(); rows = []
for i in range(meta.getNrSpectra()):
    sp = meta.getSpectrum(i)
    if sp.getMSLevel() != 2 or not sp.getPrecursors(): continue
    p = sp.getPrecursors()[0]; rows.append((i, sp.getRT(), round(p.getMZ() - p.getIsolationWindowLowerOffset(), 3), round(p.getMZ() + p.getIsolationWindowUpperOffset(), 3)))
scans = pd.DataFrame(rows, columns=['idx', 'rt', 'lo', 'hi']).sort_values('rt')
wins = scans[['lo', 'hi']].drop_duplicates().sort_values('lo').reset_index(drop=True)
wi = np.searchsorted(wins.lo.values, wmz, side='right') - 1
ok = (wi >= 0) & (wmz <= wins.hi.values[np.clip(wi, 0, None)])
print(RAW, 'MS2 scans', len(scans), 'isolation windows', len(wins), 'tol ppm', TOL, 'windows without an isolation window', int((~ok).sum()), flush=True)

L = np.zeros((n, NOFF), np.float32); M = np.zeros((n, NOFF, K), np.uint8); I = np.zeros((n, NOFF, K), np.float32); P = np.zeros((n, K), np.float32)

def ptail(p, k):
    m = p.shape[0]; dp = np.zeros((m, K + 1)); dp[:, 0] = 1
    for f in range(K):
        q = np.nan_to_num(p[:, f]); nd = dp * (1 - q)[:, None]; nd[:, 1:] += dp[:, :-1] * q[:, None]; dp = nd
    cs = np.cumsum(dp[:, ::-1], 1)[:, ::-1]
    return np.minimum(-np.log10(np.maximum(cs[np.arange(m), k], 1e-300)), LCAP)

t0 = time.time()
for g in range(len(wins)):
    sc = scans[(scans.lo == wins.lo[g]) & (scans.hi == wins.hi[g])]; srt = sc.rt.values; sidx = sc.idx.values
    w = np.where(ok & (wi == g))[0]
    if len(w) == 0 or len(sc) == 0: continue
    # All peaks of this isolation window: (mz, scan position), sorted by mz.
    spk = [np.asarray(od.getSpectrum(int(i)).get_peaks()[0], float) for i in sidx]; sint = None
    allmz = np.concatenate(spk); allj = np.repeat(np.arange(len(spk)), [len(x) for x in spk]); o = np.argsort(allmz, kind='stable'); allmz, allj = allmz[o], allj[o]
    # Empirical null per (population window, fragment): distinct null scans with a peak within tolerance.
    F = FM[wp[w]]; valid = ~np.isnan(F); Fz = np.where(valid, F, 0.0); tolmz = Fz * TOL * 1e-6
    a0 = np.searchsorted(srt, apex[w] - NULLW); a1 = np.searchsorted(srt, apex[w] - EXCL); b0 = np.searchsorted(srt, apex[w] + EXCL, 'right'); b1 = np.searchsorted(srt, apex[w] + NULLW, 'right')
    nnull = (a1 - a0) + (b1 - b0)
    lo_ = np.searchsorted(allmz, (Fz - tolmz).ravel()); hi_ = np.searchsorted(allmz, (Fz + tolmz).ravel(), 'right'); ln = hi_ - lo_
    q = np.repeat(np.arange(len(ln)), ln); pos = np.arange(ln.sum()) - np.repeat(np.cumsum(ln) - ln, ln) + np.repeat(lo_, ln)
    jj = allj[pos]; ww_ = q // K
    inn = ((jj >= a0[ww_]) & (jj < a1[ww_])) | ((jj >= b0[ww_]) & (jj < b1[ww_]))
    key = np.unique(q[inn] * 100000 + jj[inn]) // 100000  # distinct (query, scan)
    cnt = np.bincount(key, minlength=len(ln)).reshape(F.shape)
    pf = np.clip((cnt + 0.5) / (nnull[:, None] + 1.0), PMIN, PMAX); pf = np.where(valid, pf, np.nan); P[w] = np.nan_to_num(pf)
    # Observed scans around the apex (offsets -3..3), as in native_features.py.
    j0 = np.clip(np.searchsorted(srt, apex[w]), 1, len(srt) - 1); j0 -= (apex[w] - srt[j0 - 1]) < (srt[j0] - apex[w])
    pw = np.repeat(np.arange(len(w)), NOFF); po = np.tile(np.arange(NOFF), len(w)); pj = np.repeat(j0, NOFF) + po - 3
    keep = (pj >= 0) & (pj < len(srt)); pw, po, pj = pw[keep], po[keep], pj[keep]
    order = np.argsort(pj, kind='stable'); pw, po, pj = pw[order], po[order], pj[order]
    bounds = np.r_[0, np.flatnonzero(np.diff(pj)) + 1, len(pj)]
    for a, b in zip(bounds[:-1], bounds[1:]):
        mz, it = od.getSpectrum(int(sidx[pj[a]])).get_peaks(); mz = np.asarray(mz, float); it = np.asarray(it, float)
        lw, oo = pw[a:b], po[a:b]; Fq = Fz[lw]; vq = valid[lw]; tq = tolmz[lw]
        if len(mz) == 0: mat = np.zeros_like(Fq, bool); inten = np.zeros_like(Fq)
        else:
            j = np.clip(np.searchsorted(mz, Fq), 1, len(mz) - 1); jn = np.where(np.abs(mz[j - 1] - Fq) < np.abs(mz[j] - Fq), j - 1, j)
            mat = (np.abs(mz[jn] - Fq) <= tq) & vq; inten = np.where(mat, it[jn], 0.0)
        L[w[lw], oo] = ptail(pf[lw], mat.sum(1)); M[w[lw], oo] = mat; I[w[lw], oo] = inten
    if g % 30 == 0: print(RAW, f'isolation window {g}/{len(wins)} ({time.time() - t0:.0f} s)', flush=True)

valid = ~np.isnan(FM[wp])
out = np.zeros((n, 5), np.float32)
out[:, 0] = L[:, 2:5].max(1); out[:, 1] = L.sum(1); out[:, 2] = (L >= 2).sum(1)
cm = M[:, 2:5].sum(1) >= 2; qq = np.where(valid, 3 * P ** 2 - 2 * P ** 3, np.nan)
for a in range(0, n, 200000):
    out[a:a + 200000, 3] = ptail(np.clip(qq[a:a + 200000], 1e-9, 0.999), cm[a:a + 200000].sum(1))
ia = np.sqrt(I[:, 2:5].sum(1)); il = np.sqrt(FL[wp])
out[:, 4] = (ia * il).sum(1) / np.maximum(np.linalg.norm(ia, axis=1) * np.linalg.norm(il, axis=1), 1e-12)
np.save(f'{V}/{RAW}/natE.npy', out); np.save(f'{V}/{RAW}/natE_pmean.npy', P.mean(1))
dec = pop.is_decoy.values == 1
print(RAW, 'done', f'{time.time() - t0:.0f} s', 'medians target', np.median(out[~dec], 0).round(3).tolist(), 'decoy', np.median(out[dec], 0).round(3).tolist(),
      'mean p_f target/decoy', float(P[~dec][valid[~dec]].mean()), float(P[dec][valid[dec]].mean()), flush=True)
