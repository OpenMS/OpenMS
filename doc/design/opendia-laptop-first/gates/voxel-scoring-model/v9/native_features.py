"""PLAN_V9 features N: native-sampling fragment evidence with a spectrum-density chance model (label-free, model-free,
same rule for targets and decoys). Usage: native_features.py <raw01|raw02>.
For each population window (precursor p, apex t): the MS2 scans of the isolation window holding p's m/z (Astral 2 Th,
revisit about 1.65 s), the scan nearest t and 3 on each side (offsets -3..3). Fragments: p's 12-fragment library
(library12_pop.pqp; decoys carry their own decoy m/z), calibrated with the run's ms2 offset; tolerance = the run's
calibrated ms2 half width (cal.json ms2_mz.tolerance_half_width_ppm).
Per scan and fragment f: matched = a peak within tolerance (intensity = nearest peak's); chance rate from the local
density: lam_f = (peaks within +-5 Th of f outside the tolerance) * (2 tol_f) / (10 Th - 2 tol_f), p_f = 1 - exp(-lam_f),
clipped to [1e-3, 0.999]. Per scan: L = -log10 P(X >= k), X ~ Poisson-binomial(p_f over the precursor's fragments),
k = matched count (L capped at 30).
Features (5 columns, nat.npy [n, 5] float32, row order of population.tsv):
  0 nat_apex_logp    max L over offsets -1..1
  1 nat_sum_logp     sum of L over offsets -3..3
  2 nat_sig_scans    scans (offsets -3..3) with L >= 2
  3 nat_persist_logp -log10 P(Y >= k2), k2 = fragments matched in >= 2 of the 3 central scans, Y ~ Poisson-binomial with
                     q_f = 3 pbar^2 - 2 pbar^3 (pbar = mean p_f over the 3 central scans)
  4 nat_cos          cosine of sqrt(summed central-3 intensities) vs sqrt(library intensities) over the 12 fragments
Windows whose precursor m/z is in no isolation window get zeros (none expected; counted in the log)."""
import json, sqlite3, sys, time, numpy as np, pandas as pd, pyopenms as oms
RAW = sys.argv[1]
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; V = B + '/astral-voxel-pilot-20261110/v2'
CAL = {'raw01': B + '/astral-design-raw01-20261107/window300k/cal.json', 'raw02': B + '/astral-entrapment-raw02-20261108/cal/cal.json'}[RAW]
MZML = B + f'/dd-plots-20261111/{RAW}.mzML'
cal = json.load(open(CAL))['ms2_mz']; OFF, TOL = float(cal['offset_ppm']), float(cal['tolerance_half_width_ppm'])
K, NOFF, HALF_TH, PMIN, PMAX, LCAP = 12, 7, 5.0, 1e-3, 0.999, 30.0
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

# Scan index (MS2 only) per isolation window.
od = oms.OnDiscMSExperiment(); od.openFile(MZML, False); meta = od.getMetaData(); rows = []
for i in range(meta.getNrSpectra()):
    sp = meta.getSpectrum(i)
    if sp.getMSLevel() != 2 or not sp.getPrecursors(): continue
    p = sp.getPrecursors()[0]; rows.append((i, sp.getRT(), round(p.getMZ() - p.getIsolationWindowLowerOffset(), 3), round(p.getMZ() + p.getIsolationWindowUpperOffset(), 3)))
scans = pd.DataFrame(rows, columns=['idx', 'rt', 'lo', 'hi']).sort_values('rt')
wins = scans[['lo', 'hi']].drop_duplicates().sort_values('lo').reset_index(drop=True)
print(RAW, 'MS2 scans', len(scans), 'isolation windows', len(wins), 'tol ppm', TOL, 'offset ppm', OFF, flush=True)
# Assign each population window to the isolation window holding its precursor m/z.
wi = np.searchsorted(wins.lo.values, wmz, side='right') - 1
ok = (wi >= 0) & (wmz <= wins.hi.values[np.clip(wi, 0, None)])
print(RAW, 'windows without an isolation window', int((~ok).sum()), flush=True)

L = np.zeros((n, NOFF), np.float32); M = np.zeros((n, NOFF, K), np.uint8); I = np.zeros((n, NOFF, K), np.float32); P = np.zeros((n, NOFF, K), np.float32)

def ptail(p, k):
    """-log10 P(X >= k), X Poisson-binomial with rows of p (nan = absent fragment)."""
    m = p.shape[0]; dp = np.zeros((m, K + 1)); dp[:, 0] = 1
    for f in range(K):
        q = np.nan_to_num(p[:, f]); nd = dp * (1 - q)[:, None]; nd[:, 1:] += dp[:, :-1] * q[:, None]; dp = nd
    cs = np.cumsum(dp[:, ::-1], 1)[:, ::-1]  # cs[:, j] = P(X >= j)
    return np.minimum(-np.log10(np.maximum(cs[np.arange(m), k], 1e-300)), LCAP)

t0 = time.time()
for g in range(len(wins)):
    sc = scans[(scans.lo == wins.lo[g]) & (scans.hi == wins.hi[g])]; srt = sc.rt.values; sidx = sc.idx.values
    w = np.where(ok & (wi == g))[0]
    if len(w) == 0 or len(sc) == 0: continue
    j0 = np.clip(np.searchsorted(srt, apex[w]), 1, len(srt) - 1); j0 -= (apex[w] - srt[j0 - 1]) < (srt[j0] - apex[w])
    pw = np.repeat(w, NOFF); po = np.tile(np.arange(NOFF), len(w)); pj = np.repeat(j0, NOFF) + po - 3
    keep = (pj >= 0) & (pj < len(srt)); pw, po, pj = pw[keep], po[keep], pj[keep]
    order = np.argsort(pj, kind='stable'); pw, po, pj = pw[order], po[order], pj[order]
    bounds = np.r_[0, np.flatnonzero(np.diff(pj)) + 1, len(pj)]
    for a, b in zip(bounds[:-1], bounds[1:]):
        mz, it = od.getSpectrum(int(sidx[pj[a]])).get_peaks(); mz = np.asarray(mz, float); it = np.asarray(it, float)
        ww, oo = pw[a:b], po[a:b]; F = FM[wp[ww]]; valid = ~np.isnan(F); Fz = np.where(valid, F, 0.0)
        tolmz = Fz * TOL * 1e-6
        if len(mz) == 0:
            lam = np.zeros_like(Fz); mat = np.zeros_like(Fz, bool); inten = np.zeros_like(Fz)
        else:
            j = np.clip(np.searchsorted(mz, Fz), 1, len(mz) - 1); jn = np.where(np.abs(mz[j - 1] - Fz) < np.abs(mz[j] - Fz), j - 1, j)
            mat = (np.abs(mz[jn] - Fz) <= tolmz) & valid; inten = np.where(mat, it[jn], 0.0)
            ntol = np.searchsorted(mz, Fz + tolmz, 'right') - np.searchsorted(mz, Fz - tolmz, 'left')
            nloc = np.searchsorted(mz, Fz + HALF_TH, 'right') - np.searchsorted(mz, Fz - HALF_TH, 'left') - ntol
            lam = nloc * (2 * tolmz) / (2 * HALF_TH - 2 * tolmz)
        p = np.clip(1 - np.exp(-lam), PMIN, PMAX); p = np.where(valid, p, np.nan)
        L[ww, oo] = ptail(p, mat.sum(1)); M[ww, oo] = mat; I[ww, oo] = inten; P[ww, oo] = np.nan_to_num(p)
    if g % 30 == 0: print(RAW, f'isolation window {g}/{len(wins)} ({time.time() - t0:.0f} s)', flush=True)

valid = ~np.isnan(FM[wp])
out = np.zeros((n, 5), np.float32)
out[:, 0] = L[:, 2:5].max(1); out[:, 1] = L.sum(1); out[:, 2] = (L >= 2).sum(1)
cm = M[:, 2:5].sum(1) >= 2; pbar = P[:, 2:5].mean(1); q = np.where(valid, 3 * pbar ** 2 - 2 * pbar ** 3, np.nan)
for a in range(0, n, 200000):
    out[a:a + 200000, 3] = ptail(np.clip(q[a:a + 200000], 1e-9, 0.999), cm[a:a + 200000].sum(1))
ia = np.sqrt(I[:, 2:5].sum(1)); il = np.sqrt(FL[wp])
out[:, 4] = (ia * il).sum(1) / np.maximum(np.linalg.norm(ia, axis=1) * np.linalg.norm(il, axis=1), 1e-12)
np.save(f'{V}/{RAW}/nat.npy', out)
dec = pop.is_decoy.values == 1
print(RAW, 'done', f'{time.time() - t0:.0f} s', 'medians target', np.median(out[~dec], 0).round(3).tolist(), 'decoy', np.median(out[dec], 0).round(3).tolist(), flush=True)
