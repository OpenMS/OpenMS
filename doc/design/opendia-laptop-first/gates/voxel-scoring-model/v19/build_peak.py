"""Native-sampling peak geometry for every population window (label-free, model-free, same rule for targets and
decoys). Replaces the three search-run summaries that measure a 30-s-gap signal run instead of the peak
(transition_support, log_matched_scans, log_peak_width; median 46-66 s against DIA-NN FWHM 3-4 s).
Per window (precursor p, search apex t): the MS2 scans of p's isolation window, scan nearest t and 4 on each side
(Astral 2 Th, revisit ~1.65 s). Fragments: p's 6 search-library fragments (library.pqp, decoys their own m/z),
calibrated with the run's ms2 offset, tolerance = calibrated half width (as native_features.py).
Trace = per-scan sum of matched fragment intensities. Native apex = max of the trace within +-2 scans of t.
Columns (peak.npy [n, 3] float32, population.tsv row order):
  0 nat_support     fragments matched in >= 2 of the 3 scans centred on the native apex (0..6)  -> replaces col 0
  1 log1p(nat_pts)  contiguous scans above half maximum around the native apex                -> replaces col 1
  2 log1p(nat_fwhm) FWHM in s (linear half-maximum crossings; trace edge if never below half)  -> replaces col 2
Usage: build_peak.py <raw01|raw02>"""
import json, sqlite3, sys, time, numpy as np, pandas as pd, pyopenms as oms
from peak_core import local_peak, support
RAW = sys.argv[1]
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; V = B + '/astral-voxel-pilot-20261110/v2'
CAL = {'raw01': B + '/astral-design-raw01-20261107/window300k/cal.json', 'raw02': B + '/astral-entrapment-raw02-20261108/cal/cal.json'}[RAW]
MZML = B + f'/dd-plots-20261111/{RAW}.mzML'
cal = json.load(open(CAL))['ms2_mz']; OFF, TOL = float(cal['offset_ppm']), float(cal['tolerance_half_width_ppm'])
K, HALFW = 6, 4; NOFF = 2 * HALFW + 1
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
pmz = pd.read_sql('select ID precursor_id, PRECURSOR_MZ pmz from PRECURSOR', c).set_index('precursor_id').pmz
tr = pd.read_sql("""select m.PRECURSOR_ID precursor_id, t.PRODUCT_MZ mz from TRANSITION_PRECURSOR_MAPPING m
  join TRANSITION t on t.ID=m.TRANSITION_ID""", c).sort_values(['precursor_id', 'mz'])
tr['k'] = tr.groupby('precursor_id').cumcount(); assert tr.k.max() < K
pop = pd.read_csv(f'{V}/{RAW}/population.tsv', sep='\t'); n = len(pop)
pids = np.unique(pop.precursor_id.values); pix = pd.Series(np.arange(len(pids)), index=pids)
FM = np.full((len(pids), K), np.nan)
t_ = tr[tr.precursor_id.isin(pix.index)]; FM[pix[t_.precursor_id].values, t_.k.values] = t_.mz.values * (1 + OFF * 1e-6)
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
print(RAW, 'MS2 scans', len(scans), 'isolation windows', len(wins), 'windows without isolation window', int((~ok).sum()), flush=True)
M = np.zeros((n, NOFF, K), bool); I = np.zeros((n, NOFF), np.float64); T = np.full((n, NOFF), np.nan)
t0 = time.time()
for g in range(len(wins)):
    sc = scans[(scans.lo == wins.lo[g]) & (scans.hi == wins.hi[g])]; srt = sc.rt.values; sidx = sc.idx.values
    w = np.where(ok & (wi == g))[0]
    if len(w) == 0 or len(sc) < 2: continue
    j0 = np.clip(np.searchsorted(srt, apex[w]), 1, len(srt) - 1); j0 -= (apex[w] - srt[j0 - 1]) < (srt[j0] - apex[w])
    pw = np.repeat(w, NOFF); po = np.tile(np.arange(NOFF), len(w)); pj = np.repeat(j0, NOFF) + po - HALFW
    keep = (pj >= 0) & (pj < len(srt)); pw, po, pj = pw[keep], po[keep], pj[keep]
    T[pw, po] = srt[pj]
    order = np.argsort(pj, kind='stable'); pw, po, pj = pw[order], po[order], pj[order]
    bounds = np.r_[0, np.flatnonzero(np.diff(pj)) + 1, len(pj)]
    for a, b in zip(bounds[:-1], bounds[1:]):
        mz, it = od.getSpectrum(int(sidx[pj[a]])).get_peaks(); mz = np.asarray(mz, float); it = np.asarray(it, float)
        if len(mz) == 0: continue
        ww, oo = pw[a:b], po[a:b]; F = FM[wp[ww]]; valid = ~np.isnan(F); Fz = np.where(valid, F, 0.0); tolmz = Fz * TOL * 1e-6
        j = np.clip(np.searchsorted(mz, Fz), 1, len(mz) - 1); jn = np.where(np.abs(mz[j - 1] - Fz) < np.abs(mz[j] - Fz), j - 1, j)
        mat = (np.abs(mz[jn] - Fz) <= tolmz) & valid
        M[ww, oo] = mat; I[ww, oo] = np.where(mat, it[jn], 0.0).sum(1)
    if g % 50 == 0: print(RAW, f'isolation window {g}/{len(wins)} ({time.time() - t0:.0f} s)', flush=True)
out = np.zeros((n, 3), np.float32)
for r in range(n):
    m = ~np.isnan(T[r])
    if m.sum() < 2: continue
    a, fw, pts = local_peak(T[r, m], I[r, m], apex[r], search=2)
    if a < 0: continue
    out[r] = (support(M[r][m], a), np.log1p(pts), np.log1p(fw))
assert np.isfinite(out).all()
np.save(f'{V}/{RAW}/peak.npy', out)
dec = pop.is_decoy.values == 1
print(RAW, 'done', f'{time.time() - t0:.0f} s', 'medians target', np.median(out[~dec], 0).round(3).tolist(), 'decoy', np.median(out[dec], 0).round(3).tolist(), flush=True)
