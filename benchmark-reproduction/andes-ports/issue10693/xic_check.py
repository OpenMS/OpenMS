"""Is there MS1 signal where PIP-ECHO found no candidate? For precursor x run pairs of one run (the alignment
reference, where the predicted RT is the run's own RT), extract the M0 and M+1 traces (+-10 ppm) within the local
RT window around the predicted RT (theoretical m/z of the precursor), smooth both over 9 MS1 scans, and call the
signal present when, within +-15 scans of the smoothed M0 apex, >= 10 scans carry both isotopes and the smoothed traces
correlate (Pearson >= 0.6). The test is calibrated on the accepted transfers.
Usage: xic_check.py run.mzML cand_pairs.tsv.gz acceptor_map n_per_status seed [mass_shift_Da]"""
import sys
import numpy as np, pandas as pd
import pyopenms as oms
mzml, pairs, amap, n, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
# optional null: shift the precursor mass by this many Da (divided by the charge), like PIP-ECHO's decoy donors
shift_da = float(sys.argv[6]) if len(sys.argv) > 6 else 0.0
g = pd.read_csv(pairs, sep='\t')
g = g[g.acceptor_map == amap]
sample = pd.concat([x.sample(min(n, len(x)), random_state=seed) for _, x in g.groupby('status')])
opts = oms.PeakFileOptions(); opts.setMSLevels([1])
f = oms.MzMLFile(); f.setOptions(opts); exp = oms.MSExperiment(); f.load(mzml, exp)
rts = np.array([s.getRT() for s in exp.getSpectra()])
peaks = [s.get_peaks() for s in exp.getSpectra()]
C13 = 1.0033548378
def trace(mz, lo, hi):
    i0, i1 = np.searchsorted(rts, lo), np.searchsorted(rts, hi, 'right')
    out = []
    for i in range(i0, i1):
        m, it = peaks[i]
        a, b = np.searchsorted(m, mz * (1 - 10e-6)), np.searchsorted(m, mz * (1 + 10e-6), 'right')
        out.append(it[a:b].sum() if b > a else 0.0)
    return np.array(out)
def theo(key):
    seq, z = key.rsplit('_', 1)
    return oms.AASequence.fromString(seq).getMZ(int(z))
sm = lambda t: np.convolve(t, np.ones(9) / 9, mode='same')
res = []
for r in sample.itertuples():
    lo, hi = r.rt - r.gate, r.rt + r.gate
    mz = theo(r.key) + shift_da / r.charge
    t0 = trace(mz, lo, hi); t1 = trace(mz + C13 / r.charge, lo, hi)
    corr, nboth, present = np.nan, 0, False
    if len(t0) >= 9 and t0.max() > 0:
        s0, s1 = sm(t0), sm(t1)
        k = int(s0.argmax()); w = slice(max(0, k - 15), k + 16)
        nboth = int(((t0[w] > 0) & (t1[w] > 0)).sum())
        if np.std(s0[w]) > 0 and np.std(s1[w]) > 0: corr = float(np.corrcoef(s0[w], s1[w])[0, 1])
        present = nboth >= 10 and corr >= 0.6
    both = np.array([nboth])
    res.append((r.status, r.species if 'species' in sample else '', t0.max() if len(t0) else 0.0, int(both.sum()), corr, present))
x = pd.DataFrame(res, columns=['status', 'species', 'apex', 'scans_m0_m1', 'corr', 'present'])
print(f'run map {amap}, {len(x)} sampled pairs')
print(x.groupby('status').agg(n=('present', 'size'), signal_present_pct=('present', lambda s: round(100 * s.mean(), 1)),
      median_apex=('apex', 'median'), apex_below_1e4_pct=('apex', lambda s: round(100 * (s < 1e4).mean(), 1))).to_string())
pres = x[x.present]
print('\napex intensity where signal is present (quartiles):')
print(pres.groupby('status').apex.describe(percentiles=[.25, .5, .75])[['count', '25%', '50%', '75%']].round(0).to_string())
