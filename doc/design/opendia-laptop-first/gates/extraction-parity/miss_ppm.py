"""For a sample of DIA-NN nonzero points that our extraction has as zero, find the raw MS2 peak at that scan near the
fragment m/z and report its ppm offset from our calibrated query m/z (our tolerance is +-5.017 ppm around it)."""
import numpy as np, pandas as pd, pyopenms as oms
X = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/extraction-parity-20261112'
m = pd.read_parquet(X + '/parity_points_raw01.parquet')
miss = m[(m.value > 0) & (m.intensity == 0)].sample(3000, random_state=1)
hit = m[(m.value > 0) & (m.intensity > 0)].sample(1000, random_state=1)
exp = oms.OnDiscMSExperiment(); exp.openFile('/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/dd-plots-20261111/raw01.mzML')
meta = exp.getMetaData()
rts = np.array([s.getRT() for s in meta.getSpectra()]); lv = np.array([s.getMSLevel() for s in meta.getSpectra()])
prec = np.array([s.getPrecursors()[0].getMZ() if s.getPrecursors() else 0 for s in meta.getSpectra()])
def lookup(df, tag):
    rows = []
    for r in df.itertuples():
        # the scan of the isolation window holding this precursor at this RT
        idx = np.where((np.abs(rts - r.rt) < 1e-3) & (lv == 2))[0]
        best = None
        for i in idx:
            s = exp.getSpectrum(int(i)); mz, it = s.get_peaks()
            if len(mz) == 0: continue
            sel = np.abs(mz - r.fragment_mz) / r.fragment_mz * 1e6 < 25
            for a, b in zip(mz[sel], it[sel]):
                rel = abs(b - r.value) / r.value
                if best is None or rel < best[2]:
                    best = ((a - r.fragment_mz) / r.fragment_mz * 1e6, b, rel, len(idx))
        rows.append((tag,) + (best if best else (np.nan, np.nan, np.nan, len(idx))))
    return pd.DataFrame(rows, columns=['tag', 'ppm_vs_query', 'raw_int', 'rel_int_diff', 'scans_at_rt'])
res = pd.concat([lookup(miss, 'miss'), lookup(hit, 'hit')])
res.to_csv(X + '/miss_ppm_raw01.tsv', sep='\t', index=False)
with open(X + '/MISS_PPM_raw01.txt', 'w') as f:
    for tag, s in res.groupby('tag'):
        ok = s[s.rel_int_diff < 1e-3]
        line = (f'{tag}: n={len(s)} raw peak with DIA-NN intensity found {len(ok) / len(s):.3f}; ppm vs our calibrated query '
                f'{ok.ppm_vs_query.quantile([.01, .1, .5, .9, .99]).round(2).to_dict()}; |ppm|>5.017 {(ok.ppm_vs_query.abs() > 5.017).mean():.3f}; '
                f'scans at that RT {s.scans_at_rt.value_counts().to_dict()}')
        print(line); f.write(line + '\n')
    ok = res[(res.tag == 'miss') & (res.rel_int_diff < 1e-3)]
    h = np.histogram(ok.ppm_vs_query, bins=[-25, -15, -11, -8, -5.017, 0, 5.017, 8, 11, 15, 25])
    line = 'miss ppm histogram ' + str(dict(zip([f'{a:g}..{b:g}' for a, b in zip(h[1][:-1], h[1][1:])], h[0].tolist())))
    print(line); f.write(line + '\n')
