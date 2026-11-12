"""Raw ppm (vs our calibrated query) of DIA-NN-recorded points within 2.5 s of the apex of q<=1% IDs, by intensity decile."""
import numpy as np, pandas as pd, pyopenms as oms
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; X = B + '/extraction-parity-20261112'
m = pd.read_parquet(X + '/parity_points_raw01.parquet')
rep = pd.read_parquet(X + '/diann/report_raw01.parquet', columns=['Precursor.Id', 'Q.Value', 'RT']).set_index('Precursor.Id')
m['d'] = (m.rt - rep.RT.reindex(m.pr).values * 60).abs(); m['q'] = rep['Q.Value'].reindex(m.pr).values
a = m[(m.value > 0) & (m.q <= 0.01) & (m.d <= 2.5)].sample(6000, random_state=3)
exp = oms.OnDiscMSExperiment(); exp.openFile(B + '/dd-plots-20261111/raw01.mzML')
meta = exp.getMetaData(); rts = np.array([s.getRT() for s in meta.getSpectra()]); lv = np.array([s.getMSLevel() for s in meta.getSpectra()])
pp = []
for r in a.itertuples():
    i = np.where((np.abs(rts - r.rt) < 1e-3) & (lv == 2))[0]; s = exp.getSpectrum(int(i[0])); mz, it = s.get_peaks()
    sel = np.where(np.abs(mz - r.fragment_mz) / r.fragment_mz * 1e6 < 25)[0]
    j = sel[np.argmin(np.abs(it[sel] - r.value))]; pp.append((mz[j] - r.fragment_mz) / r.fragment_mz * 1e6)
a['ppm'] = pp; a['dec'] = pd.qcut(a.value, 10, labels=False)
def rs(x): return 1.4826 * np.median(np.abs(x - np.median(x)))
t = a.groupby('dec').agg(int_med=('value', 'median'), n=('ppm', 'size'), median_ppm=('ppm', 'median'), robust_sigma=('ppm', rs),
                         beyond_5ppm=('ppm', lambda x: (x.abs() > 5.017).mean()))
a['mzb'] = pd.cut(a.fragment_mz, [0, 400, 600, 800, 1000, 2100])
t2 = a.groupby('mzb', observed=True).agg(n=('ppm', 'size'), median_ppm=('ppm', 'median'), robust_sigma=('ppm', rs), beyond_5ppm=('ppm', lambda x: (x.abs() > 5.017).mean()))
with open(X + '/PPM_VS_INTENSITY_raw01.txt', 'w') as f:
    for s in [t.round(3).to_string(), t2.round(3).to_string(), f'overall robust sigma {rs(a.ppm):.2f} ppm, beyond 5.017 {(a.ppm.abs() > 5.017).mean():.3f}']:
        print(s); f.write(s + '\n')
