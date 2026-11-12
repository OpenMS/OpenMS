"""Where real signal sits in m/z: raw ppm (vs our calibrated query) of the DIA-NN-recorded peak at the 3 scans around the
DIA-NN apex of q<=1% IDs, versus points >= 30 s from the apex (background). Also: which library fragments DIA-NN leaves
out of its export."""
import numpy as np, pandas as pd, pyopenms as oms, sqlite3
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; X = B + '/extraction-parity-20261112'
m = pd.read_parquet(X + '/parity_points_raw01.parquet')
rep = pd.read_parquet(X + '/diann/report_raw01.parquet', columns=['Precursor.Id', 'Q.Value', 'RT']).set_index('Precursor.Id')
m['d'] = (m.rt - rep.RT.reindex(m.pr).values * 60).abs(); m['q'] = rep['Q.Value'].reindex(m.pr).values
apex = m[(m.value > 0) & (m.q <= 0.01) & (m.d <= 2.5)].sample(2500, random_state=2)
tail = m[(m.value > 0) & (m.q <= 0.01) & (m.d >= 30)].sample(2500, random_state=2)
exp = oms.OnDiscMSExperiment(); exp.openFile(B + '/dd-plots-20261111/raw01.mzML')
meta = exp.getMetaData()
rts = np.array([s.getRT() for s in meta.getSpectra()]); lv = np.array([s.getMSLevel() for s in meta.getSpectra()])
def ppm_of(df):
    out = []
    for r in df.itertuples():
        i = np.where((np.abs(rts - r.rt) < 1e-3) & (lv == 2))[0]
        s = exp.getSpectrum(int(i[0])); mz, it = s.get_peaks()
        sel = np.where(np.abs(mz - r.fragment_mz) / r.fragment_mz * 1e6 < 25)[0]
        j = sel[np.argmin(np.abs(it[sel] - r.value))] if len(sel) else None
        out.append((mz[j] - r.fragment_mz) / r.fragment_mz * 1e6 if j is not None else np.nan)
    return np.array(out)
bins = [-8.5, -7, -6, -5.017, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5.017, 6, 7, 8.5]
with open(X + '/SIGNAL_PPM_raw01.txt', 'w') as f:
    for name, df in [('apex+-2.5s', apex), ('tail>=30s', tail)]:
        p = ppm_of(df)
        h = np.histogram(p, bins=bins)[0] / np.isfinite(p).sum()
        line = (f'{name}: n={len(p)} share |ppm|>5.017 {np.mean(np.abs(p) > 5.017):.3f}; median {np.nanmedian(p):.2f}; '
                f'MAD-sigma {1.4826 * np.nanmedian(np.abs(p - np.nanmedian(p))):.2f}; hist ' +
                ' '.join(f'[{a:g},{b:g}):{v:.3f}' for a, b, v in zip(bins[:-1], bins[1:], h)))
        print(line); f.write(line + '\n')

    # fragments DIA-NN leaves out of the export: by ion type/number and m/z
    c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
    idm = pd.read_csv(X + '/idmap.tsv', sep='\t')
    tr = pd.read_sql("""select tp.PRECURSOR_ID precursor_id, t.PRODUCT_MZ product_mz, t.ANNOTATION annotation, t.LIBRARY_INTENSITY li
      from TRANSITION t join TRANSITION_PRECURSOR_MAPPING tp on tp.TRANSITION_ID=t.ID""", c).merge(idm, on='precursor_id')
    dn = pd.read_parquet(X + '/diann/report_raw01_xic/raw01.xic.parquet', columns=['pr', 'feature']).drop_duplicates()
    k = set(zip(dn.pr, dn.feature)); tr['kept'] = [(p, a) in k for p, a in zip(tr['Precursor.Id'], tr.annotation)]
    tr['ion'] = tr.annotation.str[0]; tr['num'] = tr.annotation.str.extract(r'^[a-z](\d+)')[0].astype(int)
    tr['z'] = tr.annotation.str.extract(r'\^(\d+)')[0].astype(int)
    for col in ['ion', 'z']:
        line = f'kept by {col}: ' + str(tr.groupby(col).kept.mean().round(3).to_dict()); print(line); f.write(line + '\n')
    tr['mz'] = pd.cut(tr.product_mz, [0, 250, 300, 350, 400, 2100])
    t = tr.pivot_table(index='num', columns='mz', values='kept', aggfunc='mean', observed=True).round(2).head(8)
    print(t.to_string()); f.write(t.to_string() + '\n')
