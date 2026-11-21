"""IM calibration headroom (diagnostic): current linear 1/K0 line vs per-charge lines, isotonic, and the best monotone
map fitted on DIA-NN's own IDs (floor), all on DIA-NN IDs. Measurement = DIA-NN 1/K0 (OpenDIA writes none per ID)."""
import json, glob, sqlite3, numpy as np, pandas as pd
from sklearn.isotonic import IsotonicRegression
C = '/nfs/research/juan/singj/opendia'
LIB = f'{C}/entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp'
c = sqlite3.connect(f'file:{LIB}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID pid, p.DECOY d, p.CHARGE z, p.PRECURSOR_MZ mz, p.LIBRARY_DRIFT_TIME lim, pe.MODIFIED_SEQUENCE ms,
  min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=pe.ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c).set_index('pid')
tg = pq[pq.d == 0]; sk = (tg.ms + '/' + tg.z.astype(str)).reset_index().drop_duplicates(0).set_index(0).pid
for run in ('D', 'H'):
    cal = json.load(open(glob.glob(f'{C}/entlib-mt-20261120/runs/{run}_mt_s0/*.calibration.json')[0]))['im']
    d = pd.read_parquet(f'{C}/entlib-diann5-20261119/{run}/report.parquet', columns=['Modified.Sequence', 'Precursor.Charge', 'Q.Value', 'Decoy', 'IM'])
    d = d[(d['Q.Value'] <= 0.01) & (d.Decoy == 0)].copy(); d['pid'] = sk.reindex(d['Modified.Sequence'] + '/' + d['Precursor.Charge'].astype(str)).values
    d = d.dropna(subset=['pid']); d['pid'] = d.pid.astype(np.int64); d = d[tg.e.reindex(d.pid).values == 0].drop_duplicates('pid').set_index('pid')
    lim = tg.lim.reindex(d.index).values; z = tg.z.reindex(d.index).values; y = d.IM.values; ok = np.isfinite(lim) & np.isfinite(y); lim, z, y = lim[ok], z[ok], y[ok]
    p95 = lambda r: np.quantile(np.abs(r), .95)
    lin = y - (cal['intercept_vssc'] + cal['slope'] * lim)
    A = np.vstack([lim, np.ones_like(lim)]).T; refit = y - A @ np.linalg.lstsq(A, y, rcond=None)[0]
    pc = np.zeros_like(y)
    for zz in np.unique(z):
        m = z == zz; Az = np.vstack([lim[m], np.ones(m.sum())]).T; pc[m] = y[m] - Az @ np.linalg.lstsq(Az, y[m], rcond=None)[0]
    o = np.argsort(lim); iso = y - IsotonicRegression(out_of_bounds='clip').fit(lim[o], y[o]).predict(lim)
    q2 = y - np.polyval(np.polyfit(lim, y, 2), lim)
    print(f'{run}: n {len(y)}, |1/K0 residual| p95: run line {p95(lin):.4f} | line refit on DIA-NN {p95(refit):.4f} | per-charge lines {p95(pc):.4f} | quadratic {p95(q2):.4f} | isotonic floor {p95(iso):.4f} | current tolerance {0.6 * cal["tolerance_half_width_vssc"]:.4f}')
