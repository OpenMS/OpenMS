"""Write run-calibrated RT (seconds) into a per-run copy of merged.pqp, so OSW extracts at the same centres the reference
analysis uses (the OpenDIACalibrate mapping, via the search's rt_apex - apex_rt_error_seconds). That calibrated RT is a
monotone function of library RT (checked: spread 0.001 s at equal library RT), so it is applied as a 1-D
interpolation to every precursor, targets and decoys alike (decoys carry their target's library RT).
Then OSW runs with calibration off, so its RT window is centred on this value."""
import os, shutil, sqlite3, numpy as np, pandas as pd
RUN = os.environ['OPENDIA_RUN']
HERE = os.path.dirname(os.path.abspath(__file__)); B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
f = B + '/astral-design-raw01-20261107/pipeline/search.features.tsv' if RUN == '01' else B + '/astral-entrapment-raw02-20261108/search.features.tsv'
lib = pd.read_sql('select ID precursor_id, LIBRARY_RT lrt from PRECURSOR where DECOY=0', sqlite3.connect(B + '/astral-entrapment-20261107/lib/library.pqp'))
s = pd.read_csv(f, sep='\t', usecols=['precursor_id', 'rt_apex', 'apex_rt_error_seconds'])
cal = (s.rt_apex - s.apex_rt_error_seconds).groupby(s.precursor_id).first().rename('cal')
d = lib.merge(cal, left_on='precursor_id', right_index=True).groupby('lrt').cal.median().sort_index()
x, y = d.index.values, d.values
out = f'{HERE}/merged_cal{RUN}.pqp'; shutil.copy(f'{HERE}/merged.pqp', out)
c = sqlite3.connect(out)
p = pd.read_sql('select ID, LIBRARY_RT from PRECURSOR', c)
v = np.interp(p.LIBRARY_RT.values, x, y)
lo, hi = p.LIBRARY_RT.values < x[0], p.LIBRARY_RT.values > x[-1]
v[lo] = y[0] + (p.LIBRARY_RT.values[lo] - x[0]) * (y[1] - y[0]) / max(x[1] - x[0], 1e-9)
v[hi] = y[-1] + (p.LIBRARY_RT.values[hi] - x[-1]) * (y[-1] - y[-2]) / max(x[-1] - x[-2], 1e-9)
c.executemany('update PRECURSOR set LIBRARY_RT=? where ID=?', list(zip(v.tolist(), p.ID.tolist())))
c.commit()
print(RUN, 'precursors', len(p), 'calibrated RT range', round(v.min(), 1), round(v.max(), 1), 'extrapolated', int(lo.sum() + hi.sum()))
