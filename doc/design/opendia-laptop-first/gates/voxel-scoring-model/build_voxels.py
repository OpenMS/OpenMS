"""Voxel pilot step 2: per-window voxel arrays from the raw 01 .xipm export (export01/pilot.xipm).
For every training window (train_windows.tsv) and evaluation window (eval_windows.tsv), RT grid of 32 bins over the
candidate apex +- 15 s (the base tensor's rt_half_width), fragments in library order (intensity descending, as the
base tensor), and these additive arrays (so synthetic mixing in train_pilot.py is a plain sum):
  I   [6, 7, 32]  fragment intensity in 7 m/z sub-bins over the calibrated offset +- 17.5 ppm (5 ppm each)
  P   [6, 32]     sum of intensity x (ppm - calibrated offset) over the same points
  M   [4, 32]     MS1 M+0..M+3 intensity within +- 10 ppm of the calibrated offset (MS1 offset taken as the MS2 one)
  lib [6]         library intensities / max
Deviation from today's tensor (recorded in DESIGN.md before training): m/z is cut at the export's +- 20 ppm, not
+- 0.02 Da, so option a is an approximation of today's tensor. All cells share it."""
import json, sqlite3, sys, numpy as np, pandas as pd, pyarrow.parquet as pq, zlib, pynumpress
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
O = B + '/astral-voxel-pilot-20261110'
RT_BINS, HALF, SUB, SUB_PPM, MS1_PPM = 32, 15.0, 7, 5.0, 10.0
cal = json.load(open(B + '/astral-design-raw01-20261107/window300k/cal.json'))
OFF = float(cal['ms2_mz']['offset_ppm'])
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
tr = pd.read_sql("""select m.PRECURSOR_ID precursor_id, t.PRODUCT_MZ mz, t.LIBRARY_INTENSITY li from TRANSITION_PRECURSOR_MAPPING m
  join TRANSITION t on t.ID=m.TRANSITION_ID""", c)
win = pd.concat([pd.read_csv(O + '/train_windows.tsv', sep='\t', usecols=['sample_id', 'precursor_id', 'rt_apex']),
                 pd.read_csv(O + '/eval_windows.tsv', sep='\t', usecols=['sample_id', 'precursor_id', 'rt_apex'])]).drop_duplicates('sample_id')
need = set(win.precursor_id)
tr = tr[tr.precursor_id.isin(need)].sort_values(['precursor_id', 'li', 'mz'], ascending=[True, False, True])
frag = {pid: (g.mz.values, g.li.values / max(g.li.values.max(), 1e-12)) for pid, g in tr.groupby('precursor_id')}
by_pid = {pid: g for pid, g in win.groupby('precursor_id')}
ids = win.sample_id.values; pos = {s: i for i, s in enumerate(ids)}; n = len(ids)
I = np.zeros((n, 6, SUB, RT_BINS), np.float32); P = np.zeros((n, 6, RT_BINS), np.float32)
M = np.zeros((n, 4, RT_BINS), np.float32); LIB = np.zeros((n, 6), np.float32); seen = np.zeros(n, bool)

def dec(blob, kind):
    if blob is None or len(blob) == 0: return np.zeros(0)
    raw = zlib.decompress(blob)
    if kind == 5 and len(raw) == 12:
        # One-point linear array: 8-byte big-endian fixed point, then the value as int32 little-endian.
        # pynumpress.decode_linear rejects this case.
        fp = np.frombuffer(raw[:8], '>f8')[0]
        return np.array([np.frombuffer(raw[8:12], '<i4')[0] / fp])
    return np.asarray(pynumpress.decode_linear(raw) if kind == 5 else pynumpress.decode_slof(raw), float)

f = pq.ParquetFile(O + '/export01/pilot.xipm')
cols = ['MS_LEVEL', 'PRECURSOR_ID', 'ANNOTATION', 'TARGET_MZ', 'MZ_DATA', 'RT_DATA', 'INTENSITY_DATA', 'MZ_COMPRESSION', 'RT_COMPRESSION', 'INTENSITY_COMPRESSION']
rows = 0
for rg in range(f.num_row_groups):
    t = f.read_row_group(rg, columns=cols).to_pydict()
    for k in range(len(t['MS_LEVEL'])):
        pid = int(t['PRECURSOR_ID'][k]); rows += 1
        if pid not in by_pid: continue
        mz = dec(t['MZ_DATA'][k], t['MZ_COMPRESSION'][k]); rt = dec(t['RT_DATA'][k], t['RT_COMPRESSION'][k])
        it = dec(t['INTENSITY_DATA'][k], t['INTENSITY_COMPRESSION'][k])
        if len(mz) == 0: continue
        target = t['TARGET_MZ'][k]
        # MS2 targets in the export are already calibrated (library m/z shifted by the calibrated offset), MS1 targets
        # are library m/z. So the offset is removed only for MS1.
        ppm = (mz - target) / target * 1e6 - (OFF if t['MS_LEVEL'][k] == 1 else 0.0)
        if t['MS_LEVEL'][k] == 1:
            ann = t['ANNOTATION'][k] or ''  # 'Precursor_iK'
            if '_i' not in ann: continue
            iso = int(ann.rsplit('_i', 1)[1])
            if iso > 3: continue
            keep = np.abs(ppm) <= MS1_PPM
            for _, w in by_pid[pid].iterrows():
                j = pos[w.sample_id]; d = rt - w.rt_apex; m = keep & (np.abs(d) < HALF)
                b = np.clip(((d[m] + HALF) / (2 * HALF) * RT_BINS).astype(int), 0, RT_BINS - 1)
                np.add.at(M[j, iso], b, it[m])
            continue
        mzs, li = frag[pid]
        lib_target = target / (1 + OFF * 1e-6)
        fi = int(np.argmin(np.abs(mzs - lib_target)))
        if abs(mzs[fi] - lib_target) / lib_target * 1e6 > 0.5: continue
        sb = np.floor((ppm + SUB * SUB_PPM / 2) / SUB_PPM).astype(int); keep = (sb >= 0) & (sb < SUB)
        for _, w in by_pid[pid].iterrows():
            j = pos[w.sample_id]; LIB[j] = li; seen[j] = True
            d = rt - w.rt_apex; m = keep & (np.abs(d) < HALF)
            b = np.clip(((d[m] + HALF) / (2 * HALF) * RT_BINS).astype(int), 0, RT_BINS - 1)
            np.add.at(I[j, fi], (sb[m], b), it[m]); np.add.at(P[j, fi], b, it[m] * ppm[m])
    if rg % 50 == 0: print('row group', rg, '/', f.num_row_groups, flush=True)
np.savez_compressed(O + '/voxels01.npz', sample_id=ids, I=I, P=P, M=M, LIB=LIB, seen=seen)
print('peak maps', rows, 'windows', n, 'with fragments', int(seen.sum()), 'nonzero I', float((I.sum((1, 2, 3)) > 0).mean()),
      'nonzero MS1', float((M.sum((1, 2)) > 0).mean()), 'offset_ppm', OFF)
