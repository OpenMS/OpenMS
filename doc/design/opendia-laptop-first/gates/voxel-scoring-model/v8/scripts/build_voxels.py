"""Voxel v2 step 2: per-window arrays for the v2 population from an .xipm export.
Usage: build_voxels.py <raw01|raw02>. Same geometry as the pilot build_voxels.py (32 RT bins over apex +- 15 s,
7 m/z sub-bins of 5 ppm over the calibrated offset +- 17.5 ppm, MS1 M+0..M+3 within +- 10 ppm), with the two pilot
fixes (one-point numpress arrays; MS2 targets in the export are already calibrated). Stored as float16 in a memmap
(I is summed over m/z sub-bins in v2 because option c is not carried forward; P keeps the ppm sum for option b).
Outputs <raw>/xic.f16 [n,6,32], <raw>/ppm.f16 [n,6,32], <raw>/ms1.f16 [n,4,32], <raw>/lib.f16 [n,6], <raw>/index.tsv."""
import json, sqlite3, sys, zlib, numpy as np, pandas as pd, pyarrow.parquet as pq, pynumpress
RAW = sys.argv[1]
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
O = B + '/rtlin-20261111'  # rtlin rebuild root
CAL = {'raw01': B + '/rtlin-20261111/raw01/cal/cal.json', 'raw02': B + '/rtlin-20261111/raw02/cal/cal.json'}[RAW]
RT_BINS, HALF, SUB, SUB_PPM, MS1_PPM = 32, 15.0, 7, 5.0, 10.0
OFF = float(json.load(open(CAL))['ms2_mz']['offset_ppm'])
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
pop = pd.read_csv(f'{O}/{RAW}/population.tsv', sep='\t')
need = set(pop.precursor_id)
tr = pd.read_sql("""select m.PRECURSOR_ID precursor_id, t.PRODUCT_MZ mz, t.LIBRARY_INTENSITY li from TRANSITION_PRECURSOR_MAPPING m
  join TRANSITION t on t.ID=m.TRANSITION_ID""", c)
tr = tr[tr.precursor_id.isin(need)].sort_values(['precursor_id', 'li', 'mz'], ascending=[True, False, True])
frag = {pid: (g.mz.values, (g.li.values / max(g.li.values.max(), 1e-12)).astype(np.float32)) for pid, g in tr.groupby('precursor_id')}
n = len(pop)
rows_of = {pid: (g.index.values, g.rt_apex.values) for pid, g in pop.groupby('precursor_id')}
mm = lambda name, shape: np.lib.format.open_memmap(f'{O}/{RAW}/{name}.npy', mode='w+', dtype=np.float32, shape=shape)
XIC, PPM, MS1, LIB = mm('xic', (n, 6, RT_BINS)), mm('ppm', (n, 6, RT_BINS)), mm('ms1', (n, 4, RT_BINS)), mm('lib', (n, 6))

def dec(blob, kind):
    if blob is None or len(blob) == 0: return np.zeros(0)
    raw = zlib.decompress(blob)
    if kind == 5 and len(raw) == 12:
        fp = np.frombuffer(raw[:8], '>f8')[0]
        return np.array([np.frombuffer(raw[8:12], '<i4')[0] / fp])
    return np.asarray(pynumpress.decode_linear(raw) if kind == 5 else pynumpress.decode_slof(raw), float)

def bins(rt, apexes):
    """[windows, points] RT bin index (or -1 outside +-HALF) for each window apex."""
    d = rt[None, :] - apexes[:, None]
    b = ((d + HALF) / (2 * HALF) * RT_BINS).astype(int)
    return np.where(np.abs(d) < HALF, np.clip(b, 0, RT_BINS - 1), -1)

f = pq.ParquetFile(f'{O}/export{RAW[-2:]}/v2.xipm')
cols = ['MS_LEVEL', 'PRECURSOR_ID', 'ANNOTATION', 'TARGET_MZ', 'MZ_DATA', 'RT_DATA', 'INTENSITY_DATA', 'MZ_COMPRESSION', 'RT_COMPRESSION', 'INTENSITY_COMPRESSION']
maps = 0
for rg in range(f.num_row_groups):
    t = f.read_row_group(rg, columns=cols).to_pydict()
    for k in range(len(t['MS_LEVEL'])):
        pid = int(t['PRECURSOR_ID'][k])
        if pid not in rows_of: continue
        mz = dec(t['MZ_DATA'][k], t['MZ_COMPRESSION'][k]); rt = dec(t['RT_DATA'][k], t['RT_COMPRESSION'][k])
        it = dec(t['INTENSITY_DATA'][k], t['INTENSITY_COMPRESSION'][k]).astype(np.float32)
        if len(mz) == 0: continue
        maps += 1
        idx, apexes = rows_of[pid]; target = t['TARGET_MZ'][k]
        if t['MS_LEVEL'][k] == 1:
            ann = t['ANNOTATION'][k] or ''
            if '_i' not in ann: continue
            iso = int(ann.rsplit('_i', 1)[1])
            if iso > 3: continue
            keep = np.abs((mz - target) / target * 1e6 - OFF) <= MS1_PPM
            bb = bins(rt[keep], apexes)
            for w, j in enumerate(idx):
                m = bb[w] >= 0
                MS1[j, iso] += np.bincount(bb[w][m], it[keep][m], RT_BINS).astype(np.float32)
            continue
        mzs, li = frag[pid]
        lib_target = target / (1 + OFF * 1e-6)
        fi = int(np.argmin(np.abs(mzs - lib_target)))
        if abs(mzs[fi] - lib_target) / lib_target * 1e6 > 0.5: continue
        ppm = (mz - target) / target * 1e6
        keep = np.abs(ppm) < SUB * SUB_PPM / 2
        bb = bins(rt[keep], apexes); ik, pk = it[keep], ppm[keep]
        for w, j in enumerate(idx):
            LIB[j] = li; m = bb[w] >= 0
            XIC[j, fi] += np.bincount(bb[w][m], ik[m], RT_BINS).astype(np.float32)
            PPM[j, fi] += np.bincount(bb[w][m], ik[m] * pk[m], RT_BINS).astype(np.float32)
    if rg % 200 == 0: print('row group', rg, '/', f.num_row_groups, flush=True)
for a in (XIC, PPM, MS1, LIB): a.flush()
print(RAW, 'peak maps', maps, 'windows', n, 'nonzero xic', float((XIC.reshape(n, -1).max(1) > 0).mean()),
      'nonzero ms1', float((MS1.reshape(n, -1).max(1) > 0).mean()), 'lib set', float((LIB.max(1) > 0).mean()))
