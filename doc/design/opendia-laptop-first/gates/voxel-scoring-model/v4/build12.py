"""PLAN_V4 phase 2: 12-fragment per-window arrays from the 12-fragment .xipm export (astral-lib12-20261111/export<r>/v12.xipm).
Same geometry and m/z cut as build_voxels.py (32 RT bins over the window apex +- 15 s, |ppm to calibrated target| < 17.5).
Fragments are in predicted-intensity order (DIALibGen top 12). Usage: build12.py <raw01|raw02>.
Outputs in v2/<raw>/: xic12.npy [n,12,32] float32, lib12.npy [n,12] (predicted rel. intensity, 0 = no fragment),
orig12.npy [n,12] uint8 (1 = same annotation as one of the original 6-fragment library transitions).
Check printed at the end: for windows whose 6 original fragments are all among the 12, the summed original-fragment XIC
must equal the summed xic.npy (same raw points), up to float error."""
import sys, sqlite3, zlib, numpy as np, pandas as pd, pyarrow.parquet as pq, pynumpress
RAW = sys.argv[1]; r2 = RAW[-2:]
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
O = B + '/astral-voxel-pilot-20261110/v2'; L = B + '/astral-lib12-20261111'
RT_BINS, HALF, MZ_HALF, K = 32, 15.0, 17.5, 12
pop = pd.read_csv(f'{O}/{RAW}/population.tsv', sep='\t', usecols=['precursor_id', 'rt_apex'])
n = len(pop)
rows_of = {pid: (g.index.values, g.rt_apex.values) for pid, g in pop.groupby('precursor_id')}
c12 = sqlite3.connect(f'file:{L}/library12_pop.pqp?mode=ro', uri=True)
t12 = pd.read_sql("""select t.ID tid, m.PRECURSOR_ID pid, t.ANNOTATION ann, t.LIBRARY_INTENSITY li from TRANSITION t
  join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID order by m.PRECURSOR_ID, t.ID""", c12)
c6 = sqlite3.connect(f'file:{B}/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
t6 = pd.read_sql("""select m.PRECURSOR_ID pid, t.ANNOTATION ann from TRANSITION t join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID""", c6)
orig = set(zip(t6.pid, t6.ann)); n6 = t6.groupby('pid').size()
t12['rank'] = t12.groupby('pid').cumcount()
t12['orig'] = [(p, a) in orig for p, a in zip(t12.pid, t12.ann)]
slot = dict(zip(t12.tid, t12['rank']))
mm = lambda name, shape, dt: np.lib.format.open_memmap(f'{O}/{RAW}/{name}.npy', mode='w+', dtype=dt, shape=shape)
X, LI, OR = mm('xic12', (n, K, RT_BINS), np.float32), mm('lib12', (n, K), np.float32), mm('orig12', (n, K), np.uint8)
for pid, g in t12.groupby('pid'):
    if pid not in rows_of: continue
    idx = rows_of[pid][0]; k = len(g)
    li = g.li.values / max(g.li.values.max(), 1e-12)
    LI[idx, :k] = li.astype(np.float32); OR[idx, :k] = g.orig.values.astype(np.uint8)

def dec(blob, kind):
    if blob is None or len(blob) == 0: return np.zeros(0)
    raw = zlib.decompress(blob)
    if kind == 5 and len(raw) == 12:
        fp = np.frombuffer(raw[:8], '>f8')[0]
        return np.array([np.frombuffer(raw[8:12], '<i4')[0] / fp])
    return np.asarray(pynumpress.decode_linear(raw) if kind == 5 else pynumpress.decode_slof(raw), float)

def bins(rt, apexes):
    d = rt[None, :] - apexes[:, None]
    b = ((d + HALF) / (2 * HALF) * RT_BINS).astype(int)
    return np.where(np.abs(d) < HALF, np.clip(b, 0, RT_BINS - 1), -1)

f = pq.ParquetFile(f'{L}/export{r2}/v12.xipm')
cols = ['MS_LEVEL', 'PRECURSOR_ID', 'TRANSITION_ID', 'TARGET_MZ', 'MZ_DATA', 'RT_DATA', 'INTENSITY_DATA', 'MZ_COMPRESSION',
        'RT_COMPRESSION', 'INTENSITY_COMPRESSION']
maps = 0
for rg in range(f.num_row_groups):
    t = f.read_row_group(rg, columns=cols).to_pydict()
    for k in range(len(t['MS_LEVEL'])):
        if t['MS_LEVEL'][k] != 2: continue
        pid = int(t['PRECURSOR_ID'][k])
        if pid not in rows_of: continue
        s = slot.get(int(t['TRANSITION_ID'][k]))
        if s is None: continue
        mz = dec(t['MZ_DATA'][k], t['MZ_COMPRESSION'][k])
        if len(mz) == 0: continue
        rt = dec(t['RT_DATA'][k], t['RT_COMPRESSION'][k]); it = dec(t['INTENSITY_DATA'][k], t['INTENSITY_COMPRESSION'][k])
        target = t['TARGET_MZ'][k]; keep = np.abs((mz - target) / target * 1e6) < MZ_HALF
        idx, apexes = rows_of[pid]; bb = bins(rt[keep], apexes); ik = it[keep]; maps += 1
        for w, j in enumerate(idx):
            m = bb[w] >= 0
            X[j, s] += np.bincount(bb[w][m], ik[m], RT_BINS).astype(np.float32)
    if rg % 200 == 0: print('row group', rg, '/', f.num_row_groups, flush=True)
for a in (X, LI, OR): a.flush()
# Check against the 6-fragment arrays.
xic = np.load(f'{O}/{RAW}/xic.npy', mmap_mode='r')
full = (np.asarray(OR).sum(1) == n6.reindex(pop.precursor_id.values).fillna(-1).values)
sel = np.flatnonzero(full)[:: max(1, int(full.sum()) // 50000)]
a = (np.asarray(X[sel]) * np.asarray(OR[sel])[:, :, None]).sum((1, 2)); b = np.asarray(xic[sel]).sum((1, 2))
rel = np.abs(a - b) / np.maximum(b, 1)
print(RAW, 'ms2 maps', maps, 'windows', n, 'windows with all 6 originals among the 12', round(float(full.mean()), 4),
      'original-fragment XIC check on', len(sel), ': median rel diff', float(np.median(rel)), 'p99', float(np.quantile(rel, 0.99)),
      'fragments per window', np.bincount((np.asarray(LI) > 0).sum(1), minlength=13).tolist(), flush=True)
