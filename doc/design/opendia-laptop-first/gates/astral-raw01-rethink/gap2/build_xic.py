"""Decode the OSW XICs of raw 01 (same library content, same raw) into dense arrays for a model-free gap analysis.
OSW precursor/transition IDs differ from our PQP, so XIC precursors are mapped to ours by (modified sequence, charge,
decoy) and fragments by annotation within precursor. Output (scratch, indexed by OUR precursor_id order in meta):
frag [N,6,T] (fragments sorted by our library intensity desc), lib [N,6], ms1 [N,3,T] (i0..i2 interpolated onto the
MS2 grid), rt [N,T] (NaN padded), npts [N]."""
import os, zlib, sqlite3, numpy as np, pandas as pd, pyarrow.parquet as pq, pyopenms as po
from multiprocessing import Pool
RUN = os.environ.get('OPENDIA_RUN', '01')
X = f'/home/singjc/Documents/datasets/proteobench-astral-entrapment/100K_sampled_fasta/opendia_results/LFQ_Astral_DIA_15min_50ng_Human_{RUN}_extracted.chrom.xic'
LIB = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-entrapment-20261107/lib/library.pqp'
OUT = os.environ['JCODE_SCRATCH_DIR'] + '/xic' + RUN
T = 128
def dec(b, lin):
    c = po.MSNumpressCoder(); cfg = po.NumpressConfig(); out = []
    cfg.np_compression = po.MSNumpressCoder.NumpressCompression.LINEAR if lin else po.MSNumpressCoder.NumpressCompression.SLOF
    c.decodeNPRaw(zlib.decompress(b), out, cfg); return np.asarray(out, np.float32)
def work(rg):
    t = pq.ParquetFile(X).read_row_group(rg, columns=['PRECURSOR_ID', 'MS_LEVEL', 'ANNOTATION', 'RT_DATA', 'INTENSITY_DATA']).to_pandas()
    return [(r.PRECURSOR_ID, r.MS_LEVEL, r.ANNOTATION, dec(r.RT_DATA, True), dec(r.INTENSITY_DATA, False))
            for r in t.itertuples(index=False) if not (r.MS_LEVEL == 1 and r.ANNOTATION == 'Precursor_i3')]
if __name__ == '__main__':
    pf = pq.ParquetFile(X)
    c = sqlite3.connect('file:' + LIB + '?mode=ro', uri=True)
    p = pd.read_sql('select pr.ID pid, pr.CHARGE ch, pr.DECOY decoy, pe.MODIFIED_SEQUENCE seq from PRECURSOR pr join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=pr.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID', c)
    x = pf.read(columns=['PRECURSOR_ID', 'MODIFIED_SEQUENCE', 'PRECURSOR_CHARGE', 'PRECURSOR_DECOY', 'MS_LEVEL']).to_pandas()
    x = x[x.MS_LEVEL == 1].drop_duplicates('PRECURSOR_ID')
    j = x.merge(p, left_on=['MODIFIED_SEQUENCE', 'PRECURSOR_CHARGE', 'PRECURSOR_DECOY'], right_on=['seq', 'ch', 'decoy'])
    o2u = dict(zip(j.PRECURSOR_ID.values, j.pid.values))
    tr = pd.read_sql('select m.PRECURSOR_ID pid, t.LIBRARY_INTENSITY li, t.ANNOTATION ann, t.ID tid from TRANSITION t join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID', c)
    tr = tr.sort_values(['pid', 'li', 'tid'], ascending=[True, False, True]); tr['slot'] = tr.groupby('pid').cumcount()
    smap = {(a, b): (s, l) for a, b, s, l in zip(tr.pid.values, tr.ann.values, tr.slot.values, tr.li.values)}
    meta = np.sort(j.pid.unique()); idx = {q: i for i, q in enumerate(meta)}; N = len(meta)
    frag = np.zeros((N, 6, T), np.float32); lib = np.zeros((N, 6), np.float32); rt = np.full((N, T), np.nan, np.float32)
    npts = np.zeros(N, np.int16); ms1raw = {}; miss = 0
    with Pool(14) as pool:
        for k, res in enumerate(pool.imap_unordered(work, range(pf.num_row_groups), chunksize=8)):
            for opid, lvl, ann, a, v in res:
                u = o2u.get(opid)
                if u is None: continue
                i = idx[u]
                if lvl == 2:
                    sl = smap.get((u, ann))
                    if sl is None or sl[0] >= 6: miss += 1; continue
                    n = min(len(a), T); frag[i, sl[0], :n] = v[:n]; lib[i, sl[0]] = sl[1]
                    if npts[i] == 0: rt[i, :n] = a[:n]; npts[i] = n
                else:
                    ms1raw[(i, int(ann[-1]))] = (a, v)
            if k % 1000 == 0: print('row groups', k, flush=True)
    print('unmapped MS2 transitions', miss)
    ms1 = np.zeros((N, 3, T), np.float32)
    for (i, iso), (a, v) in ms1raw.items():
        n = npts[i]
        if n: ms1[i, iso, :n] = np.interp(rt[i, :n], a, v, left=0, right=0)
    dd = p.set_index('pid').decoy.reindex(meta).values
    np.save(OUT + '/frag.npy', frag); np.save(OUT + '/lib.npy', lib); np.save(OUT + '/rt.npy', rt); np.save(OUT + '/ms1.npy', ms1)
    np.save(OUT + '/npts.npy', npts); pd.DataFrame({'precursor_id': meta, 'decoy': dd}).to_parquet(OUT + '/meta.parquet')
    print('N', N, 'npts pct', np.percentile(npts, [1, 50, 99]), 'frag slots filled', (lib > 0).sum(1).mean())
