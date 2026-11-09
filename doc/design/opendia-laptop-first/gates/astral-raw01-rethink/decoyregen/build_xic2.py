"""Decode the regenerated-decoy XICs (decoyregen/raw{RUN}.chrom.xic, extracted with decoyregen/merged.pqp) into dense
arrays, like ../gap2/build_xic.py, but mapping precursors by ID. OSW was run on merged.pqp itself, so the IDs are
shared. Fragments are mapped by annotation within the precursor, sorted by library intensity (top 6).
Output: $JCODE_SCRATCH_DIR/xicR{RUN}/ frag, lib, rt, ms1, npts, meta.parquet (precursor_id, decoy, family, pair)."""
import os, sys, zlib, sqlite3, numpy as np, pandas as pd, pyarrow.parquet as pq, pyopenms as po
from multiprocessing import Pool
RUN = os.environ.get('OPENDIA_RUN', '02')
HERE = os.path.dirname(os.path.abspath(__file__))
X = f'{HERE}/raw{RUN}.chrom.xic'; LIB = f'{HERE}/merged.pqp'
OUT = os.environ['JCODE_SCRATCH_DIR'] + '/xicR' + RUN; os.makedirs(OUT, exist_ok=True)
T = 144

def dec(b, lin):
    c = po.MSNumpressCoder(); cfg = po.NumpressConfig(); out = []
    cfg.np_compression = po.MSNumpressCoder.NumpressCompression.LINEAR if lin else po.MSNumpressCoder.NumpressCompression.SLOF
    c.decodeNPRaw(zlib.decompress(b), out, cfg); return np.asarray(out, np.float32)

def work(rg):
    t = pq.ParquetFile(X).read_row_group(rg, columns=['PRECURSOR_ID', 'MS_LEVEL', 'ANNOTATION', 'RT_DATA', 'INTENSITY_DATA']).to_pandas()
    return [(int(r.PRECURSOR_ID), r.MS_LEVEL, r.ANNOTATION, dec(r.RT_DATA, True), dec(r.INTENSITY_DATA, False))
            for r in t.itertuples(index=False) if not (r.MS_LEVEL == 1 and r.ANNOTATION == 'Precursor_i3')]

if __name__ == '__main__':
    pf = pq.ParquetFile(X)
    c = sqlite3.connect('file:' + LIB + '?mode=ro', uri=True)
    p = pd.read_sql('select ID pid, TRAML_ID tid, DECOY decoy from PRECURSOR', c)
    fam = np.where(p.decoy == 0, 'T', np.where(p.tid.str.startswith('DECOYP_'), 'P', np.where(p.tid.str.startswith('DECOYS_'), 'S', 'O')))
    p['family'] = fam
    # pair = merged-PQP ID of the target. A decoy TRAML_ID is <tag><targets.pqp ID> (checked: DECOYS_374935 is the
    # shuffle of targets.pqp ID 374935). Map targets.pqp ID -> target TRAML_ID -> merged ID, not assuming equal IDs.
    tp = pd.read_sql('select ID tpid, TRAML_ID tid from PRECURSOR', sqlite3.connect('file:' + HERE + '/targets.pqp?mode=ro', uri=True))
    traml2merged = dict(zip(p.tid[p.decoy == 0], p.pid[p.decoy == 0]))
    tp2merged = {a: traml2merged.get(b, -1) for a, b in zip(tp.tpid, tp.tid)}
    suf = pd.to_numeric(p.tid.str.split('_').str[-1], errors='coerce')
    p['pair'] = np.where(p.decoy == 0, p.pid, [tp2merged.get(int(s), -1) if np.isfinite(s) else -1 for s in suf]).astype('int64')
    print('decoys without a mapped target:', int(((p.decoy == 1) & (p.pair < 0)).sum()))
    tr = pd.read_sql('select m.PRECURSOR_ID pid, t.LIBRARY_INTENSITY li, t.ANNOTATION ann, t.ID tid from TRANSITION t join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID', c)
    tr = tr.sort_values(['pid', 'li', 'tid'], ascending=[True, False, True]); tr['slot'] = tr.groupby('pid').cumcount()
    smap = {(a, b): (s, l) for a, b, s, l in zip(tr.pid.values, tr.ann.values, tr.slot.values, tr.li.values)}
    ids = pf.read(columns=['PRECURSOR_ID', 'MS_LEVEL']).to_pandas()
    meta = np.sort(ids[ids.MS_LEVEL == 1].PRECURSOR_ID.unique().astype(np.int64))
    idx = {q: i for i, q in enumerate(meta)}; N = len(meta)
    frag = np.zeros((N, 6, T), np.float32); lib = np.zeros((N, 6), np.float32); rt = np.full((N, T), np.nan, np.float32)
    npts = np.zeros(N, np.int16); ms1raw = {}; miss = 0
    with Pool(12) as pool:
        for k, res in enumerate(pool.imap_unordered(work, range(pf.num_row_groups), chunksize=8)):
            for opid, lvl, ann, a, v in res:
                i = idx.get(opid)
                if i is None: continue
                if lvl == 2:
                    sl = smap.get((opid, ann))
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
    pm = p.set_index('pid').reindex(meta)
    np.save(OUT + '/frag.npy', frag); np.save(OUT + '/lib.npy', lib); np.save(OUT + '/rt.npy', rt); np.save(OUT + '/ms1.npy', ms1)
    np.save(OUT + '/npts.npy', npts)
    pd.DataFrame({'precursor_id': meta, 'decoy': pm.decoy.values, 'family': pm.family.values, 'pair': pm.pair.values}).to_parquet(OUT + '/meta.parquet')
    print('N', N, 'by family', pd.Series(pm.family.values).value_counts().to_dict(), 'npts pct', np.percentile(npts, [1, 50, 99]))
