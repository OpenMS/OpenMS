"""PRECHECK_IMFEAT step 2: per-candidate mobility features and 2-D patches from raw points (xipm export).
Sample: faint2 (6,000 DIA-NN targets + 30,000 decoys, r~1.56 D s0, default search = -im_tolerance_scale 0.6).
For each candidate peak (peaks.tsv of the dev-tool search with the run calibration and the scaled 1/K0 window), from the
MS2 points of its 6 search fragments (xipm, +-0.08 1/K0 around the calibrated precursor mobility):
  pooled view (+-0.031 around the calibrated mobility, as the default search):
    p_cos_apex  cosine of apex +-1 scan fragment sums with sqrt library intensities
    p_local_sn  log apex consensus - log median consensus outside the peak (same scan grid)
  own-apex view: consensus mobility apex m* = mode of the intensity-weighted 1/K0 histogram (0.004 bins, +-0.06) over
  the peak's RT core (apex +-2 scans):
    b_offset    |m* - calibrated mobility|
    b_frag_sd   intensity-weighted SD across fragments of each fragment's weighted mean 1/K0 within +-0.02 of m*
    b_cos_apex  cosine as p_cos_apex but only points within +-0.01 of m*
    b_support   fragments with signal within +-0.01 of m* at apex +-1 scan
    b_contrast  log (sum within +-0.01 of m*) - log (sum in 0.02..0.05 from m*), consensus over the core
    b_local_sn  as p_local_sn on the +-0.01 band trace
  2-D patch: fragments x 33 RT bins (native cycle, centred on apex) x 8 mobility bins over +-0.031, log1p / max.
Writes imfeat.parquet (precursor_id, rta, features) and patches2d.npy + patch_keys.parquet. Diagnostic only."""
import os, sqlite3, zlib, json, glob, numpy as np, pandas as pd, pyarrow as pa, pyarrow.compute, pyarrow.parquet as pq
from multiprocessing import Pool
O = os.environ.get('FAINT2', '/home/singjc/.jcode/scratch/dia2/gap7/faint2')
NPROC = int(os.environ.get('NPROC', '4'))
CAL = os.environ.get('CAL', glob.glob('/home/singjc/.jcode/scratch/dia2/gap7/ev/ims06/*.calibration.json')[0] if glob.glob('/home/singjc/.jcode/scratch/dia2/gap7/ev/ims06/*.calibration.json') else '')
cal = json.load(open(CAL))['im']
IM_A, IM_B, IM_W = cal['intercept_vssc'], cal['slope'], 0.6 * cal['tolerance_half_width_vssc']
P = pd.read_csv(f'{O}/peaks.tsv', sep='\t'); P['rta'] = P.rt_apex.round(3)
PG = {p: g for p, g in P.groupby('precursor_id')}
c = sqlite3.connect(f'{O}/mini.pqp')
lib = pd.read_sql('select ID tid, PRODUCT_MZ mz, LIBRARY_INTENSITY li from TRANSITION', c).set_index('tid')
prec = pd.read_sql('select ID pid, LIBRARY_DRIFT_TIME lim from PRECURSOR', c).set_index('pid')
T = pq.read_table(f'{O}/maps.xipm', columns=['MS_LEVEL', 'PRECURSOR_ID', 'TRANSITION_ID', 'TARGET_ION_MOBILITY', 'MZ_DATA', 'RT_DATA',
                                              'MOBILITY_DATA', 'INTENSITY_DATA', 'RT_COMPRESSION', 'MOBILITY_COMPRESSION', 'INTENSITY_COMPRESSION'])
T = T.filter(pa.compute.equal(T.column('MS_LEVEL'), 2))
groups = pd.Series(np.arange(T.num_rows)).groupby(T.column('PRECURSOR_ID').to_numpy()).apply(list).to_dict()
print(f'MS2 rows {T.num_rows}, precursors {len(groups)}, peaks {len(P)}', flush=True)
H, NB = 16, 8
CORE_S = 4.0  # OpenDIAWorkflow core_half_width_s

def dec(b, code):
    import pyopenms as oms
    if not b: return np.zeros(0)
    out = []; cfg = oms.NumpressConfig()
    cfg.np_compression = {5: oms.MSNumpressCoder.NumpressCompression.LINEAR, 6: oms.MSNumpressCoder.NumpressCompression.SLOF}[code]
    oms.MSNumpressCoder().decodeNPRaw(zlib.decompress(b), out, cfg)
    return np.asarray(out)

def cosw(v, w):
    n = np.linalg.norm(v); return float(v @ w / n) if n > 0 else 0.0

def work(pids):
    rows, patches = [], []
    for pid in pids:
        idx = groups.get(pid)
        if idx is None: continue
        sub = T.take(idx).to_pydict()
        frags = []
        for k in range(len(idx)):
            rt = dec(sub['RT_DATA'][k], sub['RT_COMPRESSION'][k]); im = dec(sub['MOBILITY_DATA'][k], sub['MOBILITY_COMPRESSION'][k])
            it = dec(sub['INTENSITY_DATA'][k], sub['INTENSITY_COMPRESSION'][k])
            frags.append((rt, im, it, lib.li.get(sub['TRANSITION_ID'][k], 0.0)))
        nf = len(frags)
        w = np.sqrt(np.maximum([f[3] for f in frags], 0)); w = w / (np.linalg.norm(w) + 1e-12)
        mcal = IM_A + IM_B * prec.lim.get(pid, np.nan)
        scans = np.unique(np.concatenate([f[0] for f in frags])) if nf else np.array([])
        if len(scans) < 3: continue
        si = [np.searchsorted(scans, f[0]) for f in frags]
        def trace(lo, hi, centre):
            X = np.zeros((len(scans), nf))
            for k, (rt, im, it, _) in enumerate(frags):
                m = np.abs(im - centre) <= hi
                if lo > 0: m &= np.abs(im - centre) >= lo
                np.add.at(X[:, k], si[k][m], it[m])
            return X
        Xp = trace(0, IM_W, mcal); cp = Xp.sum(1)
        for _, pk in PG.get(pid, P.iloc[:0]).iterrows():
            ia = int(np.argmin(np.abs(scans - pk.rt_apex)))
            a, b = np.searchsorted(scans, pk.rt_start - 1e-6), np.searchsorted(scans, pk.rt_end + 1e-6)
            outside = np.r_[cp[:a], cp[b:]]
            r = dict(precursor_id=pid, rta=pk.rta)
            r['p_cos_apex'] = cosw(Xp[max(ia - 1, 0):ia + 2].sum(0), w)
            r['p_local_sn'] = float(np.log1p(cp[ia]) - np.log1p(np.median(outside) if len(outside) else 0))
            # own mobility apex over the core (apex +-2 scans)
            lo_rt, hi_rt = scans[max(ia - 2, 0)], scans[min(ia + 2, len(scans) - 1)]
            ims, its = [], []
            for rt, im, it, _ in frags:
                m = (rt >= lo_rt) & (rt <= hi_rt) & (np.abs(im - mcal) <= 0.06); ims.append(im[m]); its.append(it[m])
            aim, ait = np.concatenate(ims), np.concatenate(its)
            if ait.sum() <= 0:
                r.update(b_offset=0.06, b_frag_sd=0.03, b_cos_apex=0.0, b_support=0, b_contrast=0.0, b_local_sn=0.0)
            else:
                edges = np.arange(mcal - 0.06, mcal + 0.0601, 0.004); hst, _ = np.histogram(aim, edges, weights=ait)
                hst = np.convolve(hst, [0.25, 0.5, 0.25], 'same'); ms = 0.5 * (edges[:-1] + edges[1:])[np.argmax(hst)]
                r['b_offset'] = float(abs(ms - mcal))
                fm, fw = [], []
                for im, it in zip(ims, its):
                    m = np.abs(im - ms) <= 0.02
                    if it[m].sum() > 0: fm.append(np.average(im[m], weights=it[m])); fw.append(it[m].sum())
                r['b_frag_sd'] = float(np.sqrt(np.average((np.array(fm) - np.average(fm, weights=fw)) ** 2, weights=fw))) if len(fm) > 1 else 0.03
                Xb = trace(0, 0.01, ms); cb = Xb.sum(1)
                r['b_cos_apex'] = cosw(Xb[max(ia - 1, 0):ia + 2].sum(0), w)
                r['b_support'] = int((Xb[max(ia - 1, 0):ia + 2].sum(0) > 0).sum())
                inn = ait[np.abs(aim - ms) <= 0.01].sum(); side = ait[(np.abs(aim - ms) >= 0.02) & (np.abs(aim - ms) <= 0.05)].sum()
                r['b_contrast'] = float(np.log1p(inn) - np.log1p(side))
                ob = np.r_[cb[:a], cb[b:]]
                r['b_local_sn'] = float(np.log1p(cb[ia]) - np.log1p(np.median(ob) if len(ob) else 0))
            rows.append(r)
            # C++-feasible variant: per-fragment 1/K0 histogram (30 bins of 0.004 over +-0.06 around the calibrated
            # mobility) of the core window |rt - apex| <= CORE_S, computed from bins only (no band trace).
            Hh = np.zeros((nf, 30))
            for k, (rt, im, it, _) in enumerate(frags):
                mb = np.floor((im - (mcal - 0.06)) / 0.004).astype(int)
                m = (np.abs(rt - pk.rt_apex) <= CORE_S) & (mb >= 0) & (mb < 30)
                np.add.at(Hh[k], mb[m], it[m])
            cons = Hh.sum(0)
            if cons.sum() <= 0:
                r.update(c_offset=0.06, c_frag_sd=0.03, c_cos=0.0, c_support=0, c_contrast=0.0)
            else:
                sm = np.convolve(cons, [0.25, 0.5, 0.25], 'same'); j = int(np.argmax(sm)); cen = mcal - 0.06 + 0.004 * (np.arange(30) + 0.5)
                r['c_offset'] = float(abs(cen[j] - mcal))
                near = np.abs(np.arange(30) - j) <= 5; fm, fw = [], []
                for k in range(nf):
                    s = Hh[k, near].sum()
                    if s > 0: fm.append(float((Hh[k, near] * cen[near]).sum() / s)); fw.append(s)
                r['c_frag_sd'] = float(np.sqrt(np.average((np.array(fm) - np.average(fm, weights=fw)) ** 2, weights=fw))) if len(fm) > 1 else 0.03
                inn = np.abs(np.arange(30) - j) <= 2; side = (np.abs(np.arange(30) - j) >= 5) & (np.abs(np.arange(30) - j) <= 12)
                v = Hh[:, inn].sum(1); r['c_cos'] = cosw(v, w); r['c_support'] = int((v > 0).sum())
                r['c_contrast'] = float(np.log1p(cons[inn].sum()) - np.log1p(cons[side].sum()))
            # 2-D patch: fragments x 33 scans x 8 mobility bins over +-IM_W around the calibrated mobility
            Z = np.zeros((6, 2 * H + 1, NB), np.float32)
            for k, (rt, im, it, _) in enumerate(frags[:6]):
                s = si[k] - ia + H; mb = np.floor((im - (mcal - IM_W)) / (2 * IM_W) * NB).astype(int)
                m = (s >= 0) & (s <= 2 * H) & (mb >= 0) & (mb < NB)
                np.add.at(Z[k], (s[m], mb[m]), it[m])
            mx = Z.max() + 1e-6; patches.append(np.log1p(Z) / np.log1p(mx))
    return rows, patches

if __name__ == '__main__':
    pids = list(groups.keys()); chunks = [pids[i::40] for i in range(40)]
    with Pool(NPROC) as pool:
        res = pool.map(work, chunks)
    rows = [r for rr, _ in res for r in rr]; patches = [p for _, pp in res for p in pp]
    D = pd.DataFrame(rows); D.to_parquet(f'{O}/imfeat.parquet')
    np.save(f'{O}/patches2d.npy', np.stack(patches)); D[['precursor_id', 'rta']].to_parquet(f'{O}/patch_keys.parquet')
    print(f'windows {len(D)}, patches {len(patches)}', flush=True)
    print(D.describe().T.round(3).to_string())
