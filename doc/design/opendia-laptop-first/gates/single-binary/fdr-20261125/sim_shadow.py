#!/usr/bin/env python3
"""Offline simulation of a label-free shadow rule after scoring (diagnosis, dev data only).
A kept window of precursor A (target or decoy) is a shadow if another TARGET precursor B, with a higher score, has its kept
window apex within DT s of A's, precursor m/z within TOL (ppm) of A's at the same charge, and shares >= K of A's library
fragments (10 ppm). Shadows are removed from the competition (both targets and decoys, same rule), then pooled paired
TDC is recomputed on the tool's calibrated score. Reports own q<=0.01 N, E, combined bound, N@FDP_c1%, and removed counts."""
import sqlite3, sys, numpy as np, pandas as pd
from collections import defaultdict
lib = sys.argv[1]
TOL, DT, K = float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID precursor_id, p.TRAML_ID t, p.DECOY d, p.CHARGE chg, p.PRECURSOR_MZ mz,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
pq['key'] = np.where(pq.d == 1, pq.t.str.replace('^DECOY_', '', regex=True), pq.t)
P = pq.set_index('precursor_id')
frag = defaultdict(list)
for pid, mz in c.execute('select tpm.PRECURSOR_ID, tr.PRODUCT_MZ from TRANSITION tr join TRANSITION_PRECURSOR_MAPPING tpm on tpm.TRANSITION_ID=tr.ID'):
    frag[pid].append(mz)
frag = {k: np.array(v) for k, v in frag.items()}

def shared(a, b):
    fa, fb = frag[a], frag[b]
    return int((np.abs(fa[:, None] - fb[None, :]) <= fb[None, :] * 1e-5).any(axis=1).sum())

def tdc(x):
    """x: kept rows with st (calibrated score), dec, fold, key. Pooled paired TDC as in the tool (per-fold decoy median/MAD)."""
    x = x.copy()
    for fo, g in x.groupby('fold'):
        d = g[g.dec == 1].st.values
        med = np.median(d); mad = 1.4826 * np.median(np.abs(d - med))
        x.loc[g.index, 'z'] = (g.st - med) / mad
    pw = x.sort_values(['z', 'dec'], ascending=[False, False]).drop_duplicates('key')
    d = np.cumsum(pw.dec.values); t = np.cumsum(1 - pw.dec.values)
    q = np.minimum(1, (d + 1) / np.maximum(t, 1))
    return pw.assign(q=np.minimum.accumulate(q[::-1])[::-1])

def metrics(pw):
    t = pw[pw.dec == 0].sort_values(['q', 'z'], ascending=[True, False])
    acc = t[t.q <= 0.01]
    fdp = np.cumsum(t.e.values) * 2 / np.arange(1, len(t) + 1)
    ok = np.where(fdp <= 0.01)[0]
    return len(acc), int(acc.e.sum()), int(ok.max() + 1) if len(ok) else 0

tot = defaultdict(lambda: np.zeros(3))
rem = defaultdict(int)
for path in sys.argv[5:]:
    f = pd.read_csv(path, sep='\t', usecols=['precursor_id', 'is_decoy', 'fold', 'rt_apex', 'score'])
    k = f[f.score > -1e29].rename(columns={'score': 'st', 'is_decoy': 'dec'}).join(P[['key', 'e', 'mz', 'chg']], on='precursor_id')
    k = k.reset_index(drop=True)
    base = tdc(k)
    # Candidates for removal: rows above a loose depth (q <= 0.1 under the base) are the only ones that matter for q <= 0.01.
    cand = base[base.q <= 0.1].precursor_id
    cand = k[k.precursor_id.isin(cand) | (k.dec == 1)]
    cand = cand[cand.st >= base[base.q <= 0.1].z.min() * 0 + k.st.quantile(0.9)]  # top 10 % of kept rows by score
    tg = k[k.dec == 0].sort_values('rt_apex')
    trt, tmz, tz, tst, tid = tg.rt_apex.values, tg.mz.values, tg.chg.values, tg.st.values, tg.precursor_id.values
    drop = set()
    for pid, rt, mz, z, st in zip(cand.precursor_id, cand.rt_apex, cand.mz, cand.chg, cand.st):
        lo, hi = np.searchsorted(trt, rt - DT), np.searchsorted(trt, rt + DT)
        for j in range(lo, hi):
            if tid[j] == pid or tz[j] != z or tst[j] <= st: continue
            if abs(tmz[j] - mz) > mz * TOL * 1e-6: continue
            if shared(pid, tid[j]) >= K:
                drop.add(pid); break
    kk = k[~k.precursor_id.isin(drop)]
    new = tdc(kk)
    run = path.split('/')[0]
    a, b = metrics(base), metrics(new)
    dr = k[k.precursor_id.isin(drop)]
    acc_ids = set(base[(base.dec == 0) & (base.q <= 0.01)].precursor_id)
    rem['T acc'] += int(dr[(dr.dec == 0) & (dr.e == 0)].precursor_id.isin(acc_ids).sum())
    rem['E acc'] += int(dr[(dr.dec == 0) & (dr.e == 1)].precursor_id.isin(acc_ids).sum())
    rem['D'] += int((dr.dec == 1).sum()); rem['T'] += int(((dr.dec == 0) & (dr.e == 0)).sum()); rem['E'] += int(((dr.dec == 0) & (dr.e == 1)).sum())
    print(f'{run}: base N {a[0]} E {a[1]} c {2 * a[1] / a[0]:.4f} N@1% {a[2]} | shadow N {b[0]} E {b[1]} c {2 * b[1] / b[0]:.4f} N@1% {b[2]} | removed {len(drop)}', flush=True)
    tot['base'] += a; tot['shadow'] += b
print(f'TOL {TOL} ppm, DT {DT} s, K {K}; removed totals: {dict(rem)}')
for kx, (n, e, n1) in tot.items():
    print(f'  {kx}: N {int(n)} E {int(e)} combined {2 * e / n:.4f} N@FDP_c1% sum {int(n1)}')
