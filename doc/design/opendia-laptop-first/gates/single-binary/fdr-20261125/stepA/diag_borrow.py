#!/usr/bin/env python3
"""Step A: run-level borrowing test. Diagnosis only.
For each run table: accepted human IDs (target, non-entrapment, own q <= 0.01) are potential donors.
For a precursor p with kept apex rt_p: donors d != p with |rt_d - rt_p| <= DT s and |precursor m/z diff| <= W Th
(same isolation window proxy). share(p) = max over donors of #p library fragments within FTOL ppm of a donor fragment;
-1 if no donor. Groups compared (pooled over runs):
  E_acc : accepted entrapments (q <= 0.01)
  D_mt  : decoys scored at least as high as the lowest accepted target (the decoys TDC counted at q <= 0.01)
  D_band: the next decoys down, ranks 1..K beyond D_mt, K = 5 x |D_mt| (same-ish score region)
  T_acc : accepted human targets (sample of 3000 per run; the donor itself is excluded)
  E_band: entrapments in the same score band as D_band
If entrapments borrow signal, E_acc share >= 2 is more frequent than D_mt / D_band.
Usage: diag_borrow.py <lib.pqp> <W Th> <DT s> <FTOL ppm> <arm>=<glob> ..."""
import sqlite3, sys, glob, numpy as np, pandas as pd
lib, W, DT, FT = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]) * 1e-6
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
P = pd.read_sql("select ID id, PRECURSOR_MZ mz, CHARGE z from PRECURSOR", c).set_index('id')
ent = pd.read_sql("""select m.PRECURSOR_ID id, min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING m join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID
  join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by m.PRECURSOR_ID""", c).set_index('id').e
fr_cache = {}
def frags(ids):
    need = [i for i in ids if i not in fr_cache]
    for k in range(0, len(need), 20000):
        ch = need[k:k + 20000]
        q = ','.join(map(str, ch))
        for i in ch: fr_cache[i] = []
        for pid, mz in c.execute(f'select tpm.PRECURSOR_ID, tr.PRODUCT_MZ from TRANSITION tr join TRANSITION_PRECURSOR_MAPPING tpm on tpm.TRANSITION_ID=tr.ID where tpm.PRECURSOR_ID in ({q})'):
            fr_cache[pid].append(mz)
    for i in need: fr_cache[i] = np.array(fr_cache[i])
def share(pid, rt, mz, drt, dmz, did):
    lo, hi = np.searchsorted(drt, rt - DT), np.searchsorted(drt, rt + DT, side='right')
    j = np.arange(lo, hi); j = j[(np.abs(dmz[j] - mz) <= W) & (did[j] != pid)]
    if len(j) == 0: return -1, 0
    a = fr_cache[pid]
    best = 0
    for d in did[j]:
        b = fr_cache[d]
        if len(a) and len(b): best = max(best, int((np.abs(a[:, None] - b[None, :]) <= b[None, :] * FT).any(axis=1).sum()))
    return best, len(j)
rng = np.random.default_rng(0)
for spec in sys.argv[5:]:
    arm, g = spec.split('=', 1)
    files = sorted(glob.glob(g))
    out = {k: [] for k in ('E_acc', 'D_mt', 'D_band', 'E_band', 'T_acc')}
    for f in files:
        t = pd.read_csv(f, sep='\t', usecols=['precursor_id', 'is_decoy', 'rt_apex', 'score', 'q_value'])
        t = t[t.score > -1e29].sort_values('score', ascending=False).drop_duplicates('precursor_id')
        t['e'] = t.precursor_id.map(ent).fillna(0).values
        t['mz'] = P.mz.reindex(t.precursor_id).values
        acc = t[(t.is_decoy == 0) & (t.q_value <= 0.01)]
        thr = acc.score.min()
        don = acc[acc.e == 0].sort_values('rt_apex')
        dec = t[t.is_decoy == 1]
        dmt = dec[dec.score >= thr]
        dband = dec[dec.score < thr].head(5 * max(len(dmt), 20))
        lo_s = dband.score.min()
        eband = t[(t.is_decoy == 0) & (t.e == 1) & (t.score < thr) & (t.score >= lo_s)]
        tacc = don.iloc[rng.choice(len(don), size=min(3000, len(don)), replace=False)]
        groups = {'E_acc': acc[acc.e == 1], 'D_mt': dmt, 'D_band': dband, 'E_band': eband, 'T_acc': tacc}
        frags(set(don.precursor_id) | set().union(*[set(x.precursor_id) for x in groups.values()]))
        drt, dmz, did = don.rt_apex.values, don.mz.values, don.precursor_id.values
        for k, x in groups.items():
            out[k] += [share(p, r, m, drt, dmz, did) for p, r, m in zip(x.precursor_id, x.rt_apex, x.mz)]
    print(f'== {arm}: runs {len(files)} (W {W} Th, DT {DT} s, FTOL {FT * 1e6:.0f} ppm)', flush=True)
    for k, v in out.items():
        v = np.array(v) if len(v) else np.zeros((0, 2), int)
        s, n = v[:, 0], v[:, 1]
        print(f'  {k:7s} n {len(s):6d}: no donor {np.mean(s < 0):.3f} | share >=1 {np.mean(s >= 1):.3f} >=2 {np.mean(s >= 2):.3f} '
              f'>=3 {np.mean(s >= 3):.3f} >=4 {np.mean(s >= 4):.3f} | median donors {np.median(n):.0f}', flush=True)
