#!/usr/bin/env python3
"""Step B preview (offline, diagnosis only): label-free shared-fragment competition after scoring.
Per run table: rows in the top region (baseline q <= QR) are candidates. A candidate p (target or decoy) is removed when a
TARGET d with a higher score, a different peptide sequence (siblings exempt), kept apex within DT s of p's, precursor m/z
within W Th, shares >= K of p's library fragments (FTOL ppm). The same rule for targets and decoys; no entrapment labels.
Then q is recomputed with paired TDC on the tool score (winner per own-decoy pair, q = (D + 1) / T, monotone).
Reports per arm, pooled: own q <= 0.01 N, E, lower, combined, median N@FDP_c 1 %, removed counts by class.
Baseline uses the same recomputed TDC, so baseline vs rule is like for like.
Usage: sim_compete.py <lib.pqp> <W Th> <DT s> <FTOL ppm> <K> <QR> <arm>=<glob> ..."""
import sqlite3, sys, glob, numpy as np, pandas as pd
lib, W, DT, FT, K, QR = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]) * 1e-6, int(sys.argv[5]), float(sys.argv[6])
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
P = pd.read_sql("""select p.ID id, p.TRAML_ID t, p.DECOY d, p.PRECURSOR_MZ mz, pe.MODIFIED_SEQUENCE seq,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
P['key'] = np.where(P.d == 1, P.t.str.replace('^DECOY_', '', regex=True), P.t)
Pi = P.set_index('id')
T = P[P.d == 0]; r = (T.e == 1).sum() / (T.e == 0).sum(); kf = 1 + 1 / r
fr = {}
def load(ids):
    need = [i for i in ids if i not in fr]
    for k in range(0, len(need), 20000):
        ch = need[k:k + 20000]
        for i in ch: fr[i] = []
        for pid, mz in c.execute(f"select tpm.PRECURSOR_ID, tr.PRODUCT_MZ from TRANSITION tr join TRANSITION_PRECURSOR_MAPPING tpm on tpm.TRANSITION_ID=tr.ID where tpm.PRECURSOR_ID in ({','.join(map(str, ch))})"):
            fr[pid].append(mz)
    for i in need: fr[i] = np.array(fr[i])
def tdc(t):
    w = t.sort_values('score', ascending=False).drop_duplicates('key')
    d = np.cumsum(w.is_decoy.values); n = np.cumsum(1 - w.is_decoy.values)
    q = np.minimum(1, (d + 1) / np.maximum(n, 1)); q = np.minimum.accumulate(q[::-1])[::-1]
    return w.assign(q=q)
def metrics(w):
    tt = w[w.is_decoy == 0]
    acc = tt[tt.q <= 0.01]; N, E = len(acc), int(acc.e.sum())
    e = tt.sort_values(['q', 'score'], ascending=[True, False]).e.values
    fdp = np.cumsum(e) * kf / np.arange(1, len(e) + 1); ok = np.where(fdp <= 0.01)[0]
    return N, E, int(ok.max() + 1) if len(ok) else 0
print(f'r {r:.4f}; W {W} Th, DT {DT} s, FTOL {FT * 1e6:.0f} ppm, K {K}, candidate region baseline q <= {QR}', flush=True)
for spec in sys.argv[7:]:
    arm, g = spec.split('=', 1)
    res = []; rem = dict(T=0, E=0, D=0, T_acc=0, E_acc=0)
    for f in sorted(glob.glob(g)):
        t = pd.read_csv(f, sep='\t', usecols=['precursor_id', 'is_decoy', 'rt_apex', 'score', 'q_value'])
        t = t[t.score > -1e29].sort_values('score', ascending=False).drop_duplicates('precursor_id')
        t = t.join(Pi[['key', 'e', 'mz', 'seq']], on='precursor_id')
        t['e'] = np.where(t.is_decoy == 1, 0, t.e)
        base = tdc(t)
        thr = base[base.q <= QR].score.min()
        top = t[t.score >= thr]
        don = top[top.is_decoy == 0].sort_values('rt_apex')
        load(set(top.precursor_id))
        drt, dmz, dsc, did, dseq = don.rt_apex.values, don.mz.values, don.score.values, don.precursor_id.values, don.seq.values
        drop = []
        for p, rt, mz, sc, sq in zip(top.precursor_id, top.rt_apex, top.mz, top.score, top.seq):
            lo, hi = np.searchsorted(drt, rt - DT), np.searchsorted(drt, rt + DT, side='right')
            j = np.arange(lo, hi)
            j = j[(dsc[j] > sc) & (np.abs(dmz[j] - mz) <= W) & (dseq[j] != sq) & (did[j] != p)]
            a = fr[p]
            if len(a) == 0: continue
            for d in did[j]:
                b = fr[d]
                if len(b) and int((np.abs(a[:, None] - b[None, :]) <= b[None, :] * FT).any(axis=1).sum()) >= K:
                    drop.append(p); break
        dr = t[t.precursor_id.isin(drop)]
        acc_ids = set(base[(base.is_decoy == 0) & (base.q <= 0.01)].precursor_id)
        rem['D'] += int((dr.is_decoy == 1).sum()); rem['T'] += int(((dr.is_decoy == 0) & (dr.e == 0)).sum()); rem['E'] += int((dr.e == 1).sum())
        rem['T_acc'] += int(dr[(dr.is_decoy == 0) & (dr.e == 0)].precursor_id.isin(acc_ids).sum())
        rem['E_acc'] += int(dr[dr.e == 1].precursor_id.isin(acc_ids).sum())
        new = tdc(t[~t.precursor_id.isin(drop)])
        tool = t[(t.is_decoy == 0) & (t.q_value <= 0.01)]
        res.append((len(tool), int(tool.e.sum())) + metrics(base) + metrics(new))
    a = np.array(res)
    def line(lab, N, E, n1): return f'{lab}: N {N.sum()} E {E.sum()} lower {E.sum() / N.sum():.4f} combined {E.sum() * kf / N.sum():.4f} median N@FDP_c1% {np.median(n1):.0f}'
    print(f'== {arm}: runs {len(a)}; tool q: N {a[:, 0].sum()} E {a[:, 1].sum()} lower {a[:, 1].sum() / a[:, 0].sum():.4f}', flush=True)
    print('  ' + line('baseline TDC', a[:, 2], a[:, 3], a[:, 4]), flush=True)
    print('  ' + line('competition ', a[:, 5], a[:, 6], a[:, 7]), flush=True)
    print(f'  removed (all runs): {rem}', flush=True)
