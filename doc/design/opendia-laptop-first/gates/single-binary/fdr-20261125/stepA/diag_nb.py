#!/usr/bin/env python3
"""Step A (diagnosis only, no tool change). Near-neighbour exposure by class and per-class own-q bounds.
A neighbour of precursor p is a HUMAN (non-entrapment) library target t != p that could give p borrowed signal:
  def iso : same charge, precursor m/z within 10 ppm (near-isobar; any library RT/IM)
  def win : precursor m/z within +-W Th (co-isolation), |library RT diff| <= RTT (normalised scale), |IM diff| <= 0.03
            (IM ignored when the library has none)
nb(p) = max over neighbours of #p library fragments within FTOL ppm of a neighbour fragment (-1: no neighbour).
Part 1: class exposure, 40,000 sampled per class (human, entrapment, decoy of human, decoy of entrapment), seed 0.
Part 2: own q <= 0.01 of each run table, entrapments split by nb class (0-1, 2, 3+ shared incl. none in 0-1);
class bound E_c (1 + 1/r_c) / N with r_c = r x (class share among all entrapments, from the sample). Pooled per arm.
Usage: diag_nb.py <lib.pqp> <W Th> <FTOL ppm> <RT tol, library units> <arm>=<glob> ..."""
import sqlite3, sys, glob, numpy as np, pandas as pd
lib, W, FT, RTT = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]) * 1e-6, float(sys.argv[4])
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
P = pd.read_sql("select ID id, DECOY d, CHARGE z, PRECURSOR_MZ mz, LIBRARY_RT rt, LIBRARY_DRIFT_TIME im from PRECURSOR", c)
ent = pd.read_sql("""select m.PRECURSOR_ID id, min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING m join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID
  join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by m.PRECURSOR_ID""", c)
P = P.merge(ent, on='id', how='left').fillna({'e': 0})
use_im = (P.im > 0).mean() > 0.5
tm = pd.read_sql("select tpm.PRECURSOR_ID id, tr.PRODUCT_MZ f from TRANSITION tr join TRANSITION_PRECURSOR_MAPPING tpm on tpm.TRANSITION_ID=tr.ID", c)
tm = tm.sort_values('id', kind='stable')
fid, fmz = tm.id.values, tm.f.values
del tm
ids_sorted = np.unique(fid)
start = np.searchsorted(fid, ids_sorted); end = np.searchsorted(fid, ids_sorted, side='right')
pos = pd.Series(np.arange(len(ids_sorted)), index=ids_sorted)
def frags(i):
    k = pos.get(i)
    return fmz[start[k]:end[k]] if k is not None else np.empty(0)
H = P[(P.d == 0) & (P.e == 0)].sort_values('mz').reset_index(drop=True)
hmz, hz, hrt, him, hid = H.mz.values, H.z.values, H.rt.values, H.im.values, H.id.values
hst, hen = start[pos.reindex(hid).values], end[pos.reindex(hid).values]
def nb(row, mode):
    if mode == 'iso':
        lo, hi = np.searchsorted(hmz, row.mz * (1 - 1e-5)), np.searchsorted(hmz, row.mz * (1 + 1e-5), side='right')
        j = np.arange(lo, hi); j = j[(hz[j] == row.z) & (hid[j] != row.id)]
    else:
        lo, hi = np.searchsorted(hmz, row.mz - W), np.searchsorted(hmz, row.mz + W, side='right')
        j = np.arange(lo, hi); m = (np.abs(hrt[j] - row.rt) <= RTT) & (hid[j] != row.id)
        if use_im: m &= np.abs(him[j] - row.im) <= 0.03
        j = j[m]
    if len(j) == 0: return -1
    a = frags(row.id)
    if len(a) == 0: return 0
    seg = [fmz[s:e] for s, e in zip(hst[j], hen[j])]
    lens = np.array([len(s) for s in seg]); b = np.concatenate(seg)
    hit = np.abs(a[:, None] - b[None, :]) <= b[None, :] * FT          # k x M
    starts = np.concatenate([[0], np.cumsum(lens)[:-1]])
    ok = lens > 0
    per = np.logical_or.reduceat(hit, starts[ok], axis=1).sum(axis=0) if ok.any() else np.zeros(1)
    return int(per.max())
rng = np.random.default_rng(0)
cls = {'human targets': (P.d == 0) & (P.e == 0), 'entrapment targets': (P.d == 0) & (P.e == 1)}
# decoy class by its own target's class: pair through TRAML key
T = pd.read_sql("select ID id, TRAML_ID t, DECOY d from PRECURSOR", c)
T['key'] = np.where(T.d == 1, T.t.str.replace('^DECOY_', '', regex=True), T.t)
tk = T[T.d == 0].merge(P[['id', 'e']], on='id').set_index('key').e
P['pe'] = P.id.map(T.set_index('id').key).map(tk)
cls['decoys of human'] = (P.d == 1) & (P.pe == 0); cls['decoys of entrapment'] = (P.d == 1) & (P.pe == 1)
r = ((P.d == 0) & (P.e == 1)).sum() / ((P.d == 0) & (P.e == 0)).sum()
print(f'library {lib}: precursors {len(P)}, human targets {len(H)}, r {r:.4f}, IM used {use_im}, W {W} Th, FTOL {FT * 1e6:.0f} ppm, RT tol {RTT}', flush=True)
share = {}
cache = {}
def nbv(df, mode):
    out = []
    for row in df.itertuples():
        k = (row.id, mode)
        if k not in cache: cache[k] = nb(row, mode)
        out.append(cache[k])
    return np.array(out)
for lab, m in cls.items():
    x = P[m]; x = x.iloc[rng.choice(len(x), size=min(40000, len(x)), replace=False)]
    parts = []
    for mode in ('iso', 'win'):
        v = nbv(x, mode)
        parts.append(f'{mode}: has {np.mean(v >= 0):.4f} >=1 {np.mean(v >= 1):.4f} >=2 {np.mean(v >= 2):.4f} >=3 {np.mean(v >= 3):.4f}')
        if lab == 'entrapment targets':
            share[mode] = np.array([np.mean(v <= 1), np.mean(v == 2), np.mean(v >= 3)])
    print(f'{lab:22s} n {len(x)} | ' + ' | '.join(parts), flush=True)
Pi = P.set_index('id')
for spec in sys.argv[5:]:
    arm, g = spec.split('=', 1)
    files = sorted(glob.glob(g))
    acc = {m: np.zeros(3, int) for m in ('iso', 'win')}; N = 0; Et = 0
    for f in files:
        t = pd.read_csv(f, sep='\t', usecols=['precursor_id', 'is_decoy', 'score', 'q_value'])
        t = t[t.score > -1e29].sort_values('score', ascending=False).drop_duplicates('precursor_id')
        a = t[(t.is_decoy == 0) & (t.q_value <= 0.01)]
        N += len(a)
        ae = Pi.loc[a.precursor_id[a.precursor_id.map(Pi.e) == 1]].reset_index()
        Et += len(ae)
        for mode in ('iso', 'win'):
            v = nbv(ae, mode)
            acc[mode] += [int((v <= 1).sum()), int((v == 2).sum()), int((v >= 3).sum())]
    k = 1 + 1 / r
    line = f'{arm}: runs {len(files)} N {N} E {Et} lower {Et / N:.4f} combined {Et * k / N:.4f}'
    for mode in ('iso', 'win'):
        b = [acc[mode][i] * (1 + 1 / (r * share[mode][i])) / N if share[mode][i] > 0 else float('nan') for i in range(3)]
        line += f' | {mode} E(0-1/2/3+) {acc[mode].tolist()} share {np.round(share[mode], 4).tolist()} bound {[round(x, 4) for x in b]}'
    print(line, flush=True)
