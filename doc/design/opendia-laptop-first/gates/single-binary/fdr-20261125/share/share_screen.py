#!/usr/bin/env python3
"""SHARE screen (offline, used data only): can a label-free shared-fragment competition step fix the diaPASEF excess?
Written before running. Selection rule for a held-out precheck (fixed here): a RULE variant goes forward only if, on
every set (HeLa D, H, N2; Astral raw01-03 pooled), the pooled combined bound is <= 1.00 % AND the median N@FDP_c 1 % is
>= 1.00x the recomputed baseline. ORACLE rows use entrapment labels and are ceilings only, never candidates.

Per run table (tool score, kept row per precursor), with paired TDC recomputed on the tool score as the baseline:
  share(p) = max number of p's library fragments (FTOL ppm) shared with one DONOR, where donors are accepted targets
  (baseline q <= 0.01, entrapments included: the rule cannot see labels) with a higher score than p, a different
  modified sequence, apex within DT s, precursor m/z within W Th, and (IM variants) library 1/K0 within DIM.
  Computed for every row with baseline q <= QR (targets and decoys alike).
Variants:
  RULE_K / RULE_K_IM      remove every row (target or decoy) with share >= K, recompute paired TDC
  ORACLE_K / ORACLE_K_IM  remove only entrapment targets with share >= K (ceiling: perfect, lossless borrowing fix)
  ORACLE_ALLE             remove every accepted-region entrapment (ceiling of any fix; sanity row)
Usage: share_screen.py <lib.pqp> <W Th> <set label> <run tables...>   env RELABEL=<audit tsv> optional"""
import sqlite3, sys, os, numpy as np, pandas as pd
lib, W, label = sys.argv[1], float(sys.argv[2]), sys.argv[3]
DT, FT, DIM, QR = 3.0, 15e-6, 0.05, 0.05
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
P = pd.read_sql("""select p.ID id, p.TRAML_ID t, p.DECOY d, p.PRECURSOR_MZ mz, p.LIBRARY_DRIFT_TIME im, pe.MODIFIED_SEQUENCE seq,
  pe.UNMODIFIED_SEQUENCE useq,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
P['key'] = np.where(P.d == 1, P.t.str.replace('^DECOY_', '', regex=True), P.t)
if os.environ.get('RELABEL'):
    rl = set(pd.read_csv(os.environ['RELABEL'], sep='\t').seq)
    m = (P.d == 0) & P.useq.isin(rl); P.loc[m, 'e'] = 0
    print(f'relabelled as human: {int(m.sum())}', flush=True)
has_im = bool((P.im > 0).any())
Pi = P.set_index('id')
T = P[P.d == 0]; r = (T.e == 1).sum() / (T.e == 0).sum(); kf = 1 + 1 / r
print(f'{label}: r {r:.4f}; W {W} Th, DT {DT} s, FTOL {FT * 1e6:.0f} ppm, DIM {DIM} (library IM {"yes" if has_im else "no"}), QR {QR}', flush=True)
fr = {}
def load(ids):
    need = [i for i in ids if i not in fr]
    for k in range(0, len(need), 20000):
        ch = need[k:k + 20000]
        for i in ch: fr[i] = []
        for pid, mz in c.execute("select tpm.PRECURSOR_ID, tr.PRODUCT_MZ from TRANSITION tr join TRANSITION_PRECURSOR_MAPPING tpm "
                                 f"on tpm.TRANSITION_ID=tr.ID where tpm.PRECURSOR_ID in ({','.join(map(str, ch))})"):
            fr[pid].append(mz)
    for i in need: fr[i] = np.sort(np.array(fr[i]))
def nshared(a, b):
    if len(a) == 0 or len(b) == 0: return 0
    j = np.clip(np.searchsorted(b, a), 1, len(b) - 1)
    near = np.minimum(np.abs(b[j] - a), np.abs(b[j - 1] - a))
    return int((near <= a * FT).sum())
def tdc(t):
    w = t.sort_values('score', ascending=False).drop_duplicates('key')
    d = np.cumsum(w.is_decoy.values); n = np.cumsum(1 - w.is_decoy.values)
    q = np.minimum(1, (d + 1) / np.maximum(n, 1)); return w.assign(q=np.minimum.accumulate(q[::-1])[::-1])
def metrics(w):
    tt = w[w.is_decoy == 0]; acc = tt[tt.q <= 0.01]
    e = tt.sort_values(['q', 'score'], ascending=[True, False]).e.values
    fdp = np.cumsum(e) * kf / np.arange(1, len(e) + 1); ok = np.where(fdp <= 0.01)[0]
    return len(acc), int(acc.e.sum()), int(ok.max() + 1) if len(ok) else 0
rows = []
for f in sys.argv[4:]:
    t = pd.read_csv(f, sep='\t', usecols=['precursor_id', 'is_decoy', 'rt_apex', 'score'])
    t = t[t.score > -1e29].sort_values('score', ascending=False).drop_duplicates('precursor_id')
    t = t.join(Pi[['key', 'e', 'mz', 'im', 'seq']], on='precursor_id'); t['e'] = np.where(t.is_decoy == 1, 0, t.e)
    base = tdc(t); res = {'BASE': metrics(base)}
    thr_acc = base[(base.is_decoy == 0) & (base.q <= 0.01)].score.min(); thr_top = base[base.q <= QR].score.min()
    don = t[(t.is_decoy == 0) & (t.score >= thr_acc)].sort_values('rt_apex')
    top = t[t.score >= thr_top]
    load(set(top.precursor_id) | set(don.precursor_id))
    drt, dmz, dim, dsc, did, dseq = (don[x].values for x in ('rt_apex', 'mz', 'im', 'score', 'precursor_id', 'seq'))
    s_all, s_im = np.zeros(len(top), int), np.zeros(len(top), int)
    for i, (p, rt, mz, im, sc, sq) in enumerate(zip(top.precursor_id.values, top.rt_apex.values, top.mz.values, top.im.values,
                                                     top.score.values, top.seq.values)):
        lo, hi = np.searchsorted(drt, rt - DT), np.searchsorted(drt, rt + DT, side='right')
        j = np.arange(lo, hi); j = j[(dsc[j] > sc) & (np.abs(dmz[j] - mz) <= W) & (dseq[j] != sq) & (did[j] != p)]
        a = fr[p]; b1 = b2 = 0
        for jj in j:
            k = nshared(a, fr[did[jj]])
            b1 = max(b1, k)
            if has_im and abs(dim[jj] - im) <= DIM: b2 = max(b2, k)
        s_all[i], s_im[i] = b1, b2
    top = top.assign(s=s_all, s_im=s_im)
    for K in (2, 3, 4, 5, 6):
        for tag, col in (('', 's'), ('_IM', 's_im')):
            if tag and not has_im: continue
            drop = set(top.precursor_id[top[col] >= K])
            res[f'RULE_{K}{tag}'] = metrics(tdc(t[~t.precursor_id.isin(drop)]))
            dropE = set(top.precursor_id[(top[col] >= K) & (top.e == 1)])
            res[f'ORACLE_{K}{tag}'] = metrics(tdc(t[~t.precursor_id.isin(dropE)]))
    res['ORACLE_ALLE'] = metrics(tdc(t[~t.precursor_id.isin(set(top.precursor_id[top.e == 1]))]))
    acc = top[(top.is_decoy == 0) & (top.score >= thr_acc)]
    print(f'  {os.path.basename(os.path.dirname(f))}: accepted {len(acc)} (E {int(acc.e.sum())}); share>=4 accepted human '
          f'{int(((acc.s >= 4) & (acc.e == 0)).sum())} E {int(((acc.s >= 4) & (acc.e == 1)).sum())}; '
          f'IM-gated human {int(((acc.s_im >= 4) & (acc.e == 0)).sum())} E {int(((acc.s_im >= 4) & (acc.e == 1)).sum())}', flush=True)
    for k, (N, E, N1) in res.items(): rows.append(dict(run=f, v=k, N=N, E=E, N1=N1))
R = pd.DataFrame(rows); B = R[R.v == 'BASE'].set_index('run').N1
for v, g in R.groupby('v', sort=False):
    g = g.set_index('run'); rel = g.N1 / B
    print(f'{label}\t{v:14s} N {g.N.sum():7d} E {g.E.sum():5d} lower {g.E.sum() / g.N.sum():.4f} combined {kf * g.E.sum() / g.N.sum():.4f} '
          f'| median N@FDP_c1% {g.N1.median():.0f} ({g.N1.median() / B.median():.3f}x) min rel {rel.min():.3f}', flush=True)
