#!/usr/bin/env python3
"""Calibration search, step 2 (dev only: Astral raw01/02 x 5 seeds). Label-free rules applied to the tool's exported
features.tsv; entrapment labels (instr) used only to evaluate. Candidate list fixed in this file before running.
Precursor-paired TDC as the tool/eval_replay.py (winner per target/decoy key by score, q = (D+1)/T).
Rules:
  R0            the tool score (baseline)
  GATE_t        symmetric demotion: rows (target or decoy) whose network-only q exceeds t are ranked below all others
  GRP_<X>       group-wise TDC: winners split by X above the 90th percentile of decoy winners with q <= 0.2 (label-free
                cut), q computed separately per group, accept q_g <= 0.01; ranking by (q_g, score)
  GRP_ANY       two groups: any of d_lib, d_max, i_nclaim above its decoy-90 % cut vs none
Usage: cal_sim.py <lib.pqp> <label> <features.tsv> ..."""
import sqlite3, sys, numpy as np, pandas as pd
lib, label = sys.argv[1], sys.argv[2]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
P = pd.read_sql("""select p.ID id, p.TRAML_ID t, p.DECOY d,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
P['key'] = np.where(P.d == 1, P.t.str.replace('^DECOY_', '', regex=True), P.t)
P = P.set_index('id')
r = (P[P.d == 0].e == 1).sum() / (P[P.d == 0].e == 0).sum(); k = 1 + 1 / r
GX = ['d_lib', 'd_max', 'i_nclaim', 'i_explained', 'f12_extra_share', 'e_nclaim']
def qcalc(dec):
    d = np.cumsum(dec); n = np.cumsum(1 - dec)
    q = np.minimum(1, (d + 1) / np.maximum(n, 1)); return np.minimum.accumulate(q[::-1])[::-1]
def metrics(w, qcol, scol):
    t = w[w.is_decoy == 0]; a = t[t[qcol] <= 0.01]
    e = t.sort_values([qcol, scol], ascending=[True, False]).e.values
    fdp = np.cumsum(e) * k / np.arange(1, len(e) + 1); ok = np.where(fdp <= 0.01)[0]
    return len(a), int(a.e.sum()), int(ok.max() + 1) if len(ok) else 0
rows = []
for fp in sys.argv[3:]:
    hdr = pd.read_csv(fp, sep='\t', nrows=0).columns
    f = pd.read_csv(fp, sep='\t', usecols=['precursor_id', 'is_decoy', 'score', 'network_q_value'] + [x for x in GX if x in hdr])
    kf = f[f.score > -1e29].join(P[['key', 'e']], on='precursor_id'); kf['e'] = np.where(kf.is_decoy == 1, 0, kf.e)
    w = kf.sort_values('score', ascending=False).drop_duplicates('key').reset_index(drop=True)
    w['q'] = qcalc(w.is_decoy.values)
    res = {'R0': metrics(w, 'q', 'score')}
    for t in (0.01, 0.02, 0.05, 0.1):
        s = w.score - 1e6 * (w.network_q_value > t)
        z = w.assign(s=s).sort_values('s', ascending=False); z['qq'] = qcalc(z.is_decoy.values)
        res[f'GATE_{t}'] = metrics(z, 'qq', 's')
    cuts = {}
    for x in GX:
        if x not in w: continue
        dv = w[(w.is_decoy == 1) & (w.q <= 0.2)][x].replace([np.inf, -np.inf], np.nan).dropna()
        cuts[x] = np.quantile(dv, 0.9) if len(dv) else np.inf
    def grp(g):
        z = w.assign(g=g).copy(); z['qq'] = 1.0
        for gv, idx in z.groupby('g').groups.items():
            sub = z.loc[idx].sort_values('score', ascending=False)
            z.loc[sub.index, 'qq'] = qcalc(sub.is_decoy.values)
        return metrics(z, 'qq', 'score')
    for x, cu in cuts.items():
        if np.isfinite(cu): res[f'GRP_{x}'] = grp((w[x] > cu).values)
    if all(x in cuts for x in ('d_lib', 'd_max', 'i_nclaim')):
        res['GRP_ANY'] = grp(((w.d_lib > cuts['d_lib']) | (w.d_max > cuts['d_max']) | (w.i_nclaim > cuts['i_nclaim'])).values)
    for name, (N, E, N1) in res.items():
        rows.append(dict(run=fp.split('/')[-2], rule=name, N=N, E=E, N1=N1))
    print(fp.split('/')[-2], 'cuts', {x: round(float(v), 4) for x, v in cuts.items()}, flush=True)
R = pd.DataFrame(rows)
base = R[R.rule == 'R0'].set_index('run').N1
print(f'== {label}: pooled over {R.run.nunique()} runs (r = {r:.4f}); combined = E(1+1/r)/N at own q <= 0.01')
for name, g in R.groupby('rule', sort=False):
    rel = (g.set_index('run').N1 / base)
    print(f'{label}\t{name:20s} N {g.N.sum():6d} E {g.E.sum():4d} combined {k * g.E.sum() / g.N.sum():.4f} worst {(k * g.E / g.N).max():.4f} '
          f'| N@FDP_c1% sum {g.N1.sum():6d} ({g.N1.sum() / base.sum():.3f}x) median {g.N1.median():.0f} runs>=R0 {int((rel >= 1).sum())}/{len(rel)}')
