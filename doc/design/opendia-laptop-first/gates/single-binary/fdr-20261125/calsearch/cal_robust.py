#!/usr/bin/env python3
"""Calibration search, step 3 (dev only): robustness of group-wise TDC (GRP) on d_max, d_lib, i_nclaim and ANY to the
label-free cut quantile (0.75, 0.8, 0.9, 0.95 of decoy winners with q <= 0.2), and a 3-group version (0.75 and 0.95).
Per-run combined bound and N@FDP_c 1 % relative to R0. Usage: cal_robust.py <lib.pqp> <label> <features.tsv> ..."""
import sqlite3, sys, numpy as np, pandas as pd
lib, label = sys.argv[1], sys.argv[2]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
P = pd.read_sql("""select p.ID id, p.TRAML_ID t, p.DECOY d,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
P['key'] = np.where(P.d == 1, P.t.str.replace('^DECOY_', '', regex=True), P.t)
import os
if os.environ.get('RELABEL'):  # step A label audit: entrapment-labelled exact human substrings count as human
    seqs = dict(c.execute("select p.ID, pe.UNMODIFIED_SEQUENCE from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID").fetchall())
    rl = set(pd.read_csv(os.environ['RELABEL'], sep='\t').seq)
    mm = (P.d == 0) & P.id.map(seqs).isin(rl); P.loc[mm, 'e'] = 0
    print(f'relabelled as human: {int(mm.sum())}', flush=True)
P = P.set_index('id')
r = (P[P.d == 0].e == 1).sum() / (P[P.d == 0].e == 0).sum(); k = 1 + 1 / r
GX = ['d_lib', 'd_max', 'i_nclaim']
def qcalc(dec):
    d = np.cumsum(dec); n = np.cumsum(1 - dec)
    q = np.minimum(1, (d + 1) / np.maximum(n, 1)); return np.minimum.accumulate(q[::-1])[::-1]
def metrics(w, qcol):
    t = w[w.is_decoy == 0]; a = t[t[qcol] <= 0.01]
    e = t.sort_values([qcol, 'score'], ascending=[True, False]).e.values
    fdp = np.cumsum(e) * k / np.arange(1, len(e) + 1); ok = np.where(fdp <= 0.01)[0]
    return len(a), int(a.e.sum()), int(ok.max() + 1) if len(ok) else 0
def grp(w, g):
    qq = np.ones(len(w))
    for gv in np.unique(g):
        m = np.where(g == gv)[0]; o = m[np.argsort(-w.score.values[m], kind='stable')]
        qq[o] = qcalc(w.is_decoy.values[o])
    return metrics(w.assign(qq=qq), 'qq')
rows = []
for fp in sys.argv[3:]:
    f = pd.read_csv(fp, sep='\t', usecols=['precursor_id', 'is_decoy', 'score'] + GX)
    kf = f[f.score > -1e29].join(P[['key', 'e']], on='precursor_id'); kf['e'] = np.where(kf.is_decoy == 1, 0, kf.e)
    w = kf.sort_values('score', ascending=False).drop_duplicates('key').reset_index(drop=True)
    w['q'] = qcalc(w.is_decoy.values)
    res = {'R0': metrics(w, 'q')}
    ref = w[(w.is_decoy == 1) & (w.q <= 0.2)]
    for qu in (0.75, 0.8, 0.9, 0.95):
        cu = {x: np.quantile(ref[x], qu) for x in GX}
        for x in GX: res[f'{x}@{qu}'] = grp(w, (w[x] > cu[x]).values.astype(int))
        res[f'ANY@{qu}'] = grp(w, np.any([(w[x] > cu[x]).values for x in GX], axis=0).astype(int))
    for x in GX:
        lo, hi = np.quantile(ref[x], 0.75), np.quantile(ref[x], 0.95)
        res[f'{x}@3g'] = grp(w, (w[x] > lo).values.astype(int) + (w[x] > hi).values.astype(int))
    for name, (N, E, N1) in res.items(): rows.append(dict(run=fp, rule=name, N=N, E=E, N1=N1))
R = pd.DataFrame(rows); base = R[R.rule == 'R0'].set_index('run')
print(f'== {label} ({R.run.nunique()} runs, r = {r:.4f})')
for name, g in R.groupby('rule', sort=False):
    g = g.set_index('run'); rel = g.N1 / base.N1
    print(f'{label}\t{name:16s} N {g.N.sum():6d} E {g.E.sum():4d} combined {k * g.E.sum() / g.N.sum():.4f} '
          f'runs>1% {int((k * g.E / g.N > 0.01).sum())} | N1 sum {g.N1.sum():6d} {g.N1.sum() / base.N1.sum():.3f}x median {g.N1.median():.0f} '
          f'rel min {rel.min():.3f} runs>=R0 {int((rel >= 1).sum())}/{len(rel)}')
