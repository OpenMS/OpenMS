#!/usr/bin/env python3
"""Calibration search, step 1 (dev only, Astral raw01/02 x 5 seeds, A runs or G0 runs). Diagnosis, labels for evaluation.
Paired precursor TDC on the tool score (as eval_replay.py). For each label-free feature X, bins by decoy-winner quantiles.
Per bin: accepted targets N, entrapments E, accepted decoys D (winners above the q<=0.01 cut). False-target estimate from
entrapments 2E (r = 1) against the TDC estimate D. Also the null-region ratio T/D among winners with q > 0.5.
A bin where 2E >> D at q<=0.01 while T/D in the null region is about 1 is a decoy/false-target non-exchangeability that a
label-free stratum reweighting cannot see; a bin where both are high could be corrected label-free.
Usage: cal_explore.py <lib.pqp> <features.tsv> ..."""
import sqlite3, sys, numpy as np, pandas as pd
lib = sys.argv[1]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
P = pd.read_sql("""select p.ID id, p.TRAML_ID t, p.DECOY d,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
P['key'] = np.where(P.d == 1, P.t.str.replace('^DECOY_', '', regex=True), P.t)
P = P.set_index('id')
FEATS = ['i_explained', 'i_nclaim', 'i_libcos_resid', 'i_dlibcos', 'e_explained', 'e_nclaim', 's_n', 's_best', 's_best_d',
         'f12_extra_share', 'd_max', 'd_lib', 'transition_support', 'x_ms1_coelution', 'mt_llr', 'abs_library_rt_error',
         'network_score']
W = []
for fp in sys.argv[2:]:
    hdr = pd.read_csv(fp, sep='\t', nrows=0).columns
    f = pd.read_csv(fp, sep='\t', usecols=['precursor_id', 'is_decoy', 'score'] + [x for x in FEATS if x in hdr])
    k = f[f.score > -1e29].join(P[['key', 'e']], on='precursor_id'); k['e'] = np.where(k.is_decoy == 1, 0, k.e)
    w = k.sort_values('score', ascending=False).drop_duplicates('key')
    d = np.cumsum(w.is_decoy.values); n = np.cumsum(1 - w.is_decoy.values)
    q = np.minimum(1, (d + 1) / np.maximum(n, 1)); w = w.assign(q=np.minimum.accumulate(q[::-1])[::-1], run=fp)
    W.append(w)
W = pd.concat(W)
acc = W[W.q <= 0.01]; nul = W[W.q > 0.5]
print(f'winners {len(W)}; at q<=0.01: T {int((acc.is_decoy == 0).sum())} E {int(acc.e.sum())} D {int(acc.is_decoy.sum())}; '
      f'2E/D {2 * acc.e.sum() / max(acc.is_decoy.sum(), 1):.2f}')
for x in FEATS:
    if x not in W: continue
    v = W[x].replace([np.inf, -np.inf], np.nan)
    if v.notna().sum() == 0: continue
    dec = v[(W.is_decoy == 1) & (W.q <= 0.2)].dropna()
    edges = np.unique(np.quantile(dec, [0, .25, .5, .75, .9, 1])) if len(dec) > 20 else np.unique(np.quantile(v.dropna(), [0, .5, 1]))
    b = np.clip(np.searchsorted(edges, v.values, side='right') - 1, 0, len(edges) - 2)
    out = []
    for i in range(len(edges) - 1):
        a = acc[b[W.q.values <= 0.01] == i]; z = nul[b[W.q.values > 0.5] == i]
        T, E, D = int((a.is_decoy == 0).sum()), int(a.e.sum()), int(a.is_decoy.sum())
        zT, zD = int((z.is_decoy == 0).sum()), int(z.is_decoy.sum())
        out.append(f'[{edges[i]:.3g},{edges[i + 1]:.3g}] T {T} E {E} D {D} 2E/D {2 * E / max(D, 1):.2f} null T/D {zT / max(zD, 1):.2f}')
    print(f'{x}: ' + ' | '.join(out))
