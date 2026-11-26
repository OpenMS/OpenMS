#!/usr/bin/env python3
"""PROMO audit (diagnosis): what the GBM adds. Same seed, network-only run (NET, entlib-gbt *_cur_*) vs default (G0, *_gbt_*).
Per run: winners of paired TDC on each arm's tool score. Targets accepted at q <= 0.01 split into
  kept     accepted in both arms
  promoted accepted by G0 only
  dropped  accepted by NET only
For each class: N, entrapments E, lower E/N and combined E(1+1/r)/N. The counted decoys (decoy winners above the arm's
q <= 0.01 score cut) split the same way (decoy counted in both / G0 only / NET only).
Under exchangeability, promoted E(1+1/r) should be about promoted D (the decoys the GBM newly lets in).
Usage: promo_audit.py <lib.pqp> <label> <net glob> <g0 glob>   (files matched by sorted order = seed)  env RELABEL optional"""
import sqlite3, sys, os, glob, numpy as np, pandas as pd
lib, label, gn, gg = sys.argv[1:5]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
P = pd.read_sql("""select p.ID id, p.TRAML_ID t, p.DECOY d, pe.UNMODIFIED_SEQUENCE useq,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
P['key'] = np.where(P.d == 1, P.t.str.replace('^DECOY_', '', regex=True), P.t)
if os.environ.get('RELABEL'):
    rl = set(pd.read_csv(os.environ['RELABEL'], sep='\t').seq); P.loc[(P.d == 0) & P.useq.isin(rl), 'e'] = 0
T = P[P.d == 0]; r = (T.e == 1).sum() / (T.e == 0).sum(); k = 1 + 1 / r
Pi = P.set_index('id')
def win(f):
    t = pd.read_csv(f, sep='\t', usecols=['precursor_id', 'is_decoy', 'score'])
    t = t[t.score > -1e29].sort_values('score', ascending=False).drop_duplicates('precursor_id')
    t = t.join(Pi[['key', 'e']], on='precursor_id'); t['e'] = np.where(t.is_decoy == 1, 0, t.e)
    w = t.drop_duplicates('key')
    d = np.cumsum(w.is_decoy.values); n = np.cumsum(1 - w.is_decoy.values)
    q = np.minimum(1, (d + 1) / np.maximum(n, 1)); w = w.assign(q=np.minimum.accumulate(q[::-1])[::-1])
    return w.set_index('key')
tot = {}
fn, fg = sorted(glob.glob(gn)), sorted(glob.glob(gg))
assert len(fn) == len(fg) and fn, (len(fn), len(fg))
for a, b in zip(fn, fg):
    wn, wg = win(a), win(b)
    acc_n = set(wn.index[(wn.is_decoy == 0) & (wn.q <= 0.01)]); acc_g = set(wg.index[(wg.is_decoy == 0) & (wg.q <= 0.01)])
    dec_n = set(wn.index[(wn.is_decoy == 1) & (wn.q <= 0.01)]); dec_g = set(wg.index[(wg.is_decoy == 1) & (wg.q <= 0.01)])
    e = pd.concat([wn.e, wg.e]).groupby(level=0).max()
    for cls, ks in (('kept', acc_n & acc_g), ('promoted', acc_g - acc_n), ('dropped', acc_n - acc_g)):
        x = tot.setdefault(cls, [0, 0]); x[0] += len(ks); x[1] += int(e.reindex(list(ks)).fillna(0).sum())
    for cls, ks in (('dec_both', dec_n & dec_g), ('dec_G0only', dec_g - dec_n), ('dec_NETonly', dec_n - dec_g)):
        tot[cls] = tot.get(cls, 0) + len(ks)
print(f'== {label}: runs {len(fn)}, r {r:.4f}')
for cls in ('kept', 'promoted', 'dropped'):
    N, E = tot[cls]
    print(f'  {cls:9s} N {N:7d} E {E:5d} lower {E / max(N, 1):.4f} combined {k * E / max(N, 1):.4f}  (false targets implied by E: {k * E:.0f})')
print(f'  counted decoys: both {tot["dec_both"]}, G0 only {tot["dec_G0only"]}, NET only {tot["dec_NETonly"]}')
print(f'  net change: false targets by E {k * (tot["promoted"][1] - tot["dropped"][1]):+.0f} vs counted decoys {tot["dec_G0only"] - tot["dec_NETonly"]:+d}')
