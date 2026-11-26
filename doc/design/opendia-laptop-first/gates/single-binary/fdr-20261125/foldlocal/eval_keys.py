#!/usr/bin/env python3
"""eval_replay.py with the library replaced by a key table (id, traml, decoy, entrapment) that was exported from it read-only.
Usage: eval_keys.py <keys_e.tsv> <variant> <kept.tsv> ...  (gbm_replay output: row, precursor_id, is_decoy, score; kept rows only)
Metrics are the same as devAB/eval_replay.py: precursor-level paired TDC at q<=1 %. N = accepted targets, E = accepted entrapments
(labels for evaluation only), lower bound E/N, combined bound 2E/N, and N@FDP_c 1 %."""
import sys, numpy as np, pandas as pd
P = pd.read_csv(sys.argv[1], sep='\t', header=None, names=['id', 't', 'd', 'e']).set_index('id')
P['key'] = np.where(P.d == 1, P.t.str.replace('^DECOY_', '', regex=True), P.t)
var = sys.argv[2]
def tdc(w):
    w = w.sort_values('score', ascending=False).drop_duplicates('key')
    d = np.cumsum(w.is_decoy.values); n = np.cumsum(1 - w.is_decoy.values)
    q = np.minimum(1, (d + 1) / np.maximum(n, 1)); return w.assign(q=np.minimum.accumulate(q[::-1])[::-1])
def met(w):
    t = w[w.is_decoy == 0]; a = t[t.q <= 0.01]
    e = t.sort_values(['q', 'score'], ascending=[True, False]).e.values
    fdp = np.cumsum(e) * 2 / np.arange(1, len(e) + 1); ok = np.where(fdp <= 0.01)[0]
    return np.array([len(a), int(a.e.sum()), int(ok.max() + 1) if len(ok) else 0])
tot = np.zeros(3, int)
for fp in sys.argv[3:]:
    f = pd.read_csv(fp, sep='\t')
    k = f[f.score > -1e29].join(P[['key', 'e']], on='precursor_id'); k['e'] = np.where(k.is_decoy == 1, 0, k.e)
    m = met(tdc(k)); tot += m
    run = fp.split('/')[-1]
    print(f'{var}\t{run}\tN {m[0]}\tE {m[1]}\tlower {m[1] / m[0]:.4f}\tcombined {2 * m[1] / m[0]:.4f}\tN@FDPc1% {m[2]}')
print(f'{var}\tPOOLED\tN {tot[0]}\tE {tot[1]}\tlower {tot[1] / tot[0]:.4f}\tcombined {2 * tot[1] / tot[0]:.4f}\tN@FDPc1% {tot[2]}')
