#!/usr/bin/env python3
"""Score gbm_replay outputs on Astral dev runs (entrapment labels for evaluation only; precursor-level paired TDC as screen_ab.py).
Usage: eval_replay.py <lib.pqp> <variant> <features.tsv>:<replay.tsv | self> ...  (self: the features file's own score)
Prints per run and pooled: N (targets q<=1%), E (accepted entrapments), lower E/N, combined 2E/N, N@FDP_c 1%."""
import sqlite3, sys, numpy as np, pandas as pd
lib, var = sys.argv[1], sys.argv[2]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
P = pd.read_sql("""select p.ID id, p.TRAML_ID t, p.DECOY d,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
P['key'] = np.where(P.d == 1, P.t.str.replace('^DECOY_', '', regex=True), P.t)
P = P.set_index('id')
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
for pair in sys.argv[3:]:
    fp, rp = pair.split(':')
    f = pd.read_csv(fp, sep='\t', usecols=['precursor_id', 'is_decoy'] + (['score'] if rp == 'self' else []))
    if rp != 'self':
        r = pd.read_csv(rp, sep='\t')
        assert len(r) == len(f)
        f['score'] = r.score.values
    k = f[f.score > -1e29].join(P[['key', 'e']], on='precursor_id'); k['e'] = np.where(k.is_decoy == 1, 0, k.e)
    m = met(tdc(k)); tot += m
    run = fp.split('/')[-2]
    print(f'{var}\t{run}\tN {m[0]}\tE {m[1]}\tlower {m[1] / m[0]:.4f}\tcombined {2 * m[1] / m[0]:.4f}\tN@FDPc1% {m[2]}')
print(f'{var}\tPOOLED\tN {tot[0]}\tE {tot[1]}\tlower {tot[1] / tot[0]:.4f}\tcombined {2 * tot[1] / tot[0]:.4f}\tN@FDPc1% {tot[2]}')
