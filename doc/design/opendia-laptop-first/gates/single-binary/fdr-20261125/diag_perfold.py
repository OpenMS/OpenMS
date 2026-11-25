#!/usr/bin/env python3
"""Pooled vs per-fold TDC on the tool's kept scores. Pooled TDC aligns folds by a decoy-fitted calibration (calibrateTail),
which reuses the same decoys that TDC then counts; per-fold TDC is invariant to any per-fold monotone transform, so it
cannot suffer from that. Per score and run: N at q<=0.01, E, combined bound, N@FDP_c1%; pooled over runs at the end.
Also per fold: decoys above the pooled threshold vs entrapments (2E), to see whether folds with few decoys carry excess E."""
import sqlite3, sys, numpy as np, pandas as pd
lib = sys.argv[1]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID precursor_id, p.TRAML_ID t, p.DECOY d,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
pq['key'] = np.where(pq.d == 1, pq.t.str.replace('^DECOY_', '', regex=True), pq.t)
P = pq.set_index('precursor_id')

def q_of(x):
    """x: pair winners with column st, dec; returns q aligned with x sorted desc (ties: decoy first)."""
    x = x.sort_values(['st', 'dec'], ascending=[False, False])
    d = np.cumsum(x.dec.values); t = np.cumsum(1 - x.dec.values)
    q = np.minimum(1, (d + 1) / np.maximum(t, 1))
    return x.assign(q=np.minimum.accumulate(q[::-1])[::-1])

def metrics(acc_targets, ranked_targets):
    n = len(acc_targets); e = int(acc_targets.e.sum())
    fdp = np.cumsum(ranked_targets.e.values) * 2 / np.arange(1, len(ranked_targets) + 1)
    ok = np.where(fdp <= 0.01)[0]
    return n, e, (2 * e / n if n else 0.0), int(ok.max() + 1) if len(ok) else 0

tot = {}
for path in sys.argv[2:]:
    f = pd.read_csv(path, sep='\t', usecols=['precursor_id', 'is_decoy', 'fold', 'score', 'q_value', 'network_score', 'network_q_value'])
    f['key'] = f.precursor_id.map(P.key); f['e'] = f.precursor_id.map(P.e).fillna(0).astype(int)
    out = [path.split('/')[0]]
    for s, q in (('network_score', 'network_q_value'), ('score', 'q_value')):
        k = f[f[s] > -1e29].rename(columns={s: 'st'}).rename(columns={'is_decoy': 'dec'})
        pw = k.sort_values(['st', 'dec'], ascending=[False, False]).drop_duplicates('key')
        # pooled (the tool's own q, recomputed on the calibrated score) and per fold
        pooled = q_of(pw)
        perf = pd.concat([q_of(g) for _, g in pw.groupby('fold')])
        for lab, x in (('pooled', pooled), ('perfold', perf)):
            t = x[x.dec == 0]
            # rank for N@FDP: by q then score (pooled score is comparable; per fold use q then st)
            t = t.sort_values(['q', 'st'], ascending=[True, False])
            m = metrics(t[t.q <= 0.01], t)
            r = tot.setdefault((s, lab), np.zeros(4)); r += m
            out.append(f'{s[:3]} {lab} N {m[0]} E {m[1]} c {m[2]:.4f} N@1% {m[3]}')
        # per-fold decoys above the pooled threshold vs 2E
        thr = pooled[(pooled.dec == 0) & (pooled.q <= 0.01)].st.min()
        a = pooled[pooled.st >= thr]
        out.append(f'{s[:3]} per-fold D/2E above pooled thr ' + ' '.join(f'f{fo}:{int(g.dec.sum())}/{2 * int(((g.dec == 0) & (g.e == 1)).sum())}' for fo, g in a.groupby('fold')))
    print('  |  '.join(out), flush=True)
print('pooled over runs (sum N, sum E, mean combined, sum N@1%):')
n_runs = len(sys.argv) - 2
for (s, lab), (n, e, c_, n1) in tot.items():
    print(f'  {s:14s} {lab:8s} N {int(n)} E {int(e)} combined {2 * e / n:.4f}  N@FDP_c1% {int(n1)} (mean {n1 / n_runs:.0f})')
