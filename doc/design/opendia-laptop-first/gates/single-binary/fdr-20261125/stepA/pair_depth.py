#!/usr/bin/env python3
"""Step A: is the target-over-own-decoy asymmetry specific to entrapments, or generic to natural sequences vs their
pseudo-reverse decoys? Diagnosis only. For human pairs and entrapment pairs separately: P(target scores above its own decoy)
by best-member q band, including deep bands where targets are mostly absent (null-like). Pooled over runs.
Usage: pair_depth.py <lib.pqp> <arm>=<glob> ..."""
import sqlite3, sys, glob, numpy as np, pandas as pd
lib = sys.argv[1]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
P = pd.read_sql("""select p.ID id, p.TRAML_ID t, p.DECOY d,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
P['key'] = np.where(P.d == 1, P.t.str.replace('^DECOY_', '', regex=True), P.t)
pairs = P[P.d == 0].set_index('key')[['id', 'e']].join(P[P.d == 1].set_index('key')[['id']].rename(columns={'id': 'did'}), how='inner')
bands = [(0, 0.001), (0.001, 0.01), (0.01, 0.05), (0.05, 0.2), (0.2, 0.5), (0.5, 0.8), (0.8, 1.01)]
for spec in sys.argv[2:]:
    arm, g = spec.split('=', 1)
    acc = {(e, b): [0, 0] for e in (0, 1) for b in bands}
    both = [0, 0, 0]
    for f in sorted(glob.glob(g)):
        t = pd.read_csv(f, sep='\t', usecols=['precursor_id', 'score', 'q_value'])
        t = t[t.score > -1e29].sort_values('score', ascending=False).drop_duplicates('precursor_id').set_index('precursor_id')
        st, sd = t.score.reindex(pairs.id).values, t.score.reindex(pairs.did).values
        qb = np.fmin(t.q_value.reindex(pairs.id).values, t.q_value.reindex(pairs.did).values)
        ok = ~np.isnan(st) & ~np.isnan(sd)
        both[0] += int(ok.sum()); both[1] += int((~np.isnan(st) & np.isnan(sd)).sum()); both[2] += int((np.isnan(st) & ~np.isnan(sd)).sum())
        for e in (0, 1):
            for b in bands:
                m = ok & (pairs.e.values == e) & (qb > b[0]) & (qb <= b[1])
                acc[(e, b)][0] += int((m & (st > sd)).sum()); acc[(e, b)][1] += int(m.sum())
    print(f'== {arm}: pairs with both scored {both[0]}, only target {both[1]}, only decoy {both[2]}', flush=True)
    for e, lab in ((0, 'human'), (1, 'entrapment')):
        print(f'  {lab:10s} ' + ' | '.join(f'q({b[0]},{min(b[1], 1):g}] P {w / max(n, 1):.3f} (n {n})' for b in bands for w, n in [acc[(e, b)]]), flush=True)
