#!/usr/bin/env python3
"""FDR mechanism diagnosis from OpenDIA -features_out tables (entrapment labels read only here).

Usage: diag_fdr.py <library.pqp> <features.tsv> [...]

Per run, for the network score (pre-GBM) and the final (GBM) score, best window per precursor:
- N, entrapments E, combined bound E(1+1/r)/N at the tool's own q <= 0.01;
- the same for targets admitted by GBM but not by the network, and dropped by GBM;
- "decoy calibration": in the top T targets ranked by score, compare the decoy count from TDC to the entrapment
  estimate of false targets (E(1+1/r)). Ratio > 1 means decoys under-count false targets (optimistic q).
- Per-feature standardised mean difference between entrapment winners and decoy winners at q <= 0.05 (which features
  separate false targets from decoys: those are where a model can rank false targets above decoys)."""
import sqlite3, sys, numpy as np, pandas as pd

lib = sys.argv[1]
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID precursor_id, p.DECOY d,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID
  group by p.ID""", c)
tg = pq[pq.d == 0]
r = (tg.e == 1).sum() / (tg.e == 0).sum()
k = 1 + 1 / r
ent = tg.set_index('precursor_id').e
print(f'library targets {len(tg)}, entrapment ratio r = {r:.4f}, k = {k:.4f}')

for path in sys.argv[2:]:
    f = pd.read_csv(path, sep='\t')
    feats = [x for x in f.columns if x not in ('candidate_row', 'precursor_id', 'is_decoy', 'fold', 'rt_apex', 'score',
                                                'q_value', 'network_score', 'network_q_value')]
    f['e'] = f.precursor_id.map(ent).fillna(0).astype(int)
    print(f'\n== {path}: windows {len(f)}')
    best = {}
    for tag, s, q in [('net', 'network_score', 'network_q_value'), ('gbm', 'score', 'q_value')]:
        b = f[f[s] > -1e29].sort_values(s, ascending=False).drop_duplicates('precursor_id').copy()
        best[tag] = b
        t = b[b.is_decoy == 0]
        acc = t[t[q] <= 0.01]
        N, E = len(acc), int(acc.e.sum())
        print(f'{tag}: q<=0.01 N {N} E {E} combined {E * k / max(N, 1):.4f}')
        # Decoy count vs entrapment-estimated false count among the top-T by score (pairs competed: use q).
        for T in (500, 1000, 1500, 2000, 2500, 3000):
            top = b.sort_values(s, ascending=False)
            top_t = top[top.is_decoy == 0].head(T)
            if len(top_t) < T: break
            thr = top_t[s].min()
            D = int(((top.is_decoy == 1) & (top[s] >= thr)).sum())
            Ef = top_t.e.sum() * k
            print(f'   top {T} targets: decoys above threshold {D}, entrapment-estimated false {Ef:.0f}, ratio {Ef / max(D, 1):.2f}')
    n, g = best['net'], best['gbm']
    nt = set(n[(n.is_decoy == 0) & (n.network_q_value <= 0.01)].precursor_id)
    gt = set(g[(g.is_decoy == 0) & (g.q_value <= 0.01)].precursor_id)
    for name, ids in [('gbm-only', gt - nt), ('net-only', nt - gt), ('both', gt & nt)]:
        e = sum(ent.get(i, 0) for i in ids)
        print(f'{name}: {len(ids)} targets, E {e}, combined {e * k / max(len(ids), 1):.4f}')
    # Which features separate entrapment winners from decoy winners near the threshold (network q <= 0.05)?
    near = n[n.network_q_value <= 0.05]
    E_rows, D_rows = near[(near.is_decoy == 0) & (near.e == 1)], near[near.is_decoy == 1]
    T_rows = near[(near.is_decoy == 0) & (near.e == 0)]
    print(f'near threshold (net q<=0.05): entrapment {len(E_rows)}, decoy {len(D_rows)}, real targets {len(T_rows)}')
    if len(E_rows) > 20 and len(D_rows) > 20:
        out = []
        for x in feats:
            sd = np.sqrt(0.5 * (E_rows[x].var() + D_rows[x].var())) or 1.0
            out.append((x, (E_rows[x].mean() - D_rows[x].mean()) / sd))
        out.sort(key=lambda v: -abs(v[1]))
        print('features, entrapment minus decoy (pooled SD units):', ', '.join(f'{a} {b:+.2f}' for a, b in out[:12]))
