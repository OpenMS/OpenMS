#!/usr/bin/env python3
"""Cross-fold relative test. Pseudo-reverse decoys keep the substring structure of the target set (nested and
missed-cleavage peptides reverse into decoys that share long substrings, so shared fragments). Entrapment shuffles do
not. Folds are by unmodified family, so related decoys can sit in different folds. If the GBM memorises its hard decoy
negatives, test-fold decoys with a related hard decoy in a training fold are pushed down, entrapments are not, and TDC
under-counts false targets. Hard set = top 2 % of decoy precursors by network score (the round-1 rule).
Relative = shares a 7-mer of unmodified sequence. Diagnosis only."""
import sqlite3, sys, numpy as np, pandas as pd
from collections import defaultdict
lib = sys.argv[1]
K = 7
c = sqlite3.connect(f'file:{lib}?mode=ro', uri=True)
pq = pd.read_sql("""select p.ID precursor_id, p.DECOY d, pe.UNMODIFIED_SEQUENCE seq,
  min(case when instr(pr.PROTEIN_ACCESSION, '_p_target') > 0 then 1 else 0 end) e
  from PRECURSOR p join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID
  join PEPTIDE_PROTEIN_MAPPING pp on pp.PEPTIDE_ID=m.PEPTIDE_ID join PROTEIN pr on pr.ID=pp.PROTEIN_ID group by p.ID""", c)
P = pq.set_index('precursor_id')
tot = defaultdict(lambda: np.zeros(3))
for path in sys.argv[2:]:
    f = pd.read_csv(path, sep='\t', usecols=['precursor_id', 'is_decoy', 'fold', 'score', 'network_score'])
    n = f[f.network_score > -1e29].set_index('precursor_id')[['network_score', 'fold', 'is_decoy']]
    g = f[f.score > -1e29].set_index('precursor_id')[['score']]
    b = n.join(g, how='inner').join(P[['seq', 'e']])
    dec = b[b.is_decoy == 1]
    hard = dec.sort_values('network_score', ascending=False).head(int(np.ceil(0.02 * len(dec))))
    idx = defaultdict(set)  # k-mer -> folds of hard decoys holding it
    for s, fo in zip(hard.seq, hard.fold):
        for i in range(len(s) - K + 1):
            idx[s[i:i + K]].add(fo)
    def rel(s, fo):
        return any(any(x != fo for x in idx.get(s[i:i + K], ())) for i in range(len(s) - K + 1))
    for s_ in ('network_score', 'score'):
        dv = np.sort(dec[s_].values)
        b[s_ + '_u'] = -np.log10((len(dv) - np.searchsorted(dv, b[s_].values, side='left') + 1) / (len(dv) + 1))
    b['dsh'] = b.score_u - b.network_score_u
    tail = b[(b.network_score_u >= 1.0) & ((b.is_decoy == 1) | (b.e == 1))].copy()
    tail['rel'] = [rel(s, fo) for s, fo in zip(tail.seq, tail.fold)]
    tail['cls'] = np.where(tail.is_decoy == 1, 'D', 'E')
    print(f'== {path.split("/")[0]}: hard decoys {len(hard)}; network-tail (u>=1) rows {len(tail)}')
    for (cl, r), x in tail.groupby(['cls', 'rel']):
        hi = x[x.network_score_u >= 2]
        print(f'   {cl} cross-fold hard relative={r!s:5s}: n {len(x):6d} mean GBM shift {x.dsh.mean():+.3f} | net u>=2: n {len(hi):5d} shift {hi.dsh.mean():+.3f}, '
              f'GBM u>=2.5 kept {int((hi.score_u >= 2.5).sum())}')
        t = tot[(cl, r)]
        t += [len(hi), hi.dsh.sum(), int((hi.score_u >= 2.5).sum())]
print('pooled, network u>=2:')
for k, (nn, ss, kk) in sorted(tot.items()):
    print(f'   {k}: n {int(nn)} mean shift {ss / max(nn, 1):+.3f} GBM u>=2.5 {int(kk)} ({kk / max(nn, 1):.3f})')
