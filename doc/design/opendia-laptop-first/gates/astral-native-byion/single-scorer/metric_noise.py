"""Diagnostic (reported only): how noisy is N@1% (entrapment FDP_c) on the design raws?
Paired bootstrap over target precursors (in-window best window per arm): resample precursors with replacement, recompute
N at FDP_c 1% and 2% for base, R1 and R2, and the paired difference to base. 1000 replicates, seed 7."""
import sqlite3, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
O = B + '/astral-single-scorer-20261109'; R = 148983 / 151239
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c)
BASE = {'raw01': B + '/astral-design-raw01-20261107/stage2-byion-design/byion', 'raw02': B + '/astral-raw02-byion-design-20261109/byion'}
def n_at(score, e, t):
    o = np.argsort(-score, kind='stable'); e = e[o]
    fdp = np.cumsum(e) * (1 + 1 / R) / np.arange(1, len(e) + 1); ok = np.where(fdp <= t)[0]
    return int(ok.max() + 1) if len(ok) else 0
rng = np.random.default_rng(7)
for raw in ('raw01', 'raw02'):
    arms = {'base': BASE[raw], 'R1': f'{O}/{raw}/R1', 'R2': f'{O}/{raw}/R2'}
    tab = None
    for a, d in arms.items():
        cr = pd.read_csv(d + '/crossfit.tsv', sep='\t', usecols=['precursor_id', 'is_decoy', 'precursor_best', 'in_rt_window', 'stage2_score'])
        w = cr[(cr.precursor_best == 1) & (cr.is_decoy == 0) & (cr.in_rt_window == 1)][['precursor_id', 'stage2_score']].rename(columns={'stage2_score': a})
        tab = w if tab is None else tab.merge(w, on='precursor_id', how='outer')
    tab = tab.merge(ent, on='precursor_id').fillna(-np.inf)
    e = tab.e.values.astype(int); n = len(tab)
    point = {a: [n_at(tab[a].values, e, t) for t in (0.01, 0.02)] for a in arms}
    boot = {a: [] for a in arms}
    for _ in range(1000):
        idx = rng.integers(0, n, n)
        for a in arms: boot[a].append([n_at(tab[a].values[idx], e[idx], t) for t in (0.01, 0.02)])
    boot = {a: np.array(v) for a, v in boot.items()}
    print(f'{raw}: target precursors {n}, entrapment {e.sum()}')
    for a in arms:
        lo1, hi1 = np.percentile(boot[a][:, 0], [5, 95]); lo2, hi2 = np.percentile(boot[a][:, 1], [5, 95])
        line = f'  {a:4s} N@1% {point[a][0]} (90% {lo1:.0f}-{hi1:.0f}) | N@2% {point[a][1]} (90% {lo2:.0f}-{hi2:.0f})'
        if a != 'base':
            d1 = boot[a][:, 0] - boot['base'][:, 0]; d2 = boot[a][:, 1] - boot['base'][:, 1]
            line += f' | diff to base @1%: 90% {np.percentile(d1, 5):.0f}..{np.percentile(d1, 95):.0f}, @2%: {np.percentile(d2, 5):.0f}..{np.percentile(d2, 95):.0f}'
        print(line, flush=True)
