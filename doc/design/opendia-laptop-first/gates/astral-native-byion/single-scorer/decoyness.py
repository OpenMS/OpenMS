"""Diagnostic (reported only): does a scorer separate decoy winners from entrapment target winners (both false)?
Per arm: AUC(entrapment winner score > decoy winner score) within each fold, using in-window precursor winners.
AUC 0.5 = decoys model false targets (TDC valid); AUC > 0.5 = the scorer learned decoy-specific features (anti-conservative q).
The LDA base uses stage2_score from the native byion runs."""
import sqlite3, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c)
RUNS = {'raw01': {'LDA': (B + '/astral-design-raw01-20261107/stage2-byion-design/byion/crossfit.tsv', 'stage2_score')},
        'raw02': {'LDA': (B + '/astral-raw02-byion-design-20261109/byion/crossfit.tsv', 'stage2_score')}}
for raw in RUNS:
    for arm in ('C0', 'N1', 'L1', 'A1', 'T1', 'M1'):
        RUNS[raw][arm] = (f'{raw}/{arm}/crossfit.tsv', 'logit')
for raw, arms in RUNS.items():
    for arm, (path, col) in arms.items():
        try:
            cr = pd.read_csv(path, sep='\t')
        except FileNotFoundError:
            continue
        w = cr[cr.in_rt_window == 1].sort_values([col, 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').merge(ent, on='precursor_id', how='left')
        aucs = []
        for f, g in w.groupby('fold'):
            d = g[g.is_decoy == 1][col].values; e = g[(g.is_decoy == 0) & (g.e == 1)][col].values
            aucs.append(mannwhitneyu(e, d).statistic / (len(e) * len(d)))
        print(f'{raw} {arm:4s}: AUC(entrapment winner > decoy winner) per fold ' + ' '.join(f'{a:.3f}' for a in aucs), flush=True)
