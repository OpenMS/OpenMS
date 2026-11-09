"""Diagnostic (reported only): which inputs separate false target winners (entrapment) from decoy winners?
Per heuristic feature (native byion features.tsv, cnn_z = grid2d CNN logit), AUC(entrapment > decoy) over the N1 TOPAZ
in-window precursor winners and over the LDA winners. |AUC - 0.5| >= 0.02 marks a decoy-separating input."""
import sqlite3, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c)
BASE = {'raw01': B + '/astral-design-raw01-20261107/stage2-byion-design/byion', 'raw02': B + '/astral-raw02-byion-design-20261109/byion'}
def auc(e, d): return mannwhitneyu(e, d).statistic / (len(e) * len(d))
for raw, base in BASE.items():
    f = pd.read_csv(base + '/features.tsv', sep='\t')
    cf = pd.read_csv(base + '/crossfit.tsv', sep='\t', usecols=['sample_id', 'in_rt_window'])
    n1 = pd.read_csv(f'{raw}/N1/crossfit.tsv', sep='\t', usecols=['sample_id', 'logit'])
    f = f.merge(cf, on='sample_id').merge(n1, on='sample_id')
    f = f[f.in_rt_window == 1].merge(ent, on='precursor_id', how='left')
    feats = [x for x in f.columns if x not in ('sample_id', 'precursor_id', 'is_decoy', 'fold', 'in_rt_window', 'e')]
    for winner_by in ('logit', 'stage2_score'):
        w = f.sort_values([winner_by, 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id')
        e = w[(w.is_decoy == 0) & (w.e == 1)]; d = w[w.is_decoy == 1]
        row = ' '.join(f'{x}={auc(e[x].values, d[x].values):.3f}' for x in feats)
        print(f'{raw} winners by {"N1" if winner_by == "logit" else "LDA"} (entrapment {len(e)}, decoy {len(d)}): {row}', flush=True)
