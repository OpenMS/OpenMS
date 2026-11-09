"""Direct test of 'the decoys cause the separation' (reported only). Same input-level AUC(entrapment winner > decoy winner)
for the native base library vs the k11 commonness-matched library (same targets, different decoys), stage-2 LDA winners.
Also: decoy vs ALL target winners for comparison (mostly absent targets)."""
import sqlite3, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
def ent(lib):
    c = sqlite3.connect('file:' + lib + '?mode=ro', uri=True)
    return pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
      from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c)
def auc(e, d): return mannwhitneyu(e, d).statistic / (len(e) * len(d))
RUNS = [('raw01 base', B + '/astral-design-raw01-20261107/stage2-byion-design/byion', B + '/astral-entrapment-20261107/lib/library.pqp'),
        ('raw01 k11 ', B + '/astral-commatched-decoys-20261109/raw01/crossfit', B + '/astral-commatched-decoys-20261109/library.pqp'),
        ('raw02 base', B + '/astral-raw02-byion-design-20261109/byion', B + '/astral-entrapment-20261107/lib/library.pqp'),
        ('raw02 k11 ', B + '/astral-commatched-decoys-20261109/raw02/crossfit', B + '/astral-commatched-decoys-20261109/library.pqp')]
FEATS = ['log_matched_scans', 'transition_support', 'boundary_library_correlation', 'cnn_z', 'byion_b', 'stage2_score']
for name, d, lib in RUNS:
    f = pd.read_csv(d + '/features.tsv', sep='\t').merge(pd.read_csv(d + '/crossfit.tsv', sep='\t', usecols=['sample_id', 'in_rt_window']), on='sample_id')
    f = f[f.in_rt_window == 1].merge(ent(lib), on='precursor_id', how='left')
    w = f.sort_values(['stage2_score', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id')
    e = w[(w.is_decoy == 0) & (w.e == 1)]; t = w[w.is_decoy == 0]; dd = w[w.is_decoy == 1]
    print(f'{name} entrapment-vs-decoy: ' + ' '.join(f'{x}={auc(e[x].values, dd[x].values):.3f}' for x in FEATS) +
          f' | all-target-vs-decoy stage2={auc(t.stage2_score.values, dd.stage2_score.values):.3f} log_matched_scans={auc(t.log_matched_scans.values, dd.log_matched_scans.values):.3f}', flush=True)
