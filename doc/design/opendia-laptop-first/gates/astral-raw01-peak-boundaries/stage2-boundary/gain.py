"""Where do the search+boundary IDs differ from search? Gained vs lost at q<=0.01 by signal level (PLAN.md reading:
the learned segmenter is only worth porting if the gain is at low SNR, where the rule picker fails)."""
import sys
import numpy as np, pandas as pd, sqlite3
sys.path.insert(0, '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-design-raw01-20261107/cnn-topaz-design')
from evaluate import paired_pooled

D = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-design-raw01-20261107'
LIB = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-entrapment-20261107/lib/library.pqp'
c = sqlite3.connect('file:' + LIB + '?mode=ro', uri=True)
p = pd.read_sql('select ID precursor_id, TRAML_ID, DECOY from PRECURSOR', c)
p['pair_key'] = np.where(p.DECOY == 1, p.TRAML_ID.str.replace('^DECOY_', '', regex=True), p.TRAML_ID)
acc = {}
for name, path in [('search', D + '/stage2-search-design/search'), ('boundary', D + '/stage2-boundary-design/boundary')]:
    df = pd.read_csv(path + '/crossfit.tsv', sep='\t').merge(p[['precursor_id', 'pair_key']], on='precursor_id', how='left')
    best, _ = paired_pooled(df, 'stage2_score')
    acc[name] = best[(best.is_decoy == 0) & (best.q <= 0.01)][['precursor_id', 'sample_id']]
s = pd.read_csv(D + '/pipeline/search.features.tsv', sep='\t', usecols=['candidate_row', 'apex_intensity', 'peak_width_seconds']).rename(columns={'candidate_row': 'sample_id'})
a, b = set(acc['search'].precursor_id), set(acc['boundary'].precursor_id)
print(f'accepted search {len(a)}, boundary {len(b)}, both {len(a & b)}, gained {len(b - a)}, lost {len(a - b)}')
both = pd.concat([acc['search'], acc['boundary']]).drop_duplicates('precursor_id').merge(s, on='sample_id')
both['grp'] = np.where(both.precursor_id.isin(a & b), 'both', np.where(both.precursor_id.isin(b - a), 'gained', 'lost'))
both['log_apex'] = np.log10(both.apex_intensity.clip(lower=1))
print(both.groupby('grp')[['log_apex', 'peak_width_seconds']].median().round(2).assign(n=both.groupby('grp').size()))
