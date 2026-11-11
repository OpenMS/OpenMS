"""Voxel v2 step 1 (no runs): scoring population, folds and training set for the full-population voxel arms.
Label-free: uses only the search candidate table (search.features.tsv), the calibration and the library.
Population: every candidate window of every target-decoy pair whose best transition_support over both members is
>= 3. On raw 01 that is 94,109 pairs, 186,702 precursors and 1,270,510 windows (of 2.8M). Reason: compute (the raw-point
export and the voxels of 2.8M windows do not fit this laptop). Effect on base: all base winners accepted at FDP_c 2 %
have transition_support >= 3 (checked), so the base N at FDP_c <= 2 % is computed on the same population.
Folds: FNV-1a 64 of the parent target's unmodified sequence mod 3 (as make_train_set.py).
Training groups: the make_train_set.py rules, restricted to the population, seed 13."""
import json, os, sqlite3, sys, numpy as np, pandas as pd
RAW = sys.argv[1] if len(sys.argv) > 1 else 'raw01'
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
O = B + '/rtlin-20261111/' + RAW  # rtlin rebuild root
os.makedirs(O, exist_ok=True)
# Search table and calibration of each raw (the base pipeline inputs).
P = B + '/rtlin-20261111/' + RAW + '/pipeline'
CAL = B + '/rtlin-20261111/' + RAW + '/cal'
# Base stage-2 features (search 4, boundary 5, byion 2) of each raw, same row order as the search table.
BASEF = B + '/rtlin-20261111/' + RAW + '/byion/features.tsv'
FEATS = ['transition_support', 'log_matched_scans', 'log_peak_width', 'abs_library_rt_error', 'boundary_log_fwhm', 'boundary_log_width',
         'boundary_coelution', 'boundary_library_correlation', 'boundary_apex_agreement', 'byion_b', 'byion_y']
rng = np.random.default_rng(13)
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
pre = pd.read_sql("""select p.ID precursor_id, p.TRAML_ID traml, p.DECOY decoy, pe.UNMODIFIED_SEQUENCE seq from PRECURSOR p
  join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID""", c)
target_seq = pre[pre.decoy == 0].set_index('traml').seq
pre['pair_key'] = np.where(pre.decoy == 1, pre.traml.str.replace('^DECOY_', '', regex=True), pre.traml)
pre['family'] = pre.pair_key.map(target_seq).fillna(pre.seq)
def fnv(s):
    h = 0xcbf29ce484222325
    for b in s.encode(): h = ((h ^ b) * 0x100000001b3) & 0xFFFFFFFFFFFFFFFF
    return h
pre['fold'] = pre.family.map(lambda s: fnv(s) % 3)
cand = pd.read_csv(P + '/search.features.tsv', sep='\t', usecols=['candidate_row', 'precursor_id', 'is_decoy', 'rt_apex', 'transition_support',
                                                                  'matched_scans', 'has_library_rt', 'apex_rt_error_seconds'])
cand = cand.rename(columns={'candidate_row': 'sample_id'}).merge(pre[['precursor_id', 'pair_key', 'fold', 'family']], on='precursor_id')
pair_max = cand.groupby('pair_key').transition_support.max()
cand = cand[cand.pair_key.isin(pair_max[pair_max >= 3].index)].copy()
cand = cand.sort_values(['precursor_id', 'transition_support', 'matched_scans', 'sample_id'], ascending=[True, False, False, True])
cand['rank'] = cand.groupby('precursor_id').cumcount()
popt = cand.sort_values('sample_id')[['sample_id', 'precursor_id', 'is_decoy', 'rt_apex', 'fold', 'pair_key']].reset_index(drop=True)
popt.to_csv(O + '/population.tsv', sep='\t', index=False)
bf = pd.read_csv(BASEF, sep='\t', usecols=['sample_id', 'precursor_id'] + FEATS).set_index('sample_id').loc[popt.sample_id]
assert (bf.precursor_id.values == popt.precursor_id.values).all()
np.save(O + '/heuristics.npy', bf[FEATS].values.astype(np.float32))

sigma = float(json.load(open(CAL + '/cal.json'))['rt']['sigma_s'])
anchors = pd.read_csv(CAL + '/anchors.tsv', sep='\t'); anchors = anchors[anchors.used == 1]
rows = []
a = cand.merge(anchors[['precursor_id', 'run_rt']], on='precursor_id'); a['d'] = (a.rt_apex - a.run_rt).abs()
a = a[a.d <= 5.0].sort_values('d').drop_duplicates('precursor_id')
rows.append(a.assign(group='anchor', label=1, weight=1.0))
part = pre.set_index('precursor_id')
partner = ('DECOY_' + part.loc[a.precursor_id, 'traml']).map(pre.set_index('traml').precursor_id).dropna().astype(int)
ad = cand[cand.precursor_id.isin(partner.values) & (cand['rank'] == 0)]
rows.append(ad.assign(group='anchor_decoy', label=0, weight=1.0))
dec = cand[(cand.is_decoy == 1) & ~cand.sample_id.isin(ad.sample_id)]
rows.append(dec.iloc[rng.choice(len(dec), 6000, replace=False)].assign(group='decoy', label=0, weight=1.0))
far = (cand.has_library_rt == 1) & (cand.apex_rt_error_seconds.abs() > 3 * sigma)
dl = cand[(cand.is_decoy == 0) & ~cand.precursor_id.isin(set(anchors.precursor_id)) & ((cand['rank'] > 0) | far)]
rows.append(dl.iloc[rng.choice(len(dl), 6000, replace=False)].assign(group='decoy_like', label=0, weight=0.5))
train = pd.concat(rows, ignore_index=True)
train[['sample_id', 'precursor_id', 'is_decoy', 'rt_apex', 'group', 'label', 'weight', 'fold', 'family']].to_csv(O + '/train_windows.tsv', sep='\t', index=False)
with open(O + '/export_precursors.txt', 'w') as f:
    f.write(f'# voxel v2 {RAW} population: pairs with max transition_support >= 3 (make_population.py)\n')
    for pid in sorted(cand.precursor_id.unique()): f.write(f'{pid}\n')
print('pairs', cand.pair_key.nunique(), 'precursors', cand.precursor_id.nunique(), 'windows', len(cand), 'sigma', sigma,
      train.groupby('group').size().to_dict(), 'anchors used', len(anchors))
