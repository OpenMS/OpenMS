"""Voxel pilot step 1 (no runs): the full raw 01 window set the pilot exports and trains on, following DESIGN.md.
Inputs are calibration anchors, the label-free search candidate table and the library only; entrapment status is
read only to write the separate evaluation file. Seed 13.

Training groups (one row per candidate window; candidate_row = window sample_id):
  anchor       label 1  the anchor's candidate window nearest its calibrated run RT (within 5 s)
  anchor_decoy label 0  the paired decoy's best window by transition support, then matched scans
  decoy        label 0  random decoy windows (any rank)
  decoy_like   label 0  (weight 0.5) target windows that are not their precursor's best by transition support, or whose
                        apex is > 3 sigma (calibration RT residual) from the library RT; precursor not an anchor
Evaluation (never trained on): the make_sample.py sample (present / entrap / decoy / low_target).
Folds: FNV-1a 64 of the parent target's unmodified sequence mod 3, so a target and its decoy share a fold."""
import json, sqlite3, numpy as np, pandas as pd
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
O = B + '/astral-voxel-pilot-20261110'
P = B + '/astral-design-raw01-20261107/pipeline'
CAL = B + '/astral-design-raw01-20261107/window300k'
rng = np.random.default_rng(13)
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
pre = pd.read_sql("""select p.ID precursor_id, p.TRAML_ID traml, p.DECOY decoy, pe.UNMODIFIED_SEQUENCE seq from PRECURSOR p
  join PRECURSOR_PEPTIDE_MAPPING m on m.PRECURSOR_ID=p.ID join PEPTIDE pe on pe.ID=m.PEPTIDE_ID""", c)
target_seq = pre[pre.decoy == 0].set_index('traml').seq
pre['parent'] = np.where(pre.decoy == 1, pre.traml.str.replace('^DECOY_', '', regex=True), pre.traml)
pre['family'] = pre.parent.map(target_seq).fillna(pre.seq)
def fnv(s):
    h = 0xcbf29ce484222325
    for b in s.encode(): h = ((h ^ b) * 0x100000001b3) & 0xFFFFFFFFFFFFFFFF
    return h
pre['fold'] = pre.family.map(lambda s: fnv(s) % 3)
pair = pre.set_index('traml').precursor_id
pre['partner'] = np.where(pre.decoy == 0, ('DECOY_' + pre.traml).map(pair), pre.parent.map(pair))

cand = pd.read_csv(P + '/search.features.tsv', sep='\t', usecols=['candidate_row', 'precursor_id', 'is_decoy', 'rt_apex', 'transition_support',
                                                                  'matched_scans', 'has_library_rt', 'apex_rt_error_seconds'])
cand = cand.rename(columns={'candidate_row': 'sample_id'})
# Best window per precursor by the search's own label-free statistics.
cand = cand.sort_values(['precursor_id', 'transition_support', 'matched_scans', 'sample_id'], ascending=[True, False, False, True])
cand['rank'] = cand.groupby('precursor_id').cumcount()
cal = json.load(open(CAL + '/cal.json'))
sigma = float(cal['rt'].get('sigma_s', cal['rt'].get('residual_sigma_s', np.nan)))
anchors = pd.read_csv(CAL + '/anchors.tsv', sep='\t')
anchors = anchors[anchors.used == 1]

rows = []
# Anchors: the candidate window nearest the anchor run RT.
a = cand.merge(anchors[['precursor_id', 'run_rt']], on='precursor_id')
a['d'] = (a.rt_apex - a.run_rt).abs()
a = a[a.d <= 5.0].sort_values('d').drop_duplicates('precursor_id')
rows.append(a.assign(group='anchor', label=1, weight=1.0))
# Paired decoys of the anchors: their best window.
partners = pre.set_index('precursor_id').partner.reindex(a.precursor_id).dropna().astype(int)
ad = cand[cand.precursor_id.isin(partners.values) & (cand['rank'] == 0)]
rows.append(ad.assign(group='anchor_decoy', label=0, weight=1.0))
# Random decoy windows.
dec = cand[(cand.is_decoy == 1) & ~cand.sample_id.isin(ad.sample_id)]
rows.append(dec.iloc[rng.choice(len(dec), 6000, replace=False)].assign(group='decoy', label=0, weight=1.0))
# Decoy-like target windows.
anchor_ids = set(anchors.precursor_id)
far = (cand.has_library_rt == 1) & (cand.apex_rt_error_seconds.abs() > 3 * sigma) if np.isfinite(sigma) else False
dl = cand[(cand.is_decoy == 0) & ~cand.precursor_id.isin(anchor_ids) & ((cand['rank'] > 0) | far)]
rows.append(dl.iloc[rng.choice(len(dl), 6000, replace=False)].assign(group='decoy_like', label=0, weight=0.5))
train = pd.concat(rows, ignore_index=True)
# Evaluation sample (make_sample.py): base winners; never a training row.
ev = pd.read_csv(O + '/sample.tsv', sep='\t')[['precursor_id', 'sample_id']]
ev = ev.merge(pd.read_csv(O + '/sample_eval_groups.tsv', sep='\t').drop_duplicates('precursor_id'), on='precursor_id')
# No overlap filter: evaluation is crossfit (fold k scored by the model trained on the other folds), so an anchor
# window may also be an evaluation window without being seen by the model that scores it.
train = train.merge(pre[['precursor_id', 'fold', 'family']], on='precursor_id')
ev = ev.merge(pre[['precursor_id', 'fold', 'family']], on='precursor_id').merge(cand[['sample_id', 'rt_apex']], on='sample_id')
train[['sample_id', 'precursor_id', 'is_decoy', 'rt_apex', 'group', 'label', 'weight', 'fold', 'family']].to_csv(O + '/train_windows.tsv', sep='\t', index=False)
ev[['sample_id', 'precursor_id', 'rt_apex', 'fold', 'family']].to_csv(O + '/eval_windows.tsv', sep='\t', index=False)
ev[['sample_id', 'precursor_id', 'group', 'e']].to_csv(O + '/eval_labels.tsv', sep='\t', index=False)
with open(O + '/export_precursors.txt', 'w') as f:
    f.write('# voxel pilot raw 01 export set (make_train_set.py, seed 13)\n')
    for pid in sorted(set(train.precursor_id) | set(ev.precursor_id)): f.write(f'{pid}\n')
print('sigma_s', sigma, train.groupby('group').size().to_dict(), 'eval', len(ev),
      'export precursors', len(set(train.precursor_id) | set(ev.precursor_id)), 'folds', train.fold.value_counts().to_dict())
