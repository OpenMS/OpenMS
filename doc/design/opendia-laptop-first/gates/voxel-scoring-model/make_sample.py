"""Voxel pilot step 0 (no runs): choose the raw 01 precursor sample whose raw points the pilot will export with OpenDIA
-out_xipm. Groups, from the native base (CNN + LDA) raw 01 outputs, in-window precursor winners by stage2_score:
  present   : base target winners with the tool's decoy-based q <= 0.01 (no entrapment information; training positives)
  entrap    : entrapment target winners (absent peptides; evaluation only, never a training label)
  decoy     : decoy winners (training negatives)
  low_target: target winners with q > 0.05 (mixed, mostly absent; unlabeled)
Sizes: all present; 3000 each of the others, drawn with seed 11. Training labels use only q and decoy status; the
entrapment group (and which present precursors are entrapment) is written to a separate evaluation file."""
import sqlite3, numpy as np, pandas as pd, os
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
O = B + '/astral-voxel-pilot-20261110'; os.makedirs(O, exist_ok=True)
R = 148983 / 151239
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
ent = pd.read_sql("""select pp.PRECURSOR_ID precursor_id, min(case when pr.PROTEIN_ACCESSION like '%_p_target%' then 1 else 0 end) e
  from PRECURSOR_PEPTIDE_MAPPING pp join PEPTIDE_PROTEIN_MAPPING m on m.PEPTIDE_ID=pp.PEPTIDE_ID join PROTEIN pr on pr.ID=m.PROTEIN_ID group by pp.PRECURSOR_ID""", c)
base = B + '/astral-design-raw01-20261107/stage2-byion-design/byion'
cr = pd.read_csv(base + '/crossfit.tsv', sep='\t', usecols=['sample_id', 'precursor_id', 'is_decoy', 'fold', 'precursor_best', 'q_value', 'in_rt_window', 'stage2_score'])
w = cr[(cr.precursor_best == 1) & (cr.in_rt_window == 1)].merge(ent, on='precursor_id', how='left')
t = w[w.is_decoy == 0]
rng = np.random.default_rng(11)
def pick(df, k): return df.iloc[rng.choice(len(df), size=min(k, len(df)), replace=False)]
groups = {
    'present': t[t.q_value <= 0.01],
    'entrap': pick(t[(t.e == 1) & (t.q_value > 0.01)], 3000),
    'decoy': pick(w[w.is_decoy == 1], 3000),
    'low_target': pick(t[(t.q_value > 0.05) & (t.e == 0)], 3000),
}
rows = pd.concat([g.assign(group=k) for k, g in groups.items()], ignore_index=True)
rows[['precursor_id', 'sample_id', 'fold', 'is_decoy', 'q_value', 'stage2_score', 'group']].assign(
    train_label=np.where(rows.group == 'decoy', 0, np.where(rows.group == 'present', 1, -1))
).drop(columns='group').to_csv(O + '/sample.tsv', sep='\t', index=False)
# Evaluation-only file: group and entrapment status (present precursors that are entrapment are false positives of base).
rows[['precursor_id', 'group', 'e']].to_csv(O + '/sample_eval_groups.tsv', sep='\t', index=False)
with open(O + '/xipm_precursors.txt', 'w') as f:
    f.write('# voxel pilot raw 01 sample (make_sample.py, seed 11)\n')
    for pid in rows.precursor_id.drop_duplicates(): f.write(f'{pid}\n')
print({k: len(g) for k, g in groups.items()}, 'present that are entrapment', int(groups['present'].e.sum()),
      'unique precursors', rows.precursor_id.nunique())
