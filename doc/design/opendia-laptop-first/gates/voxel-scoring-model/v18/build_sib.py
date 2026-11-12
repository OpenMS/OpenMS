"""V18 step 1: charge-sibling features S (3 columns) for every window. Label-free: uses the library peptide mapping,
window apex RTs, the existing held-out D11F full scores, W (lin_D11F) winners and the decoy flag for the per-fold
normalisation (as the evaluator). No entrapment or DIA-NN labels. Usage: build_sib.py <raw>. Writes v2/<raw>/featS.npy."""
import sys, sqlite3, numpy as np, pandas as pd
from sib_core import fold_z, sib_features
RAW = sys.argv[1]
O = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2'
L = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-entrapment-20261107/lib/library.pqp'
c = sqlite3.connect(f'file:{L}?mode=ro', uri=True)
m = pd.read_sql("select PRECURSOR_ID precursor_id, PEPTIDE_ID peptide_id from PRECURSOR_PEPTIDE_MAPPING", c)
assert m.precursor_id.is_unique
pop = pd.read_csv(f'{O}/{RAW}/population.tsv', sep='\t')
pep = pop.precursor_id.map(m.set_index('precursor_id').peptide_id)
assert pep.notna().all()
# Siblings share the peptide, hence the family and the fold: features of held-out rows use held-out scores only.
assert (pd.DataFrame({'p': pep.values, 'f': pop.fold.values}).groupby('p').f.nunique() == 1).all()
full = np.load(f'{O}/{RAW}/scores_D11F.npy'); lin = np.load(f'{O}/{RAW}/lin_D11F.npy')
z = fold_z(full, lin, pop.precursor_id.values, pop.sample_id.values, pop.is_decoy.values, pop.fold.values)
S = sib_features(pop.precursor_id.values, pep.values.astype(np.int64), pop.rt_apex.values, z)
assert S.shape == (len(pop), 3) and np.isfinite(S).all()
np.save(f'{O}/{RAW}/featS.npy', S.astype(np.float32))
for name, d in (('target', pop.is_decoy.values == 0), ('decoy', pop.is_decoy.values == 1)):
    print(RAW, name, 'windows', int(d.sum()), 'has_sib', round(float(S[d, 0].mean()), 4), 'has_near', round(float(S[d, 1].mean()), 4),
          'sib_z_near mean (near)', round(float(S[d & (S[:, 1] > 0), 2].mean()), 4), flush=True)
