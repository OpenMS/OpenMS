"""Structural asymmetry check: do decoys get fewer chances than false (entrapment) targets? Design diagnostic.

Under TDC, a decoy must be an exchangeable copy of a false target. Per precursor, compare across human targets,
entrapment targets and decoys:
  - library fragment slots filled (of 6) and their library-intensity entropy
  - XIC points (npts) and extracted RT span
  - number of peak candidates (of up to 5)
  - fraction of fragment traces that are all-zero
  - precursor m/z and charge distribution
Any difference between entrapment and decoys here means the competition is not between exchangeable copies.
"""
import os, numpy as np, pandas as pd
RUN = os.environ.get('OPENDIA_RUN', '01')
S = os.environ['JCODE_SCRATCH_DIR'] + '/xic' + RUN
meta = pd.read_parquet(S + '/meta.parquet'); lib = np.load(S + '/lib.npy'); npts = np.load(S + '/npts.npy')
rt = np.load(S + '/rt.npy', mmap_mode='r'); frag = np.load(S + '/frag.npy', mmap_mode='r')
om = pd.read_parquet(S + '/oswlike_meta.parquet')
ent = om.groupby('precursor_id').ent.first().reindex(meta.precursor_id).fillna(0).astype(bool).values
ncand = om.groupby('precursor_id').size().reindex(meta.precursor_id).fillna(0).values
dec = meta.decoy.values.astype(bool)
grp = np.where(dec, 'decoy', np.where(ent, 'entrapment', 'human'))
slots = (lib > 0).sum(1)
p = lib / np.maximum(lib.sum(1, keepdims=True), 1e-12)
H = -(np.where(p > 0, p * np.log(np.where(p > 0, p, 1)), 0)).sum(1)
span = np.nanmax(rt, 1) - np.nanmin(rt, 1)
zero = np.zeros(len(meta))
for b in range(0, len(meta), 50000):
    f = np.asarray(frag[b:b + 50000]); zero[b:b + 50000] = ((f.max(2) <= 0) & (lib[b:b + 50000] > 0)).sum(1) / np.maximum(slots[b:b + 50000], 1)
d = pd.DataFrame({'grp': grp, 'slots': slots, 'lib_entropy': H, 'npts': npts, 'rt_span': span, 'n_cand': ncand, 'zero_trace_frac': zero})
print(f'raw {RUN}: precursors per group', d.grp.value_counts().to_dict())
print(d.groupby('grp').agg(['mean']).round(3).T.to_string())
print('share with <6 slots:', d.groupby('grp').slots.apply(lambda s: (s < 6).mean()).round(4).to_dict())
print('share with <5 candidates:', d.groupby('grp').n_cand.apply(lambda s: (s < 5).mean()).round(4).to_dict())
