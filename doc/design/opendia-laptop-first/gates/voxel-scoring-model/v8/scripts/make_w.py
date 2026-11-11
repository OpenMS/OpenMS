"""PLAN_V4B W rule: per precursor, keep the window with the highest linear part (z.w_z of the same trained model) and
give it the model's full score; other windows of the precursor get -1e30. Same for targets and decoys, no labels.
Usage: make_w.py <raw01|raw02> <src arm> <lin file> <out arm>. Ties broken by sample_id (as the evaluator)."""
import sys, numpy as np, pandas as pd
RAW, SRC, LIN, OUT = sys.argv[1:5]
O = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/rtlin-20261111'  # rtlin rebuild root
pop = pd.read_csv(f'{O}/{RAW}/population.tsv', sep='\t', usecols=['sample_id', 'precursor_id'])
s = np.load(f'{O}/{RAW}/scores_{SRC}.npy'); lin = np.load(f'{O}/{RAW}/{LIN}.npy')
order = pop.assign(l=lin, r=np.arange(len(pop))).sort_values(['l', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').r.values
out = np.full(len(pop), -1e30, np.float32); out[order] = s[order]
np.save(f'{O}/{RAW}/scores_{OUT}.npy', out)
print(RAW, OUT, 'from', SRC, LIN, 'precursors', len(order), 'kept window = full-score argmax for', round(float(np.isin(order,
      pop.assign(s=s, r=np.arange(len(pop))).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id').r.values).mean()), 4))
