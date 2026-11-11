"""Diagnostic (not a gated arm; recorded in PLAN_V4.notes): retrain D11 exactly as train_v4.py with alpha = 0 and no MS1
(same seeds, so it should reproduce scores_D11.npy), and split each window's score s = z.w_z + tanh(g(z)).w_t into its
linear part and its tanh part. Reports, per raw: correlation with scores_D11; the window-choice rate against the DIA-NN RT
(reference only) for s, the linear part alone, and the tanh part alone; the share of score variance from each part.
Usage: diag_readout.py <raw01|raw02>. Writes <raw>/scores_D11lin.npy and <raw>/scores_D11tanh.npy."""
import sys, inspect, numpy as np, pandas as pd, torch
src = open('/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2/v4/train_v4.py').read()
exec(src.replace('for a in ARMS: run(a)', ''))
lin = np.zeros(len(pop), np.float32); tnh = np.zeros(len(pop), np.float32); full = np.zeros(len(pop), np.float32)
for k in range(3):
    tr_all = np.where(fold != k)[0]
    mu, sd = H[tr_all].mean(0), H[tr_all].std(0); sd[sd == 0] = 1; Z = (H - mu) / sd
    t = train_rows[(train_rows.fold != k) & train_rows.group.isin(['anchor', 'anchor_decoy', 'decoy'])]
    rows = row_of[t.sample_id].values; y = t.label.values.astype(np.float32)
    valid = np.array([fnv('v' + str(s)) % 100 < 15 for s in t.family.values])
    sc, ep, vo = fit(rows, y, Z, False, True, 100 + k, valid, 100, 0.0)
    for r in (1, 2):
        win, q, dec, wtab = pooled_q(tr_all, sc(tr_all))
        pos = win[(dec == 0) & (q <= (0.15 if r == 1 else 0.05))]; negs = wtab[wtab.dec == 1].r.values
        rows = np.concatenate([pos, negs]); y = np.r_[np.ones(len(pos)), np.zeros(len(negs))].astype(np.float32)
        valid = np.array([fnv('v' + str(s)) % 100 < 15 for s in pair_key[rows]])
        sc, ep, vo = fit(rows, y, Z, False, True, 200 + 10 * k + r, valid, 100, 0.0)
    cv = inspect.getclosurevars(sc).nonlocals; net, w = cv['net'], cv['w']
    te = np.where(fold == k)[0]; h = embed(net, te, Z, False); d = Z.shape[1]
    lin[te] = (h[:, :d] @ w[:d]).numpy(); tnh[te] = (h[:, d:] @ w[d:]).numpy(); full[te] = sc(te)
    print(RAW, 'fold', k, 'w_z norm', float(w[:d].norm()), 'w_t norm', float(w[d:].norm()), flush=True)
np.save(f'{O}/{RAW}/scores_D11lin.npy', lin); np.save(f'{O}/{RAW}/scores_D11tanh.npy', tnh)
ref = np.load(f'{O}/{RAW}/scores_D11.npy')
print(RAW, 'reproduces scores_D11: max abs diff', float(np.abs(full - ref).max()), 'corr', float(np.corrcoef(full, ref)[0, 1]))
print(RAW, 'variance share (all windows): linear', float(lin.var() / full.var()), 'tanh', float(tnh.var() / full.var()),
      'corr(lin, tanh)', float(np.corrcoef(lin, tnh)[0, 1]))
st = pd.read_csv(f'{O}/v3seq/misses_D11_{RAW}.tsv', sep='\t')
idx = pop.reset_index().merge(st[['precursor_id', 'rt_s']], on='precursor_id')
idx['right'] = (idx.rt_apex - idx.rt_s).abs() <= 10
idx = idx[idx.groupby('precursor_id').right.transform('any')]
for name, s in (('D11', full), ('linear part', lin), ('tanh part', tnh), ('lin + tanh / 10', lin + tnh / 10)):
    x = idx.assign(s=s[idx['index'].values]).sort_values(['s', 'sample_id'], ascending=[False, True]).drop_duplicates('precursor_id')
    print(RAW, f'window choice right ({name}):', int(x.right.sum()), 'of', len(x), flush=True)
