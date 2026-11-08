"""Fragment-count ablation, raw 01 (design). Same search on the 6-transition library and its top-4 reduction.
Compares (a) the DIA-NN-confirmed funnel and (b) a supervised oracle on the search-level features only (no CNN/stage-2
available for the top-4 run) at decoy/target ratio 1% and 5%. A drop from 6 to 4 indicates how much fragment count
carries; it is suggestive (not proof) for the gain from 6 -> 12."""
import sqlite3, sys, numpy as np, pandas as pd, torch, torch.nn as nn
torch.set_num_threads(5)
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'; D = B + '/astral-design-raw01-20261107'
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
p = pd.read_sql('select ID precursor_id, TRAML_ID, DECOY from PRECURSOR', c)
p['pair_key'] = np.where(p.DECOY == 1, p.TRAML_ID.str.replace('^DECOY_', '', regex=True), p.TRAML_ID)
ref = pd.read_parquet(B + '/astral-diann-gap-20261108/ref_exact.parquet'); ref = ref[ref.Run.str.endswith('_01')]
ref = ref[['precursor_id']].assign(diann_rt=ref.RT.values * 60)
COLS = ['candidate_row', 'precursor_id', 'is_decoy', 'rt_apex', 'peak_width_seconds', 'area', 'apex_intensity', 'transition_support', 'matched_scans', 'apex_rt_error_seconds']
def load(path):
    s = pd.read_csv(path, sep='\t', usecols=COLS).merge(p[['precursor_id', 'pair_key']], on='precursor_id')
    s['fold'] = pd.util.hash_pandas_object(s.pair_key, index=False).values % 3
    s = s.merge(ref, on='precursor_id', how='left')
    s['ts_frac'] = s.transition_support / s.groupby('precursor_id').transition_support.transform('max').clip(lower=1)
    X = np.c_[s.transition_support, np.log1p(s.matched_scans), np.log1p(s.peak_width_seconds), s.apex_rt_error_seconds.abs(),
              np.log1p(s.apex_intensity), np.log1p(s.area)].astype(np.float32)
    return s, X
def oracle(s, X, seed=0):
    torch.manual_seed(seed); rng = np.random.default_rng(seed)
    right = s.diann_rt.notna() & ((s.rt_apex - s.diann_rt).abs() <= 10)
    pos = s[right].sort_values('transition_support', ascending=False).drop_duplicates('precursor_id').index.values
    dec = np.where(s.is_decoy.values == 1)[0]
    mu, sd = X.mean(0), X.std(0) + 1e-9; Z = (X - mu) / sd
    sc = np.zeros(len(s), np.float32); fold = s.fold.values
    for k in range(3):
        pi = pos[fold[pos] != k]; ni = dec[fold[dec] != k]; ni = rng.choice(ni, min(len(ni), 400000), replace=False)
        Xt = torch.tensor(np.r_[Z[pi], Z[ni]]); yt = torch.tensor(np.r_[np.ones(len(pi)), np.zeros(len(ni))], dtype=torch.float32)
        net = nn.Sequential(nn.Linear(X.shape[1], 32), nn.GELU(), nn.Linear(32, 16), nn.GELU(), nn.Linear(16, 1))
        opt = torch.optim.AdamW(net.parameters(), lr=3e-3, weight_decay=1e-3); lf = nn.BCEWithLogitsLoss()
        for _ in range(5):
            perm = torch.randperm(len(Xt))
            for i in range(0, len(Xt), 4096):
                b = perm[i:i + 4096]; l = lf(net(Xt[b]).squeeze(1), yt[b]); opt.zero_grad(); l.backward(); opt.step()
        with torch.no_grad():
            sc[fold == k] = net(torch.tensor(Z[fold == k])).squeeze(1).numpy()
    s = s.assign(sc=sc)
    w = s.sort_values('sc', ascending=False).drop_duplicates('precursor_id')
    w = w.sort_values(['sc', 'is_decoy'], ascending=[False, False]).drop_duplicates('pair_key').sort_values('sc', ascending=False)
    d = w.is_decoy.values; T = np.cumsum(1 - d); Dc = np.cumsum(d)
    q = np.minimum.accumulate((np.minimum(1, (Dc + 1) / np.maximum(T, 1)))[::-1])[::-1]
    out = []
    for qq in (0.01, 0.05):
        a = w[(d == 0) & (q <= qq)]; out.append((len(a), int(((a.rt_apex - a.diann_rt).abs() <= 10).sum())))
    return out
for name, path in [('6 transitions', D + '/pipeline/search.features.tsv'), ('top-4 transitions', D + '/rethink/frag4/search.features.tsv')]:
    s, X = load(path)
    cand = s[s.diann_rt.notna()]
    near = cand.assign(n=(cand.rt_apex - cand.diann_rt).abs() <= 10).groupby('precursor_id').n.any()
    print(f'{name}: candidates {len(s)} | confirmed with any candidate {cand.precursor_id.nunique()} | within 10 s {int(near.sum())}', flush=True)
    for seed in (0, 1):
        print(f'   oracle(search features) seed {seed}: q<=0.01 / 0.05 (accepted, confirmed-right) = {oracle(s, X, seed)}', flush=True)
