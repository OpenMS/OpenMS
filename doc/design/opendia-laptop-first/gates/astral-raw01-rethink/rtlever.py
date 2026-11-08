"""How much would a sharper RT prediction buy? (raw 01 design diagnostic; DIA-NN RT used only to SIMULATE a predictor.)

Restricted to the 2754 library pairs whose target DIA-NN confirms (target + its paired decoy). For each pair a simulated
predicted RT = DIA-NN RT + N(0, sigma) is shared by target and decoy (decoys share the library RT). Candidates are
limited to +-3 sigma, scored by an RT-free oracle (cross-fitted MLP on all current non-RT features, positives = confirmed
right peaks, negatives = decoy candidates) plus the Gaussian-vs-uniform RT log-likelihood ratio. Winner per precursor,
paired TDC within these pairs, then count confirmed targets accepted at q<=0.01 with the winner within 10 s of DIA-NN RT.
Absolute counts are optimistic (oracle labels, DIA-NN-selected pairs); the trend across sigma is the result.
"""
import sqlite3, numpy as np, pandas as pd, torch, torch.nn as nn
torch.manual_seed(0); torch.set_num_threads(5); rng = np.random.default_rng(0)
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
D = B + '/astral-design-raw01-20261107'
f = pd.read_csv(D + '/stage2-boundary-design/boundary/features.tsv', sep='\t')
x = pd.read_csv(D + '/stage2-design/features.tsv', sep='\t', usecols=['sample_id', 'fragment_fraction', 'fragment_coelution', 'library_correlation', 'ms2_abs_ppm', 'ms1_abs_ppm', 'ms1_coelution', 'ms1_isotope'])
cr = pd.read_csv(D + '/stage2-boundary-design/boundary/crossfit.tsv', sep='\t', usecols=['sample_id', 'in_rt_window'])
s = pd.read_csv(D + '/pipeline/search.features.tsv', sep='\t', usecols=['candidate_row', 'rt_apex']).rename(columns={'candidate_row': 'sample_id'})
df = f.merge(x, on='sample_id').merge(cr, on='sample_id').merge(s, on='sample_id')
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
p = pd.read_sql('select ID precursor_id, TRAML_ID, DECOY from PRECURSOR', c)
p['pair_key'] = np.where(p.DECOY == 1, p.TRAML_ID.str.replace('^DECOY_', '', regex=True), p.TRAML_ID)
df = df.merge(p[['precursor_id', 'pair_key']], on='precursor_id')
ref = pd.read_parquet(B + '/astral-diann-gap-20261108/ref_exact.parquet')
ref = ref[ref.Run.str.endswith('_01')][['precursor_id', 'RT']]
ref['diann_rt'] = ref.RT * 60.0
ref = ref.merge(p[['precursor_id', 'pair_key']], on='precursor_id')[['pair_key', 'diann_rt']]
sub = df.merge(ref, on='pair_key')
print('pairs', ref.pair_key.nunique(), 'candidates in those pairs', len(sub), 'target', int((sub.is_decoy == 0).sum()))
FEAT = ['cnn_z', 'transition_support', 'log_matched_scans', 'log_peak_width'] + [c for c in f.columns if c.startswith('boundary_')] + \
       ['fragment_fraction', 'fragment_coelution', 'library_correlation', 'ms2_abs_ppm', 'ms1_abs_ppm', 'ms1_coelution', 'ms1_isotope']
right = (sub.is_decoy == 0) & ((sub.rt_apex - sub.diann_rt).abs() <= 10)
pos = sub[right].sort_values('cnn_z', ascending=False).drop_duplicates('precursor_id').index
# negatives: decoy candidates from ALL decoys (not only these pairs) to keep it general, subsampled
negall = df[df.is_decoy == 1].sample(200000, random_state=0)
Xs = sub[FEAT].to_numpy(np.float64); Xn = negall[FEAT].to_numpy(np.float64)
mu, sd = Xn.mean(0), Xn.std(0) + 1e-9
Xs = (Xs - mu) / sd; Xn = (Xn - mu) / sd
score = np.full(len(sub), np.nan)
for k in sorted(df.fold.unique()):
    pi = np.array([i for i in pos if sub.fold[i] != k]); ni = negall.fold.values != k
    Xtr = np.r_[Xs[sub.index.get_indexer(pi)], Xn[ni]]; y = np.r_[np.ones(len(pi)), np.zeros(ni.sum())]
    net = nn.Sequential(nn.Linear(len(FEAT), 32), nn.GELU(), nn.Linear(32, 16), nn.GELU(), nn.Linear(16, 1))
    opt = torch.optim.AdamW(net.parameters(), lr=3e-3, weight_decay=1e-3)
    Xt = torch.tensor(Xtr, dtype=torch.float32); yt = torch.tensor(y, dtype=torch.float32)
    lf = nn.BCEWithLogitsLoss()
    for _ in range(6):
        perm = torch.randperm(len(Xt))
        for i in range(0, len(Xt), 4096):
            b = perm[i:i + 4096]; l = lf(net(Xt[b]).squeeze(1), yt[b]); opt.zero_grad(); l.backward(); opt.step()
    m = sub.fold.values == k
    with torch.no_grad():
        score[m] = net(torch.tensor(Xs[m], dtype=torch.float32)).squeeze(1).numpy()
sub = sub.assign(s=score)
# prior logit shift: pos/neg ratio in training is ~ (2449*2/3)/(133k); keep scores as-is (common shift cancels in ranking)
for sigma in (40.0, 20.0, 10.0, 5.0, 2.5):
    for rep in range(2):
        pred = ref.assign(pred=ref.diann_rt + rng.normal(0, sigma, len(ref)))[['pair_key', 'pred']]
        t = sub.merge(pred, on='pair_key')
        err = (t.rt_apex - t.pred).abs()
        t = t[err <= 3 * sigma].assign(fs=lambda d: d.s - 0.5 * ((d.rt_apex - d.pred) / sigma) ** 2 + np.log(6 / np.sqrt(2 * np.pi)))
        w = t.sort_values('fs', ascending=False).drop_duplicates('precursor_id')
        w = w.sort_values(['fs', 'is_decoy'], ascending=[False, False]).drop_duplicates('pair_key').sort_values('fs', ascending=False)
        d = w.is_decoy.values; T = np.cumsum(1 - d); Dc = np.cumsum(d)
        q = np.minimum.accumulate((np.minimum(1, (Dc + 1) / np.maximum(T, 1)))[::-1])[::-1]
        acc = w[(d == 0) & (q <= 0.01)]
        ok = int(((acc.rt_apex - acc.diann_rt).abs() <= 10).sum())
        acc5 = w[(d == 0) & (q <= 0.05)]
        print(f'sigma {sigma:5.1f} s rep {rep}: accepted q<=0.01 {len(acc)} (right RT {ok}) | q<=0.05 {len(acc5)}', flush=True)
