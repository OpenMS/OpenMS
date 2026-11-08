"""RT-precision lever in the FULL search universe (raw 01 design diagnostic; DIA-NN used only to simulate a predictor
and to define the oracle's positives, never as a pipeline signal).

Oracle score s = cross-fitted MLP on all current non-RT features (positives: DIA-NN-confirmed right peaks; negatives:
decoy candidates). Final score = s + Gaussian RT log-likelihood with width sigma, candidates restricted to +-3 sigma.
Predicted RT: confirmed targets get DIA-NN RT + N(0, sigma) (a simulated predictor of that precision); every other
target and every decoy keeps the calibrated library RT (decoy = null, any centre is equivalent). Paired TDC over all
287K pairs; report all accepted targets and the confirmed ones among them. sigma = 24 s with the calibrated library RT
for everyone is the current situation.
"""
import sqlite3, numpy as np, pandas as pd, torch, torch.nn as nn
torch.manual_seed(0); torch.set_num_threads(5); rng = np.random.default_rng(0)
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
D = B + '/astral-design-raw01-20261107'
f = pd.read_csv(D + '/stage2-boundary-design/boundary/features.tsv', sep='\t')
x = pd.read_csv(D + '/stage2-design/features.tsv', sep='\t', usecols=['sample_id', 'fragment_fraction', 'fragment_coelution', 'library_correlation', 'ms2_abs_ppm', 'ms1_abs_ppm', 'ms1_coelution', 'ms1_isotope'])
cr = pd.read_csv(D + '/stage2-boundary-design/boundary/crossfit.tsv', sep='\t', usecols=['sample_id', 'in_rt_window'])
s = pd.read_csv(D + '/pipeline/search.features.tsv', sep='\t', usecols=['candidate_row', 'rt_apex', 'apex_rt_error_seconds']).rename(columns={'candidate_row': 'sample_id'})
df = f.merge(x, on='sample_id').merge(cr, on='sample_id').merge(s, on='sample_id')
df['cal_rt'] = df.rt_apex - df.apex_rt_error_seconds
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
p = pd.read_sql('select ID precursor_id, TRAML_ID, DECOY from PRECURSOR', c)
p['pair_key'] = np.where(p.DECOY == 1, p.TRAML_ID.str.replace('^DECOY_', '', regex=True), p.TRAML_ID)
df = df.merge(p[['precursor_id', 'pair_key']], on='precursor_id')
ref = pd.read_parquet(B + '/astral-diann-gap-20261108/ref_exact.parquet')
ref = ref[ref.Run.str.endswith('_01')][['precursor_id', 'RT']]
ref['diann_rt'] = ref.RT * 60.0
df = df.merge(ref[['precursor_id', 'diann_rt']], on='precursor_id', how='left')
conf = df.diann_rt.notna()
right = conf & ((df.rt_apex - df.diann_rt).abs() <= 10)
FEAT = ['cnn_z', 'transition_support', 'log_matched_scans', 'log_peak_width'] + [c for c in f.columns if c.startswith('boundary_')] + \
       ['fragment_fraction', 'fragment_coelution', 'library_correlation', 'ms2_abs_ppm', 'ms1_abs_ppm', 'ms1_coelution', 'ms1_isotope']
X = df[FEAT].to_numpy(np.float32); mu, sd = X.mean(0), X.std(0) + 1e-9; X = (X - mu) / sd
pos = df[right].sort_values('cnn_z', ascending=False).drop_duplicates('precursor_id').index.values
dec = np.where(df.is_decoy.values == 1)[0]
fold = df.fold.values
score = np.zeros(len(df), np.float32)
for k in np.unique(fold):
    pi = pos[fold[pos] != k]; ni = dec[fold[dec] != k]; ni = rng.choice(ni, min(len(ni), 400000), replace=False)
    Xt = torch.tensor(np.r_[X[pi], X[ni]]); yt = torch.tensor(np.r_[np.ones(len(pi)), np.zeros(len(ni))], dtype=torch.float32)
    net = nn.Sequential(nn.Linear(len(FEAT), 32), nn.GELU(), nn.Linear(32, 16), nn.GELU(), nn.Linear(16, 1))
    opt = torch.optim.AdamW(net.parameters(), lr=3e-3, weight_decay=1e-3); lf = nn.BCEWithLogitsLoss()
    for _ in range(6):
        perm = torch.randperm(len(Xt))
        for i in range(0, len(Xt), 4096):
            b = perm[i:i + 4096]; l = lf(net(Xt[b]).squeeze(1), yt[b]); opt.zero_grad(); l.backward(); opt.step()
    m = fold == k
    with torch.no_grad():
        score[m] = net(torch.tensor(X[m])).squeeze(1).numpy()
df['s'] = score

def tdc(t, col):
    w = t.sort_values(col, ascending=False).drop_duplicates('precursor_id')
    w = w.sort_values([col, 'is_decoy'], ascending=[False, False]).drop_duplicates('pair_key').sort_values(col, ascending=False)
    d = w.is_decoy.values; T = np.cumsum(1 - d); Dc = np.cumsum(d)
    q = np.minimum.accumulate((np.minimum(1, (Dc + 1) / np.maximum(T, 1)))[::-1])[::-1]
    return w.assign(q=q)

def report(name, w):
    out = []
    for qq in (0.01, 0.05):
        a = w[(w.is_decoy == 0) & (w.q <= qq)]
        out.append(f'q<={qq}: all {len(a)} confirmed-right {int(((a.rt_apex - a.diann_rt).abs() <= 10).sum())}')
    print(f'{name:40s}', ' | '.join(out), flush=True)

report('oracle, no RT term (all candidates)', tdc(df, 's'))
pk = df.drop_duplicates('pair_key')[['pair_key']]
for sigma in (24.0, 12.0, 6.0, 3.0):
    pr = df[['precursor_id', 'cal_rt', 'diann_rt']].drop_duplicates('precursor_id')
    if sigma == 24.0:
        pr['pred'] = pr.cal_rt
    else:
        pr['pred'] = np.where(pr.diann_rt.notna(), pr.diann_rt + rng.normal(0, sigma, len(pr)), pr.cal_rt)
    t = df.merge(pr[['precursor_id', 'pred']], on='precursor_id')
    z = (t.rt_apex - t.pred) / sigma
    t = t[z.abs() <= 3].assign(fs=lambda d: d.s - 0.5 * ((d.rt_apex - d.pred) / sigma) ** 2)
    report(f'oracle + RT sigma {sigma:4.1f} s' + (' (current)' if sigma == 24 else ' (confirmed only)'), tdc(t, 'fs'))
