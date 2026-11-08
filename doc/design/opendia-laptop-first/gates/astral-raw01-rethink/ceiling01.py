"""Raw 01 ceiling funnel and an oracle separability bound (design diagnostic; DIA-NN is reference only, never a training
signal for the pipeline). Positives here = DIA-NN 1% precursors that exactly exist in our library, at DIA-NN RT +-10 s.

A. Funnel: candidate near DIA-NN RT, winner there, accepted at q<=0.01 by the tool's stage-2 (search+boundary).
B. Oracle: supervised cross-fit (pair-grouped folds) of DIA-NN-confirmed right-peak candidates vs decoy winners, with
   (i) LDA and (ii) MLP on all current features. Report how many confirmed targets rank above the decoy level that
   would give 1% / 5% FDR. If the oracle cannot recover most of the ceiling, information (not learning) is the limit.
"""
import sqlite3, numpy as np, pandas as pd, torch, torch.nn as nn
torch.manual_seed(0); torch.set_num_threads(5); np.random.seed(0)
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
D = B + '/astral-design-raw01-20261107'
f = pd.read_csv(D + '/stage2-boundary-design/boundary/features.tsv', sep='\t')
x = pd.read_csv(D + '/stage2-design/features.tsv', sep='\t', usecols=['sample_id', 'fragment_fraction', 'fragment_coelution', 'library_correlation', 'ms2_abs_ppm', 'ms1_abs_ppm', 'ms1_coelution', 'ms1_isotope'])
cr = pd.read_csv(D + '/stage2-boundary-design/boundary/crossfit.tsv', sep='\t', usecols=['sample_id', 'in_rt_window', 'precursor_best', 'q_value'])
s = pd.read_csv(D + '/pipeline/search.features.tsv', sep='\t', usecols=['candidate_row', 'rt_apex', 'apex_intensity', 'area']).rename(columns={'candidate_row': 'sample_id'})
df = f.merge(x, on='sample_id').merge(cr, on='sample_id').merge(s, on='sample_id')
df['q'] = pd.to_numeric(df.q_value, errors='coerce')
c = sqlite3.connect('file:' + B + '/astral-entrapment-20261107/lib/library.pqp?mode=ro', uri=True)
p = pd.read_sql('select ID precursor_id, TRAML_ID, DECOY from PRECURSOR', c)
p['pair_key'] = np.where(p.DECOY == 1, p.TRAML_ID.str.replace('^DECOY_', '', regex=True), p.TRAML_ID)
df = df.merge(p[['precursor_id', 'pair_key']], on='precursor_id')
ref = pd.read_parquet(B + '/astral-diann-gap-20261108/ref_exact.parquet')
ref = ref[ref.Run.str.endswith('_01')][['precursor_id', 'RT', 'entrap']].rename(columns={'RT': 'diann_rt'})
ref['diann_rt'] *= 60.0
print('DIA-NN confirmed in-library precursors raw 01:', len(ref), 'entrap', int(ref.entrap.sum()))

# A. funnel
cand = df.merge(ref, on='precursor_id')
cand['near'] = (cand.rt_apex - cand.diann_rt).abs() <= 10
anyc = ref.precursor_id.isin(df.precursor_id)
near = cand.groupby('precursor_id').near.any()
win = cand[cand.precursor_best == 1]
print('funnel: any candidate', int(anyc.sum()), '| candidate within 10 s', int(near.sum()),
      '| winner within 10 s', int(win.near.sum()), '| accepted q<=0.01', int((win.near & (win.q <= 0.01)).sum()),
      '| q<=0.05', int((win.near & (win.q <= 0.05)).sum()))
# intensity profile of confirmed right-peak winners, accepted or not
rw = win[win.near].copy()
rw['tert'] = pd.qcut(np.log(rw.apex_intensity), 3, labels=['low', 'mid', 'high'])
print('recall at q<=0.01 by intensity tertile:', rw.groupby('tert', observed=True).apply(lambda g: round((g.q <= 0.01).mean(), 2)).to_dict())

# B. oracle separability (supervised, diagnostic only)
FEAT = ['cnn_z', 'transition_support', 'log_matched_scans', 'log_peak_width', 'abs_library_rt_error'] + \
       [c for c in f.columns if c.startswith('boundary_')] + \
       ['fragment_fraction', 'fragment_coelution', 'library_correlation', 'ms2_abs_ppm', 'ms1_abs_ppm', 'ms1_coelution', 'ms1_isotope']
pos_ids = set(cand[cand.near].sample_id)
dw = df[(df.is_decoy == 1) & (df.precursor_best == 1)]
pos = df[df.sample_id.isin(pos_ids)].sort_values('cnn_z', ascending=False).drop_duplicates('precursor_id')
data = pd.concat([pos.assign(y=1), dw.assign(y=0)])
X = data[FEAT].to_numpy(np.float64); y = data.y.values; fold = data.fold.values
mu, sd = X.mean(0), X.std(0) + 1e-9; X = (X - mu) / sd
def lda(Xa, ya):
    a, b = Xa[ya == 1], Xa[ya == 0]
    S = (np.cov(a.T) * (len(a) - 1) + np.cov(b.T) * (len(b) - 1)) / (len(Xa) - 2) + 1e-3 * np.eye(Xa.shape[1])
    w = np.linalg.solve(S, a.mean(0) - b.mean(0)); return lambda Z: Z @ w
def mlp(Xa, ya):
    net = nn.Sequential(nn.Linear(Xa.shape[1], 32), nn.GELU(), nn.Linear(32, 16), nn.GELU(), nn.Linear(16, 1))
    opt = torch.optim.AdamW(net.parameters(), lr=3e-3, weight_decay=1e-3)
    Xt = torch.tensor(Xa, dtype=torch.float32); yt = torch.tensor(ya, dtype=torch.float32)
    lf = nn.BCEWithLogitsLoss(pos_weight=torch.tensor((ya == 0).sum() / (ya == 1).sum()))
    for _ in range(8):
        perm = torch.randperm(len(Xt))
        for i in range(0, len(Xt), 2048):
            b = perm[i:i + 2048]; l = lf(net(Xt[b]).squeeze(1), yt[b]); opt.zero_grad(); l.backward(); opt.step()
    return lambda Z: net(torch.tensor(Z, dtype=torch.float32)).squeeze(1).detach().numpy()
n_dec_total = len(dw)
for name, fitter in [('LDA', lda), ('MLP', mlp)]:
    sc = np.full(len(data), np.nan)
    for k in np.unique(fold):
        m = fold == k
        sc[m] = fitter(X[~m], y[~m])(X[m])
    # decoy-score thresholds: keep top decoys equal to 1% / 5% of the confirmed set (~ FDR-like level)
    dsc = np.sort(sc[y == 0])[::-1]; tsc = sc[y == 1]
    out = []
    for frac in (0.001, 0.01, 0.05):
        # threshold such that #decoys above / #confirmed above ~ frac
        order = np.sort(np.r_[tsc, dsc])[::-1]
        best = 0
        for thr in order[::50]:
            T = (tsc >= thr).sum(); Dd = (dsc >= thr).sum()
            if T and Dd / T <= frac: best = max(best, T)
        out.append(best)
    print(f'oracle {name}: confirmed recovered at decoy/target ratio 0.1% / 1% / 5% =', out, 'of', int(y.sum()))
