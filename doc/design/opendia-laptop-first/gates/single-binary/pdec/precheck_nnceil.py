"""PRECHECK_NNCEIL: can a learned trace model add information over the 36 features for faint peaks?
Sample (faint_mini.py): 6,000 DIA-NN targets + 30,000 decoys, r~1.56 D s0, search valley 0.8 / max 8 (= ev/v08m8 arm).
Input per candidate window: the search's own trace matrix (6 fragments x native scans, mobility-pooled within the
calibrated window), 33 scans centred on the apex; channels = log1p(I)/max, I/max (per window), sqrt library intensity
(broadcast). A small 1D CNN (two conv layers, max+mean pool) is trained per workflow fold (cross-fitted).
Arms, TPR of each DIA-NN target's right window vs decoy precursors' best windows at fixed decoy FPR:
  G36   GBM on the 36 workflow features, DIA-NN-labelled (the PRECHECK_FAINT ceiling, recomputed here)
  CNN   CNN alone, DIA-NN-labelled
  G+C   GBM on 36 features + out-of-fold CNN logit (stacked, same folds)
  CNNlf CNN trained label-free: positives = tool q<=0.01 target winners, negatives = decoy windows
  G+Clf GBM (DIA-NN labels) on 36 + label-free CNN logit: does the label-free trace model carry the information?
Labels for evaluation only; DIA-NN-labelled arms are ceilings, not methods."""
import glob, sqlite3, time, numpy as np, pandas as pd, torch, torch.nn as nn
from sklearn.ensemble import HistGradientBoostingClassifier
torch.set_num_threads(16); torch.manual_seed(0); rng = np.random.default_rng(0)
G = '/home/singjc/.jcode/scratch/dia2/gap7'; O = f'{G}/faint'
S = pd.read_csv(f'{O}/sample.tsv', sep='\t').set_index('pid')
P = pd.read_csv(f'{O}/peaks.tsv', sep='\t'); P['rta'] = P.rt_apex.round(3)
c = sqlite3.connect(f'{O}/mini.pqp')
lib = pd.read_sql('select m.PRECURSOR_ID pid, t.PRODUCT_MZ mz, t.LIBRARY_INTENSITY li from TRANSITION t join TRANSITION_PRECURSOR_MAPPING m on m.TRANSITION_ID=t.ID', c)
libg = {p: g for p, g in lib.groupby('pid')}
t0 = time.time()
T = pd.read_csv(f'{O}/traces.tsv', sep='\t', dtype={'precursor_id': np.int32, 'scan': np.int32, 'slot': np.int8})
H = 16; W = 2 * H + 1
X, keys = [], []
for pid, t in T.groupby('precursor_id', sort=False):
    m = t.pivot_table(index='rt', columns='slot', values='intensity', aggfunc='max').fillna(0.0)
    fm = t.drop_duplicates('slot').set_index('slot').fragment_mz.reindex(m.columns).values
    lg = libg[pid]; li = np.array([lg.li.values[np.argmin(np.abs(lg.mz.values - x))] for x in fm])
    w = np.sqrt(np.maximum(li, 0)); w = w / (w.max() + 1e-12)
    rt = m.index.values; A = m.values.T.astype(np.float32)  # fragments x scans
    if A.shape[0] < 6: A = np.vstack([A, np.zeros((6 - A.shape[0], A.shape[1]), np.float32)]); w = np.r_[w, np.zeros(6 - len(w))]
    A = np.pad(A, ((0, 0), (H, H)))
    for ra in P.rta.values[P.precursor_id.values == pid]:
        i = int(np.argmin(np.abs(rt - ra))) + H
        patch = A[:, i - H:i + H + 1]; mx = patch.max() + 1e-6
        X.append(np.concatenate([np.log1p(patch) / np.log1p(mx), patch / mx, np.repeat(w[:, None], W, 1)], 0).astype(np.float32))
        keys.append((pid, ra))
X = np.stack(X); K = pd.DataFrame(keys, columns=['precursor_id', 'rta'])
print(f'patches {X.shape} in {time.time() - t0:.0f} s', flush=True)
ids = set(S.index); parts = []
for ch in pd.read_csv(glob.glob(f'{G}/ev/v08m8/*.features.tsv')[0], sep='\t', chunksize=1_000_000):
    parts.append(ch[ch.precursor_id.isin(ids)])
F = pd.concat(parts); F['rta'] = F.rt_apex.round(3)
K['k'] = np.arange(len(K)); F = F.merge(K.drop_duplicates(['precursor_id', 'rta']), on=['precursor_id', 'rta'], how='inner')
X = X[F.k.values]
meta = ['candidate_row', 'precursor_id', 'is_decoy', 'fold', 'rt_apex', 'score', 'q_value', 'rta', 'k']
X36 = [c for c in F.columns if c not in meta]
dec = F.is_decoy.values == 1; isd = F.precursor_id.map(S.cls).values == 'dnn'
pos = isd & (np.abs(F.rt_apex.values - F.precursor_id.map(S.rt_s).values) <= 6)
fold = F.fold.values
# label-free positives: tool q<=0.01 target winners (best kept window per precursor)
F['_s'] = np.where(F.score > -1e29, F.score, -1e30)
win = F.sort_values('_s', ascending=False).drop_duplicates('precursor_id')
lf_rows = set(win.index[(win.is_decoy == 0) & (win.q_value <= 0.01)])
pos_lf = F.index.isin(lf_rows)
print(f'windows {len(F)}, features {len(X36)}, DIA-NN positives {pos.sum()}, label-free positives {pos_lf.sum()}, decoy windows {dec.sum()}', flush=True)

class Net(nn.Module):
    def __init__(self, cin=18, h=32):
        super().__init__()
        self.c = nn.Sequential(nn.Conv1d(cin, h, 5, padding=2), nn.ReLU(), nn.Conv1d(h, h, 5, padding=2), nn.ReLU(), nn.Conv1d(h, h, 3, padding=1), nn.ReLU())
        self.head = nn.Sequential(nn.Linear(3 * h, 32), nn.ReLU(), nn.Linear(32, 1))
    def forward(self, x):
        z = self.c(x); centre = z[:, :, H - 2:H + 3].mean(2)
        return self.head(torch.cat([z.max(2).values, z.mean(2), centre], 1)).squeeze(1)

Xt = torch.from_numpy(X.reshape(len(X), 18, W))
def cnn(posm, seed=0, epochs=8, steps=300):
    out = np.zeros(len(F), np.float32)
    for f in np.unique(fold):
        torch.manual_seed(seed + int(f)); g = np.random.default_rng(seed + int(f))
        tr = fold != f; tp = np.where(tr & posm)[0]; tn = np.where(tr & dec)[0]
        va_p, va_n = tp[g.random(len(tp)) < 0.1], tn[g.random(len(tn)) < 0.1]
        tp = np.setdiff1d(tp, va_p); tn = np.setdiff1d(tn, va_n)
        net = Net(); opt = torch.optim.AdamW(net.parameters(), 1e-3, weight_decay=1e-3); lossf = nn.BCEWithLogitsLoss()
        best, best_state = 1e9, None
        for ep in range(epochs):
            net.train()
            for _ in range(steps):
                b = np.r_[g.choice(tp, 128), g.choice(tn, 128)]; y = torch.tensor(np.r_[np.ones(128), np.zeros(128)], dtype=torch.float32)
                opt.zero_grad(); l = lossf(net(Xt[b]), y); l.backward(); opt.step()
            net.eval()
            with torch.no_grad():
                vb = np.r_[va_p, g.choice(va_n, min(len(va_n), 20000), replace=False)]
                vy = torch.tensor(np.r_[np.ones(len(va_p)), np.zeros(len(vb) - len(va_p))], dtype=torch.float32)
                vl = float(lossf(net(Xt[vb]), vy))
            if vl < best: best, best_state = vl, {k: v.clone() for k, v in net.state_dict().items()}
        net.load_state_dict(best_state); net.eval()
        te = np.where(fold == f)[0]
        with torch.no_grad():
            for s in range(0, len(te), 8192): out[te[s:s + 8192]] = net(Xt[te[s:s + 8192]]).numpy()
        print(f'  fold {f}: best val loss {best:.4f}', flush=True)
    return out

def gbm(cols_or_mat, posm, seed=0):
    Xv = cols_or_mat if isinstance(cols_or_mat, np.ndarray) else F[cols_or_mat].values.astype(np.float32)
    s = np.zeros(len(F))
    for f in np.unique(fold):
        tr = (fold != f) & (posm | dec); y = posm[tr].astype(int); wts = np.where(y == 1, (y == 0).sum() / max(y.sum(), 1), 1.0)
        m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=31, random_state=seed).fit(Xv[tr], y, sample_weight=wts)
        s[fold == f] = m.decision_function(Xv[fold == f])
    return s

qq = S.qq; dec_ids = S.index[S.cls == 'decoy']
def report(name, s):
    F['_r'] = s; best = F.groupby('precursor_id')._r.max(); rb = F[pos].groupby('precursor_id')._r.max()
    neg = best.reindex(dec_ids).dropna().values; out = []
    for fpr in (1e-3, 3e-3, 1e-2):
        thr = np.quantile(neg, 1 - fpr)
        out.append(f'FPR {fpr:g}: ' + ' '.join(f'{k} {(rb[rb.index.map(qq) == k] > thr).mean():.3f}' for k in ['Q1', 'Q2', 'Q3', 'Q4']) + f' all {(rb > thr).mean():.3f}')
    print(f'{name:8s} ' + ' | '.join(out), flush=True)
    return rb, neg

report('tool', F._s.values)
g36 = gbm(X36, pos); report('G36', g36)
c_dn = cnn(pos); report('CNN', c_dn)
M = np.c_[F[X36].values.astype(np.float32), c_dn]; gc = gbm(M, pos); report('G+C', gc)
c_lf = cnn(pos_lf, seed=100); report('CNNlf', c_lf)
M2 = np.c_[F[X36].values.astype(np.float32), c_lf]; gclf = gbm(M2, pos); report('G+Clf', gclf)
np.savez(f'{O}/nnceil_scores.npz', g36=g36, cnn=c_dn, gc=gc, cnn_lf=c_lf, gclf=gclf, pid=F.precursor_id.values, rta=F.rta.values)
