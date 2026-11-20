"""PRECHECK_IMFEAT step 3: do mobility-aware hand features or a 2-D (RT x 1/K0) trace CNN add information over the 36
default features for faint peaks? Sample faint2 (6,000 DIA-NN targets + 30,000 decoys, r~1.56 D s0, default search
with -im_tolerance_scale 0.6, i.e. ev/ims06 features). DIA-NN labels for evaluation only; GBM arms are ceilings.
Arms (GBM, DIA-NN-labelled, cross-fitted on workflow folds, 2 seeds):
  36 | 36+faint2 (p_cos_apex, p_local_sn) | 36+band6 (b_*) | 36+faint2+band6 | 36+CNN2d_lf | 36+all
Label-free CNNs: positives = tool q<=0.01 target winners, negatives = decoy windows (as the C++ scorer seeds).
  CNN1d: 2-D patch summed over mobility; CNN2d: full 6 x 33 x 8 patch.
Measure: TPR of the DIA-NN targets' right window vs decoy precursors' best windows at decoy FPR 1e-3 / 3e-3 / 1e-2.
Proceed rule (fixed before running): a component goes to step 2 if 36+component gains >= +0.020 TPR at FPR 3e-3 over
36 in both GBM seeds."""
import os, glob, numpy as np, pandas as pd, torch, torch.nn as nn
from sklearn.ensemble import HistGradientBoostingClassifier
NT = int(os.environ.get('NTHREADS', '4')); torch.set_num_threads(NT)
G = '/home/singjc/.jcode/scratch/dia2/gap7'; O = os.environ.get('FAINT2', f'{G}/faint2')
FEAT = os.environ.get('FEATS', glob.glob(f'{G}/ev/ims06/*.features.tsv')[0] if glob.glob(f'{G}/ev/ims06/*.features.tsv') else '')
S = pd.read_csv(f'{O}/sample.tsv', sep='\t').set_index('pid')
D = pd.read_parquet(f'{O}/imfeat.parquet')
Xp = np.load(f'{O}/patches2d.npy'); K = pd.read_parquet(f'{O}/patch_keys.parquet'); K['k'] = np.arange(len(K))
ids = set(S.index); parts = []
for ch in pd.read_csv(FEAT, sep='\t', chunksize=1_000_000):
    parts.append(ch[ch.precursor_id.isin(ids)])
F = pd.concat(parts); F['rta'] = F.rt_apex.round(3)
n0 = len(F)
F = F.merge(D.drop_duplicates(['precursor_id', 'rta']), on=['precursor_id', 'rta']).merge(K.drop_duplicates(['precursor_id', 'rta']), on=['precursor_id', 'rta'])
print(f'workflow windows for sample {n0}, joined {len(F)}', flush=True)
FAINT = ['p_cos_apex', 'p_local_sn']; BAND = ['b_offset', 'b_frag_sd', 'b_cos_apex', 'b_support', 'b_contrast', 'b_local_sn']
meta = ['candidate_row', 'precursor_id', 'is_decoy', 'fold', 'rt_apex', 'score', 'q_value', 'rta', 'k']
X36 = [c for c in F.columns if c not in meta + FAINT + BAND]
dec = F.is_decoy.values == 1
pos = (F.precursor_id.map(S.cls).values == 'dnn') & (np.abs(F.rt_apex.values - F.precursor_id.map(S.rt_s).values) <= 6)
fold = F.fold.values; qq = S.qq; dec_ids = S.index[S.cls == 'decoy']
F['_s'] = np.where(F.score > -1e29, F.score, -1e30)
win = F.sort_values('_s', ascending=False).drop_duplicates('precursor_id')
pos_lf = F.index.isin(set(win.index[(win.is_decoy == 0) & (win.q_value <= 0.01)]))
print(f'features {len(X36)}, DIA-NN right windows {pos.sum()} ({F[pos].precursor_id.nunique()} precursors), label-free positives {pos_lf.sum()}', flush=True)
P2 = torch.from_numpy(Xp[F.k.values])                    # n x 6 x 33 x 8
P1 = P2.sum(3) / 8.0                                     # n x 6 x 33

class Net1(nn.Module):
    def __init__(s, h=32):
        super().__init__(); s.c = nn.Sequential(nn.Conv1d(6, h, 5, padding=2), nn.ReLU(), nn.Conv1d(h, h, 5, padding=2), nn.ReLU())
        s.o = nn.Sequential(nn.Linear(3 * h, 32), nn.ReLU(), nn.Linear(32, 1))
    def forward(s, x):
        z = s.c(x); return s.o(torch.cat([z.max(2).values, z.mean(2), z[:, :, 14:19].mean(2)], 1)).squeeze(1)
class Net2(nn.Module):
    def __init__(s, h=32):
        super().__init__(); s.c = nn.Sequential(nn.Conv2d(6, h, (5, 3), padding=(2, 1)), nn.ReLU(), nn.Conv2d(h, h, (5, 3), padding=(2, 1)), nn.ReLU())
        s.o = nn.Sequential(nn.Linear(3 * h, 32), nn.ReLU(), nn.Linear(32, 1))
    def forward(s, x):
        z = s.c(x); return s.o(torch.cat([z.amax((2, 3)), z.mean((2, 3)), z[:, :, 14:19, 2:6].mean((2, 3))], 1)).squeeze(1)

def cnn(Net, X, posm, seed=100, epochs=8, steps=300):
    out = np.zeros(len(F), np.float32)
    for f in np.unique(fold):
        torch.manual_seed(seed + int(f)); g = np.random.default_rng(seed + int(f))
        tr = fold != f; tp = np.where(tr & posm)[0]; tn = np.where(tr & dec)[0]
        vp, vn = tp[g.random(len(tp)) < 0.1], tn[g.random(len(tn)) < 0.1]; tp, tn = np.setdiff1d(tp, vp), np.setdiff1d(tn, vn)
        net = Net(); opt = torch.optim.AdamW(net.parameters(), 1e-3, weight_decay=1e-3); L = nn.BCEWithLogitsLoss(); best, st = 1e9, None
        for ep in range(epochs):
            net.train()
            for _ in range(steps):
                b = np.r_[g.choice(tp, 128), g.choice(tn, 128)]; y = torch.tensor(np.r_[np.ones(128), np.zeros(128)], dtype=torch.float32)
                opt.zero_grad(); l = L(net(X[b]), y); l.backward(); opt.step()
            net.eval()
            with torch.no_grad():
                vb = np.r_[vp, g.choice(vn, min(len(vn), 20000), replace=False)]
                vl = float(L(net(X[vb]), torch.tensor(np.r_[np.ones(len(vp)), np.zeros(len(vb) - len(vp))], dtype=torch.float32)))
            if vl < best: best, st = vl, {k: v.clone() for k, v in net.state_dict().items()}
        net.load_state_dict(st); net.eval(); te = np.where(fold == f)[0]
        with torch.no_grad():
            for s in range(0, len(te), 8192): out[te[s:s + 8192]] = net(X[te[s:s + 8192]]).numpy()
        print(f'  {Net.__name__} fold {f}: val loss {best:.4f}', flush=True)
    return out

def gbm(M, seed):
    s = np.zeros(len(F))
    for f in np.unique(fold):
        tr = (fold != f) & (pos | dec); y = pos[tr].astype(int); w = np.where(y == 1, (y == 0).sum() / max(y.sum(), 1), 1.0)
        s[fold == f] = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=31, random_state=seed).fit(M[tr], y, sample_weight=w).decision_function(M[fold == f])
    return s
def tpr(s):
    F['_r'] = s; best = F.groupby('precursor_id')._r.max(); rb = F[pos].groupby('precursor_id')._r.max(); neg = best.reindex(dec_ids).dropna().values
    return [(rb > np.quantile(neg, 1 - f)).mean() for f in (1e-3, 3e-3, 1e-2)], (rb[rb.index.map(qq) == 'Q1'] > np.quantile(neg, 1 - 3e-3)).mean()

c1 = cnn(Net1, P1, pos_lf); c2 = cnn(Net2, P2, pos_lf)
F['cnn1d'] = c1; F['cnn2d'] = c2
a, q = tpr(F._s.values); print(f'{"tool score":22s} TPR 1e-3/3e-3/1e-2 {a[0]:.3f}/{a[1]:.3f}/{a[2]:.3f}  Q1@3e-3 {q:.3f}', flush=True)
for name, s in (('CNN1d alone (lf)', c1), ('CNN2d alone (lf)', c2)):
    a, q = tpr(s); print(f'{name:22s} TPR 1e-3/3e-3/1e-2 {a[0]:.3f}/{a[1]:.3f}/{a[2]:.3f}  Q1@3e-3 {q:.3f}', flush=True)
arms = [('36', X36), ('36+faint2', X36 + FAINT), ('36+band6', X36 + BAND), ('36+faint2+band6', X36 + FAINT + BAND),
        ('36+CNN1d', X36 + ['cnn1d']), ('36+CNN2d', X36 + ['cnn2d']), ('36+all', X36 + FAINT + BAND + ['cnn2d'])]
res = {}
for name, cols in arms:
    M = F[cols].values.astype(np.float32); r = [tpr(gbm(M, sd)) for sd in (0, 1)]; res[name] = r
    print(f'{name:22s} ' + ' | '.join(f'seed{sd} {a[0]:.3f}/{a[1]:.3f}/{a[2]:.3f} Q1 {q:.3f}' for sd, (a, q) in enumerate(r)), flush=True)
print('\nProceed rule (>= +0.020 TPR at 3e-3 over 36 in both seeds):')
for name in res:
    if name == '36': continue
    d = [res[name][sd][0][1] - res['36'][sd][0][1] for sd in (0, 1)]
    print(f'  {name:20s} delta {d[0]:+.3f} / {d[1]:+.3f} -> {"PROCEED" if min(d) >= 0.020 else "no"}')
F[['precursor_id', 'rta', 'cnn1d', 'cnn2d']].to_parquet(f'{O}/cnn_scores.parquet')
