"""Information-content test on raw 01 (diagnostic; DIA-NN labels used ONLY to bound what is learnable, never in OpenDIA).
Same candidates as regimes2 'picker' (top-3 local maxima within +-60 s of calibrated library RT, OSW XICs).
Inputs per candidate: 6 fragments + 3 MS1 isotopes, 33 scans centred on the apex (+-26 s), per-sample max-normalised
(log1p of raw scale kept as one scalar). Compare, on identical candidates and folds:
  (a) hand-feature MLP oracle (11 features from regimes2)
  (b) trace CNN oracle (learned representation from the traces only)
  (c) CNN + hand features
Positives: DIA-NN-confirmed right peaks; negatives: decoy groups (subsampled). Cross-fitted by pair-hash fold; winner per
precursor; paired TDC; report q<=0.01/0.05 with entrapment FDP_c. If (b) or (c) >> (a), the raw traces carry
information our summaries discard and representation learning is the lever; the remaining problem is then labels."""
import sys; sys.path.insert(0, '.')
from regimes_lib import *
import torch, torch.nn as nn, time
torch.set_num_threads(12); W = 16
res, ap = run(frag, ms1, lib, npts, rt, cal_a := meta.cal_rt.fillna(-1e9).values.astype(np.float32), ok, WIN, 1)
have = np.isfinite(res[:, :, 0]); ii, kk = np.where(have)
X = res[ii, kk].astype(np.float32); apex_rt = ap[ii, kk]
m = meta.iloc[ii].reset_index(drop=True)
prec = m.precursor_id.values; dec = m.decoy.values.astype(bool); fold = m.fold.values; pair = m.pair.values; ent = m.ent.fillna(0).values
right = (~dec) & m.dn_rt.notna().values & (np.abs(apex_rt - m.dn_rt.values) <= 10)

@nb.njit(parallel=True, cache=True)
def windows(frag, ms1, rt, npts, ii, apex_rt, W):
    n_ = len(ii); out = np.zeros((n_, 9, 2 * W + 1), np.float32); sc = np.zeros(n_, np.float32)
    for g in nb.prange(n_):
        i = ii[g]; n = npts[i]
        a = 0
        for t in range(n):
            if rt[i, t] == apex_rt[g]: a = t; break
        for d in range(-W, W + 1):
            t = a + d
            if 0 <= t < n:
                for f in range(6): out[g, f, d + W] = frag[i, f, t]
                for f in range(3): out[g, 6 + f, d + W] = ms1[i, f, t]
        mf = out[g, :6].max(); m1 = out[g, 6:].max()
        sc[g] = np.log1p(mf)
        if mf > 0: out[g, :6] /= mf
        if m1 > 0: out[g, 6:] /= m1
    return out, sc
t0 = time.time()
Win, scale = windows(frag, ms1, rt, npts, ii.astype(np.int64), apex_rt, W)
print('windows', Win.shape, f'{time.time() - t0:.0f} s', flush=True)
Xh = np.c_[X, scale]; Xh = np.nan_to_num((Xh - np.nanmean(Xh, 0)) / (np.nanstd(Xh, 0) + 1e-9)).astype(np.float32)

class Net(nn.Module):
    def __init__(s, use_trace, use_hand):
        super().__init__(); s.t, s.h = use_trace, use_hand
        s.conv = nn.Sequential(nn.Conv1d(9, 32, 5, padding=2), nn.GELU(), nn.Conv1d(32, 32, 5, padding=2), nn.GELU(),
                               nn.Conv1d(32, 32, 3, padding=1), nn.GELU())
        s.pair = nn.Sequential(nn.Conv2d(1, 8, (1, 5), padding=(0, 2)), nn.GELU())
        d = (64 + 36 if use_trace else 0) + (Xh.shape[1] if use_hand else 0)
        s.head = nn.Sequential(nn.Linear(d, 64), nn.GELU(), nn.Linear(64, 16), nn.GELU(), nn.Linear(16, 1))
    def forward(s, w, h):
        z = []
        if s.t:
            c = s.conv(w); z += [c.mean(2), c.amax(2)]
            f = w[:, :6]; f = f - f.mean(2, keepdim=True); f = f / (f.norm(dim=2, keepdim=True) + 1e-6)
            z.append((f @ f.transpose(1, 2)).flatten(1))  # 6x6 fragment correlation matrix (co-elution structure)
        if s.h: z.append(h)
        return s.head(torch.cat(z, 1)).squeeze(1)

def fit_score(use_trace, use_hand, seed=0, epochs=4):
    torch.manual_seed(seed); rng = np.random.default_rng(seed); out = np.zeros(len(Xh), np.float32)
    for k in range(3):
        pos = np.where(right & (fold != k))[0]; neg = np.where(dec & (fold != k))[0]
        neg = rng.choice(neg, min(len(neg), 300000), replace=False)
        idx = np.r_[pos, neg]; y = np.r_[np.ones(len(pos)), np.zeros(len(neg))].astype(np.float32)
        net = Net(use_trace, use_hand); opt = torch.optim.AdamW(net.parameters(), lr=2e-3, weight_decay=1e-3)
        lf = nn.BCEWithLogitsLoss(pos_weight=torch.tensor(20.0))
        for ep in range(epochs):
            perm = rng.permutation(len(idx))
            for b in range(0, len(idx), 1024):
                j = idx[perm[b:b + 1024]]
                loss = lf(net(torch.from_numpy(Win[j]), torch.from_numpy(Xh[j])), torch.from_numpy(y[perm[b:b + 1024]]))
                opt.zero_grad(); loss.backward(); opt.step()
        net.eval(); te = np.where(fold == k)[0]
        with torch.no_grad():
            for b in range(0, len(te), 20000):
                j = te[b:b + 20000]; out[j] = net(torch.from_numpy(Win[j]), torch.from_numpy(Xh[j])).numpy()
    return out

def tdc(score):
    allw = winners(np.arange(len(score)), prec, score)
    w = pd.DataFrame({'pair': pair[allw], 's': score[allw], 'dec': dec[allw], 'right': right[allw], 'ent': ent[allw]})
    w = w.sort_values(['s', 'dec'], ascending=[False, False]).drop_duplicates('pair').sort_values('s', ascending=False)
    d = w.dec.values; Tn = np.cumsum(~d); Dn = np.cumsum(d)
    q = np.minimum.accumulate(np.minimum(1, (Dn + 1) / np.maximum(Tn, 1))[::-1])[::-1]
    r = []
    for qq in (0.01, 0.05):
        a = w[(~d) & (q <= qq)]
        r.append(f'{len(a)} (conf-right {int(a.right.sum())}, FDP_c {a.ent.sum() * (1 + 1 / R_ENT) / max(1, len(a)):.3f})')
    return ' | '.join(r)
R_ENT = 148983 / 151239
print('confirmed-right candidates', int(right.sum()), '| decoy groups', int(dec.sum()), flush=True)
for name, t_, h_ in [('(a) hand MLP', False, True), ('(b) trace CNN', True, False), ('(c) CNN + hand', True, True)]:
    for seed in (0, 1):
        t0 = time.time(); sc = fit_score(t_, h_, seed)
        print(f'{name:16s} seed {seed}: q<=0.01 / 0.05 = {tdc(sc)}  [{time.time() - t0:.0f} s]', flush=True)
