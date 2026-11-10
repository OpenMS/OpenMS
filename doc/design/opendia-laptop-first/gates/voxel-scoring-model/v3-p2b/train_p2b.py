"""Voxel v3 phase 2b (PLAN_V3_P2b.md): one network, decoy-whitened (deep LDA) objective.
Usage: train_p2b.py <raw01|raw02> <arm> [...]; arms D11, D11n (11 features), DV (voxel + 11 features).
Output: v2/<raw>/scores_<arm>.npy (row order of population.tsv). Never reads entrapment labels."""
import sys, time, numpy as np, pandas as pd, torch, torch.nn as nn
sys.path.insert(0, '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110')
from train_pilot import Branch, Coelution, cross, dual, rowcos, CORE  # v2 voxel body (pilot option b)
RAW, ARMS = sys.argv[1], sys.argv[2:]
O = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2'
torch.set_num_threads(10)
LAM, NPOS, NNEG = 1e-2, 256, 1024
pop = pd.read_csv(f'{O}/{RAW}/population.tsv', sep='\t')
H = np.load(f'{O}/{RAW}/heuristics.npy').astype(np.float64)
X = {k: np.load(f'{O}/{RAW}/{k}.npy', mmap_mode='r') for k in ('xic', 'ppm', 'ms1', 'lib')}
train_rows = pd.read_csv(f'{O}/{RAW}/train_windows.tsv', sep='\t')
row_of = pd.Series(np.arange(len(pop)), index=pop.sample_id.values)
fold = pop.fold.values; is_decoy = pop.is_decoy.values; pair_key = pop.pair_key.values; pid = pop.precursor_id.values

def fnv(s):
    h = 0xcbf29ce484222325
    for b in s.encode(): h = ((h ^ b) * 0x100000001b3) & 0xFFFFFFFFFFFFFFFF
    return h

def pooled_q(rows, s):
    """As train_v2.pooled_q: precursor winners, per-fold decoy-winner median/MAD, best of each pair, q = (D+1)/T."""
    d = pd.DataFrame({'r': rows, 's': s, 'pid': pid[rows], 'f': fold[rows], 'dec': is_decoy[rows], 'key': pair_key[rows]})
    w = d.sort_values(['s', 'r'], ascending=[False, True]).drop_duplicates('pid')
    for f_, g in w[w.dec == 1].groupby('f'):
        med = g.s.median(); mad = 1.4826 * (g.s - med).abs().median()
        w.loc[w.f == f_, 's'] = (w.loc[w.f == f_, 's'] - med) / mad
    best = w.sort_values(['s', 'dec'], ascending=[False, False]).drop_duplicates('key').sort_values('s', ascending=False)
    T = np.cumsum(1 - best.dec.values); D = np.cumsum(best.dec.values)
    q = np.minimum.accumulate(np.minimum(1, (D + 1) / np.maximum(T, 1))[::-1])[::-1]
    return best.r.values, q, best.dec.values, w

class VoxelBody(nn.Module):
    """v2 Scorer without its head: 176 values (fragment branch, fragment coelution, MS1 part, library fit, library)."""
    def __init__(self, w=32):
        super().__init__()
        self.xic, self.ms1 = Branch(18, w), Branch(6, w)
        self.xco, self.mco = Coelution(6), Coelution(3)
        self.dim = 2 * w + 18 + (2 * w + 18 + 4) + 2 + 6
    def forward(self, xic, ppm, ms1, lib):
        has = (ms1.flatten(1).amax(1, keepdim=True) > 0).float()
        ms1p = torch.cat([self.ms1(dual(ms1)), self.mco(ms1), cross(ms1, xic)], 1) * has
        core = xic[:, :, CORE[0]:CORE[1]].sum(2)
        fit = torch.stack([rowcos(torch.sqrt(core), torch.sqrt(lib)), rowcos(core, lib)], 1)
        return torch.cat([self.xic(torch.cat([dual(xic), ppm], 1)), self.xco(xic), ms1p, fit, lib], 1)

class Net(nn.Module):
    """h = [z, tanh(g(u))] (skip) or tanh(g(u)) (no skip); u = z or [voxel body, z]."""
    def __init__(self, voxel, skip, n=11):
        super().__init__()
        self.body = VoxelBody() if voxel else None
        din = (self.body.dim if voxel else 0) + n
        self.g = nn.Sequential(nn.Linear(din, 64), nn.ReLU(), nn.Linear(64, 16))
        self.skip = skip
    def forward(self, z, vox=None):
        u = torch.cat([self.body(*vox), z], 1) if self.body is not None else z
        t = torch.tanh(self.g(u))
        return torch.cat([z, t], 1) if self.skip else t

def load(rows, Z, voxel):
    """Tensors for rows in the given order (memmaps are read in sorted order)."""
    z = torch.from_numpy(Z[rows].astype(np.float32))
    if not voxel: return z, None
    o = np.argsort(rows); sr = rows[o]; inv = np.argsort(o)
    xic = np.asarray(X['xic'][sr], np.float32)[inv]; p = np.asarray(X['ppm'][sr], np.float32)[inv]
    ppm = np.clip(np.where(xic > 0, p / np.maximum(xic, 1e-12), 0) / 5.0, -4, 4).astype(np.float32)
    ms1 = np.asarray(X['ms1'][sr, :3], np.float32)[inv]; lib = np.asarray(X['lib'][sr], np.float32)[inv]
    return z, tuple(torch.from_numpy(a) for a in (xic, ppm, ms1, lib))

def embed(net, rows, Z, voxel, bs=4096):
    net.eval(); out = []
    with torch.no_grad():
        for b in range(0, len(rows), bs):
            z, v = load(rows[b:b + bs], Z, voxel); out.append(net(z, v))
    return torch.cat(out).double()

def readout(hp, hn):
    """w = (S_- + lam I)^-1 (mu_+ - mu_-), and the decoy-whitened objective dmu' w."""
    mp, mn = hp.mean(0), hn.mean(0); c = hn - mn
    S = c.T @ c / max(len(hn) - 1, 1) + LAM * torch.eye(hn.shape[1], dtype=hn.dtype)
    dmu = mp - mn; w = torch.linalg.solve(S, dmu)
    return w, (dmu * w).sum()

def fit(rows, y, Z, voxel, skip, seed, valid, steps):
    torch.manual_seed(seed); rng = np.random.default_rng(seed)
    net = Net(voxel, skip); opt = torch.optim.Adam(net.parameters(), 1e-3)
    tp, tn = rows[(~valid) & (y == 1)], rows[(~valid) & (y == 0)]
    vp, vn = rows[valid & (y == 1)], rows[valid & (y == 0)]
    best, stale, state = -np.inf, 0, None
    for ep in range(30):
        net.train()
        for _ in range(steps):
            bp = rng.choice(tp, NPOS, replace=len(tp) < NPOS); bn = rng.choice(tn, NNEG, replace=len(tn) < NNEG)
            r = np.concatenate([bp, bn]); z, v = load(r, Z, voxel); h = net(z, v).double()
            _, obj = readout(h[:NPOS], h[NPOS:])
            opt.zero_grad(); (-obj).backward(); opt.step()
        hv = embed(net, np.concatenate([vp, vn]), Z, voxel)
        _, vo = readout(hv[:len(vp)], hv[len(vp):]); vo = float(vo)
        if vo > best + 1e-4: best, stale, state = vo, 0, {a: t.clone() for a, t in net.state_dict().items()}
        else:
            stale += 1
            if stale >= 5: break
    net.load_state_dict(state)
    ht = embed(net, np.concatenate([tp, tn]), Z, voxel); w, _ = readout(ht[:len(tp)], ht[len(tp):])
    def score(r): return (embed(net, r, Z, voxel) @ w).numpy().astype(np.float32)
    return score, ep + 1, best

def run(arm):
    voxel, skip = arm == 'DV', arm != 'D11n'; steps = 50 if voxel else 100
    scores = np.zeros(len(pop), np.float32); t0 = time.time()
    for k in range(3):
        tr_all = np.where(fold != k)[0]
        mu, sd = H[tr_all].mean(0), H[tr_all].std(0); sd[sd == 0] = 1; Z = (H - mu) / sd
        t = train_rows[(train_rows.fold != k) & train_rows.group.isin(['anchor', 'anchor_decoy', 'decoy'])]
        rows = row_of[t.sample_id].values; y = t.label.values.astype(np.float32)
        valid = np.array([fnv('v' + str(s)) % 100 < 15 for s in t.family.values])
        sc, ep, vo = fit(rows, y, Z, voxel, skip, 100 + k, valid, steps); log = [f'r0 {len(rows)} {ep} ep obj {vo:.2f}']
        for r in (1, 2):
            win, q, dec, wtab = pooled_q(tr_all, sc(tr_all))
            pos = win[(dec == 0) & (q <= (0.15 if r == 1 else 0.05))]; negs = wtab[wtab.dec == 1].r.values
            rows = np.concatenate([pos, negs]); y = np.r_[np.ones(len(pos)), np.zeros(len(negs))].astype(np.float32)
            valid = np.array([fnv('v' + str(s)) % 100 < 15 for s in pair_key[rows]])
            sc, ep, vo = fit(rows, y, Z, voxel, skip, 200 + 10 * k + r, valid, steps)
            log.append(f'r{r} pos {len(pos)} neg {len(negs)} {ep} ep obj {vo:.2f}')
        te = np.where(fold == k)[0]; scores[te] = sc(te)
        print(f'{RAW} {arm} fold {k}: ' + ' | '.join(log), flush=True)
    np.save(f'{O}/{RAW}/scores_{arm}.npy', scores); print(f'{RAW} {arm} done in {time.time() - t0:.0f} s', flush=True)

for a in ARMS: run(a)
