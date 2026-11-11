"""Voxel v4 phase 1 (PLAN_V4.md): D11 (PLAN_V3_P2b decoy-whitened single network on 11 features) with a within-precursor
ranking term (R) and/or 6 MS1 isotope features (M). Usage: train_v4.py <raw01|raw02> <arm> [...]; arms D11R, D11M, D11MR.
Output: v2/<raw>/scores_<arm>.npy (row order of population.tsv). Never reads entrapment labels."""
import os, sys, time, inspect, numpy as np, pandas as pd, torch, torch.nn as nn
sys.path.insert(0, '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110')
from train_pilot import Branch, Coelution, cross, dual, rowcos, CORE  # v2 voxel body (pilot option b)
RAW, ARMS = sys.argv[1], sys.argv[2:]
O = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2'
# TRAIN_THREADS (default 10) lowers the thread count when the CPU is shared (PLAN_V7.notes).
torch.set_num_threads(int(os.environ.get('TRAIN_THREADS', '10')))
LAM, NPOS, NNEG = 1e-2, 256, 1024
pop = pd.read_csv(f'{O}/{RAW}/population.tsv', sep='\t')
H = np.load(f'{O}/{RAW}/heuristics.npy').astype(np.float64)
X = {k: np.load(f'{O}/{RAW}/{k}.npy', mmap_mode='r') for k in ('xic', 'ppm', 'ms1', 'lib')}
if any('I' in a[3:] for a in ARMS):
    X['xic_corrected'] = np.load(f'{O}/{RAW}/xic_corrected.npy', mmap_mode='r')
train_rows = pd.read_csv(f'{O}/{RAW}/train_windows.tsv', sep='\t')
row_of = pd.Series(np.arange(len(pop)), index=pop.sample_id.values)
fold = pop.fold.values; is_decoy = pop.is_decoy.values; pair_key = pop.pair_key.values; pid = pop.precursor_id.values
M1 = np.load(f'{O}/{RAW}/ms1feat.npy').astype(np.float64)
F12 = np.load(f'{O}/{RAW}/feat12.npy').astype(np.float64) if any('F' in a for a in ARMS) else None
FPC = np.load(f'{O}/{RAW}/featPC.npy').astype(np.float64) if any(('P' in a[3:]) or ('C' in a[3:]) for a in ARMS) else None
# PLAN_V7: cross-run features of the other run's matching window (xrun_features.py, label-free, model-free).
FX = np.load(f'{O}/{RAW}/featX.npy').astype(np.float64) if any('X' in a[3:] for a in ARMS) else None
# PLAN_V9: native-sampling fragment evidence with a spectrum-density chance model (native_features.py, label-free).
FN = np.load(f'{O}/{RAW}/nat.npy').astype(np.float64) if any('N' in a[3:] for a in ARMS) else None
# PLAN_V10: the same evidence with an empirical per-m/z local null (native_features_emp.py, label-free).
FE = np.load(f'{O}/{RAW}/natE.npy').astype(np.float64) if any('E' in a[3:] for a in ARMS) else None
# PLAN_V11: the same evidence on fragments unique within the isolation window and RT tolerance (native_features_uniq.py).
FU = np.load(f'{O}/{RAW}/natU.npy').astype(np.float64) if any('U' in a[3:] for a in ARMS) else None
PC_COLS = {'P': [0, 1, 2], 'C': [3, 4]}
ALPHA, MAXSIB, SEP = 1.0, 8, 10.0
# Siblings: the windows of the same precursor with apex > SEP s from the positive's apex.
rt = pop.rt_apex.values
_order = np.argsort(pid, kind='stable'); _p = pid[_order]; _start = np.r_[0, np.flatnonzero(np.diff(_p)) + 1]; _end = np.r_[_start[1:], len(_p)]
GROUP = {int(_p[a]): _order[a:b] for a, b in zip(_start, _end)}
def siblings(r, rng):
    g = GROUP[int(pid[r])]; g = g[np.abs(rt[g] - rt[r]) > SEP]
    return g if len(g) <= MAXSIB else rng.choice(g, MAXSIB, replace=False)

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
    """Dual-view body: raw and corrected six-fragment profiles, MS1, mass error and library intensities."""
    def __init__(self, w=32):
        super().__init__()
        self.xic, self.ms1 = Branch(30, w), Branch(6, w)
        self.xco, self.xco_corr, self.mco = Coelution(6), Coelution(6), Coelution(3)
        self.dim = 2 * w + 18 + 18 + (2 * w + 18 + 8) + 4 + 6
    def forward(self, xic, xic_corr, ppm, ms1, lib):
        has = (ms1.flatten(1).amax(1, keepdim=True) > 0).float()
        ms1p = torch.cat([self.ms1(dual(ms1)), self.mco(ms1), cross(ms1, xic), cross(ms1, xic_corr)], 1) * has
        core = xic[:, :, CORE[0]:CORE[1]].sum(2)
        corrected_core = xic_corr[:, :, CORE[0]:CORE[1]].sum(2)
        fit = torch.stack([rowcos(torch.sqrt(core), torch.sqrt(lib)), rowcos(core, lib),
                           rowcos(torch.sqrt(corrected_core), torch.sqrt(lib)), rowcos(corrected_core, lib)], 1)
        return torch.cat([self.xic(torch.cat([dual(xic), dual(xic_corr), ppm], 1)),
                          self.xco(xic), self.xco_corr(xic_corr), ms1p, fit, lib], 1)

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
    xic_corr = np.asarray(X['xic_corrected'][sr], np.float32)[inv] if voxel else None
    ppm = np.clip(np.where(xic > 0, p / np.maximum(xic, 1e-12), 0) / 5.0, -4, 4).astype(np.float32)
    ms1 = np.asarray(X['ms1'][sr, :3], np.float32)[inv]; lib = np.asarray(X['lib'][sr], np.float32)[inv]
    if voxel:
        return z, tuple(torch.from_numpy(a) for a in (xic, xic_corr, ppm, ms1, lib))
    return z, None

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

def rank_loss(net, w_scale, bp, z_all_fn, rng):
    """Mean over positives with siblings of logsumexp(s_pos, s_sib) - s_pos, with s = h.w / sqrt(obj) (decoy-SD units)."""
    w, scale = w_scale; idx, sib = [], []
    for i, r in enumerate(bp):
        s_ = siblings(r, rng)
        if len(s_): idx.append(i); sib.append(s_)
    if not idx: return None, 0
    flat = np.concatenate(sib); seg = np.repeat(np.arange(len(idx)), [len(s) for s in sib])
    hs = net(z_all_fn(flat)).double(); hp = net(z_all_fn(bp[idx])).double()
    sp = hp @ w / scale; ss = hs @ w / scale
    m = torch.full((len(idx), MAXSIB + 1), -torch.inf, dtype=sp.dtype); m[:, 0] = sp
    pos_in = np.concatenate([np.arange(len(s)) for s in sib]) + 1
    m = m.index_put((torch.from_numpy(seg), torch.from_numpy(pos_in)), ss)
    return (torch.logsumexp(m, 1) - sp).mean(), len(idx)

def fit(rows, y, Z, voxel, skip, seed, valid, steps, alpha=0.0):
    torch.manual_seed(seed); rng = np.random.default_rng(seed)
    net = Net(voxel, skip, Z.shape[1]); opt = torch.optim.Adam(net.parameters(), 1e-3)
    zf = lambda r: torch.from_numpy(Z[r].astype(np.float32))
    tp, tn = rows[(~valid) & (y == 1)], rows[(~valid) & (y == 0)]
    vp, vn = rows[valid & (y == 1)], rows[valid & (y == 0)]
    best, stale, state = -np.inf, 0, None
    for ep in range(30):
        net.train()
        for _ in range(steps):
            bp = rng.choice(tp, NPOS, replace=len(tp) < NPOS); bn = rng.choice(tn, NNEG, replace=len(tn) < NNEG)
            r = np.concatenate([bp, bn]); z, v = load(r, Z, voxel); h = net(z, v).double()
            w, obj = readout(h[:NPOS], h[NPOS:]); loss = -obj
            if alpha > 0:
                lr_, _ = rank_loss(net, (w, torch.sqrt(obj.clamp_min(1e-9))), bp, zf, rng)
                if lr_ is not None: loss = loss + alpha * obj.detach() * lr_
            opt.zero_grad(); loss.backward(); opt.step()
        hv = embed(net, np.concatenate([vp, vn]), Z, voxel)
        wv, vo = readout(hv[:len(vp)], hv[len(vp):]); vo = float(vo)
        crit = vo
        if alpha > 0:
            net.eval()
            with torch.no_grad():
                lv, _ = rank_loss(net, (wv, np.sqrt(max(vo, 1e-9))), vp, zf, np.random.default_rng(seed + 7))
            crit = np.log(max(vo, 1e-9)) - alpha * (float(lv) if lv is not None else 0.0)
        if crit > best + 1e-4: best, stale, state = crit, 0, {a: t.clone() for a, t in net.state_dict().items()}
        else:
            stale += 1
            if stale >= 5: break
    net.load_state_dict(state)
    ht = embed(net, np.concatenate([tp, tn]), Z, voxel); w, _ = readout(ht[:len(tp)], ht[len(tp):])
    def score(r): return (embed(net, r, Z, voxel) @ w).numpy().astype(np.float32)
    return score, ep + 1, best

def run(arm):
    voxel, skip = 'I' in arm[3:], True; steps = 100; alpha = ALPHA if arm.endswith('R') else 0.0
    Hx = np.c_[H, M1] if 'M' in arm[3:] else H
    if 'F' in arm[3:]: Hx = np.c_[Hx, F12]
    for key, cols in PC_COLS.items():
        if key in arm[3:]: Hx = np.c_[Hx, FPC[:, cols]]
    if 'X' in arm[3:]: Hx = np.c_[Hx, FX]
    if 'N' in arm[3:]: Hx = np.c_[Hx, FN]
    if 'E' in arm[3:]: Hx = np.c_[Hx, FE]
    if 'U' in arm[3:]: Hx = np.c_[Hx, FU]
    lin = np.zeros(len(pop), np.float32)
    scores = np.zeros(len(pop), np.float32); t0 = time.time()
    for k in range(3):
        tr_all = np.where(fold != k)[0]
        mu, sd = Hx[tr_all].mean(0), Hx[tr_all].std(0); sd[sd == 0] = 1; Z = (Hx - mu) / sd
        t = train_rows[(train_rows.fold != k) & train_rows.group.isin(['anchor', 'anchor_decoy', 'decoy'])]
        rows = row_of[t.sample_id].values; y = t.label.values.astype(np.float32)
        valid = np.array([fnv('v' + str(s)) % 100 < 15 for s in t.family.values])
        sc, ep, vo = fit(rows, y, Z, voxel, skip, 100 + k, valid, steps, alpha); log = [f'r0 {len(rows)} {ep} ep crit {vo:.3f}']
        for r in (1, 2):
            win, q, dec, wtab = pooled_q(tr_all, sc(tr_all))
            pos = win[(dec == 0) & (q <= (0.15 if r == 1 else 0.05))]; negs = wtab[wtab.dec == 1].r.values
            rows = np.concatenate([pos, negs]); y = np.r_[np.ones(len(pos)), np.zeros(len(negs))].astype(np.float32)
            valid = np.array([fnv('v' + str(s)) % 100 < 15 for s in pair_key[rows]])
            sc, ep, vo = fit(rows, y, Z, voxel, skip, 200 + 10 * k + r, valid, steps, alpha)
            log.append(f'r{r} pos {len(pos)} neg {len(negs)} {ep} ep crit {vo:.3f}')
        te = np.where(fold == k)[0]; scores[te] = sc(te)
        # Linear part z.w_z of the final readout (PLAN_V4B W rule; the full score is unchanged).
        cv = inspect.getclosurevars(sc).nonlocals; d = Z.shape[1]
        lin[te] = (embed(cv['net'], te, Z, voxel)[:, :d] @ cv['w'][:d]).numpy()
        print(f'{RAW} {arm} fold {k}: ' + ' | '.join(log) + f' ({time.time() - t0:.0f} s)', flush=True)
    np.save(f'{O}/{RAW}/scores_{arm}.npy', scores); np.save(f'{O}/{RAW}/lin_{arm}.npy', lin)
    print(f'{RAW} {arm} done in {time.time() - t0:.0f} s', flush=True)

if __name__ == '__main__':
    for a in ARMS: run(a)
