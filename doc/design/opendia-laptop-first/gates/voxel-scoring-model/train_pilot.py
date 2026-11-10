"""Voxel pilot step 3: train and score the 8 pilot cells (inputs a-d x training sets T0/T1) under PLAN_PILOT.md.
Reads voxels01.npz (build_voxels.py), train_windows.tsv and eval_windows.tsv. Never reads eval_labels.tsv
(evaluate_pilot.py does). Writes scores_<cell>.tsv with one crossfit logit per evaluation window.

Model: TOPAZ body as in OpenDIACNN.cpp (Branch: dual view -> conv 1x1, 3, 3 -> max+mean pool -> proj; Coelution head:
10 statistics + learned pair similarity; MS1 branch + MS1 coelution + MS1-MS2 cross-modal; library fit; MLP head).
Only the input stage changes:
  a  fragment XIC [6,32] (m/z sub-bins summed), MS1 M+0..2
  b  a + ppm-offset map [6,32] (intensity-weighted, relative to the calibrated offset, /5 ppm, 0 without signal) as
     extra fragment-branch channels
  c  a with the 7 m/z sub-bins: per fragment conv over (m/z, view) -> 4 channels, then the RT branch
  d  b + MS1 M+0..3 (M-1 is not exported by -out_xipm; recorded in DESIGN.md)
Training sets (DESIGN.md table): T0 anchors + paired decoys + random decoys; T1 adds decoy-like windows (weight 0.5),
synthetic false targets and synthetic weak positives, built only from training-fold windows. Class weights are
normalised so positives and negatives carry the same total weight. 3 family folds; inside the training folds 15 % of
families (by hash) are held out for early stopping (patience 5, max 40 epochs, Adam 1e-3, batch 64, width 32)."""
import sys, time, numpy as np, pandas as pd, torch, torch.nn as nn, torch.nn.functional as F
B = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107'
O = B + '/astral-voxel-pilot-20261110'
torch.set_num_threads(8)
EPS, BETA, LAG, CORE = 1e-6, 10.0, 2, (12, 20)
V = np.load(O + '/voxels01.npz')
idx = {s: i for i, s in enumerate(V['sample_id'])}
I_all, P_all, M_all, LIB_all = V['I'], V['P'], V['M'], V['LIB']
train = pd.read_csv(O + '/train_windows.tsv', sep='\t'); ev = pd.read_csv(O + '/eval_windows.tsv', sep='\t')
def fam_hash(s):
    h = 0xcbf29ce484222325
    for b in s.encode(): h = ((h ^ b) * 0x100000001b3) & 0xFFFFFFFFFFFFFFFF
    return h
train['valid'] = train.family.map(lambda s: (fam_hash('v' + s) % 100) < 15)

def roll0(x, k):
    """Shift along the last (RT) axis with zero fill."""
    if k == 0: return x
    y = np.zeros_like(x)
    if k > 0: y[..., k:] = x[..., :-k]
    else: y[..., :k] = x[..., -k:]
    return y

def synthetics(tr, rng):
    """T1 synthetic rows from training windows only: false targets (partial / shifted / wrong-ratio co-elution of an
    anchor added onto a decoy or decoy-like window, label 0) and weak positives (scaled anchor added onto a decoy-like
    window, label 1). Returns arrays and labels."""
    anc = tr[tr.group == 'anchor'].sample_id.map(idx).values
    bg = tr[tr.group.isin(['decoy', 'decoy_like'])].sample_id.map(idx).values
    bgl = tr[tr.group == 'decoy_like'].sample_id.map(idx).values
    out = {'I': [], 'P': [], 'M': [], 'L': [], 'y': []}
    for _ in range(3000):  # synthetic false targets
        b, a = rng.choice(bg), rng.choice(anc)
        Ib, Pb, Mb = I_all[b].copy(), P_all[b].copy(), M_all[b].copy()
        k = rng.integers(2, 5)  # 2-4 of 6 fragment channels get foreign signal
        dst = rng.choice(6, k, replace=False); src = rng.choice(6, k, replace=False)
        scale_b = max(I_all[b].sum((1, 2)).max(), 1.0); scale_a = max(I_all[a].sum((1, 2)).max(), 1.0)
        for d_, s_ in zip(dst, src):
            u = np.exp(rng.uniform(np.log(0.3), np.log(3.0))) * scale_b / scale_a  # wrong ratio
            sh = int(rng.integers(-3, 4))  # partial or shifted co-elution
            Ib[d_] += u * roll0(I_all[a][s_], sh); Pb[d_] += u * roll0(P_all[a][s_], sh)
        out['I'].append(Ib); out['P'].append(Pb); out['M'].append(Mb); out['L'].append(LIB_all[b]); out['y'].append(0.0)
    for _ in range(len(anc)):  # synthetic weak positives, one per training anchor
        b, a = rng.choice(bgl), rng.choice(anc); s = rng.uniform(0.05, 0.5)
        out['I'].append(I_all[b] + s * I_all[a]); out['P'].append(P_all[b] + s * P_all[a]); out['M'].append(M_all[b] + s * M_all[a])
        out['L'].append(LIB_all[a]); out['y'].append(1.0)
    return {k: np.asarray(v, np.float32) for k, v in out.items()}

def make_inputs(I, P, M, L, opt):
    """numpy arrays -> dict of tensors for one input option."""
    xic = I.sum(2)
    d = {'xic': torch.from_numpy(xic), 'lib': torch.from_numpy(L)}
    d['ms1'] = torch.from_numpy(M[:, :4] if opt == 'd' else M[:, :3])
    if opt in ('b', 'd'):
        ppm = np.where(xic > 0, P / np.maximum(xic, 1e-12), 0.0) / 5.0
        d['ppm'] = torch.from_numpy(np.clip(ppm, -4, 4).astype(np.float32))
    if opt == 'c': d['sub'] = torch.from_numpy(I)
    return d

def dual(x): return torch.cat([x / (x.amax(2, keepdim=True) + EPS), torch.log1p(x.clamp_min(0))], 1)
def maxn(x): return x / (x.amax(2, keepdim=True) + EPS)
def zs(x):
    c = x - x.mean(2, keepdim=True)
    return c / (torch.sqrt(c.pow(2).sum(2, keepdim=True) / (x.size(2) - 1)) + EPS)
def unit(z): return z / torch.sqrt(z.pow(2).sum(2, keepdim=True) + EPS)
def apex(n):
    w = torch.softmax(n * BETA, 2)
    return (w * torch.linspace(0, 1, n.size(2)).view(1, 1, -1)).sum(2)
def ustd(x): return torch.sqrt((x - x.mean(1, keepdim=True)).pow(2).sum(1) / (x.size(1) - 1) + 1e-12)
def rowcos(a, b): return (a * b).sum(1) / (torch.sqrt(a.pow(2).sum(1) + 1e-12) * torch.sqrt(b.pow(2).sum(1) + 1e-12))

class Branch(nn.Module):
    def __init__(self, cin, w):
        super().__init__()
        self.c0, self.c1, self.c2 = nn.Conv1d(cin, w, 1), nn.Conv1d(w, 2 * w, 3, padding=1), nn.Conv1d(2 * w, 2 * w, 3, padding=1)
        self.proj = nn.Linear(4 * w, 2 * w)
    def forward(self, x):
        h = F.relu(self.c2(F.relu(self.c1(F.relu(self.c0(x))))))
        return self.proj(torch.cat([h.amax(2), h.mean(2)], 1))

class Coelution(nn.Module):
    def __init__(self, c):
        super().__init__()
        iu = torch.triu_indices(c, c, 1); self.register_buffer('upper', iu[0] * c + iu[1]); k = c * (c - 1) // 2
        self.sim = nn.Sequential(nn.Linear(2 * k, 32), nn.ReLU(), nn.Linear(32, 8))
    def pc(self, a, b): return torch.bmm(unit(a), unit(b).transpose(1, 2)).flatten(1).index_select(1, self.upper)
    def forward(self, x):
        L = x.size(2); n = maxn(x); z = zs(n); same = self.pc(z, z); lag = [same]
        for l in range(1, LAG + 1):
            lag += [self.pc(z[:, :, l:], z[:, :, :L - l]), self.pc(z[:, :, :L - l], z[:, :, l:])]
        lb = torch.stack(lag, 2).amax(2); ap = apex(n); w = torch.softmax(n * BETA, 2); ent = -(w * torch.log(w + EPS)).sum(2)
        peak = torch.sigmoid((n.amax(2) / (n.mean(2) + EPS) - 2.0) * 2.0).mean(1)
        s = torch.stack([same.mean(1), same.amin(1), lb.mean(1), lb.amin(1), ustd(ap), ap.amax(1) - ap.amin(1),
                         (ap - ap.mean(1, keepdim=True)).abs().mean(1), ent.mean(1), ent.amin(1), peak], 1)
        return torch.cat([s, self.sim(torch.cat([same, lb], 1))], 1)

def cross(m1, m2):
    n1, n2 = maxn(m1), maxn(m2); c = torch.bmm(unit(zs(n1)), unit(zs(n2)).transpose(1, 2)); a1, a2 = apex(n1), apex(n2)
    return torch.stack([c.mean((1, 2)), c.flatten(1).amax(1), (a1.mean(1) - a2.mean(1)).abs(), ustd(torch.cat([a1, a2], 1))], 1)

class Scorer(nn.Module):
    def __init__(self, opt, w=32):
        super().__init__()
        self.opt = opt; iso = 4 if opt == 'd' else 3
        if opt == 'c':
            self.mz = nn.Conv2d(2, 4, (7, 1)); cin = 6 * 4
        else:
            cin = 12 + (6 if opt in ('b', 'd') else 0)
        self.xic, self.ms1 = Branch(cin, w), Branch(2 * iso, w)
        self.xco, self.mco = Coelution(6), Coelution(iso)
        f = 2 * w + 18 + (2 * w + 18 + 4) + 2 + 6
        self.head = nn.Sequential(nn.Linear(f, 4 * w), nn.ReLU(), nn.Dropout(0.1), nn.Linear(4 * w, 2 * w), nn.ReLU(), nn.Dropout(0.1), nn.Linear(2 * w, 1))
    def forward(self, d):
        xic = d['xic']; n = xic.size(0)
        if self.opt == 'c':
            s = d['sub'].reshape(n * 6, 7, 32)
            v = torch.stack([s / (s.amax((1, 2), keepdim=True) + EPS), torch.log1p(s.clamp_min(0))], 1)  # [n*6,2,7,32]
            xin = F.relu(self.mz(v)).reshape(n, 6 * 4, 32)
        else:
            xin = dual(xic)
            if 'ppm' in d: xin = torch.cat([xin, d['ppm']], 1)
        m = d['ms1']; has = (m.flatten(1).amax(1, keepdim=True) > 0).float()
        ms1 = torch.cat([self.ms1(dual(m)), self.mco(m), cross(m, xic)], 1) * has
        core = xic[:, :, CORE[0]:CORE[1]].sum(2)
        fit = torch.stack([rowcos(torch.sqrt(core), torch.sqrt(d['lib'])), rowcos(core, d['lib'])], 1)
        return self.head(torch.cat([self.xic(xin), self.xco(xic), ms1, fit, d['lib']], 1)).squeeze(1)

def subset(d, i): return {k: v[i] for k, v in d.items()}

def fit_predict(opt, tset, k, seed=1):
    rng = np.random.default_rng(seed * 100 + k); torch.manual_seed(seed * 100 + k)
    tr = train[train.fold != k]
    if tset == 'T0': tr = tr[tr.group.isin(['anchor', 'anchor_decoy', 'decoy'])]
    fit, val = tr[~tr.valid], tr[tr.valid]
    def real(df):
        j = df.sample_id.map(idx).values
        return I_all[j], P_all[j], M_all[j], LIB_all[j], df.label.values.astype(np.float32), df.weight.values.astype(np.float32)
    parts = [real(fit)]
    if tset == 'T1':
        s = synthetics(fit, rng); parts.append((s['I'], s['P'], s['M'], s['L'], s['y'], np.ones(len(s['y']), np.float32)))
    cat = [np.concatenate([p[c] for p in parts]) for c in range(6)]
    X = make_inputs(*cat[:4], opt); y = torch.from_numpy(cat[4]); w = cat[5].copy()
    w[cat[4] == 1] *= w[cat[4] == 0].sum() / w[cat[4] == 1].sum(); w = torch.from_numpy(w / w.mean())
    vr = real(val); Xv = make_inputs(*vr[:4], opt); yv = torch.from_numpy(vr[4]); wv = vr[5].copy()
    wv[vr[4] == 1] *= wv[vr[4] == 0].sum() / max(wv[vr[4] == 1].sum(), 1e-9); wv = torch.from_numpy(wv / wv.mean())
    model = Scorer(opt); optm = torch.optim.Adam(model.parameters(), 1e-3)
    best, stale, state = np.inf, 0, None
    for epoch in range(40):
        model.train(); perm = torch.randperm(len(y))
        for b in range(0, len(y), 64):
            i = perm[b:b + 64]; optm.zero_grad()
            loss = (F.binary_cross_entropy_with_logits(model(subset(X, i)), y[i], reduction='none') * w[i]).mean()
            loss.backward(); optm.step()
        model.eval()
        with torch.no_grad():
            vl = float((F.binary_cross_entropy_with_logits(model(Xv), yv, reduction='none') * wv).mean())
        if vl < best - 1e-4: best, stale, state = vl, 0, {a: t.clone() for a, t in model.state_dict().items()}
        else:
            stale += 1
            if stale >= 5: break
    model.load_state_dict(state); model.eval()
    e = ev[ev.fold == k]; j = e.sample_id.map(idx).values
    Xe = make_inputs(I_all[j], P_all[j], M_all[j], LIB_all[j], opt)
    with torch.no_grad():
        s = torch.cat([model(subset(Xe, torch.arange(b, min(b + 2048, len(j))))) for b in range(0, len(j), 2048)]).numpy()
    return e.assign(score=s), epoch + 1, best, len(y)

if __name__ == '__main__':
    cells = sys.argv[1:] or [f'{o}_{t}' for t in ('T0', 'T1') for o in 'abcd']
    for cell in cells:
        opt, tset = cell.split('_'); t0 = time.time(); outs = []
        for k in range(3):
            o, ep, vl, n = fit_predict(opt, tset, k); outs.append(o)
            print(f'{cell} fold {k}: {n} training rows, {ep} epochs, val loss {vl:.4f}', flush=True)
        pd.concat(outs)[['sample_id', 'precursor_id', 'fold', 'score']].to_csv(f'{O}/scores_{cell}.tsv', sep='\t', index=False)
        print(f'{cell} done in {time.time() - t0:.0f} s', flush=True)
