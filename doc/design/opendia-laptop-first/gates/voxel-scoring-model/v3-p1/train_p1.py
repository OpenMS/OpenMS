"""Voxel v3 phase 1 (PLAN_V3_P1.md): heuristic-only controls with the v2 schedule.
Usage: train_p1.py <raw01|raw02> <arm> [...]; arms L11, N11, N11a. Output: <raw>/scores_<arm>.npy (row order of population.tsv)."""
import sys, time, numpy as np, pandas as pd, torch, torch.nn as nn, torch.nn.functional as F
RAW, ARMS = sys.argv[1], sys.argv[2:]
O = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2'
torch.set_num_threads(10)
pop = pd.read_csv(f'{O}/{RAW}/population.tsv', sep='\t')
H = np.load(f'{O}/{RAW}/heuristics.npy').astype(np.float64)
train_rows = pd.read_csv(f'{O}/{RAW}/train_windows.tsv', sep='\t')
row_of = pd.Series(np.arange(len(pop)), index=pop.sample_id.values)
fold = pop.fold.values; is_decoy = pop.is_decoy.values; pair_key = pop.pair_key.values; pid = pop.precursor_id.values

def fnv(s):
    h = 0xcbf29ce484222325
    for b in s.encode(): h = ((h ^ b) * 0x100000001b3) & 0xFFFFFFFFFFFFFFFF
    return h

def pooled_q(rows, s):
    """As train_v2.pooled_q."""
    d = pd.DataFrame({'r': rows, 's': s, 'pid': pid[rows], 'f': fold[rows], 'dec': is_decoy[rows], 'key': pair_key[rows]})
    w = d.sort_values(['s', 'r'], ascending=[False, True]).drop_duplicates('pid')
    for f_, g in w[w.dec == 1].groupby('f'):
        med = g.s.median(); mad = 1.4826 * (g.s - med).abs().median()
        w.loc[w.f == f_, 's'] = (w.loc[w.f == f_, 's'] - med) / mad
    best = w.sort_values(['s', 'dec'], ascending=[False, False]).drop_duplicates('key').sort_values('s', ascending=False)
    T = np.cumsum(1 - best.dec.values); D = np.cumsum(best.dec.values)
    q = np.minimum.accumulate(np.minimum(1, (D + 1) / np.maximum(T, 1))[::-1])[::-1]
    return best.r.values, q, best.dec.values, w

class Head(nn.Module):
    def __init__(self, n, w=32):
        super().__init__()
        self.net = nn.Linear(n, 1) if ARM_LINEAR else nn.Sequential(nn.Linear(n, 4 * w), nn.ReLU(), nn.Dropout(0.1), nn.Linear(4 * w, 2 * w), nn.ReLU(), nn.Dropout(0.1), nn.Linear(2 * w, 1))
    def forward(self, x): return self.net(x).squeeze(1)

def fit_lda(Z, rows, y, w, cov='pooled'):
    """cov: 'pooled' (base LDA, dominated by the decoys), 'balanced' (S+ + S-, the Fisher-loss form), 'decoy' (S- only)."""
    pos, neg = rows[y == 1], rows[y == 0]
    mp, mn = Z[pos].mean(0), Z[neg].mean(0)
    Sp, Sn = (Z[pos] - mp).T @ (Z[pos] - mp), (Z[neg] - mn).T @ (Z[neg] - mn)
    if cov == 'pooled': S = (Sp + Sn) / (len(pos) + len(neg) - 2)
    elif cov == 'balanced': S = Sp / len(pos) + Sn / len(neg)
    else: S = Sn / len(neg)
    S = S + 1e-3 * np.eye(Z.shape[1])
    wv = np.linalg.solve(S, mp - mn)
    return lambda r: Z[r] @ wv

def fit_net(Z, rows, y, w, seed, valid):
    """As train_v2.fit: class-balanced weights, Adam 1e-3, batch 64, early stopping (patience 5, max 40)."""
    torch.manual_seed(seed); rng = np.random.default_rng(seed)
    w = w.copy(); w[y == 1] *= w[y == 0].sum() / max(w[y == 1].sum(), 1e-9); w /= w.mean()
    Zt = torch.from_numpy(Z.astype(np.float32))
    tr, va = np.where(~valid)[0], np.where(valid)[0]
    yv, wv = torch.from_numpy(y.astype(np.float32)), torch.from_numpy(w.astype(np.float32))
    m = Head(Z.shape[1]); opt = torch.optim.Adam(m.parameters(), 1e-3); best, stale, state = np.inf, 0, None
    for ep in range(40):
        m.train(); perm = rng.permutation(tr)
        for b in range(0, len(perm), 64):
            i = perm[b:b + 64]
            loss = (F.binary_cross_entropy_with_logits(m(Zt[rows[i]]), yv[i], reduction='none') * wv[i]).mean()
            opt.zero_grad(); loss.backward(); opt.step()
        m.eval()
        with torch.no_grad(): vl = float((F.binary_cross_entropy_with_logits(m(Zt[rows[va]]), yv[va], reduction='none') * wv[va]).mean())
        if vl < best - 1e-4: best, stale, state = vl, 0, {a: t.clone() for a, t in m.state_dict().items()}
        else:
            stale += 1
            if stale >= 5: break
    m.load_state_dict(state); m.eval()
    def score(r):
        out = np.empty(len(r), np.float32)
        with torch.no_grad():
            for b in range(0, len(r), 65536): out[b:b + 65536] = m(Zt[r[b:b + 65536]]).numpy()
        return out
    return score

ARM_LINEAR = False

def fit_fisher(Z, rows, y, seed, valid, linear):
    """Fisher loss -(mu+ - mu-)^2 / (var+ + var-) on balanced batches of 1024 + 1024; early stopping on the validation loss."""
    torch.manual_seed(seed); rng = np.random.default_rng(seed)
    global ARM_LINEAR
    ARM_LINEAR = linear; m = Head(Z.shape[1]); opt = torch.optim.Adam(m.parameters(), 1e-3)
    Zt = torch.from_numpy(Z.astype(np.float32))
    tp, tn = rows[(~valid) & (y == 1)], rows[(~valid) & (y == 0)]
    vp, vn = rows[valid & (y == 1)], rows[valid & (y == 0)]
    # Signed Fisher criterion (deviation D1 in PLAN_V3_P2a.deviations: the squared form in the plan is sign-free, so the score can come out inverted).
    def floss(sp, sn): return -(sp.mean() - sn.mean()) / torch.sqrt(sp.var() + sn.var() + 1e-6)
    best, stale, state = np.inf, 0, None
    for ep in range(40):
        m.train()
        for _ in range(300):
            bp = rng.choice(tp, 1024, replace=len(tp) < 1024); bn = rng.choice(tn, 1024, replace=len(tn) < 1024)
            loss = floss(m(Zt[bp]), m(Zt[bn])); opt.zero_grad(); loss.backward(); opt.step()
        m.eval()
        with torch.no_grad(): vl = float(floss(m(Zt[vp]), m(Zt[vn])))
        if vl < best - 1e-4: best, stale, state = vl, 0, {a: t.clone() for a, t in m.state_dict().items()}
        else:
            stale += 1
            if stale >= 5: break
    m.load_state_dict(state); m.eval()
    def score(r):
        out = np.empty(len(r), np.float32)
        with torch.no_grad():
            for b in range(0, len(r), 65536): out[b:b + 65536] = m(Zt[r[b:b + 65536]]).numpy()
        return out
    return score

def run(arm):
    global ARM_LINEAR
    ARM_LINEAR = arm in ('G11', 'G11a')  # phase 1b / 2a: linear head
    scores = np.zeros(len(pop), np.float32); t0 = time.time()
    for k in range(3):
        tr_all = np.where(fold != k)[0]
        mu, sd = H[tr_all].mean(0), H[tr_all].std(0); sd[sd == 0] = 1; Z = (H - mu) / sd
        t = train_rows[(train_rows.fold != k) & train_rows.group.isin(['anchor', 'anchor_decoy', 'decoy'])]
        rows = row_of[t.sample_id].values; y = t.label.values.astype(np.float32); w = t.weight.values.astype(np.float64)
        valid = np.array([fnv('v' + str(s)) % 100 < 15 for s in t.family.values])
        if arm in ('L11', 'L11b', 'L11d'): fit = lambda r, y_, w_, sd_, v_: fit_lda(Z, r, y_, w_, {'L11': 'pooled', 'L11b': 'balanced', 'L11d': 'decoy'}[arm])
        elif arm in ('FL11', 'FM11'): fit = lambda r, y_, w_, sd_, v_: fit_fisher(Z, r, y_, sd_, v_, arm == 'FL11')
        else: fit = lambda r, y_, w_, sd_, v_: fit_net(Z, r, y_, w_, sd_, v_)
        sc = fit(rows, y, w, 100 + k, valid); log = [f'r0 {len(rows)}']
        for r in (1, 2):
            win, q, dec, wtab = pooled_q(tr_all, sc(tr_all))
            pos = win[(dec == 0) & (q <= (0.15 if r == 1 else 0.05))]; negs = wtab[wtab.dec == 1].r.values
            if arm in ('N11', 'G11'):
                negs = np.random.default_rng(1000 + 10 * k + r).choice(negs, min(len(negs), max(6000, 4 * len(pos))), replace=False)
            rows = np.concatenate([pos, negs]); y = np.r_[np.ones(len(pos)), np.zeros(len(negs))].astype(np.float32)
            valid = np.array([fnv('v' + str(s)) % 100 < 15 for s in pair_key[rows]])
            sc = fit(rows, y, np.ones(len(rows)), 200 + 10 * k + r, valid); log.append(f'r{r} pos {len(pos)} neg {len(negs)}')
        te = np.where(fold == k)[0]; scores[te] = sc(te)
        print(f'{RAW} {arm} fold {k}: ' + ' | '.join(log), flush=True)
    np.save(f'{O}/{RAW}/scores_{arm}.npy', scores); print(f'{RAW} {arm} done in {time.time() - t0:.0f} s', flush=True)

for a in ARMS: run(a)
