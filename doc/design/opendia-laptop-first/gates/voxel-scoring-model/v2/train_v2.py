"""Voxel v2 step 3: train and score the v2 arms on the full v2 population of one raw (PLAN_V2.md).
Usage: train_v2.py <raw01|raw02> <arm> [<arm> ...]. Arms (input option b, the pilot ppm-channel input):
  V0  voxel model, T0 training set, no rounds (pilot control on the full population)
  F0  V0 + the 11 heuristic features in the head (search 4, boundary 5, byion 2), standardised per training fold
  S0  V0 + semi-supervised rounds
  FS  F0 + semi-supervised rounds
Semi-supervised rounds (PyProphet / OpenDIACNN -stage2 schedule, inside the training folds only): round 0 trains on T0
(calibration anchors, their paired decoys, 6000 random decoys). Rounds 1 and 2 re-score every window of the training
folds with the round's model, take per-precursor winners, run pooled paired TDC on the training folds, and train a new
model from scratch on positives = target winners at q <= 0.15 (round 1) / 0.05 (round 2), negatives = all decoy
winners of the training folds (subsampled to 4x positives, at least 6000). The held-out fold is scored only by the
final model of its fold. Never reads entrapment labels.
Output: <raw>/scores_<arm>.npy, one logit per population row (row order of population.tsv)."""
import sys, time, numpy as np, pandas as pd, torch, torch.nn as nn, torch.nn.functional as F
sys.path.insert(0, '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110')
from train_pilot import Branch, Coelution, cross, dual, rowcos, CORE  # same TOPAZ body as the pilot
RAW, ARMS = sys.argv[1], sys.argv[2:]
O = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2'
torch.set_num_threads(10)
FEATS = ['transition_support', 'log_matched_scans', 'log_peak_width', 'abs_library_rt_error', 'boundary_log_fwhm', 'boundary_log_width',
         'boundary_coelution', 'boundary_library_correlation', 'boundary_apex_agreement', 'byion_b', 'byion_y']
pop = pd.read_csv(f'{O}/{RAW}/population.tsv', sep='\t')
X = {k: np.load(f'{O}/{RAW}/{k}.npy', mmap_mode='r') for k in ('xic', 'ppm', 'ms1', 'lib')}
H = np.load(f'{O}/{RAW}/heuristics.npy')  # [n, 11] in FEATS order, row order of population.tsv
train_rows = pd.read_csv(f'{O}/{RAW}/train_windows.tsv', sep='\t')
row_of = pd.Series(np.arange(len(pop)), index=pop.sample_id.values)
fold = pop.fold.values; is_decoy = pop.is_decoy.values; pair_key = pop.pair_key.values; pid = pop.precursor_id.values

class Scorer(nn.Module):
    """Pilot option b (dual view + ppm channel per fragment, MS1 M+0..2) with an optional heuristic block in the head."""
    def __init__(self, n_feat, w=32):
        super().__init__()
        self.xic, self.ms1 = Branch(18, w), Branch(6, w)
        self.xco, self.mco = Coelution(6), Coelution(3)
        self.n_feat = n_feat
        f = 2 * w + 18 + (2 * w + 18 + 4) + 2 + 6 + n_feat
        self.head = nn.Sequential(nn.Linear(f, 4 * w), nn.ReLU(), nn.Dropout(0.1), nn.Linear(4 * w, 2 * w), nn.ReLU(), nn.Dropout(0.1), nn.Linear(2 * w, 1))
    def forward(self, xic, ppm, ms1, lib, h):
        m = ms1; has = (m.flatten(1).amax(1, keepdim=True) > 0).float()
        ms1p = torch.cat([self.ms1(dual(m)), self.mco(m), cross(m, xic)], 1) * has
        core = xic[:, :, CORE[0]:CORE[1]].sum(2)
        fit = torch.stack([rowcos(torch.sqrt(core), torch.sqrt(lib)), rowcos(core, lib)], 1)
        parts = [self.xic(torch.cat([dual(xic), ppm], 1)), self.xco(xic), ms1p, fit, lib]
        if self.n_feat: parts.append(h)
        return self.head(torch.cat(parts, 1)).squeeze(1)

def batch(rows, mu, sd):
    rows = np.sort(rows)
    xic = torch.from_numpy(np.asarray(X['xic'][rows], np.float32)); p = np.asarray(X['ppm'][rows], np.float32)
    ppm = torch.from_numpy(np.clip(np.where(xic.numpy() > 0, p / np.maximum(xic.numpy(), 1e-12), 0) / 5.0, -4, 4).astype(np.float32))
    ms1 = torch.from_numpy(np.asarray(X['ms1'][rows, :3], np.float32)); lib = torch.from_numpy(np.asarray(X['lib'][rows], np.float32))
    h = torch.from_numpy(((H[rows] - mu) / sd).astype(np.float32))
    return rows, (xic, ppm, ms1, lib, h)

def score_rows(model, rows, mu, sd, bs=4096):
    model.eval(); out = np.empty(len(rows), np.float32); order = np.argsort(rows)
    with torch.no_grad():
        for b in range(0, len(rows), bs):
            r, t = batch(rows[order[b:b + bs]], mu, sd)
            out[order[b:b + bs]] = model(*t).numpy()
    return out

def fit(rows, y, w, n_feat, mu, sd, seed, valid_mask):
    """Train from scratch; early stopping on the weighted loss of the validation families (patience 5, max 40)."""
    torch.manual_seed(seed); rng = np.random.default_rng(seed)
    w = w.copy(); w[y == 1] *= w[y == 0].sum() / max(w[y == 1].sum(), 1e-9); w /= w.mean()
    tr, va = np.where(~valid_mask)[0], np.where(valid_mask)[0]
    Rva, Tva = batch(rows[va], mu, sd); yva = torch.from_numpy(y[va][np.argsort(rows[va])]); wva = torch.from_numpy(w[va][np.argsort(rows[va])].astype(np.float32))
    model = Scorer(n_feat); opt = torch.optim.Adam(model.parameters(), 1e-3); best, stale, state = np.inf, 0, None
    for epoch in range(40):
        model.train(); perm = rng.permutation(tr)
        for b in range(0, len(perm), 64):
            i = perm[b:b + 64]; r, t = batch(rows[i], mu, sd); o = np.argsort(rows[i])
            loss = (F.binary_cross_entropy_with_logits(model(*t), torch.from_numpy(y[i][o]), reduction='none') * torch.from_numpy(w[i][o].astype(np.float32))).mean()
            opt.zero_grad(); loss.backward(); opt.step()
        model.eval()
        with torch.no_grad():
            vl = float((F.binary_cross_entropy_with_logits(model(*Tva), yva, reduction='none') * wva).mean())
        if vl < best - 1e-4: best, stale, state = vl, 0, {a: t.clone() for a, t in model.state_dict().items()}
        else:
            stale += 1
            if stale >= 5: break
    model.load_state_dict(state)
    return model, epoch + 1, best

def fnv(s):
    h = 0xcbf29ce484222325
    for b in s.encode(): h = ((h ^ b) * 0x100000001b3) & 0xFFFFFFFFFFFFFFFF
    return h

def pooled_q(rows, s):
    """Pooled paired TDC over rows (OpenDIACNN -pooled_tdc): precursor winners, per-fold decoy-winner median/MAD,
    best of each pair, q = (D+1)/T. Returns winner rows, their q and their target flag."""
    d = pd.DataFrame({'r': rows, 's': s, 'pid': pid[rows], 'f': fold[rows], 'dec': is_decoy[rows], 'key': pair_key[rows]})
    w = d.sort_values(['s', 'r'], ascending=[False, True]).drop_duplicates('pid')
    for f_, g in w[w.dec == 1].groupby('f'):
        med = g.s.median(); mad = 1.4826 * (g.s - med).abs().median()
        w.loc[w.f == f_, 's'] = (w.loc[w.f == f_, 's'] - med) / mad
    best = w.sort_values(['s', 'dec'], ascending=[False, False]).drop_duplicates('key').sort_values('s', ascending=False)
    T = np.cumsum(1 - best.dec.values); D = np.cumsum(best.dec.values)
    q = np.minimum.accumulate(np.minimum(1, (D + 1) / np.maximum(T, 1))[::-1])[::-1]
    return best.r.values, q, best.dec.values, w

def run_arm(arm):
    n_feat = 11 if arm in ('F0', 'FS') else 0
    rounds = 2 if arm in ('S0', 'FS') else 0
    scores = np.zeros(len(pop), np.float32); t0 = time.time()
    for k in range(3):
        tr_mask = fold != k; tr_rows_all = np.where(tr_mask)[0]
        mu, sd = H[tr_rows_all].mean(0), H[tr_rows_all].std(0); sd[sd == 0] = 1
        t0rows = train_rows[(train_rows.fold != k) & train_rows.group.isin(['anchor', 'anchor_decoy', 'decoy'])]
        rows = row_of[t0rows.sample_id].values; y = t0rows.label.values.astype(np.float32); w = t0rows.weight.values.astype(np.float64)
        fam = t0rows.family.values
        model, ep, vl = fit(rows, y, w, n_feat, mu, sd, 100 + k, np.array([fnv('v' + str(s)) % 100 < 15 for s in fam]))
        log = [f'r0 {len(rows)} rows {ep} ep vl {vl:.4f}']
        for r in range(1, rounds + 1):
            s = score_rows(model, tr_rows_all, mu, sd)
            win, q, dec, wtab = pooled_q(tr_rows_all, s)
            cut = 0.15 if r == 1 else 0.05
            pos = win[(dec == 0) & (q <= cut)]
            negs = wtab[wtab.dec == 1].r.values
            rng = np.random.default_rng(1000 + 10 * k + r)
            negs = rng.choice(negs, min(len(negs), max(6000, 4 * len(pos))), replace=False)
            rows = np.concatenate([pos, negs]); y = np.r_[np.ones(len(pos)), np.zeros(len(negs))].astype(np.float32)
            w = np.ones(len(rows)); keys = pair_key[rows]
            model, ep, vl = fit(rows, y, w, n_feat, mu, sd, 200 + 10 * k + r, np.array([fnv('v' + str(s)) % 100 < 15 for s in keys]))
            log.append(f'r{r} pos {len(pos)} neg {len(negs)} {ep} ep vl {vl:.4f}')
        te = np.where(~tr_mask)[0]; scores[te] = score_rows(model, te, mu, sd)
        print(f'{RAW} {arm} fold {k}: ' + ' | '.join(log), flush=True)
    np.save(f'{O}/{RAW}/scores_{arm}.npy', scores)
    print(f'{RAW} {arm} done in {time.time() - t0:.0f} s', flush=True)

for a in ARMS: run_arm(a)
