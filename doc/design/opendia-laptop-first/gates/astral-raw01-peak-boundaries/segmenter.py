"""Prototype step B: small 1D boundary segmenter on augmented anchor windows, cross-fitted by anchor fold.
Compared with the rule picker on the held-out fold only (anchors never seen in training). Diagnostic, raw 01 design data."""
import sys, time
import numpy as np, torch, torch.nn as nn, torch.nn.functional as F

torch.manual_seed(0)
torch.set_num_threads(4)
d = np.load(sys.argv[1])
X = torch.tensor(d['x'])          # N x 6 x 32
Y = torch.tensor(d['y'])          # N x 3 x 32 (inside, start, end)
fold = d['fold']
neg = d['negative'].astype(bool)
snr = d['log_snr']
scan_mask = d['scan_mask']
print('bins with a scan per window, median', int(np.median(scan_mask.sum(1))), 'of 32')


def norm(x):
    # Per-window max normalisation then log1p: scale invariant, like -topaz_faithful.
    m = x.amax(dim=(1, 2), keepdim=True).clamp_min(1e-6)
    return torch.log1p(9.0 * x / m) / np.log(10.0)


class Segmenter(nn.Module):
    """Fragment-shared conv (order invariant) -> mean/max pool over fragments -> dilated RT convs -> 3 per-bin logits."""
    def __init__(self, w=12):
        super().__init__()
        self.frag = nn.Sequential(nn.Conv1d(1, w, 3, padding=1), nn.GELU(), nn.Conv1d(w, w, 3, padding=1), nn.GELU())
        self.rt = nn.ModuleList([nn.Conv1d(2 * w + 1, 2 * w, 3, padding=1, dilation=1),
                                 nn.Conv1d(2 * w, 2 * w, 3, padding=2, dilation=2),
                                 nn.Conv1d(2 * w, 2 * w, 3, padding=4, dilation=4),
                                 nn.Conv1d(2 * w, 2 * w, 3, padding=8, dilation=8)])
        self.head = nn.Conv1d(2 * w, 3, 1)

    def forward(self, x):
        n, f, t = x.shape
        h = self.frag(x.reshape(n * f, 1, t)).reshape(n, f, -1, t)
        z = torch.cat([h.mean(1), h.amax(1), x.mean(1, keepdim=True)], 1)
        for i, conv in enumerate(self.rt):
            z = F.gelu(conv(z)) + (z if i > 0 else 0)
        return self.head(z)


def decode(logits):
    """Peak present if max inside prob > 0.5; start = argmax start prob, end = argmax end prob at or after start."""
    p = torch.sigmoid(logits)
    present = p[:, 0].amax(1) > 0.5
    start = p[:, 1].argmax(1)
    cum = torch.arange(32).unsqueeze(0) >= start.unsqueeze(1)
    end = torch.where(cum, p[:, 2], torch.full_like(p[:, 2], -1)).argmax(1)
    return present, start, end


def rule(x, masks, lo=10, hi=22):
    """Rule picker (OpenDIAPeakPicker defaults) on the bins that received a scan (TraceCollector ms2_mask), as stage 2
    does; boundaries map back to bin indices. hint = strongest consensus scan bin within +-6 bins of the centre."""
    out = []
    for w, mask in zip(x.numpy(), masks):
        idx = np.where(mask)[0]
        if len(idx) < 3:
            out.append((-1, -1)); continue
        s = w[:, idx].sum(0).astype(np.float64)
        t = s.copy()
        t[1:-1] = 0.25 * s[:-2] + 0.5 * s[1:-1] + 0.25 * s[2:]
        cand = np.where((idx >= lo) & (idx <= hi))[0]
        hint = int(cand[np.argmax(s[cand])]) if len(cand) else int(np.argmax(s))
        n = len(idx) - 1
        apex = max(range(max(0, hint - 2), min(n, hint + 2) + 1), key=lambda i: t[i])
        if t[apex] <= 0:
            out.append((-1, -1)); continue
        l = apex
        while l > 0 and not t[l - 1] > t[l]:
            l -= 1
            if t[l] <= 0.05 * t[apex]: break
        r = apex
        while r < n and not t[r + 1] > t[r]:
            r += 1
            if t[r] <= 0.05 * t[apex]: break
        out.append((idx[l], idx[r]))
    return np.array(out)


def iou(a0, a1, b0, b1):
    inter = np.maximum(0, np.minimum(a1, b1) - np.maximum(a0, b0) + 1)
    return inter / (np.maximum(a1, b1) - np.minimum(a0, b0) + 1)


bce = nn.BCEWithLogitsLoss(reduction='none')
results = []
params = None
t_train = t_inf = 0.0
for k in sorted(set(fold.tolist())):
    tr, te = np.where(fold != k)[0], np.where(fold == k)[0]
    model = Segmenter()
    params = sum(p.numel() for p in model.parameters())
    opt = torch.optim.AdamW(model.parameters(), lr=3e-3, weight_decay=0.01)
    xt, yt = norm(X[tr]), Y[tr]
    # Boundary heatmaps are 1-hot; weight positives to balance the 1-in-32 rate.
    wpos = torch.tensor([1.0, 8.0, 8.0]).view(1, 3, 1)
    start = time.time()
    epochs, bs = 12, 256
    steps = epochs * (len(tr) // bs)
    sched = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=3e-3, total_steps=steps, pct_start=0.05)
    for ep in range(epochs):
        perm = torch.randperm(len(tr))
        for i in range(0, len(tr) - bs + 1, bs):
            b = perm[i:i + bs]
            out = model(xt[b])
            target = yt[b]
            loss = (bce(out, target) * (1 + (wpos - 1) * target)).mean()
            opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(model.parameters(), 5.0); opt.step(); sched.step()
    t_train += time.time() - start
    model.eval()
    with torch.no_grad():
        start = time.time()
        logits = model(norm(X[te]))
        t_inf += time.time() - start
    present, s, e = (v.numpy() for v in decode(logits))
    ts, tee = Y[te, 1].argmax(1).numpy(), Y[te, 2].argmax(1).numpy()
    rb = rule(X[te], d['scan_mask'][te])
    for j, i in enumerate(te):
        results.append((k, neg[i], snr[i], present[j], s[j], e[j], rb[j, 0], rb[j, 1], ts[j], tee[j]))
r = np.array(results, dtype=float)
k, ng, sn, pr, s, e, rl, rr, ts, te_ = r.T
ng = ng.astype(bool)
print(f'Segmenter {params} params | train {t_train:.0f} s for 3 folds | inference {1e6 * t_inf / len(r):.1f} us/window (batch, 4 threads)')
print(f'held-out negatives: model says no peak {np.mean(~pr[ng].astype(bool)):.2%} (rule picker always returns a peak)')
pos = ~ng
mi = iou(s, e, ts, te_)
ri = np.where(rl >= 0, iou(rl, rr, ts, te_), 0)
print(f'held-out positives n={pos.sum()}: IoU model {mi[pos].mean():.3f} vs rule {ri[pos].mean():.3f} | detected by model {pr[pos].mean():.2%}')
for lo, hi in [(-9, 0.0), (0.0, 0.5), (0.5, 1.0), (1.0, 9)]:
    m = pos & (sn >= lo) & (sn < hi)
    print(f'  log10 SNR [{lo:+.1f},{hi:+.1f}) n {m.sum():5d}: IoU model {mi[m].mean():.3f} rule {ri[m].mean():.3f} | '
          f'start<=1 bin model {np.mean(np.abs(s[m] - ts[m]) <= 1):.0%} rule {np.mean(np.abs(rl[m] - ts[m]) <= 1):.0%} | '
          f'end<=1 bin model {np.mean(np.abs(e[m] - te_[m]) <= 1):.0%} rule {np.mean(np.abs(rr[m] - te_[m]) <= 1):.0%}')
