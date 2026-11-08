"""Baseline for step B: the rule picker (numpy port of OpenDIAPeakPicker, same defaults) on the augmented 6x32 windows.
A learned boundary model is only worth adding if it beats this on held-out folds, especially on weak/interfered windows."""
import sys
import numpy as np

d = np.load(sys.argv[1])
x, y, neg, ls = d['x'], d['y'], d['negative'].astype(bool), d['log_snr']


def pick(w, hint=16, apex_search=2, floor=0.05):
    s = w.sum(0).astype(np.float64)
    t = s.copy()
    t[1:-1] = 0.25 * s[:-2] + 0.5 * s[1:-1] + 0.25 * s[2:]
    lo, hi = max(0, hint - apex_search), min(31, hint + apex_search)
    apex = hint
    for i in range(lo, hi + 1):
        if t[i] > t[apex]:
            apex = i
    if t[apex] <= 0:
        return None
    left = apex
    while left > 0:
        if t[left - 1] > t[left]:
            break
        left -= 1
        if t[left] <= floor * t[apex]:
            break
    right = apex
    while right < 31:
        if t[right + 1] > t[right]:
            break
        right += 1
        if t[right] <= floor * t[apex]:
            break
    return left, right, apex


def iou(a, b):
    inter = max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)
    return inter / (max(a[1], b[1]) - min(a[0], b[0]) + 1)


rows = []
for i in np.where(~neg)[0]:
    true = (int(np.argmax(y[i, 1])), int(np.argmax(y[i, 2])))
    # Hint: the strongest consensus bin (no label information), as a candidate apex from the search would be.
    hint = int(np.argmax(x[i].sum(0)))
    p = pick(x[i], hint)
    if p is None:
        rows.append((ls[i], 0.0, 99, 99, 0))
        continue
    rows.append((ls[i], iou(true, p[:2]), abs(p[0] - true[0]), abs(p[1] - true[1]), int(true[0] <= p[2] <= true[1])))
r = np.array(rows)
print(f'rule picker on {len(r)} augmented positive windows (bins of 0.94 s):')
print(f'  mean IoU {r[:, 1].mean():.3f} | median |start err| {np.median(r[:, 2]):.0f} bins, |end err| {np.median(r[:, 3]):.0f} | apex inside true peak {r[:, 4].mean():.2%}')
for lo, hi in [(-9, 0.0), (0.0, 0.5), (0.5, 1.0), (1.0, 9)]:
    k = (r[:, 0] >= lo) & (r[:, 0] < hi)
    print(f'  log10 SNR [{lo:+.1f},{hi:+.1f}): n {k.sum():5d} IoU {r[k, 1].mean():.3f} start<=1 bin {np.mean(r[k, 2] <= 1):.2%} end<=1 bin {np.mean(r[k, 3] <= 1):.2%} apex-in {r[k, 4].mean():.2%}')
