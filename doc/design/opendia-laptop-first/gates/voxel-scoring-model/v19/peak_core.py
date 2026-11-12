"""Native-sampling peak measurement helpers (pure, label-free)."""
import numpy as np

def local_peak(rt, y, t0, search=2):
    """Peak of trace y (native scans, rt ascending) near t0.
    Apex: the largest y among the scans within `search` scans of the scan nearest t0.
    FWHM: half-maximum crossings walking out from the apex, linearly interpolated; if a side never drops below half
    maximum inside the trace, that side ends at the trace edge.
    Returns (apex_index, fwhm_s, points_above_half) with points counted contiguously around the apex."""
    rt = np.asarray(rt, float); y = np.asarray(y, float); n = len(y)
    if n == 0: return -1, 0.0, 0
    j0 = int(np.argmin(np.abs(rt - t0)))
    lo, hi = max(0, j0 - search), min(n, j0 + search + 1)
    a = lo + int(np.argmax(y[lo:hi]))
    h = y[a] / 2.0
    if y[a] <= 0: return a, 0.0, 0
    i = a
    while i > 0 and y[i - 1] > h: i -= 1
    left = rt[0] if i == 0 and y[0] > h else (rt[i - 1] + (h - y[i - 1]) * (rt[i] - rt[i - 1]) / (y[i] - y[i - 1]) if i > 0 else rt[0])
    k = a
    while k < n - 1 and y[k + 1] > h: k += 1
    right = rt[-1] if k == n - 1 else rt[k] + (y[k] - h) * (rt[k + 1] - rt[k]) / (y[k] - y[k + 1])
    return a, float(right - left), int(k - i + 1)

def support(M, a):
    """Fragments matched in >= 2 of the 3 scans centred on apex index a. M: [scans, fragments] bool."""
    lo, hi = max(0, a - 1), min(M.shape[0], a + 2)
    return int((M[lo:hi].sum(0) >= 2).sum())
