"""PLAN_V5 features (label-free), per population window. Output v2/<raw>/featPC.npy [n, 5] float32, columns FEATS_PC:
  p_n_tight      original fragments (of 6) with signal in the central bins (apex +- 4.7 s) and |intensity-weighted mean
                 ppm - calibrated offset| <= 3 (from xic.npy / ppm.npy)
  p_ppm_abs_med  median over fragments with signal of that |mean ppm| (17.5 when no fragment has signal)
  p_tight_share  central-bin intensity share of the tight fragments
  c_med_top3     mean of the 3 highest correlations of a fragment (xic12, scaled to unit maximum) with the per-bin median
                 consensus of the fragments with signal
  c_med_k        fragments (of 12) with that correlation >= 0.8
Usage: features_v5.py <raw01|raw02>."""
import sys, numpy as np
RAW = sys.argv[1]
O = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2'
FEATS_PC = ['p_n_tight', 'p_ppm_abs_med', 'p_tight_share', 'c_med_top3', 'c_med_k']
X = np.load(f'{O}/{RAW}/xic.npy', mmap_mode='r'); P = np.load(f'{O}/{RAW}/ppm.npy', mmap_mode='r')
X12 = np.load(f'{O}/{RAW}/xic12.npy', mmap_mode='r'); L12 = np.load(f'{O}/{RAW}/lib12.npy', mmap_mode='r')
n = len(X); out = np.zeros((n, len(FEATS_PC)), np.float32); B = 50000

def corr(a, b):
    a = a - a.mean(-1, keepdims=True); b = b - b.mean(-1, keepdims=True)
    return (a * b).sum(-1) / np.sqrt((a * a).sum(-1) * (b * b).sum(-1) + 1e-12)

for b in range(0, n, B):
    x = np.asarray(X[b:b + B], np.float64)[:, :, 12:20]; p = np.asarray(P[b:b + B], np.float64)[:, :, 12:20]
    a = x.sum(2); sig = a > 0; m = np.where(sig, p.sum(2) / np.maximum(a, 1e-12), np.nan)
    tight = sig & (np.abs(np.nan_to_num(m, nan=99.0)) <= 3)
    with np.errstate(all='ignore'):
        med = np.nanmedian(np.abs(m), 1)
    out[b:b + B, 0] = tight.sum(1); out[b:b + B, 1] = np.nan_to_num(med, nan=17.5)
    out[b:b + B, 2] = (a * tight).sum(1) / np.maximum(a.sum(1), 1e-12)
    x12 = np.asarray(X12[b:b + B], np.float64); li = np.asarray(L12[b:b + B])
    s12 = (li > 0) & (x12[:, :, 12:20].max(2) > 0)
    xn = x12 / np.maximum(x12.max(2, keepdims=True), 1e-12)
    with np.errstate(all='ignore'):
        cons = np.nan_to_num(np.nanmedian(np.where(s12[:, :, None], xn, np.nan), 1))
    rm = np.where(s12, corr(xn, cons[:, None, :]), 0.0)
    out[b:b + B, 3] = -np.sort(-rm, 1)[:, :3].mean(1); out[b:b + B, 4] = (rm >= 0.8).sum(1)
    if b % 500000 == 0: print(RAW, b, '/', n, flush=True)
out = np.nan_to_num(out)
np.save(f'{O}/{RAW}/featPC.npy', out)
print(RAW, 'medians', dict(zip(FEATS_PC, np.median(out, 0).round(3))), 'means', dict(zip(FEATS_PC, out.mean(0).round(3))))
