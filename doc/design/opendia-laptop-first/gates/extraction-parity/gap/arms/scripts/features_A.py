"""PLAN_GAP arm A: fragment mass-consistency and dominance features per population window (featA.npy [n, 4] float32).
Label-free. From the 6-fragment arrays (xic.npy [n,6,32], ppm.npy [n,6,32] = intensity x ppm to the calibrated target,
lib.npy [n,6]) over the central bins 12:20 (apex +- 4.7 s, as features12.py / deficit_features.py), and the run's
fitted fragment m/z model in cal.json (ms2_mz: signal_weights, signal_sigmas_ppm; the uniform background spans the
extraction window +- 17.5 ppm of build_voxels.py).
  area_f   central area of fragment f;  ppm_f  intensity-weighted mean ppm of f over the central bins (area_f > 0)
  a_ppm_llr   sum over fragments with signal of log p_signal(ppm_f) - log p_background(ppm_f): how much the fragments'
              masses look like the run's real-signal mass errors rather than random ions in the extraction window
  a_ppm_sd    SD of ppm_f across the fragments with signal (>= 2); otherwise the SD of the background (35 / sqrt 12)
  a_top_share largest area_f / total central area (1 when no fragment has signal)
  a_top_excess log of (top fragment's observed share) / (its library share among the 6), 0 when no signal
Usage: features_A.py <raw01|raw02> <root>."""
import json, sys, warnings, numpy as np
RAW, O = sys.argv[1], sys.argv[2]
C0, C1, BS, HALF_PPM = 12, 20, 100000, 17.5
cal = json.load(open(f'{O}/{RAW}/cal/cal.json'))['ms2_mz']
W = np.asarray(cal['signal_weights'], float); S = np.asarray(cal['signal_sigmas_ppm'], float); W = W / W.sum()
X = np.load(f'{O}/{RAW}/xic.npy', mmap_mode='r'); P = np.load(f'{O}/{RAW}/ppm.npy', mmap_mode='r')
L = np.load(f'{O}/{RAW}/lib.npy', mmap_mode='r'); n = len(X)
warnings.filterwarnings('ignore', category=RuntimeWarning)
BG_SD = 2 * HALF_PPM / np.sqrt(12)


def log_signal(x):
    d = (W[None, None, :] / (S[None, None, :] * np.sqrt(2 * np.pi))) * np.exp(-0.5 * (x[..., None] / S[None, None, :]) ** 2)
    return np.log(d.sum(-1) + 1e-300)


out = np.zeros((n, 4), np.float32)
for b in range(0, n, BS):
    x = np.asarray(X[b:b + BS], np.float64)[:, :, C0:C1]; p = np.asarray(P[b:b + BS], np.float64)[:, :, C0:C1]
    li = np.asarray(L[b:b + BS], np.float64)
    area = x.sum(2); has = area > 0
    ppm = np.where(has, p.sum(2) / np.maximum(area, 1e-300), np.nan)
    llr = np.where(has, log_signal(np.nan_to_num(ppm)) - np.log(1 / (2 * HALF_PPM)), 0.0).sum(1)
    k = has.sum(1)
    sd = np.where(k >= 2, np.nanstd(ppm, 1), BG_SD)
    tot = area.sum(1); top = area.argmax(1); r = np.arange(len(x))
    share = np.where(tot > 0, area[r, top] / np.maximum(tot, 1e-300), 1.0)
    lshare = li[r, top] / np.maximum(li.sum(1), 1e-300)
    excess = np.where((tot > 0) & (lshare > 0), np.log(share / np.maximum(lshare, 1e-300)), 0.0)
    out[b:b + BS] = np.c_[llr, sd, share, excess]
np.save(f'{O}/{RAW}/featA.npy', out)
q = np.percentile(out, [5, 50, 95], axis=0).round(3)
print(RAW, 'featA', out.shape, 'finite', bool(np.isfinite(out).all()), 'p5/p50/p95 per column', q.T.tolist())
