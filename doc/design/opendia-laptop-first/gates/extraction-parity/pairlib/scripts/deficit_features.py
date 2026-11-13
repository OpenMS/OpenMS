"""PLAN_INTERF: label-free intensity-aware fragment-deficit features per population window (featD.npy).
For each window, from the 12-fragment arrays (xic12 [n,12,32], lib12 [n,12], orig12 [n,12]), central bins 12:20:
  area_f      central area of fragment f; library fragments are those with lib12 > 0
  scale       median over library fragments with area_f > 0 of area_f / lib12_f (0 when no fragment has signal)
  deficit_f   log10(scale * lib12_f + N) - log10(area_f + N): how far a fragment falls below the abundance implied by the
              observed fragments, with the run's noise level N as a floor
  N           median nonzero central library-fragment area over 200,000 windows drawn with numpy seed 0 (one per run)
Columns: d_max (max deficit over the 12), d6_max (max over the 6 original fragments), d_lib (library-intensity share
of fragments with deficit > 1). Never reads labels. Usage: deficit_features.py <raw01|raw02> <root>."""
import sys, warnings, numpy as np
RAW, O = sys.argv[1], sys.argv[2]
C0, C1, BS = 12, 20, 100000
X = np.load(f'{O}/{RAW}/xic12.npy', mmap_mode='r'); L = np.load(f'{O}/{RAW}/lib12.npy', mmap_mode='r')
G = np.load(f'{O}/{RAW}/orig12.npy', mmap_mode='r'); n = len(X)
samp = np.sort(np.random.default_rng(0).choice(n, 200000, replace=False))
xa = np.asarray(X[samp], np.float64)[:, :, C0:C1].sum(2); la = np.asarray(L[samp]) > 0
N = float(np.median(xa[(xa > 0) & la]))
warnings.filterwarnings('ignore', category=RuntimeWarning)


def deficit(area, li, mask):
    present = (li > 0) & mask; sig = present & (area > 0)
    ratio = np.where(sig, area / np.maximum(li, 1e-12), np.nan)
    scale = np.where(sig.any(1), np.nanmedian(ratio, 1), 0.0)
    d = np.where(present, np.log10(scale[:, None] * li + N) - np.log10(area + N), -np.inf)
    dmax = np.where(present.any(1), d.max(1), 0.0)
    big = present & (d > 1)
    return dmax, (li * big).sum(1) / np.maximum((li * present).sum(1), 1e-12)


out = np.zeros((n, 3), np.float32)
for b in range(0, n, BS):
    x = np.asarray(X[b:b + BS], np.float64)[:, :, C0:C1].sum(2); li = np.asarray(L[b:b + BS], np.float64)
    og = np.asarray(G[b:b + BS]) > 0
    d_all, d_lib = deficit(x, li, np.ones_like(og)); d6, _ = deficit(x, li, og)
    out[b:b + BS] = np.nan_to_num(np.c_[d_all, d6, d_lib]).astype(np.float32)
np.save(f'{O}/{RAW}/featD.npy', out)
print(RAW, 'windows', n, 'N', round(N, 1), 'medians d_max/d6_max/d_lib', np.median(out, 0).round(3), 'means', out.mean(0).round(3), flush=True)
