"""PLAN_V4 phase 2: 12-fragment evidence features per population window (label-free), from build12.py outputs
(xic12.npy [n,12,32], lib12.npy [n,12], orig12.npy [n,12]). Central bins C0:C1 = apex +- 4.7 s (as ms1_features.py).
Usage: features12.py <raw01|raw02>. Output: v2/<raw>/feat12.npy [n, 6] float32, columns FEATS12:
  f12_support      fragments (of 12) whose central-bin maximum is above 0 and whose profile correlates > 0.5 with the
                   summed profile of the other fragments
  f12_coelution    mean of those correlations over the fragments present in the library (0 for absent fragments)
  f12_libcos       cosine between sqrt(central area) and sqrt(predicted intensity), over all 12
  f12_extra_support  as f12_support, counted over the non-original (extra) fragments only
  f12_extra_coelution mean correlation over the extra fragments
  f12_extra_share  log10 of (extra-fragment central area + 1) / (all-fragment central area + 1)"""
import sys, numpy as np
RAW = sys.argv[1]
O = '/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-voxel-pilot-20261110/v2'
FEATS12 = ['f12_support', 'f12_coelution', 'f12_libcos', 'f12_extra_support', 'f12_extra_coelution', 'f12_extra_share']
C0, C1 = 12, 20
X = np.load(f'{O}/{RAW}/xic12.npy', mmap_mode='r'); LI = np.load(f'{O}/{RAW}/lib12.npy', mmap_mode='r')
OR = np.load(f'{O}/{RAW}/orig12.npy', mmap_mode='r'); n = len(X)
out = np.zeros((n, len(FEATS12)), np.float32); B = 100000
for b in range(0, n, B):
    x = np.asarray(X[b:b + B], np.float64); li = np.asarray(LI[b:b + B], np.float64); o = np.asarray(OR[b:b + B]) > 0
    present = li > 0; extra = present & ~o
    tot = x.sum(1, keepdims=True); other = tot - x
    xc = x - x.mean(2, keepdims=True); oc = other - other.mean(2, keepdims=True)
    r = (xc * oc).sum(2) / np.sqrt((xc * xc).sum(2) * (oc * oc).sum(2) + 1e-12)
    r = np.where(present, r, 0.0)
    sig = x[:, :, C0:C1].max(2) > 0
    ok = present & sig & (r > 0.5)
    area = x[:, :, C0:C1].sum(2)
    f = np.zeros((len(x), len(FEATS12)))
    f[:, 0] = ok.sum(1)
    f[:, 1] = r.sum(1) / np.maximum(present.sum(1), 1)
    a, l = np.sqrt(area), np.sqrt(li)
    f[:, 2] = (a * l).sum(1) / (np.linalg.norm(a, axis=1) * np.linalg.norm(l, axis=1) + 1e-12)
    f[:, 3] = (ok & extra).sum(1)
    f[:, 4] = np.where(extra, r, 0).sum(1) / np.maximum(extra.sum(1), 1)
    f[:, 5] = np.log10((np.where(extra, area, 0).sum(1) + 1) / (area.sum(1) + 1))
    out[b:b + B] = np.nan_to_num(f).astype(np.float32)
    if b % 500000 == 0: print(RAW, b, '/', n, flush=True)
np.save(f'{O}/{RAW}/feat12.npy', out)
print(RAW, 'feature medians', dict(zip(FEATS12, np.median(out, 0).round(3))), 'means', dict(zip(FEATS12, out.mean(0).round(3))))
