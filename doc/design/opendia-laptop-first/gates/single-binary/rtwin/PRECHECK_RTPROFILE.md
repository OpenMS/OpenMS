# PRECHECK_RTPROFILE: does the RT error depend on gradient position? (16 November 2026)

Written before any plan for an RT-dependent extraction window. It uses the CUR seed-0 runs of PLAN_RTWIN2
(`library12_full`, self-calibrating) and their calibration JSONs. Only target/decoy status is used, plus the accepted IDs
(q ≤ 0.01).

## 1. Accepted IDs, |apex − mapped RT| by position relative to the anchor span
| zone | raw01 n / p95 / p99 (s) | raw02 | raw03 |
|---|---|---|---|
| interior of the anchor span | 1536 / 48 / 70 | 1642 / 47 / 66 | 1789 / 53 / 72 |
| first 60 s of anchors | 84 / 52 / 82 | 83 / 91 / 117 | 20 / 38 / 48 |
| last 60 s of anchors | 61 / 100 / 116 | 66 / 99 / 113 | 52 / 82 / 98 |
| before the first anchor (extrapolated) | 72 / 130 / 149 (median −65) | 12 / 122 / 122 | 6 / 36 / 39 |
| after the last anchor (extrapolated) | 19 / 149 / 150 (median +101) | 11 / 114 / 119 | 9 / 76 / 89 |

Anchor spans: raw01 508-1176 s, raw02 410-1177 s, raw03 188-1177 s. T is 151 / 124 / 107 s.

- **Inside the anchor span the error is nearly uniform:** p95 is 47-53 s on all three raws.
- **The wide tails sit at the edges:** in the last 60 s of anchors, and where the map is extrapolated (ends without
  anchors).
- There, the p99 is at T and the medians are biased (raw01 −65 s before, +101 s after). True peaks are being cut off,
  and the map is wrong there, not only noisy.
- Raw01 has the most IDs before its first anchor (72), because its anchors start at 508 s.

## 2. Anchor residuals binned by mapped RT (Gaussian per bin; K by 5-fold held-out log likelihood)
- **raw01 (426 anchors):** best K = 2, gain +7 log-likelihood units. The rms is 15.7 s (499-802 s) and 20.0 s
  (803-1178 s).
- **raw02 (546 anchors):** best K = 10, gain +34. The rms is 6.9 s at 410-576 s and 14-22 s after that.
- The interior heterogeneity is real but modest: at most a factor of about 1.4, apart from raw02's early bin.

## 3. Rejected method
I also tried a per-bin Gaussian-plus-decoy-background fit on all windows (the calibration's spread fit, binned).
- It is unusable: the fitted offsets drift by up to ±100 s, and the true-peak share is 0.01-0.08.
- The chance background dominates, and the bins have too little signal.

## Reading
- Binning the interior would save little: the middle already needs about ±50 s for 95 %, against ±107-151 s now.
- The bigger effects are:
  1. The window is too narrow, and the centre biased, where the map is extrapolated or has few anchors.
  2. The window is wider than needed in the interior.
- An RT-dependent window should therefore model **map uncertainty**: anchor density and distance to the nearest anchors,
  as well as the local residual spread.
