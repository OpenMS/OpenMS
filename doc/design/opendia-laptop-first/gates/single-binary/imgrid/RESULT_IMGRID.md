# PLAN_IMGRID result: PASS (kept, OpenMS 963543c7f0)

Plan: `PLAN_IMGRID.md` (sha256 checked). Design verdict committed before the held-out run (`VERDICT_DESIGN.md`).
Counts are own q ≤ 0.01. DIA-NN 2.2 (same library) is evaluation-only.

| | BASE 6dd893b954 (RTCENTRE) | NEW 963543c7f0 | DIA-NN 2.2 |
|---|---|---|---|
| D HeLa 2737.d, N median (range) | 14,875 (14,804-14,927) | **16,116** (15,957-16,177), +8.3 % | 28,950 |
| D agreement with DIA-NN | 0.9344 | 0.9299 | |
| H HeLa 22474.d (held out), N median (range) | 21,682 (21,600-22,079) | **24,000** (23,889-24,121), +10.7 % | 38,074 |
| H agreement with DIA-NN | 0.9530 | 0.9524 | |
| Astral raw01 | | byte-identical | |

All rules pass. The design-run agreement drop (−0.0045) sits near the 0.005 limit. On the held-out run it is −0.0006.

Cumulative today on HeLa (as benchmarked this morning → now):
- D: 11,653 → 16,116 (+38 %). That is 40 % → 56 % of DIA-NN, and 61 % → 84 % of OpenSWATH SVM (19,132).
- H: 19,725 (with the 1/K0 fix) → 24,000. That is 52 % → 63 % of DIA-NN.
