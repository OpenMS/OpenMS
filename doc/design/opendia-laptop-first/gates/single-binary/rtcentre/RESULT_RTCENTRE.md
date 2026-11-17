# PLAN_RTCENTRE result: PASS (kept, OpenMS 6dd893b954)

Plan: `PLAN_RTCENTRE.md` (sha256 checked). Design verdict was committed before the held-out run (`VERDICT_DESIGN.md`).
Counts are each tool's own q ≤ 0.01 (no entrapment on diaPASEF). DIA-NN 2.2 on the same library is used only for the
agreement and map checks.

| | BASE 0e9e4fc3a5 | NEW 6dd893b954 | DIA-NN 2.2 |
|---|---|---|---|
| D (HeLa 20200505 2737.d), N median (range) | 13,127 (13,109-13,290) | **14,875** (14,804-14,927), +13.3 % | 28,950 |
| D agreement with DIA-NN | 0.931 | 0.934 | |
| D map abs err p95 | 19.3 s | 4.7 s | |
| H (HeLa 20201207 22474.d, held out), N median (range) | 19,725 (19,671-19,775) | **21,682** (21,600-22,079), +9.9 % | 38,074 |
| H agreement with DIA-NN | 0.9535 | 0.9530 | |
| H map abs err p95 (early median err) | 92.8 s (−101 s) | 23.2 s (0.0 s) | |
| Astral raw01 | | byte-identical (2-pass calibration, unchanged path) | |

All five rules pass. On H, agreement dropped by 0.0005, inside the 0.005 margin.

- On H the old one-pass map was about 100 s late at the start of a 21-min gradient. The final map extrapolated over
  a stretch with no anchors.
- The share of DIA-NN's IDs now reached is 51 % on D (was 45 %) and 57 % on H (was 52 %).
- The remaining gap is the diaPASEF scoring (DIAG_HELA.md, cause 3).
