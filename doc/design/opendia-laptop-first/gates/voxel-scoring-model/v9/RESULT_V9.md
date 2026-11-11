# Voxel v9 result: native-sampling fragment evidence, uniform-density chance model (D11FNWT). FAIL

**Plan.** `PLAN_V9.md`, sha256 `70b19f18…`, written 2026-11-11T20:35:04Z. The sha was checked OK by `chain_v9.sh` before training and again after the run.

## Result (`EVALUATION_V9.txt`, `NULL_V9.txt`, `VS_D11FWT_V9.txt`)

| Raw | Arm | N@FDP_c 0.5 / 1 / 1.5 / 2 % | q ≤ 0.01 N, FDP_c | AUC (entrapment > decoy) |
|---|---|---|---|---|
| 01 | D11FWT | 971 / 1496 / 1616 / 1716 | 1614, 0.014 | 0.4907 |
| 01 | D11FNWT | 1229 / 1467 / 1762 / 1820 | 1774, **0.018** | 0.4962 |
| 02 | D11FWT | 1127 / 1426 / 1625 / 1736 | 1536, 0.013 | 0.4958 |
| 02 | D11FNWT | 570 / 1485 / 1608 / **1726** | 1693, **0.019** | 0.4975 |

- **Rule 1.3 fails on both raws:** FDP_c at q ≤ 0.01 is 0.018 / 0.019.
- **Rule 2 fails on raw 02:** N@2 % is 1726 against 1736.
- **Rule 3 passes:** null-winner AUC 0.4948 / 0.4947 against 0.4933 / 0.4936.

**Reported only:**
- Paired bootstrap against D11FWT, 90 % intervals: @2 % −7..+211 on raw 01 and −137..+127 on raw 02.
- Funnel: D 866 → 910 (raw 01) and 929 → 892 (raw 02); E 1413 → 1407 and 1362 → 1429.

## Why it failed (diagnostics after the run; descriptive)

**The chance model was badly mis-calibrated** (`NAT_NULL_CHECK_raw01.txt`). We took 3000 decoy windows and scored scans 40 revisits (about 66 s) away from the apex, where the precursor's peak should be absent:

| Measure | Observed | Uniform-density model expects |
|---|---|---|
| Matched fragments per scan | 1.38 | 0.08 |
| P(L ≥ 2) | 0.41 | 0.01 |

- Mass-shifted fragment lists match only 0.37 per scan. So library fragment m/z values are much more often present than a random m/z at the same local density.
- Plausible reasons, not separated here: shared low-mass ions (y1/y2, b2), and co-isolated peptides in the same 2 Th window.

**Consequence.** `nat_apex_logp` ≥ 2 held for 80 % of all decoy windows. The features mostly re-encoded the match count. Score-matched AUC for D against decoys was 0.59, 0.40, 0.33 and 0.46 for the four L-based features (`NAT_SEPARATION_raw01.txt`).

**Corrected null (diagnostic, raw 01 only; `native_features_emp.py`, `NATE_SEPARATION_raw01.txt`).** For each fragment, the chance rate is measured empirically: the share of the same isolation window's scans 8–60 s from the apex that hold a peak within tolerance.

| Measure | Uniform null | Empirical null |
|---|---|---|
| Decoy windows with apex L ≥ 2 | 0.80 | **0.033** |
| Decoy windows with apex L ≥ 3 | 0.50 | **0.005** |
| Score-matched AUC, D vs decoys (`apex_logp`) | 0.59 | **0.84** |
| Score-matched AUC, entrapment vs decoys | 0.48 | 0.55 |
| Score-matched AUC, all human vs decoys | 0.56 | 0.59 |

The empirical null is close to calibrated: it expects about 0.03 for the max of 3 scans at L ≥ 2.

## Outcome

- **Fail.** D11FWT stays the base.
- The empirical-null version is tested under a new frozen plan (PLAN_V10). It changes only the chance rate; everything else is identical.
