# Raw 02 native byion replication (design run)

`PLAN.md` was written before the run. Binary is commit `8c49a5bbb2`. Arms A (search+boundary) and B (+byion) share byte-identical crossfit columns 1–6 and 9 (`CMP.txt`).

| | FDP_c ≤ 0.5% | **1%** | 1.5% | 2% | plain q≤0.01 | plain q≤0.05 | bulk-null q≤0.01 | bulk-null q≤0.05 |
|---|---|---|---|---|---|---|---|---|
| A search+boundary | 850 | 1030 | 1115 | 1186 | 860 / 0.005 | 1607 / 0.073 | | |
| **B +byion** | 1230 | **1298 (+26.0%)** | 1331 | 1372 | 1236 / 0.007 | 1788 / 0.074 | 1115–1125 / 0.002 | 1658–1677 / **0.062–0.063** |
| OSW+Percolator (same metric, OSWLIKE.md) | | 1287 | | | | | | |

- **N criterion: PASS.** B is +26% over A, and at FDP_c 1% it is at or above OSW+Percolator.
- **Calibration criterion: FAIL, by a narrow margin.** Bulk-null q ≤ 0.05 has FDP_c 0.062–0.063 for all 3 seeds against the 0.06 bar. q ≤ 0.01 is conservative (0.002).
- Plain pooled q is anti-conservative at q ≤ 0.05 on both arms (0.073–0.074).
- Fold weights: byion_b +2.6 to +2.7 and byion_y +3.6 to +3.7, as on raw 01.
- Cost: +1:21 wall time (14:17 against 12:56); peak RSS 24.7 GB against 24.6 GB.

Reading:
- The byion counts are a robust ID gain on both design runs: +36% on raw 01 and +26% on raw 02.
- With the native stage-2 score, bulk-null weighting fixes q ≤ 0.01 on both runs (0.006 / 0.002) but only partly fixes q ≤ 0.05 on raw 02. It brings 0.074 down to 0.062, where the offline harness reached 0.040–0.055.
- Native q calibration at 5% is still open.

## h13: MLP head on the native features (`mlp_head_native.py`, `mlp_head01.txt`, `mlp_head02.txt`): FAIL
Entrapment-matched N at FDP_c 1%, then bulk-null q≤0.05 FDP_c:

| | raw 01 | raw 02 |
|---|---|---|
| tool LDA | 1029 / 0.053 | **1298** / 0.060 |
| MLP 5-seed ensemble | 1115 / 0.060 | 1220 / 0.053 |
| MLP seed 0 | 1169 / 0.046 | 1185 / 0.035 |

- The rule (N ≥ LDA on both runs) fails on raw 02. On raw 01 the MLP gains at FDP_c 0.5–1% but loses at 1.5–2% (1222 against 1321).
- The native LDA stays the head.
- **Sensitivity of the 5% calibration on raw 02.** The same LDA score passed through this script gives bulk-null q≤0.05 FDP_c 0.060. `bulknull_native.py` gives 0.062–0.063; it differs by per-fold decoy standardisation and by classifier seed. Raw 02 5% calibration therefore sits at the bar: it is neither a clear pass nor a clear fail.
- q ≤ 0.01 is conservative under every variant (0.002–0.009).

## h15: density-ratio classifier for a native bulk-null q (`ratio_classifier.py`, `ratio_classifier.txt`)

FDP_c at q ≤ 0.01 / 0.05 on the native byion score:

| classifier | raw 01 | raw 02 |
|---|---|---|
| none (plain TDC) | 0.013 / 0.065 | 0.007 / 0.074 |
| logistic | 0.014 / 0.077 | 0.011 / 0.092 |
| boosted stumps | 0.011 / 0.056 | 0.007 / 0.068 |
| boosted depth-2 | 0.010 / 0.054 | 0.007 / 0.068 |
| HGB, 31 leaves | 0.006 / 0.050 | 0.002 / 0.062 |

- **Decision under the predeclared rule:** only HGB qualifies, and it still sits above 0.06 at 5% on raw 02.
- A native port would need a full histogram gradient-boosting learner of about 31 leaves and 200 rounds. That is a large new component, for a correction that is borderline at 5% and unnecessary at 1%.
- Plain pooled TDC already meets the 1% gate rule (≤ 0.015) on both raws.
- **So no native bulk-null for now:**
  - Report plain q at 1% as the primary q.
  - Do not report a q ≤ 0.05 list as calibrated.
  - Keep bulk-null as an offline diagnostic.
- Revisit if a decoy-side fix narrows the 5% tail mismatch (CALIB).
