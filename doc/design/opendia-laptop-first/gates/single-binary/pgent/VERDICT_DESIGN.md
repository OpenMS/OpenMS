# PLAN_PGENT design verdict: FAIL on criterion 1 (absolute FDP_c), pass on criterion 2 and the relative calibration

Written 18 November 2026, after the D design runs and before any H run. Evaluation: `eval_pgent.py`, `EVAL_D.txt`.

- **Library:** `hela_entrap.pqp` (sha256 9c4e8c7d…). 366,357 target precursors, r = 1.1225.
- **Container:** `opendia-pg-a2aa8fc56c.sif` (sha256 4455ace8…).

| D, seeds 0-4 | N at q ≤ 0.01, median | entrapments at q ≤ 0.01, pooled | FDP_c at q ≤ 0.01, pooled | N@FDP_c 1 %, median (range) |
|---|---|---|---|---|
| CUR | 10,946 | 511 / 54,780 | 0.0176 | 9,701 (9,596-9,944) |
| NEW (`-pg_scores`) | 12,922 | 565 / 64,346 | **0.0166** | **12,033** (11,934-12,332) |

1. **Criterion 1, FDP_c:**
   - (a) NEW pooled FDP_c ≤ 0.012: **fail**, 0.0166.
   - (b) NEW ≤ CUR + 0.003: pass, 0.0166 against 0.0176. NEW is better calibrated than CUR.
2. **Criterion 2:** NEW median N@FDP_c 1 % ≥ 1.05 × CUR: **pass**, +24.0 %. Every NEW seed is above every CUR seed.
3. **Criterion 3 (restated):** Astral pooled 0.0073: pass.
4. **Criterion 4 (restated):** Astral +11.8 / +2.4 / +6.2 %: pass.

**Verdict by the plan: FAIL** (criterion 1a). `-pg_scores` stays opt-in, and nothing is retuned.

**What the result shows** (a reading, not a change of verdict):
- On this predicted entrapment library, both arms are above 1 % true FDR at q ≤ 0.01: CUR 1.76 %, NEW 1.66 %.
- The failure of 1a is today's default too. It is not caused by the new scores, and the scores make it slightly better.
- At equal true FDR, the new scores add 24 % more IDs on diaPASEF.
- **Why both arms miss nominal:** this is open, and it is likely the same parent-interference effect seen on Astral
  (RESULT_PAIRLIB: rho 1.3-1.6).
  - Every HeLa entrapment shares length, cleavage sites and composition with a parent that is in the library and
    usually elutes.
  - So entrapments get real, partly matching signal that decoys (pseudo-reverse) do not.
  - That makes FDP_c with these entrapments an upper estimate of the false rate for unrelated false targets.
  - This is a property of the benchmark construction. It was not foreseen in the plan, and it does not change the
    verdict.

H is run next and reported either way.
