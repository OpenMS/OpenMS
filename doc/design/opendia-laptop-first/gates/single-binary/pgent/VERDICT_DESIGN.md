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
- **Why both arms miss nominal:** this is open. The parent-interference explanation was checked and is **not** the main
  cause (`parent_check.py`, D seed 0):
  - Of the accepted entrapments, 15 % (CUR) and 25 % (NEW) have their parent's apex within 5 s.
  - 17 % and 20 % have their parent also accepted.
  - So most accepted entrapments are not riding on their parent's signal.
- The remaining candidates, none tested yet:
  1. **Library mismatch.** A fully predicted library on an Evosep 5-min run, where the RT and IM prediction errors are
     large, may favour random matches. Targets and entrapments get the same prediction quality, but decoys
     (pseudo-reverse) inherit it as well.
  2. **Decoy model.** Pseudo-reverse decoys may be easier to reject than shuffled targets on this data.
  3. **The scorer's self-training.** On a library where far fewer targets are present, the scorer may drift.
- CUR shows the same excess, so this is a property of today's pipeline on this library, not of the new scores.

H is run next and reported either way.
