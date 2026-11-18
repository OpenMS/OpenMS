# PLAN_PGENT result: FAIL (criterion 1a); `-pg_scores` stays opt-in

Written 18 November 2026. The design verdict (`VERDICT_DESIGN.md`, fe183d4a95, 12:09:37) was committed before any H run.
The H jobs were submitted at 12:09 and resubmitted at 12:17 with 96 GB, after all ten were killed for memory at 32 GB.
Only the memory changed.

| run | arm | N at q ≤ 0.01, median | FDP_c at q ≤ 0.01, pooled over 5 seeds | N@FDP_c 1 %, median (range) |
|---|---|---|---|---|
| D (design) | CUR | 10,946 | 0.0176 | 9,701 (9,596-9,944) |
| D (design) | NEW | 12,922 | 0.0166 | **12,033** (11,934-12,332), +24.0 % |
| H (held out) | CUR | 24,365 | 0.0167 | 23,136 (23,029-23,212) |
| H (held out) | NEW | 28,421 | 0.0186 | **26,924** (26,870-27,069), +16.4 % |

- **Library:** the fully predicted HeLa entrapment library (r = 1.12).
- **Not comparable** with the experimental-library numbers (16,116 / 24,000).

## Verdict
- **FAIL on criterion 1a.** NEW's FDP_c at q ≤ 0.01 is 0.0166 on D, above 0.012. On H it is 0.0186.
- Criterion 1b (NEW ≤ CUR + 0.003) passes on D (−0.0010) and would also hold on H (+0.0019).
- Criterion 2 passes: +24 % at equal true FDR on D, +16 % on H.
- `-pg_scores` stays opt-in. Nothing is retuned.

## Reading
1. **The gain at equal true FDR is real and holds on the held-out run.**
   - The new scores give 16-24 % more IDs on diaPASEF, measured with entrapment.
   - On Astral the gain is +2 to +12 %.
   - Every NEW seed beats every CUR seed on D and on H.
2. **Neither arm is calibrated on this library.**
   - At q ≤ 0.01, CUR's true FDR is 1.67-1.76 % and NEW's is 1.66-1.86 %.
   - Criterion 1a would have failed for today's default too. It measures the pipeline on this predicted library, not
     the new scores.
3. **NEW's calibration is no worse than CUR on D and slightly worse on H** (+0.19 percentage points). Seeds of one run
   are not independent (Q7), so this difference is within one run's noise on both.
4. **Cause of the shared over-shoot:** open.
   - Parent interference is a minor share: only 15-25 % of accepted entrapments sit at their parent's apex.
   - The candidates are a predicted-library effect on this short gradient, the pseudo-reverse decoy model, or
     self-training drift.
   - Settling it is a separate question about OpenDIA's q calibration on diaPASEF, not about `-pg_scores`.

## Consequence
The evidence now says `-pg_scores` gives more IDs at equal true FDR on 3 Astral raws and 2 diaPASEF runs. The blocking
issue is the pipeline's absolute calibration on diaPASEF, which CUR shares. Two options for the next step:

- **(a)** Freeze a plan with **relative** calibration criteria: NEW FDP_c ≤ CUR + δ, and N@FDP_c 1 % ≥ CUR.
  - This would have passed here.
  - It needs a new plan, and new HeLa runs as the held-out set.
- **(b)** First fix the diaPASEF calibration itself, for example the decoy method, as its own plan. Then re-test
  `-pg_scores` against the absolute bound.

## Addendum: the decoys under-count false targets (`RHO.txt`, seed 0)
rho = (entrapments × (1 + 1/r)) / decoy winners, per run-q band. rho = 1 means the decoys model false targets.

| run, arm | q 0-0.005 | q 0.005-0.01 | q 0.01-0.02 | q 0.02-0.05 | q 0.05-0.2 |
|---|---|---|---|---|---|
| D CUR | 2.20 | 1.04 | 1.41 | 1.30 | 1.19 |
| D NEW | 1.39 | 1.58 | 1.79 | 1.47 | 1.17 |
| H CUR | 1.80 | 1.31 | 1.29 | 1.12 | 1.12 |
| H NEW | 1.69 | 1.48 | 1.49 | 1.19 | 1.10 |

- A shuffled-sequence false target outscores a pseudo-reverse decoy by a factor of about 1.3-1.8 in the accepted range.
- This is the same "rho" excess as on Astral (RESULT_PAIRLIB: 1.3-1.6), where the decoys are made the same way.
- The NEW scores raise rho a little near the cut on D (1.0-1.6 → 1.4-1.8) and are similar on H. That is consistent with
  better separation, whose benefit goes slightly more to shuffled targets than to decoys.
- The AUC(entrapment > decoy) over q ≤ 0.05 is only 0.52-0.54, so the difference is small per ID but systematic.
- It points at option (b): the decoy model, which is shared by the default. Options for (b):
  - shuffle decoys matched to how the entrapments are made;
  - or a decoy-free q correction estimated from the run.
  - Either needs its own precheck and plan.
