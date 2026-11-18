# PLAN_PGSCORES result: FAIL (D criterion 2, Astral raw 02 criterion 3); code removed

Written 18 November 2026. The design verdict (`VERDICT_DESIGN.md`, commit 0bea8cd75f, 09:13:44) was committed before any
held-out run. The Astral design runs had finished by 09:09 but were evaluated only after it; they fail a criterion too.

- **Code:** 5a7c95fd68, flag `-pg_scores`, 10 features. Reverted.
- **Container:** `opendia-pg-5a7c95fd68.sif` (sha256 314ed0bc…).

## Design
**HeLa D:**

| D | N median (min-max) | agreement with DIA-NN |
|---|---|---|
| CUR | 16,116 (15,957-16,177) | 0.9299 |
| NEW | 18,792 (18,585-18,885) | 0.9104 |

- Criterion 1 passes: +16.6 %, and every NEW seed is above every CUR seed.
- **Criterion 2 fails:** 0.9104 < 0.9249.

**Astral** (entrapment, `EVAL_astral_design.txt`; CUR = PLAN_MS1 CUR runs):

| raw | CUR median N@FDP_c 1 % | NEW median N@FDP_c 1 % | NEW FDP_c at q ≤ 0.01, median / max |
|---|---|---|---|
| 01 | 1,804 | 2,017 (+11.8 %) | 0.0054 / 0.0073 (pass) |
| 02 | 1,909 | 1,954 (+2.4 %) | **0.0103 / 0.0129 (fail)** |

- **Gate 0** passed: on the labelled common set the ceiling went 14,988 → 17,304, +2,316.

## Held out (evaluated after the verdict)
| run | CUR | NEW |
|---|---|---|
| HeLa H, N median (agreement) | 24,000 (0.9524) | 28,020 (0.9390), +16.8 % |
| Astral raw 03, median N@FDP_c 1 % | 1,907 | 2,025, +6.2 % |
| Astral raw 03, FDP_c at q ≤ 0.01, median | 0.0076 | 0.0063 |

Against DIA-NN reading .raw (2,204 / 2,356 / 2,259), NEW would reach 92 / 83 / 90 % on raws 01/02/03. Today's code
reaches 81-84 %.

## Reading
- **The features carry real evidence.** At equal true FDR (entrapment N@FDP_c 1 %), NEW is higher on all three Astral
  raws and on every raw 01 and raw 03 seed. It is the largest gain of any block since block I.
- **They also make the q-values less reliable.**
  - On raw 02, two of five seeds exceed 1 % true FDP at q ≤ 0.01 (0.0125 and 0.0129).
  - On D, the added IDs agree with DIA-NN at a marginal rate of 0.79.
- So some of the new features separate decoys from targets more than they separate entrapments from targets.
- Which features do this has not been measured. Measuring it is the obvious precheck for any follow-up plan. Any such
  plan would be a new plan with its own freeze, not a retune of this one.
- **Cost:** extraction 56 → 117 s on D. Astral raw 01 extraction was 95 s on codon, peak memory unchanged.

## Outcome
Fail, as the plan states. The code is removed with `git revert`, and nothing is retuned.
