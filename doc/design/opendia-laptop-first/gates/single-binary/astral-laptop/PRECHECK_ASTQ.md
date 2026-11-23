# PRECHECK_ASTQ: which stage made OpenDIA's own q ≤ 0.01 loose on Astral? (23 Nov 2026, diagnostic)

**Why.** RESULT_ASTRAL_LAPTOP: at its own q ≤ 0.01, the current default has entrapment FDP_c above 1 % on 10 of 15 Astral runs (medians 1.17 / 0.95 / 1.32 % on raws 01 / 02 / 03, max 1.77 %). The 16 Nov code (`ab46f97ad4`, the MS1 CUR arm) was never above 1 % (0.12-0.95 %). Since then seven HeLa stages were adopted in this order:
1. evidence block;
2. 1/K0 scale 0.6;
3. 10 ppm trace tensor;
4. matched filter;
5. RT refinement;
6. sibling rescoring;
7. GBM rescoring.

On HeLa H, the GBM already loosened the lower bound (1.15 to 1.29 %, RESULT_GBT_H). This precheck finds the stage (or stages) responsible before any change is designed. It is a diagnostic. It changes no default and has no pass rule for adoption.

## Setup (fixed before any result)

- **Container:** `opendia-pd10-gbmdef.sif` (OpenMS `bcb156e036`, sha256 `7dd8a1d6…`). It is the GBM default before the cleanup and still has the stage switches.
- **Inputs:** the same as RESULT_ASTRAL_LAPTOP: arm-E search library, `library12_full.pqp` scoring library, self-calibrating, 16 threads, 32 GB.
- **Runs:** design raws 01 and 02 only, seeds 0-4. Raw 03 is not used, because it is held back for design.
- **Arms:** a cumulative ladder that removes the adopted stages newest first (`astq.sbatch`):

  | Arm | Flags | Stages present |
  |---|---|---|
  | L0 | current default: RESULT_ASTRAL_LAPTOP batch 2 outputs, OpenMS `88e5da3023` | all |
  | L1 | `-no_gbm_rescore` | no GBM |
  | L2 | + `-no_sibling_features` | no GBM, no sibling |
  | L3 | + `-no_rt_refine` | ... no RT refinement |
  | L4 | + `-no_matched_filter` | ... no matched filter |
  | L5 | + `-trace_tolerance_ppm 0 -im_tolerance_scale 1 -no_extra_scores` | close to the 16 Nov scorer |

- **Anchor A0:** `gbmdef` with no flags on raw01_s0 and raw02_s2. Both runs completed in batch 1, without the GBM guard path. They must be byte-identical to L0, so that L0 and the `gbmdef` ladder are the same code.
- **Evaluator:** `gates/single-binary/pdec/eval_pdec.py` (sha256 `bb03387d…`), unchanged, run per raw (`P = raw01`, `raw02`). The library gives r = 151,377 / 151,102.
- **Measures (per arm, pooled over the 10 runs by summing E and N from `eval_pdec.py`'s per-run rows of both raws):**
  - lower bound E/N and combined bound E(1+1/r)/N at own q ≤ 0.01, each with a 90 % Poisson interval;
  - median N @ FDP_c 1 %;
  - pair fairness P (entrapment target beats its own decoy, best member q ≤ 0.05).
  The same is computed on L0 from the existing outputs.

## Attribution rule (written before any result)

- **Anchor gate:** if A0 is not byte-identical to L0 on both runs, stop. The ladder then measures different code. Report this and do not attribute.
- **Step size:** the step for a stage is the drop in the pooled combined bound from L(k-1) to Lk.
- **Named cause:** the stage with the largest step is named the cause if all of these hold:
  - the step is at least 0.20 pp;
  - the step is at least half of the total drop from L0 to L5;
  - after the step, the pooled combined bound is at most 1.10 %.
- **Two contributors:** if no single stage meets the rule, but two adjacent steps together meet it, both are named as joint contributors.
- **Not attributable:** if neither holds, the looseness is spread over stages. The next step would then target the final q-value step, not a stage.
- **Fairness supports:** if pair fairness rises by at least 0.02 at the same step, this is recorded as supporting evidence (decoys easier than false targets after that stage).
- **Report only:** the ID cost of each step (median N @ FDP_c 1 %) is reported, not ruled on. No default changes from this precheck. A fix would get its own precheck and a frozen plan.

## Cost

There are 52 jobs at 16 threads / 32 GB: 5 arms × 10 runs, plus 2 anchors. Each takes about 7-15 min.
