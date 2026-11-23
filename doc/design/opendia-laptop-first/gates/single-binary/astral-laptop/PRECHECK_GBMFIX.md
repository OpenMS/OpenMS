# PRECHECK_GBMFIX: a GBM variant that keeps the HeLa gain and restores Astral's own q (23 Nov 2026, diagnostic)

**Why.** RESULT_ASTQ names the hard-negative GBM rescore as the cause of the loose own q ≤ 0.01 on Astral. The pooled combined bound is 1.17 % with the GBM and 0.73 % without it (raws 01 + 02). On Astral the GBM gains +9.9 % at FDP_c 1 % on raw 01, but loses 6.1 % on raw 02 and causes the low seed raw02_s4. On HeLa it is worth +8.6 % (D) and +7.1 % (H) (PLAN_GBT), but it also moved H's combined bound from 1.31 to 1.47 %.

The Astral logs point at the later rounds:
- Round 2 or 3 often has far fewer positives (about 600-750 against about 1,900 in round 1).
- A fold then keeps no tree, and the round is discarded or crashes.

On HeLa, all three rounds are healthy. The precheck tests two minimal variants of the existing GBM. There is no new model and no new feature.

## Arms (new opt-in options; defaults unchanged)

OpenMS adds `-gbm_rounds <n>` (default 3) and `-gbm_positive_q <q>` (default 0.01). With the defaults, the output must be byte-identical to `88e5da3023`. This is checked first (gate below).

| Arm | Flags | Idea |
|---|---|---|
| G0 | none (current default) | reference |
| G1 | `-gbm_rounds 1` | one relabel round: no self-reinforcement on its own positives |
| GQ | `-gbm_positive_q 0.005` | cleaner positives: fewer false targets in the training labels |
| NG | network + sibling only, no GBM: the RESULT_ASTQ L1 outputs (container `gbmdef`, `-no_gbm_rescore`) | the no-GBM floor (Astral only) |

## Runs (fixed before any result)

- **Astral:** raws 01 and 02 (design runs only; raw 03 stays held back), seeds 0-4, 16 threads, 32 GB. The inputs are as in RESULT_ASTRAL_LAPTOP. G0 is the existing batch 2 outputs, and NG is the existing L1 outputs.
- **HeLa r ≈ 7 D** (2737.d, `hela_r5.predicted_decoys.pqp`): seeds 0 and 1, 8 threads, 160 GB, as PLAN_GBT. G0 is the existing `entlib-gbt-20261122/runs/D_gbt_s{0,1}`.
- **HeLa r ≈ 7 H** is not run in this precheck. It is run only in a later frozen plan.
- **Evaluator:** `pdec/eval_pdec.py` (sha256 `bb03387d…`), unchanged.
  - Astral: r = 1.0018, pooled over the 10 runs as in PRECHECK_ASTQ.
  - HeLa D: r = 7.0539.

## Gate

The new container with no flags must be byte-identical to G0 on:
- Astral raw01_s0 and raw02_s2;
- HeLa D r ≈ 1.56 s0 (the D156 byte-check).

If not, stop.

## Decision rule (written before any result)

A variant **passes** and goes to a frozen 5-seed plan (D, then H, then Astral), if all of these hold:
1. **Astral calibration:** the pooled combined bound at own q ≤ 0.01 is ≤ 1.00 %, and no single run is above 1.50 %.
2. **Astral IDs:** median N @ FDP_c 1 % is ≥ 0.97 × G0 on each of raw 01 and raw 02.
3. **HeLa D IDs:** mean N @ FDP_c 1 % over the 2 seeds is ≥ 0.98 × G0. This keeps at least about 75 % of the GBM's +8.6 %.
4. **HeLa D calibration:** the combined bound over the 2 seeds is ≤ G0 + 0.10 pp.

- **If both variants pass:** the one with the lower Astral pooled combined bound goes to the plan. A tie within 0.05 pp goes to the higher HeLa D mean N.
- **If neither passes:** both are recorded as fail, and the next step is a different fix (for example the GBM's ranking with the network's calibrated tail). That would need its own precheck.
- **Reported, not ruled on:** the pair fairness P, wall time and peak RSS.

## Cost

There are 22 jobs:
- Astral: 2 arms × 10 runs = 20 jobs, plus 2 gate runs.
- HeLa D: 2 arms × 2 seeds = 4 jobs, plus 1 D156 gate run.

Astral jobs take about 7-15 min each, and HeLa D jobs about 20-35 min each.
