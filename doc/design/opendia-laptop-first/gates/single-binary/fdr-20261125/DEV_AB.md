# DEV_AB: evidence sharing (A) and GBM restrictions (B) on Astral dev runs

Status: dev screen only (Astral raw01/02 x seeds 0-4, the 10 diagnostic runs). No held-out data was used. Nothing here is a gate result.

## What was built (OpenMS fork `4c6ca367fd`, all opt-in, default off)

| Option | What it does |
|---|---|
| `-evidence_share` (A) | After the first scorer pass, the interference resolver runs again. This time a window may be explained only by windows with a higher first-pass score that are at q <= 0.01 (any label, same rule for targets and decoys). The four results are added as `e_explained`, `e_libcos_resid`, `e_dlibcos` and `e_nclaim` before the sibling rescoring pass. |
| `-gbm_exclude <names>` (B) | The hard-negative GBM never splits on the listed columns. |
| `-gbm_monotone <name:+1/-1>` (B) | Monotone constraints for the GBM. Leaf weights are bounded and the split gain uses the bounded weights (the basic LightGBM method). |

Checks:
- `OpenDIAGBM_test` and `OpenDIAInterference_test` pass, including new tests for exclusion, monotone +1/-1 sweeps, the claimant filter and bad sizes.
- Flag-off run on raw01_s0 is byte-identical to the diagnostic run. All three outputs match (`devAB/identity.txt`).
- The offline replay (`devAB/gbm_replay.cpp`, linked to the tool's own `gbmRescore`) reproduces the exported score on all 1,675,700 of 1,675,700 rows of raw01_s0. B and AB were screened with this replay.

## Results (pooled over 10 runs; corrected `_p_target` labels via instr)

N is targets at the tool's q <= 0.01. "Combined" is the 2E/N entrapment bound.

| Arm | N | E | lower | combined | N@FDP_c 1 % | vs G0 |
|---|---:|---:|---:|---:|---:|---:|
| G0 (default) | 21,487 | 121 | 0.56 % | 1.13 % | 20,871 | 1.000x |
| B1 exclude the 4 `i_*` features | 20,812 | 118 | 0.57 % | 1.13 % | 20,400 | 0.977x |
| B2 monotone `i_explained`, `i_nclaim` -1 | 21,617 | 114 | 0.53 % | 1.05 % | 20,595 | 0.987x |
| B3 = B2 plus exclude `i_libcos_resid`, `i_dlibcos` | 20,820 | 142 | 0.68 % | 1.36 % | 18,514 | 0.887x |
| B4 = B1 plus exclude the `f12_extra_*` features | 20,653 | 123 | 0.60 % | 1.19 % | 19,261 | 0.923x |
| **A** `-evidence_share` | 22,140 | 124 | 0.56 % | 1.12 % | **21,603** | **1.035x** |
| AB2 = A plus monotone -1 on `i_explained`, `i_nclaim`, `e_explained`, `e_nclaim` (replay) | 22,003 | 122 | 0.55 % | 1.11 % | 21,614 | 1.036x |

Per run (N@FDP_c 1 %, files `devAB/EVAL_*.txt`):
- A beats G0 on 6 of 10 runs.
- Without the weakest G0 run, raw02_s4 (1,728), A is still 1.020x.
- Seed-to-seed spread is large. Single-run G0 ranges from 1,728 to 2,270.

## Reading

- **B does not work.** No restriction lowers the bound without losing IDs. B2 comes closest (1.05 %) but loses 1.3 %. Taking borrowing features away from the GBM does not stop it rewarding borrowed evidence, because the same signal is still carried by the other features.
- **A gains IDs but does not fix FDR control.** N@FDP_c 1 % rises 3.5 % and N rises 3.0 %, but the combined bound stays at 1.12 %, unchanged from G0.
- **AB matches A.** Adding monotone constraints on top of A changes nothing material.
- **No arm meets the goal.** None reaches a combined bound <= 1 % at no ID cost on dev. A is the only arm that meets the "no ID loss" condition.

## Cost (descriptive only)

- A adds one resolver pass and keeps the resolver inputs through the first scorer pass.
- Peak RSS 13.3-14.4 GB, against about 13.3 GB for G0 on the same runs.
- Wall times of 24-79 min were set by an unrelated 13-core job on the same machine and are not comparable (`devAB/COST_A.txt`).

## Not done

- No held-out evaluation.
- No HeLa. A on diaPASEF needs new runs and a new held-out set.
- No change to defaults.
