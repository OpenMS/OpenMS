# RESULT_BLEND: the network/GBM blend fixes Astral's own q but costs about 5 % of HeLa IDs (25 Nov 2026, FAIL)

- **Plan:** `PRECHECK_BLEND.md`, frozen and pushed in docs `02caf9aa13` before the container build or any remote run. The submission is in `SUBMITTED_BLEND.md`.
- **Code:** OpenMS `b3391c92e8`.
- **Container:** `opendia-blend-20261125.sif`, sha256 `8ba71d01…`.
- **Runs:** codon `blend-20261125/runs`, 28 jobs. All exited 0.
- **Evaluator:** `pdec/eval_pdec.py` (sha256 `bb03387d…`), unchanged. Full output is in `EVAL_BLEND.txt`, and wall and RSS are in `COST_BLEND.txt`.
- **Evaluation procedure:** `blend_eval.sh` ran the gate. Its single-stream reference fetch then broke mid-transfer (tar EOF).
  - `blend_eval2.sh` refetched each reference file with a size check and ran the same layouts and evaluator.
  - No rule, metric or reference choice changed.

**Gate: passed.** With no flag, the container is byte-identical to G0 on Astral raw01_s0, HeLa D s0 and HeLa H s0 (3/3 files each). So G0 is the existing default outputs, as the plan specifies.

## Numbers (own q ≤ 0.01; pooled bounds; median N @ FDP_c 1 % over 5 seeds)

| Set | G0 lower / combined | B5 lower / combined | G0 median N@FDP_c 1 % | B5 median N@FDP_c 1 % | B5 / G0 |
|---|---|---|---:|---:|---:|
| Astral raw 01 (dev) | 0.61 / 1.22 % | 0.52 / 1.03 % | 2,136 | 2,035 | 0.953 |
| Astral raw 02 (dev) | 0.56 / 1.11 % | 0.33 / 0.66 % | 2,146 | 2,231 | 1.040 |
| **Astral raw 03 (held out)** | 0.60 / 1.20 % | 0.37 / **0.75 %** | 2,246 | 2,288 | **1.019** |
| Astral, all 15 runs | 0.59 / 1.18 % | 0.40 / **0.81 %** | 2,170 | 2,231 | 1.028 |
| **HeLa D r ≈ 7 (held out)** | 1.11 / 1.27 % | **1.03** / 1.18 % | 21,396 | 20,288 | **0.948** |
| **HeLa H r ≈ 7 (held out)** | 1.29 / 1.47 % | **1.16** / 1.33 % | 37,342 | 35,483 | **0.950** |

- **Worst single-run combined bound:**
  - Astral: B5 1.34 % (raw01_s1) against G0 1.77 % (raw01_s0);
  - HeLa D: 1.21 % against 1.32 %;
  - HeLa H: 1.36 % against 1.50 %.
- **Every HeLa seed:** each B5 seed is below each G0 seed in N @ FDP_c 1 %. D is 19,941-20,736 against 21,016-21,468. H is 35,180-35,821 against 37,213-37,428.
- **Pair fairness P:** unchanged on Astral (0.477 against 0.471) and slightly lower on HeLa (D 0.539 against 0.544, H 0.526 against 0.540).

## Decision rule

| Criterion | Bar | Observed | |
|---|---|---|---|
| 1a raw 03 combined | ≤ 1.00 % | 0.75 % | pass |
| 1b raw 03 IDs | ≥ 0.98 × G0 | 1.019 | pass |
| 2a Astral 15-run combined | ≤ 1.00 % | 0.81 % | pass |
| 2b IDs per raw | ≥ 0.98 × G0 each | raw 01 0.953, raw 02 1.040, raw 03 1.019 | **fail** (raw 01) |
| 3a HeLa D IDs | ≥ 0.98 × G0 | 0.948 | **fail** |
| 3b HeLa D lower bound | ≤ G0 + 0.10 pp | −0.08 pp | pass |
| 3c HeLa H IDs | ≥ 0.98 × G0 | 0.950 | **fail** |
| 3d HeLa H lower bound | ≤ G0 + 0.10 pp | −0.13 pp | pass |

**Verdict: FAIL.** B5 stays opt-in, and the default is unchanged. Nothing is retuned: no other w and no other combination is tried on these runs.

## Reading

- **Astral:** the blend does what the diagnosis predicted, including on the held-out raw 03.
  - The pooled combined bound falls from 1.18 to 0.81 % (raw 03: 1.20 to 0.75 %).
  - IDs at true FDP are kept or slightly raised (15-run median 2,170 to 2,231).
  - The raw 01 median falls 4.7 %. Most of that is one seed, raw01_s1, at 994 (G0 2,270). In that run GBM round 2 kept no tree in one fold (`trees 0/212/180`), and the blend inherits a fold whose GBM score is constant. This is an observation and is not tested further here.
- **HeLa:** the blend tightens own q (lower bound D 1.11 to 1.03 %, H 1.29 to 1.16 %) but gives back most of the GBM's ID gain.
  - On D, B5 is +3.0 % over the network alone (19,704, PLAN_GBT CUR) against +8.6 % for the GBM.
  - On H, B5 is +1.8 % over the network alone (34,871) against +7.1 % for the GBM.
  - So on HeLa the blend trades about 5 % of IDs for about 0.1 pp of calibration. The user's goal (fix FDR without losing IDs) is not met there.
- **Interpretation:** on HeLa, much of the GBM's extra separation is real signal that the network misses. Averaging it with the network dilutes that signal along with the Astral excess. A fix that keeps HeLa's gain has to target the near-isomer class (DIAG_FDR Result 3) directly, not the whole GBM score.
- **HeLa calibration with G0 is itself not shown valid.** The H lower bound is 1.29 % > 1 %. That is a pre-existing limit of the default, unrelated to Astral's GBM issue, and B5 does not remove it (1.16 %).

## Next (needs a decision; not run)

1. **Keep B5 opt-in** for users who need conservative q on Astral-type data. Its trade-off is documented above.
2. **A targeted fix for the near-isomer class.** This would be a decoy or feature design that gives decoys the "co-isolated isobaric neighbour" signal that the entrapments borrow. For example, add decoys built from isobaric shuffles of library targets, or add a label-free feature "a higher-scored target shares ≥ 2 fragments at this apex" to both networks. Either would need its own precheck on a new held-out set. raw 03, D and H have now been used.
