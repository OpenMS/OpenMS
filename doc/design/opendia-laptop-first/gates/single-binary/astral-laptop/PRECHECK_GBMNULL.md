# PRECHECK_GBMNULL: can a run tell, from its own targets and decoys, whether the GBM rescore is safe? (frozen 23 Nov 2026)

**Question.** The GBM rescore adds IDs on HeLa but loosens own q on Astral (RESULT_ASTQ, RESULT_GBMFIX). Can a rule that reads only one run's targets and decoys tell these cases apart? The rule may not use entrapment labels, other datasets or tuned thresholds. Entrapment is used only afterwards, to grade the rule.

**Mechanism tested.** TDC q is valid only while decoys stand in for false targets. The GBM is trained to separate q ≤ 0.01 targets from decoys. If it also learns features that separate targets from decoys as classes (rather than signal from noise), false targets move above their decoys and q becomes optimistic. In the null region, pairs where both members are almost surely false, an exchangeable scorer has a target win probability of about 1/2, and the GBM should not move it. A shift there is visible without labels.

## Statistic (computed per run, `gbmnull.py`, reads no entrapment labels)

Inputs: the run's network score (before the GBM) and GBM score per precursor (the kept, best window), and the library target-decoy pairing (TRAML_ID, `DECOY_` prefix).
- Both scores are the tool's calibrated score `-log S` against that fold's decoy winners, so `ln 2` is the decoy-winner median in every fold.
- **Null region:** pairs where both target and decoy are scored in both runs, and both members' *network* score is ≤ ln 2 (both at or below the decoy median before the GBM). The region is defined by the score the GBM is judged against, not by the GBM.
- **Flips:** b = pairs the decoy wins under the network but the target wins under the GBM; c = the reverse.
- **Test:** exact one-sided binomial (McNemar), p = P(X ≥ b | n = b + c, 1/2).
- **Rule:** skip the GBM in this run if p < 0.01; otherwise use it. The level is fixed in advance and is not tuned.
- Reported but not used by the rule: b, c, the target win fraction in the null region under each score, and the region size.

The run already holds both scores, so this can later run inside OpenDIA. Here it is computed offline from seed-matched pairs of runs that differ only by the GBM.

## Data

| Set | Network (no GBM) | GBM | Status |
|---|---|---|---|
| HeLa D r7, seeds 0-4 | `entlib-gbt-20261122/runs/D_cur_s*` | `D_gbt_s*` (same container) | outcome known (RESULT_GBT_D) |
| HeLa H r7, seeds 0-4 | `H_cur_s*` | `H_gbt_s*` | outcome known (RESULT_GBT_H) |
| Astral raw 01, seeds 0-4 | `astq-20261123/runs/raw01_L1_s*` | `astral-laptop-20261123/runs/raw01_s*` | outcome known (RESULT_ASTQ) |
| Astral raw 02, seeds 0-4 | `raw02_L1_s*` | `raw02_s*` | outcome known |
| **Astral raw 03, seeds 0-4** | **new:** `gbmnull-20261123/runs/raw03_L1_s*` (`astq.sbatch`, container `gbmdef`, `-no_gbm_rescore`) | `astral-laptop-20261123/runs/raw03_s*` | **held out**: no no-GBM run exists yet |

- **Anchor gate for raw 03:** A0 (`gbmdef`, no flags) raw03_s0 must be byte-identical to `astral-laptop-20261123/runs/raw03_s0`. Otherwise raw 03 is not evaluated.
- **Honesty note:** I have seen the entrapment outcomes on D, H, raw 01 and raw 02, but not this statistic. Those four sets check whether the rule separates known cases. Only raw 03 is out of sample.

## Grading (entrapment, `eval_pdec.py` sha256 `bb03387d…` unchanged, r per library)

The grading layout per set (D, H, A1, A2, A3) is `<set>_never_s<k>`, `<set>_always_s<k>` and `<set>_rule_s<k>`, all symlinks to the run directories. `eval_pdec.py <lib> <layout> <set>` is then run once per set. The libraries are the r7 HeLa library for D and H and the Astral arm-E library for A1-A3.

Per set there are three outputs: **never** (network runs), **always** (GBM runs) and **rule** (per seed, the GBM run if the rule says use, else the network run). Each output gets the pooled lower and combined bound at own q ≤ 0.01 and the median N @ FDP_c 1 %.

- **Admissible arm:** an arm is admissible on a set if its pooled lower bound ≤ never + 0.0020 and its combined bound ≤ never + 0.0030. These are PLAN_GBT criterion 2's bars, unchanged.
- **Oracle:** the admissible arm, of always and never, with the higher median N @ FDP_c 1 %. Never is always admissible.
- **Rule passes a set if:**
  - (1) the rule output's pooled lower bound ≤ never + 0.0020 and its combined bound ≤ never + 0.0030;
  - (2) its median N @ FDP_c 1 % ≥ 0.98 × the oracle's.
- **Verdict:** PASS if the rule passes all five sets, including raw 03. Any failed set is a FAIL and is reported per set.

## Consequences

- **PASS:** implement the statistic inside OpenDIA as a log-only diagnostic first, then as an opt-in switch (default byte-identical, byte-checked), then a frozen plan before any default change.
- **FAIL:** report it. The GBM stays as is, and nothing is retuned (no other region, level or statistic is tried on these runs).
- **Seen-data failure:** if the rule fails on D or H (it would skip a GBM that helps), that is reported as a fail too. It would mean the null region still holds true targets that the GBM legitimately rescues.
