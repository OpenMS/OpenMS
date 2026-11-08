# R-WP4: Error control

[Master plan](../README.md)

**Objective.** Report precursor q-values with independent evidence that they control error.

**Delivered:**
- Family outer-fold cross-fitting with fold-wise q-values (`c355109de2`).
- Opt-in paired target-decoy competition, where `DECOY_<X>` competes with `<X>` (`ae4687e508`).
- Recorded failure: the log(n)-window correction failed its gate (`0c757c3d2d`).

**Open gates:**
1. An entrapment or two-species check on HeLa-like data.
2. An independent-run result without a same-run teacher.
3. A best-of-n-windows decoy bias analysis. More windows per precursor raises the best decoy score.

## Cross-run design (proposal, not implemented)

Per-run calibration and score iteration (R-WP2, R-WP3) stay the first stage. The cross-run step comes after, in two separable parts.

**A. Global score model.**
- Pool round-0 seeds from every run's calibration anchors.
- Train one CNN with outer folds keyed by peptide family across runs, so a family's windows in every run share one fold. This keeps cross-fitting honest across runs.
- Score every run with the fold model that did not see that family.
- Score iterations run on the pooled population.
- This puts scores on one scale and lets weak runs borrow from strong ones.
- Cost: one training on pooled seeds, which can be capped per run, and scoring is unchanged.

**B. Global error control** (PyProphet global context).
- For each precursor, take the best score over runs within each fold, then run one experiment-wide paired target-decoy competition on those maxima.
- Report run-specific q-values (per run, as now) and experiment-wide q-values, with run-level IDs filtered to experiment-wide q ≤ 0.01.

A without B still needs per-run q-values. B can use per-run models, but then it compares scores from different models, so it needs A or per-run score normalisation.

Gates:
- Run-specific counts with A must not fall on any run, compared with per-run models.
- B must be validated by entrapment on a multi-run set. The ProteoBench Astral entrapment runs 01–03 are a candidate, once the rewrite reads Thermo `.raw` with no 1/K0.

## First entrapment result: Astral raw 02, 8 November 2026 (FAIL)

The frozen teacher-free design with `-pooled_tdc` was evaluated once against a hashed gate (`gates/astral-raw02-20261108`).
- Result: 360 targets at pooled paired q ≤ 0.01, of which 5 are entrapment. FDP_c = 0.028, against a limit of 0.015 (95% interval 0.009 to 0.065). At q ≤ 0.05: FDP_c 0.044.
- The seed task is nearly trivial: held-out seed-decoy AUC is at least 0.9995 per fold. The entrapment hits have high logits, and decoys do not model them.
- Candidate causes, each needing a new design run and gate:
  1. Decoy windows that are easier than false targets. Shuffled-decoy traces lack co-eluting interference.
  2. The pooled MAD standardisation over small decoy-winner sets.
  3. The best-of-many window choice per precursor inflating target scores beyond what paired competition corrects.
- Until one of these is settled, the CNN q on weak runs is a decoy-count estimate that this run shows to be anti-conservative.

### Post-mortem (`gates/astral-raw02-20261108/postmortem/POSTMORTEM.md`, descriptive only)

- **Tail behaviour:** entrapment winners outnumber their decoys about 1.2 to 1.7 times in the score tail. Decoys of entrapment and decoys of human behave alike.
  - The anti-conservative part is the shallow tail, q ≤ 0.03, which rests on 2 to 3 decoys.
  - At q ≤ 0.05, FDP_c is 0.044, below q. This corrects the first RESULT.md wording.
- **`OpenDIA -out_xipm` extraction review:**

  | | co-elution | library intensity correlation |
  |---|---|---|
  | accepted human IDs | 0.91 | 0.94 |
  | accepted entrapment hits | 0.29 | 0.29 |
  | top decoys | 0.20 | 0.28 |

  - The CNN rewards fragments that carry signal individually, without requiring them to co-elute or to match the library pattern.
- **Next design direction:** a second-stage score that adds co-elution, library correlation and MS1 mass error to the CNN logit (PyProphet-style). It must be tested on a new design run and gate.

### Stage-2 feasibility on raw 01 (design run, target/decoy only; `gates/astral-raw01-stage2-feasibility/`)

- **Separation below the 1% cut:** with the frozen iteration-0 scores, co-elution, library correlation and MS1 mass error still separate targets from decoys in CNN z 6 to 10. That band lies just below the current 1% cut.
  - MS1 ppm error is the strongest single feature: 1.3 ppm for targets against 7.1 ppm for decoys above z 10.
- **Hand rule:** about 690 targets against about 3 decoys above z 6, extrapolated from 60-sample bands. The CNN-only result is 321. This is a feasibility signal only.
- **Next:** a native, opt-in second-stage discriminant over CNN logit plus extraction features, designed on raw 01. Evaluating it needs raw 03 and a new gate.

### Stage-2 score: native implementation and raw 01 design (`gates/astral-raw01-stage2-design/`)

- **Code:** `OpenDIACNN -stage2`, commit `014cb0dbf1`. It is opt-in and requires `-paired_tdc -pooled_tdc`. Without `-stage2`, the HeLa crossfit output and the fold checkpoints are byte-identical.
- **Inputs to the score:** per window, the fold-decoy-standardised CNN logit, plus seven extraction features from a core window of ±4 s:
  - fragment fraction;
  - fragment co-elution;
  - library intensity correlation;
  - absolute fragment ppm and absolute MS1 ppm, each relative to the offset of the seed windows;
  - MS1-to-fragment co-elution;
  - first-isotope presence.
- **Fitting:** a PyProphet-style semi-supervised ridge LDA, cross-fitted so that each fold's weights come only from the other folds. It runs 3 rounds, with positives first at q ≤ 0.15 and then at 0.05. Pooled paired TDC follows.
- **Raw 01 design run (T/D only):**
  - The rule was written before the run. It passed with 333 targets at q ≤ 0.01, against 321 for CNN-only.
  - Weights were stable and `cnn_z` was positive in all folds.
  - At q ≤ 0.05 the count is 770 against 634. The defaults are frozen.
- **Raw 02 re-score (`gates/astral-raw02-stage2-diagnostic/`):** a diagnostic only, because raw 02 is spent and the features came from its post-mortem.
  - 518 accepted at q ≤ 0.01, 3 of them entrapment.
  - FDP_c 0.0117, with a 95% CP interval of 0.002 to 0.034. The CNN-only figures were 360, 5 and 0.028.
  - Above q 0.02 the estimate is still anti-conservative: FDP_c 0.054 at q ≤ 0.05.
- **Open item:** an independent check needs raw 03, the last held-back replicate, under a new hashed gate. It waits for the user's go-ahead.

### Stage-2 feature groups: search peak evidence and MS1 isotopes (commit `48e40f4d35`; `gates/astral-raw01-stage2-search/`, `gates/astral-raw01-stage2-isotope/`)

- **Option:** `-stage2_feature_set` takes `+`-separated groups after the CNN logit.
  - `extraction`: the default, with output byte-identical to before.
  - `search`: transition support, log matched scans, log peak width and |library RT error|, read from the OpenDIA `-features_out` table via `-stage2_search_features`.
  - `isotope`: MS1 M+0..M+3 averagine correlation, M−1 overlap and M+1/M+2 co-elution.
- **Raw 01 design runs** (target/decoy only, rules written before each run):
  - `search`: **717** targets at q ≤ 0.01 (1412 at 0.05). The extraction set gives 333 and CNN-only gives 321. The native run reproduces the offline ablation exactly. PASS; this is the frozen design candidate.
  - `search+isotope`: 650 at q ≤ 0.01. FAIL under the +2% rule, so the group stays off.
    - The isotope features separate targets from decoys above CNN z 6.
    - They are near noise in z 4–6, where the 1% cut is decided.
- **Context:** OpenSWATH+Percolator on the same library gives 1441 on raw 01, so the remaining gap is about 2×. This has not been entrapment-checked.
