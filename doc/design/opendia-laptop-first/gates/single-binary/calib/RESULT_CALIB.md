# RESULT_CALIB, 23 November 2026

**Diagnostic outcome: calibration is not established. No scoring change or revised benchmark verdict.**

Curve precheck committed and pushed as `d129bac5ea` before real-run curves. Structural audit precheck committed and pushed as `e8284fe208` before library audit execution. Evaluators are frozen by SHA256 manifests. Synthetic suite: 6 methods pass. The new evaluator reproduces **all 50 prior per-run N/E counts at own q <= 1% exactly** against EVAL_GBMNULL.txt.

## Calibration as a curve

Full outputs are `results/*.per_run.tsv` (450 rows, 5 sets x 2 arms x 5 seeds x 9 thresholds), `results/*.curve.tsv`, and `results/calibration_curves.{png,svg}`. The grid is q = 0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5%. No threshold was selected after inspecting results. Arms: never = no GBM, always = GBM. A3 is informational because the earlier anchor failed.

Selected points below show descriptive pooled lower / combined estimates in percent. Pooling sums E and N only. It does not supply independent-replicate confidence intervals.

| Set | Arm | q=0.1% | q=0.5% | q=1% | q=2% | q=5% |
|---|---|---:|---:|---:|---:|---:|
| D | never | 0.108 / 0.123 | 0.574 / 0.655 | 1.099 / 1.254 | 2.043 / 2.333 | 4.803 / 5.484 |
| D | always | 0.212 / 0.242 | 0.583 / 0.666 | 1.112 / 1.269 | 2.215 / 2.529 | 5.075 / 5.794 |
| H | never | 0.321 / 0.367 | 0.679 / 0.775 | 1.146 / 1.309 | 2.058 / 2.350 | 4.701 / 5.368 |
| H | always | 0.322 / 0.368 | 0.747 / 0.853 | 1.288 / 1.471 | 2.250 / 2.569 | 5.083 / 5.803 |
| A1 | never | undefined | 0.426 / 0.852 | 0.479 / 0.958 | 1.095 / 2.188 | 3.115 / 6.225 |
| A1 | always | undefined | 0.371 / 0.741 | 0.612 / 1.224 | 1.189 / 2.376 | 2.582 / 5.159 |
| A2 | never | undefined | 0.075 / 0.151 | 0.261 / 0.522 | 0.859 / 1.717 | 2.392 / 4.780 |
| A2 | always | 0.096 / 0.192 | 0.218 / 0.436 | 0.557 / 1.112 | 0.866 / 1.730 | 2.351 / 4.697 |
| A3 (informational) | never | undefined | 0.337 / 0.674 | 0.380 / 0.760 | 0.684 / 1.366 | 2.674 / 5.343 |
| A3 (informational) | always | undefined | 0.240 / 0.480 | 0.603 / 1.205 | 1.185 / 2.367 | 2.359 / 4.714 |

### Reading

- **HeLa reveals a problem beyond the GBM.** At q=0.1%, H's lower estimate is about 0.32% with either arm. Every seed shows it: never 0.305-0.331%, always 0.310-0.332%. This departure is already present without GBM, so merely disabling GBM cannot resolve it. These are repeated fits of the same sample, not independent confirmations.
- At q=1-2%, HeLa's lower estimates generally exceed nominal q, with larger departures under GBM on H. D's GBM departure is especially visible at q=0.1%. This motivates examining upstream scoring, selection and decoy construction, not a global empirically tuned q multiplier.
- **Astral is threshold-dependent, not uniformly loose by a constant factor.** A1's combined curve is above the identity at 1-5%, while A2 is below it at 0.5%, 2% and 5%, slightly above at 1%. Some strict-q runs have no discoveries. A combined estimate crossing the identity is not by itself proof of violated FDR control. The earlier GBMNULL operational FAIL remains unchanged.
- The lower/combined estimates are not deterministic bounds enclosing realized FDP, and the combined estimate is not an exact calibration target. No curve-wide hypothesis test or new pass/fail threshold is asserted. Seed envelopes are algorithmic variability, not confidence bands. Entrapment-generation and null-exchangeability assumptions remain essential.

## Paired estimator: prerequisite not met

See `PAIRING_PROVENANCE.md`. Composition/charge/modification matching yields 136,509 unique candidate pairs, not a complete authenticated original-generator map. The source ZIP inspected on Codon contains FASTA only. No paired real-run statistic was computed, no ambiguous matches were guessed, and no subset or unvalidated fallback was used to produce a favorable answer. The earlier claim that same composition alone establishes applicability is withdrawn.

## Asymmetry audit

See `ASYMMETRY_PLAN.md` and structural summaries under `results/`. Both locally audited libraries are SHA256-identical to the frozen Codon libraries (HeLa `76f2a6c2...`, Astral `d9e3fbe0...`). The bounded remote inventory found no `.features.tsv[.gz]` files under the GBMNULL layouts. Winner-only precursor TSVs are not substituted as candidate-count tables. Exact observed candidate-window and co-eluting sibling asymmetry therefore remains blocked.

| Library | Own-decoy pairs | Equal library sibling count | Target has more | Decoy has more | Unpaired precursors |
|---|---:|---:|---:|---:|---:|
| Astral | 300,942 | 299,891 | 1,051 (0.349%) | 0 | 1,537 |
| HeLa | 1,388,026 | 1,386,912 | 1,114 (0.080%) | 0 | 2,137 |

Every Astral unequal pair has one extra target-side sibling. HeLa has 1,113 pairs with one extra and one pair with two extra. Mean sibling multiplicity over paired members is Astral 1.48047 target / 1.47698 decoy and HeLa 1.14675 / 1.14595. These are structural opportunities, not observed-window counts or a proof of FDR inflation. All candidate-related output fields are explicitly null.

**There is a concrete symmetry defect to investigate:** missing own-decoy partners and one-sided sibling deficits are present in the prepared libraries. Whether these hypotheses enter the search, are accepted, or materially drive calibration has not been measured by this precheck. Before claiming paired-TDC validity, audit their actual workflow membership and treatment, then freeze a complete-pair handling experiment rather than tuning a q multiplier. Complete JSON audit records remain in `/home/singjc/.jcode/scratch/dia2/{astral,hela}_structural_calib_20261123.json`; committed summaries record their SHA256 digests and library provenance.

**Source-map wording correction:** `LABEL_USE.md` calls the sibling grouping sequence unmodified. The loader actually uses the modified sequence when present (`OpenDIALibrary.cpp:260`), which is what the frozen structural audit implements. This wording correction does not change either evaluator or its frozen hashes.

## Next mechanism-based work

1. Recover the original entrapment pairing/generation and sampling metadata, or use a newly generated explicitly paired benchmark under a separate frozen plan. Do not change these run scores or gate interpretations to compensate for missing metadata.
2. Obtain matched candidate exports with unchanged scoring, hashes, timing/memory and default-byte checks under a new frozen run plan. Check candidate and sibling opportunities before final keep. Use entrapment only for grading a frozen diagnostic, not for choosing its threshold.
3. Only then design a scoring/statistical intervention with the source label-use map (`LABEL_USE.md`) in view. A final-stage split alone is not an end-to-end RESET guarantee. No q rescaling, seed selection, default flip or RESET implementation was performed here.
