# Interference correction, round 6: transition refinement in the native CNN + LDA

Written before any round 6 run, after `RESULT.md` rounds 1–5 and the diagnostics. Code: `6bca00e7a9`. Raw 03 stays untouched.

## Why

The decoy-swap and input diagnostics show that targets carry more fragment signal at their m/z values than decoys, even when the target is absent (entrapment). The AUC of entrapment winners over decoy winners is 0.531 for matched scans, 0.525 for transition support and 0.512–0.519 for the CNN logit. Commonness-matched decoys remove only about half of this.

Hypothesis: much of this signal is interference, meaning fragments of other peptides that do not share one elution shape. EncyclopeDIA's automated transition refinement (Searle et al. 2018, Supplementary Figure 4) keeps only the fragments that match the median unit-area peak shape. Two ways to use it:

- **Refine features** (`-stage2_feature_set ...+refine`). These add new evidence that counts only co-eluting fragments. The CNN and the other features do not change.
- **Trace mask** (`-trace_refinement_min_r 0.75`). This zeroes fragments with r < 0.75 in the tensor before the CNN and every feature.

Both treat targets and decoys the same way and use no labels.

## Arms

Each arm is the native base command (`astral-raw02-byion-design-20261109/run.sh`: `-stage2 -stage2_feature_set search+boundary+byion`, CNN grid2d, `-ss_iterations 0`), with one change:

| Arm | Change |
|---|---|
| **R1** | `-stage2_feature_set search+boundary+byion+refine` |
| **R2** | R1 plus `-trace_refinement_min_r 0.75` |

Order of runs:
1. HeLa default equivalence. `crossfit.tsv` must be byte-identical, or the chain stops.
2. R1 on raw 01, then R1 on raw 02, then R2 on raw 01, then R2 on raw 02.

All runs are sequential, each about 15 min and 25 GB.

**Internal check for R1 (part of rule 4).** The CNN is not changed in R1, so the `crossfit.tsv` columns 1–7 and 9 (sample, logit, fold, best, window) must be byte-identical to the base run. Feature columns 1–15 of `features.tsv` must also be byte-identical to the base run.

## Decision rule

This is the `PLAN.md` rules 1–4, unchanged, applied to each arm on both raws, with the score `stage2_score`:

1. N@1% ≥ 1050 (raw 01) and ≥ 1324 (raw 02).
2. FDP_c at q ≤ 0.01 is ≤ 0.015.
3. FDP_c at q ≤ 0.05 is ≤ 0.070 (raw 01) and ≤ 0.079 (raw 02).
4. The run completes with peak RSS ≤ 30 GB and wall time ≤ 35 min. For R1, the internal check above must also hold.

Outcomes:
- If both arms pass, the arm with the higher summed N@1% is preferred.
- A passing arm is a candidate for the native default. It does not become the default without the raw 03 gate, which needs approval.

## Mechanism check (reported, not a pass condition)

The check is the AUC of entrapment winners over decoy winners, over the in-window precursor winners by `stage2_score`. It is reported for:
- each refine column;
- `transition_support` and `log_matched_scans` (R2 only; these come from the search and are not masked);
- `cnn_z`;
- `stage2_score`.

The hypothesis is supported if two things hold on both raws:
1. `refine_quantitative` has |AUC − 0.5| < |AUC(transition_support) − 0.5|.
2. For R2, `cnn_z` has an AUC nearer 0.5 than the base (raw 01: 0.512–0.519).

Also reported: N at FDP_c 0.5 / 1.5 / 2 %, q counts, LDA weights per fold, masked-channel counts, time and RSS. All q-derived counts are decoy-count estimates.

## Correction before any run (2026-11-09T23:40:37Z)

The base features.tsv has 16 columns before the new group: ids, cnn_z, search, boundary and byion. The R1 internal check therefore compares feature columns 1-16, not 1-15. Nothing else changes. The hash of the text above this section is in PLAN6.md.sha256.
