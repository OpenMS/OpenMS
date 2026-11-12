# Proposed V17: make self-training use the deployed window selector

Written: 2026-11-12T09:21:00Z. Design proposal following V16, not an executed experiment. Keep D11FWT as the accepted reference. Implementation, tests, source hashes, and input hashes must be frozen before an intervention run. This document alone is not an execution freeze.

## Hypothesis, and what this does not claim

D11F learns and refreshes its training set using each precursor's full-score-maximizing window. Deployment W instead chooses the linear-score-maximizing window and applies the full model there. In the original final self-training round, the two selectors differ by >10 seconds for about 77%/78% of eligible decoy candidates on raws 01/02. Admitted D positives almost always already use the same peak (round 2 >10 s mismatch 3/745 and 1/931), so a story about widespread wrong positive windows is not supported. However, of the deeply suppressed D (below the decoy 90th percentile, 287/230), 251/212 were never admitted as positives and have a full-vs-linear selector disagreement >10 s: training saw a different window for those precursors than deployment scores.

The hypothesis is narrower: training and q-based positive admission use a materially different null-candidate distribution from deployment. This may contribute to excluding weak, coherent targets or learning an unnecessarily unfavorable boundary for W-selected windows. The existing diagnostic does not prove either causal effect. Maximizing full scores also supplies useful hard negatives, so removing that selection pressure may worsen specificity. The unchanged entrapment and decoy gates are essential.

This is not interference subtraction, an evidence bonus, a broader target label, a second stage, an architecture search, or another version of V4's sibling-ranking loss. No DIA-NN or species label enters training.

## One intervention

Use the **existing W window-selection rule inside self-training rounds 1 and 2**, symmetrically for target and decoy precursors:

1. Evaluate the current round model's linear contribution and full score on the eligible training-fold candidate rows.
2. Select one row per precursor by the linear contribution, using the existing W tie rule.
3. Use the full score at that row for the unchanged per-fold median/MAD normalization, paired target-decoy competition, and q estimation.
4. Admit positive target rows at the existing q<=0.15 and q<=0.05 round thresholds. Use the corresponding selected decoy rows as negatives.
5. Refit the same single D11F model with the same recipe. At deployment apply the existing W and T operations unchanged.

The only treatment is the candidate-selection policy inside `pooled_q`. Full scores still determine confidence, positive admission and ranking. The linear contribution is part of the same model, not a separately fitted classifier. Do not use the frozen baseline's linear scores as a teacher for the new model.

Keep round-0 anchors, 17 input features, preprocessing, fold and validation hashing, seeds, decoy-whitened objective, network, optimizer, training schedule, early stopping, calibration, population, and evaluation unchanged. Do not combine this with dropout simulation, additional features, native extraction, residual damping, or interference correction.

## Preconditions and technical checks

- Preserve original training source and all baseline arrays. Use an isolated experimental runner and new output names.
- Implement one pure candidate-selector helper with explicit separate selection and confidence scores. Test a case where their maxima are different: the linear score must choose the row, but the full score must be retained for confidence.
- Verify target/decoy symmetry, deterministic ties, one row per precursor, pair-decoy tie priority, empty/degenerate cases, and training-fold exclusion. For any helper simplification preserve original behavior exactly in the control path.
- Verify W's exact tie rule against `make_w.py`; if training row ordering and sample_id ordering differ, document that technical discrepancy before running rather than silently changing the control.
- Reproduce the unchanged D11FWT control with the original ten Torch threads and saved execution environment. Require original full and linear scores to match within 1e-3, Spearman>=0.99999, and W row agreement exactly 1. V16 achieved zero error on both raws.
- Freeze implementation hashes, dependency versions, input hashes, tests, runner, arm name, and write time before training. A mismatch stops the run before attribution.
- Training/selection logs must be label-blind. Join D/E/DIA-NN/entrapment annotations only after outputs are complete. Include all target-library candidates, not a D-selected subset.

## Run order

One treatment on raw01, then the identical treatment on raw02 without a raw01-guided revision. No seed sweep or threshold tuning. Save checkpoints, round-specific memberships, linear/full candidate rows, final scores and STATUS/EXIT/log files. Do not treat raw02 as a pristine holdout after the prior experiments.

The exact baseline reconstruction took a few minutes for both raws in V16. Run locally only if the estimated complete job remains below the user's 25-30 minute limit, using setsid/nohup/systemd-inhibit. Otherwise use the authorized codon compute workflow. Poll for at most 15 minutes per monitor and terminate only by PID if needed.

## Existing pass rule, unchanged

Both raws must satisfy the established V10/V11 rule:

1. PLAN_V2 rules 1-4 against original population BASE: bootstrap 90% interval for delta N@FDP_c 2% is above zero, the interval for delta N@1% is not entirely below zero, FDP_c at operational q<=0.01 is <=0.015, and entrapment-over-decoy winner AUC does not exceed that BASE's recorded value. Use the canonical scripts and records for exact bootstrap implementation and reference AUC precision, not rounded prose.
2. N@FDP_c 2% exceeds D11FWT on **both** raws: 1716 and 1736 respectively.
3. Null-window winner AUC is <= the corresponding D11FWT value +0.005.

Retain the existing distinction between the all-precursor-target empirical FDP curve and paired competition used for operational q. Do not change the comparison pool to make a result pass.

## Additional required reporting, not replacement pass rules

- Full N@FDP_c 0.5%, 1%, 1.5%, 2% and operational q results, paired bootstrap comparisons, and all existing null diagnostics.
- Total accepted human/entrapment targets and paired decoy counts at q1%, with baseline-to-treatment gain/loss tables. At baseline these are human1603/1526, entrapment11/10, paired decoys15/14. Report counts as decoy-count estimates where appropriate.
- Historical D recovery **and** previously accepted E losses. Also report genuinely operational-missed D (excluding baseline q1% D87/99) and total human retention, not only a favorable historical subset.
- Report the near-cut q1-5% subgroup, deep-null subgroup, and paired-decoy losses separately, with no subgroup-specific cutoff.
- Record whether the intervention actually changes weak-window positive admission and decoy-window selection. If it passes the global gates but does not support the proposed mechanism, say so.
- Report entrapment and hard-decoy tails, not just overall AUC. A lower FDP ratio can conceal an increase in false-target counts when the denominator grows.

The user wants fewer false discoveries, not just a larger list. Passing the historical gates while increasing absolute entrapment counts must be described as a tradeoff, not a demonstrated reduction. Any stricter deployment condition requires an explicit proposal rather than silently rewriting the pass rule.

## Decision and stop condition

A pass makes this a candidate for discussion, not automatic C++ promotion. If either raw fails, retain D11FWT and stop this selector-alignment arm. Do not salvage it by mixing in previous feature variants or selecting a favorable seed. No claim is made that the 394/366 never-admitted D precursors are all recoverable.

If the intervention fails, the next *separate* question would be training coverage under real fragment missingness and sparse temporal sampling, not amplitude-only attenuation. Such an experiment first needs a label-blind, raw-evidence-consistent perturbation with all 17 features recomputed and null controls. It is deliberately not combined with V17.
