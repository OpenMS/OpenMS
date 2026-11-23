# PRECHECK_CALIB (23 November 2026)

## Scope and status

Diagnostic reanalysis of already observed data, NOT a new confirmatory test or a replacement for previous gate verdicts. No OpenDIA code or defaults change. Freeze and push this file and evaluator before computing any new real-run curve. Existing 1% results and library composition checks have already been seen.

## Pairing prerequisite

The original generator's peptide_pair_index map is required for the published paired estimator. Same composition, C-terminal residue, charge, and modification multiset identify only 136,509 unique 1:1 precursor matches among 151,377 Astral entrapment precursors. These are candidate matches, not authenticated generator pairs. The earlier assertion that composition alone establishes applicability is withdrawn. No arbitrary matching of ambiguous groups, no unvalidated doubled-count fallback, and no subsetting discoveries to make the pairing fit. If a complete authenticated map cannot be recovered, record the paired analysis as blocked and report lower/combined curves instead. A near-one library ratio alone is not an applicability criterion.

Reference: Wen et al., Nature Methods (2025), DOI 10.1038/s41592-025-02719-x. Authors' FDRBench README requires the generator's peptide pair file for precursor paired FDP estimation.

## Frozen population and grid

Reuse GBMNULL layouts in /nfs/research/juan/singj/opendia/gbmnull-20261123/lay/{hela,astral}. Sets D, H, A1, A2, A3. Arms never and always, seeds 0,1,2,3,4, all mandatory. No new seed selection. A3 remains informational because the earlier anchor failed. HeLa uses the existing r7 library, Astral the arm-E library. Score/q TSVs as written, no recompetition or q recalibration.

Evaluate own q <= alpha at alpha = 0.001, 0.0025, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.03, 0.05. For each run count accepted target precursors N (including entrapments), entrapments E, lower E/N, combined E*(1+1/r)/N. Zero discoveries give undefined FDP, not a calibration pass. Unknown IDs or missing runs are errors. Save every per-run row, median/min/max across seeds, and descriptive pooled E/N. Existing 1% counts must reproduce prior GBMNULL results.

## Interpretation and uncertainty

Seeds are repeated fits of the same spectra and entrapment library, not independent experimental replicates. Do not pool them into independent Poisson/binomial confidence intervals. Seed spread measures algorithmic variability only. Thresholds are nested, not independent tests. No pointwise interval converted into a simultaneous curve-wide claim. Lower/combined are estimator outputs under their assumptions, not deterministic confidence brackets for realized FDP. A combined estimate > alpha alone does not prove violation, and an estimate <= alpha does not establish FDR control. Report curve shape and departures descriptively, with no new pass/fail bars or retuning. Independent raw samples and validated null exchangeability are required for stronger claims.

## Target-decoy asymmetry

The companion ASYMMETRY_PLAN.md will be frozen separately before any audit statistics are computed. It defines label-free candidate-count and sibling-availability diagnostics. Biological target/decoy differences alone do not demonstrate null asymmetry. No diagnostic changes scoring or selects an arm. Entrapment can grade but must not feed a decision rule.

## Deliverables

Per-run and aggregate curve TSVs, plotted curves with identity line and seed spread, pairing provenance outcome, asymmetry audit outputs or explicit artifact limitations, RESULT_CALIB.md. Existing gate verdicts remain unchanged. Scripts and this precheck hashed before execution. Run board_sync.py --archive-stale before each docs commit.
