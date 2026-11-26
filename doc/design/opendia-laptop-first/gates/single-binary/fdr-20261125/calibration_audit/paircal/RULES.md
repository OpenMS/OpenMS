# Fixed pair-stratified calibration screen

Frozen 2026-11-26 before running this screen. Exact coordinator DM recovered from saved session at 20:16 UTC after earlier DM was not delivered into worker context.

Inputs: D/H cur (NET) and gbt (G0) full precursor tables for seeds 0,1,2,3,4 under `/home/singjc/.jcode/scratch/blend-20261125/eval/ref/entlib-gbt-20261122/runs`. Library: `/home/singjc/.jcode/scratch/blend-20261125/eval/hela_r5.predicted_decoys.pqp`.

Baseline: best valid score per precursor then best member per library target-decoy pair by G0 score, exact member ties to decoys. Plus-one TDC q-values with reverse cumulative minimum, equal score groups processed together. Baseline and candidates use identical recomputation.

A agreement: whole pair assigned agree/disagree according to whether NET and G0 pair-winner precursor IDs agree. Pair absent from NET is conservatively its own missing stratum, never silently dropped.
B net99: whole pair assigned high/low according to pair maximum NET score > 99th percentile of ALL NET-kept decoy precursor scores in that run. Missing NET pairs get missing stratum. Quantile uses numpy linear interpolation.

Within each fixed stratum compute G0 pair competition and plus-one TDC. No grid search or adjustment after results. Selection consumes only numeric precursor identity, library pairing, target/decoy flag, and scores. Protein, organism, human similarity and entrapment annotations are forbidden inputs to selection. All candidates retain whole pairs before competition.

Evaluation is a separate script after selection: original and exact-label-corrected pooled combined bound E*(1+1/r)/N at q<=0.01, five per-seed counts, 90% Poisson interval, median N@FDP_c1% versus same-recomputed baseline. Exact-label correction comes from `gates/single-binary/fdr-20261125/stepA/label_audit_hela.tsv`, evaluation only. D/H screen is exploratory and not held-out validation. No source edits or remote runs.
