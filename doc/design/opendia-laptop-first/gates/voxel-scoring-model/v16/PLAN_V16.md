# PLAN_V16: Stage D mechanism audit (diagnostic only)

Written: 2026-11-12T08:54:44Z. Freeze before new quantitative analyses. No recovery claim or model promotion is authorized by this diagnostic.

## Fixed scope

- Keep the original v2 population and D11FWT reference on raws 01/02. Preserve all existing inputs, reference arrays, scoring recipe, comparison base, and pass rules.
- Stage D is the existing DIA-NN-reference/right-window/below-empirical-FDP_c-1%-cut label, not automatically a failure of the operational q<=0.01 filter. Explicitly report this distinction and overlap. DIA-NN agreement is not ground truth.
- Raw01 is discovery; raw02 is descriptive replication of the same questions, not an untouched holdout after the many prior experiments.
- Entrapment and DIA-NN annotations are evaluation-only. They must not select model inputs, positives, negatives, features, or trained rescue rules. No intervention model training, threshold search, parameter sweep, or C++ changes.
- An exact reconstruction of the EXISTING D11F training recipe is allowed solely to save training-selection traces, model states, and held-out feature sensitivities. Use the original seeds/rounds/validation split and stop attribution if reference scores cannot be reproduced (Spearman >=0.99999, W row agreement 1, maximum absolute score difference <=1e-3). Save under this diagnostic directory, never overwrite reference arrays. Report any reconstruction discrepancy.

## Questions to answer

1. Threshold semantics and distance: which D rows already pass q<=0.01? How many are near the cut versus deep in the null tail? Preserve the original empirical FDP curve, compute pooled TDC as the existing evaluator, and never treat decoy competitors as accepted target IDs.
2. Score mechanism: compare the SAME model's linear contribution and nonlinear contribution, normalize within fold against the same W-selected decoy rows, and report where the residual improves or damages target/null separation. Attribution is descriptive, not proof that removing the residual is useful; linear-only was already worse.
3. Training coverage: reconstruct rounds if feasible, save label-blind positive/negative/validation membership before joining evaluation labels, then count D and E inclusion in the two models for which their family was eligible to train. Determine whether D-like rows were absent, removed during tightening, or included yet still poorly ranked. Report selection-conditioned limitations.
4. Measurement mechanism: quantify RT disagreement relative to DIA-NN FWHM, local peak width/fragment support, native vs binned profile coherence, and mass-error/support where present. Do not infer every D is at the right apex merely from a median RT difference of zero. Use existing raw evidence before expensive re-extraction. Identify whether signals are missing, poorly summarized, or simply weak/ambiguous.
5. Controls: D, accepted E, all human target winners, entrapment winners, and decoy winners. Use same-fold score strata or documented score matching, report coverage and tie-corrected AUC (0.5 for ties). If matching on area/support, report added matching bias and do not claim independent validation. All-human comparison is essential because DIA-NN-selected D can exaggerate apparent feature separation.
6. Prior failures: quantify why a proposed action is distinct from V4/V5 MS1/coelution/mass evidence, V7 cross-run features, V9-V11 native-count variants, V12 profile correction, V13 deconvolution, and V14/V15 amplitude-only diagnostics. No repeated approach without new measured mechanism.

## Outputs and integrity

- Scripts, tests, execution logs, hashes before execution, and tabular outputs under v4/stageD_v16/. Freeze each script revision before its associated analysis and record changes transparently.
- Test scoring/matching/round-logging helpers and verify population row alignment, fold isolation, finite values, and known D/E totals. Report missingness and controls rather than silently dropping rows.
- If informative, provide a compact diagnostic figure showing subgroups, score effects, and training selection. No hand-picked success-only panels.
- Synthesize a result with supported causes vs open hypotheses. Propose a single primary recovery experiment, not a collection of tuned arms. Preserve the existing global pass rule and explicitly audit newly gained human/entrapment/decoy counts. Any additional go/no-go condition must be presented as a proposal, not substituted for the current gate.
- Stop this task at diagnosis plus a concrete next plan. Do not claim increased IDs without a separately frozen intervention run.
