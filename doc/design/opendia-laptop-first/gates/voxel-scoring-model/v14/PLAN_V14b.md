# PLAN_V14b: Local-flank weak-profile counterfactual

written_at: 2026-11-12T00:40:00Z
status: FROZEN BEFORE V14b FEATURE GENERATION AND SCORING

## Question

Does the existing single D11F model down-rank a correct, coherent Stage-D-like peak group when the peak's positive signal component is attenuated to the observed Stage-D/accepted abundance ratio, while its measured local fragment background is retained?

This is a diagnostic of the current scoring function, not a proposed model change. The model's 11 heuristic inputs remain unchanged. Only the six F12 summaries are recomputed from the counterfactual XIC profiles. Therefore a null result applies only to the F12 path and does not rule out attenuation sensitivity in the unmodified heuristic inputs.

## Frozen evidence and scope

- Use raw01 only. Raw02 is not part of this diagnostic.
- Base model is the existing D11FWT pipeline on the V2 population. Do not change the comparison base, training objective, hyperparameters, pass rule, or C++ implementation.
- The attenuation factor is fixed at 0.25655589613105056, calculated from the already-recorded raw01 Stage-D and accepted DIA-NN `log_area` medians (4.314779 and 4.905597). DIA-NN is used only to set this diagnostic's signal scale and for evaluation, never as a training label.
- The feature-generation script reads no entrapment, DIA-NN, or OSW labels. It selects target donors by the existing baseline target-decoy winner/q procedure, using the library decoy flag only.
- Do not train an augmented or alternative model. Re-run only the existing D11F cross-fold training recipe as a deterministic baseline reconstruction because its in-memory fold scorers are needed to score the counterfactual inputs. This reconstruction must reproduce the saved D11F output before counterfactual scores are accepted.
- Do not promote this approach or change the C++ model based on this diagnostic alone.

## V14 deviation recorded

The initial V14 paired-donor idea was stopped before feature extraction. Its q<=0.001 target donor set was empty because the add-one target-decoy q procedure has a minimum attainable q of about 0.001499. At q<=0.005, only 531 of 1,327 target winners with a same-key decoy candidate had a decoy window within 5 seconds (about 40 percent), so the frozen paired-profile coverage criterion was not met. No V14 trained model or scored result was produced. V14b replaces that invalid matching assumption with a per-trace local flank baseline and must be treated as a distinct frozen diagnostic.

## V14b method

1. Select raw01 non-decoy precursor winners from the existing `scores_D11FWT.npy` using the existing pooled target-decoy q calculation. Keep donors with q<=0.005. Expected donor count is roughly 1.3k, subject to input integrity checks.
2. For every donor and every one of its 12 XIC fragment traces, estimate local baseline as the median of bins 0:8 and 24:32. The central signal region remains bins 12:20.
3. Transform every bin as `baseline + alpha * max(observed - baseline, 0)`, where alpha is the frozen abundance ratio above. This preserves observed local flank background and attenuates positive signal above that baseline. Do not add synthetic noise or use a paired decoy trace.
4. Recompute the six F12 features using the existing `features12.py` definitions. Validate the implementation against stored F12 for up to the first 256 donors before saving. Require maximum absolute feature error <=2e-5.
5. Run the unchanged D11F training recipe with arm name D11FQ to reconstruct the baseline. Score the saved counterfactual F12 summaries with the exact held-out fold scorer, while leaving all other inputs unchanged. This is an instrumentation path only. No intervention model is fit.
6. Require Spearman correlation >=0.99999 between reconstructed and saved D11F scores. Also verify that the W-selected precursor/window identities are unchanged. If either check fails, stop and report the counterfactual as invalid.
7. Convert counterfactual logits through the existing fold-specific W/T calibration mapping, using the baseline's saved D11FWT output and fold-selected decoy calibration rows. Validate the mapping on baseline rows before comparing score distributions.

## Evaluation, after scoring only

- Report donor count, feature reproduction error, baseline reconstruction correlation, W-selection agreement, and counterfactual score deltas.
- Join target/decoy/entrapment labels only in the evaluation step. Entrapment labels must not enter candidate selection, feature generation, baseline reconstruction, or scoring.
- Compare attenuated donors to the observed Stage-D raw01 score distribution and to score-matched decoys. Report target-vs-decoy AUC for the transformed donor/control panel and entrapment-vs-decoy AUC on the unchanged baseline control panel, with exact panel and matching definition documented.
- Record whether transformed target scores move into the Stage-D range and whether the target/decoy separation degrades. This is a diagnostic, not a performance pass claim.
- Stop without any augmented training if the synthetic target-vs-decoy AUC is below 0.60 or if the existing entrapment-vs-decoy control AUC exceeds 0.53. Any subsequent augmentation experiment needs its own plan and freeze.

## Integrity and execution

- Preserve the original `PLAN_V14.md` and its hash. Record V14's pre-extraction stop as a deviation, not as a result.
- Freeze this plan's SHA256, the feature-generator SHA256, and the instrumented trainer SHA256 before execution.
- Keep generated data and logs under `v4/weak_aug_v14/`. Do not overwrite reference arrays.
- No C++ promotion, no pass-rule change, and no changes to the D11FWT base.
