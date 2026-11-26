# DIAG_CALIBRATION: species-independent training and competition audit

26 November 2026. Follow-up to DIAG_NONBORROW. The user requires organism-independent FDR calibration, not a human-specific exclusion. Entrapment annotations remain evaluation-only. No default or production source was changed. Source revision is `4c6ca367fd`.

## Bottom line

Two concrete training/extraction isolation concerns are now identified. Neither is yet proved to cause the observed D/H FDP excess:

1. GBM row fitting excludes the held-out fold, but positive selection and hard-negative mining are computed globally first. An executable perturbation test demonstrates that held-out scores can alter that fold's training set and fitted predictions.
2. RT refinement uses a full fit for the search centers of its own training targets, whereas decoys never train it. Its out-of-fold residuals select the window width, not those search centers. A frozen ten-run D/H knockout is now running to measure the broad refinement effect.

Search/competition conditional-symmetry checks pass. A species-free, iid-null synthetic GBM stress does not reproduce inflation, so the isolation gap alone is not evidence of miscalibration on all data.

## 1. GBM selection depends on the held-out fold

`OpenDIAGBM.cpp:383-408` computes TDC positives and the hardest 2% of decoy precursor winners using all folds. Only at `:426-429` are held-out rows removed. Validation split/weights and the fit follow at `:437-451`. `:456` predicts the held-out fold. Later rounds select from the pooled previous round again, giving a further route for cross-fold feedback.

By contrast, the network stage in `OpenDIAScorer.cpp:628-689` computes feature scaling, seed fitting and self-training inside each outer training partition. However, using its pooled out-of-fold predictions as GBM training labels is not the same as nesting the entire two-stage learner: another fold's network prediction may come from a network trained on the current outer fold. Merely moving the GBM threshold into the fold does not remove this indirect dependency.

### Executed synthetic intervention

`calibration_audit/training/selection_isolation.cpp` uses production `pooledPairedTDC` and `HistGBM`. 15,000 numeric pairs, three folds, four features. Shift only held-out-fold seed scores by +0.5 for decoys and -0.5 for targets, keeping every feature row and the other folds' scores fixed. Compare selected training sets and raw predictions on the same 10,000 fixed feature probes. The diagnostic fits at most 40 trees rather than the production default 300 iterations.

| Selection policy | Positive-row symmetric difference | Negative-row symmetric difference | Changed probe predictions |
|---|---:|---:|---:|
| Production global-before-exclusion rule | 78 | 83 | 10,000 / 10,000 |
| Exclude held-out fold before selection | 0 | 0 | 0 / 10,000 |

The control establishes direct selection isolation for fixed seed scores only. It is not a fully nested scorer, a production fix or an FDP experiment. This numerical test was corrected during development to use the actual `max_iter` option and a mild perturbation; only the final version is archived. No biological labels or reference tables were edited.

### Full D/H network-table check, seed 0

`export_selection.py` reads PRECURSOR and PEPTIDE data, never PROTEIN/entrapment annotations. It exports full precursor-winner score tables, not q-truncated feature rows: D 2,420,158 rows, H 2,395,697 rows. Family-fold reconstruction agrees with all 1,000 checked production feature rows in each set.

| Set | Held-out fold | Unperturbed global vs local positive-row difference | Perturbation: global positive-row difference | Perturbation: global negative-row difference |
|---|---:|---:|---:|---:|
| D | 0 | 102 | 3 | 2,835 |
| D | 1 | 47 | 4 | 2,837 |
| D | 2 | 16 | 9 | 2,837 |
| H | 0 | 36 | 0 | 2,808 |
| H | 1 | 5 | 0 | 2,808 |
| H | 2 | 59 | 0 | 2,809 |

The same fixed +/-0.5 score perturbation is used. Local-before-selection changes zero positive and zero negative rows in every fold. These negative counts are winner-row set differences, not the number of all candidate windows from each hard precursor that production ultimately trains on. The positive sets correspond to first-round selection from the existing network-kept scores. No GBM was refitted on these full D/H exports, so no improved FDP is inferred.

## 2. Conditional search/competition symmetry

The separate actual-source executable in `calibration_audit/symmetry/` passed 12/12 checks: search label-swap invariance at fixed query coordinates/evidence, equal-evidence area/support, row-order invariance for the competition fixture, strict-score pair wins, decoy tie preference and the +1 q correction. It also documents missing-pair behavior and a pair-max support gate replica.

These checks do not establish end-to-end exchangeability. They do not measure whether the library decoys mimic all false-target mechanisms, learned search coordinates, real sequence overlap, or correlated candidates. Target/decoy-specific tail normalization is intentional and was inspected, not claimed to be label-swap invariant. Detailed source references and test limits are in the archived REPORT.md.

## 3. Known-null stress, fixed seeds 0-9

The rules were written before this stress ran. Production GBM defaults, 12,000 numeric pairs, three paired folds, one candidate per precursor, six features. In 80% of pairs the true target's feature 0 has mean +3. All remaining target/decoy features are independent standard normal. The initial score is feature 0. Known truth is used only to count false target calls.

| Arm | Successful seeds | Accepted targets | Known-null accepted targets | Pooled realized FDP |
|---|---:|---:|---:|---:|
| Initial score | 10/10 | 74,696 | 336 | 0.450 % |
| Production GBM | 10/10 | 75,042 | 347 | 0.462 % |

Maximum observed individual GBM FDP is 0.560%. No seed failures. This easy synthetic model is conservative and does not reproduce the real-data excess. It is neither a new proteomics decoy type nor cross-species validation. No confidence/coverage claim follows from ten synthetic seeds.

## 4. RT refinement: new measured-control plan

`OpenDIAWorkflow.cpp:598` selects only q<=0.01 targets for RT training. `:637` fits the ridge model and `:640-649` uses its full-fit predictions for all library precursors, training targets included. `:658` uses cross-fitted residuals only for the tolerance ratio, then `:609` reruns search. This is an inspected asymmetry in training exposure, not yet a measured FDP effect.

PRECHECK_RTNULL was frozen and pushed as `83eed8ac1a` before jobs 44052409-44052418. The existing `-no_rt_refine` flag is compared with fixed D/H G0 references, original five seeds each. Original corrected-label bounds and no-ID-loss bar are unchanged. The knockout removes beneficial refinement as well as the suspected exposure, so it is not a specific test of cross-fitting alone. See PRECHECK_RTNULL.md and SUBMITTED_RTNULL.md. All ten jobs were RUNNING at 20:19:47 UTC. Six PXD017703 held-out runs remain untouched.

## 5. Fixed pair-stratified calibration screen: neither rule passes

Two rules were specified before results and run on D/H seeds 0-4. A groups whole pairs by whether their network and GBM winner identities agree. B groups whole pairs by whether their maximum network score exceeds the 99th percentile of network-kept decoy scores. Missing-network pairs form an explicit separate stratum. Neither rule uses organism or entrapment annotations. A separate evaluator process reads those annotations only after all selections are complete.

| Set | Recomputed baseline bound | A agreement bound / ID ratio | B network99 bound / ID ratio |
|---|---:|---:|---:|
| D | 1.0747 % | 1.0871 % / 0.9983x | 1.0680 % / 0.9988x |
| H | 1.2765 % | 1.2615 % / 0.9948x | 1.2193 % / 0.9960x |

ID ratios are median N@FDP_c 1% relative to the identically recomputed baseline. Median ID changes are D -38/-26 and H -198/-150 for A/B. Neither reaches <=1% or preserves IDs on both sets. No rule is promoted and there is no retuning. This screen does not reproduce the tool's full fold-standardized q pipeline, so its baseline is the identical offline recomputation, even though its q<=1% N/E totals match the original tables here. Detailed outputs and helper checks are archived in `calibration_audit/paircal/`; the large compressed selected-row file remains scratch-only.

## Next implementation direction

A credible species-independent training candidate must keep the full network-to-GBM stack inside the outer family/pair split, including label selection and negative mining. It should initialize each outer GBM from that outer fold's own network model, evolve training scores only inside its training partition, and score the excluded fold without refitting on it. Independently, RT-center predictions should be compared using models excluding the query family's RT fold for both target and paired decoy, rather than using in-sample predictions for selected targets. Hyperparameter selection also needs an explicit isolation contract.

These are candidate designs, not implemented defaults or claims of universal FDR control. A passing used-data screen still needs a prospectively frozen held-out evaluation and additional species/instrument settings. If decoys remain non-exchangeable with false targets, correcting a training dependency alone need not solve calibration.

## Evidence

Small executable sources, rules, test outputs and source hashes are in `calibration_audit/`. Large D/H exports remain scratch-only and can be regenerated from read-only references. `training/watch_rtnull.sh` is follow-through code, outside the immutable remote freeze. The original input labels, binary and reference outputs were not modified.
