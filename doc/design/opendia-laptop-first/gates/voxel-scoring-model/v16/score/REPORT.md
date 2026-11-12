# Stage D score mechanism audit

## Conclusion

**Deep Stage D losses predominantly occur in rows never admitted as positives, strongly associated with full-versus-linear window disagreement, not in admitted positives trained at another apex.** The association is strong on both raws but does not establish causation: positive admission and the later score are selected by related models. A large, directly measured training/inference discrepancy also affects negative windows. This motivates one selector-consistency experiment, not residual removal or an ID-recovery claim.

## Exact reconstruction and integrity

With the approved original **10 Torch threads**, both raws reproduce every saved D11F full and linear score with **maximum absolute difference 0.0** and **W-row agreement 1.0**. Spearman is 1.0 / 0.9999999999999999 over 1,270,510 / 1,275,383 population rows. The original `train_v4.py` source was imported without editing it or calling `run()`. Original seeds, normalization, 100 steps, rounds, validation hashes, and OMP settings were retained. No intervention model was fitted.

All 18 fold/round checkpoints, positive/negative train/validation memberships, and 12 prior-model full/linear selector snapshots were saved label-blind under `original_threads/`. Their hashes and completion timestamps precede evaluation joins. Source/settings hashes were frozen before execution. All 16 frozen input hashes remain unchanged. Three original helper tests, three mismatch-helper tests, and the earlier brute-force pairwise-reversal test passed. Population/fold isolation and expected D/E totals were checked.

**Preserved failed attempt:** four-thread raw01 reconstruction failed (rho 0.93682994, max score error 1520.55664, corrected W agreement 0.95764909). Its checkpoints were never attributed or joined to evaluation labels. `AMENDMENTS.md`, `raw01/`, and `REPORT_before_original_threads.md` preserve that failure. The approved retry is separate in `original_threads/`, with successful `EXIT` and `ANALYSIS_EXIT` both 0.

## 1. Training coverage: the main new result

Each precursor is eligible for two of the three models. The table counts **distinct W-selected winners**, pooled over those two eligible models, as positives in training **or validation**. “Lost” means admitted in r0/r1 but absent from r2 in both eligible models. “Retained” means present in at least one r2 selection. Exact-row coverage is distinguished from any window of the same precursor.

| Raw/group | N | Exact row never admitted | Admitted then lost | Retained r2 | Deep-score never / lost / retained |
|---|---:|---:|---:|---:|---|
| 01 D | 866 | 397 | 60 | 409 | 268 / 18 / 1 |
| 02 D | 929 | 371 | 51 | 507 | 222 / 7 / 1 |
| 01 E | 1413 | 0 | 0 | 1413 | 0 / 0 / 0 |
| 02 E | 1362 | 0 | 0 | 1362 | 0 / 0 / 0 |

“Deep” is held-out full score below the same-fold W-decoy 90th percentile, not a fitted threshold. Thus **268/287 (93.4%) and 222/230 (96.5%) deep-D rows were never admitted**, while only one deep-D row per raw was retained. Conversely, the retained D groups have median full-score decoy percentiles 0.999646/0.999714. They remain Stage D because the empirical FDP cut is stringent, not generally because their scores are ordinary-null-like.

This is not merely a window-identity bookkeeping issue: precursor-level never/lost/retained counts are 394/61/411 and 366/54/509. Only **3/5 D precursors** were ever admitted at another row while the eventual W row was never admitted. Any positive admission more than 10 seconds from the final W apex occurs for only **5/6 distinct D precursors** across all rounds.

Training-only membership is materially different from training-or-validation membership. Exact-row never/lost/retained counts are **462/52/352 and 450/45/434 for D**, versus **137/62/1214 and 121/71/1170 for E**. Validation rows do not contribute to the final readout fit, although validation chooses the checkpoint. The full per-role and per-model tables prevent conflating these roles. Decoys are negative-only by construction, so their zero positive coverage is expected.

## 2. Full-selector training versus linear-selector inference

Training rounds r1/r2 choose per-precursor **full-score maxima**, then apply pooled q thresholds and select full-max decoys as negatives. Inference W chooses **linear-contribution maxima** but assigns the full score. We measured two distinct discrepancies, rather than attributing cross-model differences entirely to the selector.

### Same prior-round model, changing only the selector

Below are r2 selection exposures summed over the two eligible models per precursor. Every decoy full winner becomes a negative before the validation split. The selected-positive rows are a conditional subset, not an unbiased comparison population.

| Group | raw01 eligible exposures | Exact mismatch / >10s | raw02 eligible exposures | Exact mismatch / >10s |
|---|---:|---|---:|---|
| D | 1732 | 428 / 425 | 1858 | 403 / 399 |
| E | 2826 | 1 / 1 | 2724 | 0 / 0 |
| allH | 95340 | 71196 / 70674 | 96588 | 73131 / 72605 |
| ENT | 92070 | 71952 / 71405 | 93570 | 74740 / 74162 |
| DEC | 185994 | 143940 / 143026 | 188770 | 148735 / 147752 |

The negative example distribution therefore differs at an apex separation >10s for **76.9%/78.3% of decoy exposures**. This is a measured selection-distribution discrepancy, not evidence that those decoys are incorrectly labeled.

However, among D exposures actually admitted as r2 positives, full-versus-linear mismatches are only **3/745 and 3/931**, with >10s mismatches **3/745 and 1/931**. Against the final held-out W row, admitted D positive mismatches are **3/745 and 5/931**, with the same >10s counts. Therefore, **“the model trained D positives at the wrong apex” is not supported as the main explanation**. The discrepancy largely occurs among omitted candidates and negatives.

### Final held-out model, full winner versus W winner

| Group | raw01 exact / >10s mismatches | raw02 exact / >10s mismatches |
|---|---|---|
| D | 283 / 281 of 866 | 231 / 224 of 929 |
| E | 1 / 1 of 1413 | 0 / 0 of 1362 |
| allH | 36296 / 36032 of 47670 | 36535 / 36257 of 48294 |
| ENT | 36759 / 36487 of 46035 | 37079 / 36797 of 46785 |
| DEC | 73254 / 72763 of 92997 | 73788 / 73289 of 94385 |

Among never-admitted D rows, final >10s selector disagreement affects **261/397 and 216/371**. Among retained D rows it affects only **1/409 and 1/507**. Of the deep-D rows, **251/287 and 212/230 simultaneously were never admitted and have >10s disagreement**. These are strong associations, not an isolated causal effect of selection: admission, the nonlinear score and full winner are mutually dependent.

## 3. Score decomposition and tail-specific controls

Full, linear and nonlinear (= full minus linear) are evaluated on identical W rows, against identical same-fold W decoys, with tie-corrected percentiles/AUC. Controls cover D/E/allH/ENT/DEC, and every D row has same-fold controls.

| Same-fold AUC | raw01 linear / full | raw02 linear / full |
|---|---|---|
| D vs all decoys | 0.9468 / 0.7355 | 0.9549 / 0.7906 |
| D vs full-score top-1% decoys | 0.4980 / 0.5356 | 0.5574 / 0.6329 |
| E vs full-score top-1% decoys | 0.9727 / 0.9984 | 0.9803 / 0.9990 |
| allH vs all decoys | 0.5603 / 0.5134 | 0.5572 / 0.5170 |
| ENT vs all decoys | 0.5407 / 0.4907 | 0.5398 / 0.4958 |

The residual suppresses a deep-D subgroup but also counteracts linear target/decoy bias and separates accepted E from the difficult tail. The tail controls are explicitly **selected by the full score**, so these are conditional diagnostics, not alternative-model performance claims. Fixed full-percentile strata yield D nonlinear AUC 0.4603/0.5171 versus ENT 0.4619/0.4623, not a consistently D-specific residual defect.

D crosses **upward** through the decoy 99th percentile for 147/192 rows and downward for 17/22, despite >=0.10 percentile losses for 253/216 rows. Actual pairwise comparisons against the same individual decoys give D down/up inversion fractions **0.2420/0.0308 and 0.1928/0.0286**. E inversions are negligible. Full D/E/allH/ENT/DEC counts and pairwise rates are retained in the reference and pairwise tables. This mixed effect rejects blanket residual removal.

## Limitations and one proposed next experiment

- Stage D means below the **empirical FDP_c 1% cut**, not necessarily failing operational pooled q<=0.01. The coordinator's threshold audit supplies exact overlap counts. DIA-NN agreement is not ground truth; raw02 is descriptive replication, not an untouched holdout.
- Positive selection is endogenous. Never-admitted status can be a consequence of weak/ambiguous signal, not the cause of low scores. Merely being retained in an eligible model does not mean the separate held-out scorer saw that row. Broad per-feature 1st–99th percentile support bounds are descriptive only: D/E have median zero out-of-range features overall, but omitted D has median one. This does not prove multivariate support or identify a causal feature.
- **Proposed single recovery experiment:** freeze a selector-consistent version of the same scorer: in r1/r2 choose windows by the same linear W rule for both target candidates and decoy negatives, retain full scores for pooled q/readout and inference, and keep features, seeds, objective, validation, thresholds, and final W/T pipeline otherwise unchanged. Do not tune or select rows using D/ENT labels. This directly tests the measured selector/negative-distribution discrepancy, rather than fitting a label-driven rescue.
- This differs from V4R's sibling-ranking penalty: no new rank loss is proposed, only consistency between the training examples and inference selector. It is also distinct from V4/V5 feature additions, V7 cross-run evidence, V9–V11 native counts, V12 profile correction, V13 deconvolution and V14/V15 amplitude changes. Those failures are not overridden by this diagnosis.
- The proposed experiment is **not run here** and could lose positives or worsen tail calibration. Preserve the existing global pass rule against D11FWT, and audit newly gained/lost H/ENT/DEC and empirical FDP at operational q1%. No increased-ID or causal-recovery claim is supported yet.

## Output map

All outputs are under `stageD_v16/score/`. Main validated reconstruction tables are `original_threads/coverage_categories.tsv`, `coverage_by_model_round.tsv`, `selector_mismatch_pooled.tsv`, `admitted_row_mismatch.tsv`, `final_selector_mismatch.tsv`, `distinct_winner_mismatch.tsv`, `coverage_selector_overlap.tsv`, and `feature_support.tsv` with `feature_names.tsv`. Raw-level `coverage_winners.tsv.gz`, `coverage_arrays.npz`, `selector_rows.tsv.gz`, and `admitted_row_mismatches.tsv.gz` retain identities and denominators. `SUMMARY.txt`/`OVERLAP.txt` provide compact readouts. Checkpoints, membership/selector NPZ files, hashes, timestamps, tests and logs are retained per raw.

Earlier validated reference tables remain `reference_score_summary.tsv`, `reference_control_auc.tsv`, `reference_score_strata.tsv`, `reference_winners_raw*.tsv.gz`, and `pairwise_reversal_summary.tsv`/`pairwise_reversals_raw*.tsv.gz`. No reference inputs, outside documents, or main-repository commits were changed.
