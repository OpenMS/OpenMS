# Candidate-window and sibling-availability precheck

Status: implementation and fabricated-fixture validation only, 2026-11-23. No real-run statistics have been computed for this audit. No production changes. This document must be committed and pushed with the calibration precheck before any real-data evaluation. This task does not authorize committing, pushing, rerunning production, or choosing a scoring arm.

## Questions and limits

1. Does the search library offer targets and their **own generated decoys** the same sibling multiplicity?
2. In the exported candidate population, how many windows and distinct sibling precursors are available on each side of each own-decoy pair?
3. Is any difference also visible in a frozen, label-blind low-reference-score proxy region?

All-target versus all-decoy differences are descriptive. True biological signal can give targets more peaks, stronger co-elution, and more observed charge states. These differences alone do not establish failure of null exchangeability or justify disabling a feature. Library multiplicity differences are structural opportunities, not proof of FDR inflation. Candidate availability is already affected by population selection.

The optional proxy region requires **both** target and own-decoy best kept scores in the same frozen no-GBM reference run to be at most `ln(2)`, following `gbmnull.py`. This is not known-null truth. Selection on scores can induce bias. The script emits strata, not p-values, significance labels, policy decisions, or an exchangeability verdict. Missing reference scores are a separate stratum, never silently discarded. No entrapment labels, protein accessions, composition keys, or calibration FDP results enter these diagnostic rules. Own-decoy TRAML pairing is distinct from original-generator pairing required for paired entrapment FDP estimation. No composition-only pairing is used here.

Seeds of the same raw file are correlated technical perturbations, not independent experiments. Keep run, arm, raw file, and seed separate. Do not pool seeds to manufacture sample size. Do not tune the proxy threshold, RT tolerance, or diagnostic choice after seeing calibration results.

## Source/schema inspection

The docs checkout's current `src/topp/OpenDIA.cpp` is the newer front-end and does not define these standalone outputs. The actual current writer inspected is `/home/singjc/Documents/github/OpenMS-r1-sparse-experiment/src/topp/OpenDIA/OpenDIA.cpp`, checkout HEAD `3b1d90ed6786b4712b62e587d74c4a9e1a30a6c9`, `writeRun_` at line 366 and `writeFeatures_` at line 385. Its sibling block is in `source/OpenDIAWorkflow.cpp` near line 784. Its `source/OpenDIALibrary.cpp` selects `PEPTIDE.MODIFIED_SEQUENCE AS FullPeptideName` at line 162 and assigns `c.sequence = modified.empty() ? peptide : modified` at line 260. `source/OpenDIAWorkflow.cpp` line 224 assigns `p.sequence = c.sequence`. Thus the loaded sequence is modified when provided, with unmodified fallback. The exact PQP inputs must have nonempty modified sequence for this minimal audit, which rejects empty values.

Also inspected the historical source in this repository with `git show bcb156e036:<path>`. The remote `gbmnull.sbatch` identifies the container as `opendia-pd10-gbmdef.sif`, OpenMS `bcb156e036`. This is a source/container assertion from the submission script, not an independently authenticated binary digest.

Relevant files under `src/topp/OpenDIA/` at that revision:

- `OpenDIAWorkflow.cpp`, `writeRun_`: `<run>.precursors.tsv` exports only kept per-precursor windows, with `run,candidate_row,precursor_id,is_decoy,sequence,charge,rt_apex,rt_start,rt_end,area,transition_support,score,q_value`. It is **not a candidate-count table**. The combined `precursors.tsv` is further target/q filtered and must not be used for this audit.
- `OpenDIAWorkflow.cpp`, `writeFeatures_`: optional `-features_out` produces `<run>.features.tsv`, with `candidate_row,precursor_id,is_decoy,fold,rt_apex`, feature columns, `score,q_value`. Every `result.windows` row is written, including final NOT_KEPT rows. Use this table, one run per invocation. Candidate IDs need not be contiguous.
- `source/OpenDIAWorkflow.cpp`, population construction: `result.windows` is already restricted by pair-level fragment support. Counts from features are **post-population candidate-window counts**, not raw nomination, extraction, or pre-support counts. A library precursor with no feature rows has zero exported candidates, but that does not identify the upstream reason.
- `source/OpenDIAWorkflow.cpp`, `AnalysisLibrary::read`: own-decoy key strips exactly one leading `DECOY_` from decoy `TRAML_ID`. Family/fold assignment is not the sibling key.
- `source/OpenDIAWorkflow.cpp`, sibling block: siblings group by side plus `p.sequence` (modified sequence), across library precursors. `include/OpenDIAStage2.h` documents `s_n = group_precursors - 1` and `s_co` as distinct other precursor IDs within inclusive 3 seconds. This includes windows not ultimately kept. `s_best` uses first-pass scores, falling back to first-pass `smin` for unkept windows. Final exported scores cannot reconstruct those first-pass values.

The audit reads PQP `PRECURSOR`, `PRECURSOR_PEPTIDE_MAPPING`, and `PEPTIDE.MODIFIED_SEQUENCE` only. One peptide mapping per precursor and globally unique TRAML IDs are required. PQP modified-sequence strings must match the workflow's loaded sequence normalization. Exported `s_n` and `s_co`, when available, are checked against reconstructed values and mismatches stop the audit. Where those columns are absent, sequence normalization remains a provenance limitation. Library records without transitions may not enter the actual workflow library, so use the exact prepared search PQP and authenticate membership before interpreting counts.

## Artifacts confirmed by metadata inspection only

Paths were obtained from `/home/singjc/.jcode/scratch/dia2/gbmnull_grade.sh`, the submission script, and remote directory listings. No result table was read, counted, graded, or summarized.

Codon root: `/nfs/research/juan/singj/opendia`.

- Confirmed Astral search library: `gap-20261114/home/lib/library.pqp`.
- Confirmed HeLa predicted-decoy library: `entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp`.
- Confirmed Astral no-GBM example directory: `gbmnull-20261123/runs/raw03_L1_s0`. Its complete top-level file listing contains only `LFQ_Astral_DIA_15min_50ng_Human_03.precursors.tsv`, the corresponding `.calibration.json`, and combined `precursors.tsv`. **No features table is available in this inspected directory. Exact candidate-window counts and per-window sibling availability are therefore blocked for this run.**
- `gbmnull_grade.sh` references `astral-gbmfix-20261123/runs/raw03_s<seed>` and other `gbmnull-20261123/stat/...` layouts. Their feature-table availability has not been checked. No blanket availability claim is made for other runs or HeLa.
- The remote submission uses scoring library `rtwin2-20261116/lib/library12_full.pqp`. Sibling membership comes from the search library, not this scoring library. Existence/content of this scoring library was not inspected.
- Remote `gbmnull-20261123/gbmnull.py` is available as code and uses own-decoy TRAML pairing and the frozen no-GBM low-score proxy. It is not executed here.

## Implemented outputs

`audit_asymmetry.py` uses Python standard library only, read-only SQLite, TSV or gzip TSV input, and exclusive-create JSON output. It records input SHA256 digests and explicit precheck provenance. It does not alter input files.

Per library precursor: sibling multiplicity excluding itself, exported candidate-window count including NOT_KEPT, distinct siblings with any exported candidate, number of own windows with any co-eluting sibling, and maximum co-eluting distinct siblings (`null` when there are no own windows). Counts use distinct precursor IDs, not duplicate windows, for sibling availability. Groups retain target/decoy separation and exact modification state. The implementation intentionally counts all other precursors in the group, matching code, rather than assuming exactly one precursor per charge state.

Per unambiguous own-decoy pair: target and decoy values and signed target-minus-decoy differences for the first four metrics, plus frozen reference stratum. Unpaired precursor IDs are reported separately, not silently dropped. Missing reference information does not remove a pair. Duplicate candidate IDs, inconsistent sides, ambiguous library mappings, nonfinite required numeric fields, and available sibling-feature disagreements fail closed. This minimal script retains window rows in memory and compares within sibling groups, so memory and runtime should be checked before any later large-library deployment.

`audit_asymmetry.py --library-only` is the frozen fallback when features are absent. It emits only structural library sibling multiplicities and own-decoy paired differences. Every candidate-related field and difference is JSON `null`, never zero, and `candidate_population` is `unavailable_library_only`. No reference score table is permitted in this mode. The same precheck commit/push embargo applies to library-only results. This mode does not authenticate transition-supported workflow membership and cannot diagnose observed co-elution or candidate-window asymmetry.

The exact structural input manifest is frozen to these two libraries, with no substitution or path selection based on outcomes:

```text
/nfs/research/juan/singj/opendia/gap-20261114/home/lib/library.pqp
/nfs/research/juan/singj/opendia/entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp
```

After the precheck is committed and pushed, execute one library-only invocation per path and preserve both complete JSON outputs with SHA256 provenance. This is a structural diagnostic only, not a null-exchangeability test.

## Frozen execution protocol

1. Finish synthetic tests and review this plan with PRECHECK_CALIB. The authorized maintainer must commit and push the precheck before real statistics. Record full commit ID, remote branch visibility, exact raw/run/arm/seed manifest, binary and library digests, and feature-export completeness. The CLI push flag is an operator attestation, not a remote Git verification.
2. Inventory existing complete `.features.tsv[.gz]` files without reading statistical values. If absent, stop the candidate-window branch with an artifact gap, but run the frozen structural-only branch after precheck publication. Do not substitute per-precursor winner rows or the target-only combined table, and do not infer missing candidate windows. A future authorized run with `-features_out` is a separate action requiring approval and reproducibility checks, not part of this task.
3. Freeze no-GBM reference and exact same raw/seed/library matching before evaluating the proxy. Do not use the arm under diagnosis to redefine its own proxy region. Hash both inputs. The script does not infer or validate raw/arm/seed identities, so this manifest check is mandatory.
4. Execute one audit per run. Report the full library/unpaired population first, complete own-decoy paired differences second, and all three reference strata separately if requested. Never select only pairs with two observed candidates. Keep structural library opportunity separate from observed-window availability. This minimal implementation emits records, not pooled curves or hypothesis tests.
5. Any future inferential null-exchangeability test needs its own preregistered null definition, clustering unit and error control. These descriptive records cannot establish a null or a production policy.

Future real-data commands, intentionally not executed:

```bash
python3 audit_asymmetry.py --library SEARCH.pqp --library-only --output STRUCTURAL.json \
  --precheck-commit FULL_40_HEX_COMMIT --precheck-pushed
python3 audit_asymmetry.py --library SEARCH.pqp --features RUN.features.tsv \
  --reference FROZEN_NO_GBM_RUN.precursors.tsv --output AUDIT.json \
  --precheck-commit FULL_40_HEX_COMMIT --precheck-pushed
```

Without both precheck arguments, real-data CLI execution is blocked before reading inputs. `--synthetic` is exclusively for fabricated fixtures and is not authorization to inspect real data. Importing the Python function is not a security boundary. The execution embargo is procedural as well as a CLI guard.

## Synthetic validation

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover \
  -s doc/design/opendia-laptop-first/gates/single-binary/calib \
  -p 'test_asymmetry.py' -v
```

Validation result on 2026-11-23: all 5 unittest methods passed (including their malformed-input subcases), with no real-data inputs. `git diff --check` also passed. The files remain uncommitted by this task.

Fixtures are constructed in temporary directories, with no production inputs. Coverage includes multiple windows, final NOT_KEPT inclusion, same-sequence charge siblings, independently different decoy sequences, RT boundary at exactly 3 seconds, library members with zero windows, unpaired precursors, low/outside/missing reference strata, gzip input, invalid schema/side/duplicate/nonfinite/sibling values, the precheck gate, synthetic provenance and refusal to overwrite outputs.
