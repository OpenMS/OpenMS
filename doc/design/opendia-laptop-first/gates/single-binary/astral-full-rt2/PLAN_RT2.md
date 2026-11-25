# Frozen plan: Astral full-library runtime round (mem5 + rt2) identity and resource experiment

Date: 2026-11-25. Status: FROZEN PROTOCOL, SUBMISSION STILL REQUIRES PUSHED PREFLIGHT APPROVAL.
Source commit: `bdb0159ff440a46d47990c6a7fb29273cbe2edcc` (`singjc/OpenMS`, `feature/opendia-clean-rewrite`).

The source adds two output-identical rounds on top of the compact Parquet source `8c0a89f053`:

- `a11a95ec05` (mem5): extraction collectors and stage-2 blocks are released before scoring, features are moved (not copied) into the scorer, sibling columns widened in place.
- `bdb0159ff4` (rt2): vectorised scorer forward pass with unchanged per-unit summation order, thread-split decoy covariance with unchanged per-entry row order, compact calibration slot scan with unchanged skip and close rules, and the bucketed m/z index for by-ion counting (returns exactly `MZBegin`).

Local prerequisite: all 19 OpenDIA CTest gates pass, including new scorer tests for forward-pass bit identity and thread-count invariance of the readout. Local raw01 seed0 direct-Parquet runs of `a11a95ec05` and `bdb0159ff4` (and an intermediate rt1 arm) each match all three baseline output files byte-for-byte. The local library has 603k precursors, so its peak RSS does not predict the full-library peak.

## Scope and immutable inputs

WORK=/nfs/research/juan/singj/opendia/astral-rt2-20261125

B=/nfs/research/juan/singj/opendia/astral-full-20261123 (frozen attempt-3 baseline, PQP input)

C=/nfs/research/juan/singj/opendia/astral-parquet-20261124 (compact Parquet result, 15/15 PASS)

CONTAINER=/nfs/research/juan/singj/opendia/opendia-rt2-20261125.sif

The primary runs read the two Parquet caches that the compact experiment's bootstrap wrote, `C/library/search6.opendia.parquet` (SHA256 `1943673979d97b8da9b87cd6b9df0a48ba5519d72af1537bf64f582f05dafaeb`) and `C/library/scoring12.opendia.parquet` (SHA256 `425c9e02a3e9fe0a21ccfce44d63b542e1a4664a6d5a9607b0e349a0310603fa`). They are symlinked unchanged into WORK/library. Direct Parquet input is read only, with no cache writes. Never modify baseline or compact inputs, results, plans or jobs. No tuning, filtering, prediction or library change.

There is no bootstrap. Its purpose (first cache creation) is already met and recorded in the compact result.

## Mandatory root release gates

1. Commit the implementation code, pass local byte identity and all 19 tests, and push to `singjc/OpenMS`.
2. Run board_sync, freeze and push this PLAN before any remote action for this experiment.
3. The container is built from the identical binary (local, Docker and SIF SHA256 equal). Container-native `OpenDIAWorkflow_test`, `OpenDIALibraryParquet_test` and `OpenDIAScorer_test` pass.
4. Push a preflight amendment `PREFLIGHT_APPROVED.md` recording the source commit, plan commit, board_sync evidence, local evidence and exact SHA256 of the container, source provenance, raw inputs, both Parquet caches and every bundle script and plan. `preflight.sha256` lists these with absolute Codon paths for container, caches and raws. `runtime.sha256` lists exactly `run.sbatch`, `evaluate.sbatch`, `submit.sh` and `compare_resources.py` with WORK-relative paths.
5. Transfer the frozen bundle into a fresh WORK. No jobs before the amendment is pushed.

## Fixed job sequence

`submit.sh` uses `/ebi/slurm/codon/bin`, refuses existing job_ids or runs, reserves `submission.lock` atomically, verifies the full preflight once, creates the two cache symlinks and verifies them by hash. A partial submission is inspected, not cleaned or resubmitted. Never cancel or touch other jobs.

Primary array is `0-14%3`, raw=floor(i/5)+1 and seed=i%5. Every run uses the direct WORK Parquet paths, 16 cores, 256G, 24h, standard partition and unchanged algorithm options, the same as the baseline and the compact runs. Each task verifies `runtime.sha256` and `library.sha256` (the two caches only) before running and makes no sacct query. Afterany evaluation follows the primary array, even on failure. If dependency cancellation prevents evaluation, root collects a read-only evaluation manually and keeps all 15 entries.

## Evaluation and interpretation

The primary comparison uses the unchanged `compare_resources.py` against B, with the same three files (`precursors.tsv`, `LFQ_Astral_DIA_15min_50ng_Human_XX.precursors.tsv`, `LFQ_Astral_DIA_15min_50ng_Human_XX.calibration.json`). PASS requires all three SHA256 matches, candidate actual exit 0, baseline exit 0 and exact candidate task accounting COMPLETED 0:0. All absent, failed, cancelled or different runs stay in the fixed denominator of 15. No replacement seeds or selective exclusion.

A second run of the same evaluator against C gives the matched compact resources (`results/comparison_compact.json`). It is a resource reference only. Its identity verdict must agree with the primary one because C is itself byte-identical to B.

Report matched wall seconds and peak RSS from GNU time against both B and C, with node and accounting provenance, and the per-stage seconds and RSS from the OpenDIA logs. The expected direction is a lower peak RSS than C's 64.5 GiB median (scoring-time peak) and a shorter scoring stage. These expectations are not acceptance criteria. Resource comparisons are descriptive only: nodes are assigned by Slurm, and filesystem caches are not controlled. No FDR or quantification claim.

## Local bundle checks

Run `python3 -m unittest discover -s <this-directory> -p 'test_*.py' -v` and `bash -n` on all three shell files. These are local synthetic checks only.
