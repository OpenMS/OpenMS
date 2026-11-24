# Frozen plan: Astral direct-Parquet resource and identity experiment

Date: 2026-11-24. Status: FROZEN PROTOCOL, SUBMISSION STILL REQUIRES PUSHED PREFLIGHT APPROVAL.
Source commit: `8c0a89f053284c9aebb620d36e351b14a642d15c` (`singjc/OpenMS`, `feature/opendia-clean-rewrite`).

Local prerequisite: all 19 OpenDIA CTest gates pass. Raw01 seed0 direct-Parquet and exact-release automatic-cache arms each match all three baseline output files byte-for-byte. Evidence and the retained missing-input failure are recorded in [local validation](../compact-library-20261124/RESULT_LOCAL.md). The exact release-binary arm includes a system suspend, so its GNU wall time is not a benchmark. Container construction and remote execution have not started at protocol freeze.

## Scope and immutable inputs

WORK=/nfs/research/juan/singj/opendia/astral-parquet-20261124

B=/nfs/research/juan/singj/opendia/astral-full-20261123

CONTAINER=/nfs/research/juan/singj/opendia/opendia-compact-parquet-20261124.sif

Use existing `B/library/search6.pqp` (SHA256
`200dea366ea358f9f97674715c070a6bd4ade8a839d2310fbf214fbf845847c7`)
and `B/library/scoring12.pqp` (SHA256
`95a19dec2d169d6b19f973b27eda9763df26fcad69d1a3847115cb76186c49a0`).
Symlink unchanged into WORK/library. Never modify baseline inputs, results, plans,
or unrelated active baseline jobs. No tuning, filtering, prediction or library
modification. Cache creation must write adjacent to the WORK symlink, not its
resolved baseline target. Confirm this in local validation before release.

## Mandatory root release gates

1. Commit implementation code, pass local byte identity and all 19 relevant tests.
2. Root performs board_sync, fills source commit above, freezes and pushes this PLAN.
3. Only after code validation does root build the named container.
4. Push a preflight amendment `PREFLIGHT_APPROVED.md` recording source commit,
   plan commit, board_sync evidence, local test/identity evidence and exact SHA256
   of container, source provenance, raw inputs and every bundle script/plan.
   Supply `preflight.sha256` with those immutable files and baseline PQPs using
   their baseline paths (WORK symlinks do not exist until submit). Include the
   approval amendment itself and `runtime.sha256` in the full manifest. The
   small runtime manifest lists exactly `bootstrap.sbatch`, `run.sbatch`,
   `evaluate.sbatch`, `submit.sh`, and `compare_resources.py`, using WORK-relative
   paths. It must not include PQPs, raw files, container images, or the full
   preflight manifest. Record amendment push commit.
5. Transfer the frozen bundle into fresh WORK. No jobs before the amendment is
   pushed. This plan commit intentionally supplies neither approval nor either hash manifest.

The root owns release, provenance, transfer and eventual submission. This protocol
is frozen before container construction and before any full-library candidate result.
No remote workload is authorized by this plan alone.

## Fixed job sequence

`submit.sh` uses `/ebi/slurm/codon/bin`, refuses existing job_ids/runs/bootstrap,
and reserves `submission.lock` atomically. A partial submission must be inspected,
not automatically cleaned or resubmitted. Never cancel or touch baseline jobs.

First run `bootstrap.sbatch`: full raw01 seed0 OpenDIA with PQP inputs, triggering
first automatic cache creation at WORK/library/search6.opendia.parquet and
scoring12.opendia.parquet. Preexisting cache files are forbidden. Record cold
bootstrap GNU time wall/RSS, actual executable exit, node information and Slurm
accounting separately. Cold here means no preexisting application cache, not an
OS page-cache eviction guarantee. Exactly three outputs must byte-match
B/runs/raw01_s0 before bootstrap succeeds.

Primary array is `0-14%3`, raw=floor(i/5)+1 and seed=i%5. All 15 use the direct
WORK Parquet paths, never PQP. Bootstrap and primary both retain baseline
16 cores, 256G, 24h, standard partition and unchanged algorithm options.
Primary depends afterok on bootstrap and verifies bootstrap actual exit, exact
bytes, small runtime hashes and cache-only hashes. The full preflight manifest
is verified once at submit and once at bootstrap, never per primary task.
`bootstrap/library.sha256` lists only the two generated Parquet caches, not PQPs.
Primary performs no immediate sacct query: accounting propagation lag must not
turn successful bootstrap into failed primary runs. Bootstrap and every primary
still require COMPLETED/0:0 accounting in the final evaluation. Missing accounting
there remains a failure, never a success inferred from afterok or sentinels.
Afterany evaluation follows primary even on failure.
If scheduler dependency cancellation prevents evaluation, root must collect a
read-only evaluation manually, preserving all 15 failed/absent entries.

## Evaluation and interpretation

Compare only these three files against the same raw and seed baseline:
`precursors.tsv`, `LFQ_Astral_DIA_15min_50ng_Human_XX.precursors.tsv`, and
`LFQ_Astral_DIA_15min_50ng_Human_XX.calibration.json`.
PASS requires all three SHA256 matches, candidate actual exit 0, baseline exit 0,
and exact candidate Slurm task accounting COMPLETED with ExitCode 0:0.
A completion sentinel alone is never evidence of success. Accounting includes
step records for MaxRSS inspection, while success uses the exact task record.
All absent, failed, cancelled or different runs remain in the fixed denominator 15.
No replacement seeds or selective exclusion. Bootstrap is separate, not run 16.

Report matched candidate/baseline GNU time wall seconds and peak RSS with node
and scheduler accounting provenance. Resource comparisons are descriptive only:
hardware, filesystem and cache variability prevent causal speedup claims. Do not
pool bootstrap with primary timings. Primary speed/memory/identity is sufficient.
No new FDR or quantification claim. Existing eval_full.py label correction is
separate and applies to future analyses only. Do not rerun historical results.

## Local bundle checks

Run `python3 -m unittest discover -s <this-directory> -p 'test_*.py' -v` and
`bash -n` on all four shell files. These are local synthetic checks, not evidence
that code identity, container validation, or remote execution gates have passed.
