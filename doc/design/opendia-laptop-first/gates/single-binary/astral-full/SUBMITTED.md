# Submission status, 23 November 2026, 15:38 UTC

Protocol frozen/pushed in89ae34a351, scheduler-only corrections a6ffde1c29/9cb53f57e0, CUDA runtime retry frozen/pushed f209e5d7aa. Nine synthetic acceptance tests pass. No OpenMS code/default changes.

| Phase | Slurm ID | Observed status |
|---|---:|---|
| Full target prediction, A100,16CPU,384GiB |41518555|RUNNING on codon-gpu-012|
| Predicted own-decoy preparation and search6/scoring12 |41518556|PENDING, dependency|
| Three raws x five seeds, max3 concurrent |41518557|PENDING, dependency|
| ID/calibration curves/quant/resource evaluation |41518558|PENDING, dependency|

The repaired job logs both `[RT] provider CUDA device 0 sessions 4` and `[MS2] provider CUDA device 0 sessions 4`. This verifies GPU initialization, not completed prediction, a valid library or benchmark success. Digest eligibility is8,649,692 precursors and164,893,268 transitions before prediction, not final search counts.

Attempt1 jobs41517072-41517075 and every original failed artifact are preserved remotely under `attempt1/`, with accounting and previous frozen scripts. Its CUDA failure remains part of the report. The retry uses the same scientific recipe and pinned libcudart12.8.90, SHA256 c3a75b33af334a3486d197dbd1584a2985183ba4688d237a2be5f2f679329920. No CPU fallback.

Remote work directory `/nfs/research/juan/singj/opendia/astral-full-20261123`. A local bounded-SSH watcher is running as background task303078jlhs and alerts on failed/cancelled jobs or evaluation completion. Job dependencies continue independently of the local client. Complete requested benchmark endpoints are still pending. No ID, calibration, quantification or laptop-memory claim has been made from this submission.
