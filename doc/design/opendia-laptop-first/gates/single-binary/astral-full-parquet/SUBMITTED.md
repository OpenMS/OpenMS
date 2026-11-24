# Frozen Parquet experiment submission, 24 November 2026

Status: submitted, no candidate result yet.

- Source pushed: `8c0a89f053284c9aebb620d36e351b14a642d15c`.
- Frozen protocol pushed: `045aae4160719d1a56b127d5ab093bd3efb71705`.
- Hashed preflight approval pushed: `97a803dc71dda2df9cdb7d999791ffd3d7fe20af`, before any remote workload submission.
- New SIF transferred to its unique approved path and SHA256 verified: `035cd7004896d636dc83800230876b4a504fbc247433ea180c02e7ddaffb9caf`.
- Eleven exact approved files transferred into fresh `/nfs/research/juan/singj/opendia/astral-parquet-20261124`.
- Full preflight passed for every listed container, script, provenance, raw and baseline PQP file before first `sbatch`. Both work-library PQP symlinks point to unchanged baseline files. Submission guard retained `submission.lock` and job IDs.

| Stage | Slurm ID | Frozen dependency and resources |
|---|---|---|
| Uncached PQP-to-Parquet bootstrap, raw01 seed0 | `42109936` | 16 CPUs, 256G, standard, 24h |
| Direct-Parquet primary array | `42109937` | `afterok:42109936`, array `0-14%3`, same resources |
| Read-only identity/resource evaluation | `42109938` | `afterany:42109937`, 2 CPUs, 32G, standard, 4h |

The guarded transfer/submission command completed at 18:21:32 UTC. Initial read-only check at 18:21:46 found bootstrap RUNNING on `hl-codon-cn-045`, with primary and evaluation waiting on dependencies. This is scheduler state only, not evidence of analysis completion or identity. Bootstrap rechecks full preflight before OpenDIA. Cold means no preexisting application caches, not evicted filesystem caches.

The preserved `submission-evidence/transfer_submit.log` records exact preflight output and job IDs. `monitor_codon.py` reads only those three jobs and preserves local state snapshots every five minutes. It cannot submit or cancel jobs. No baseline inputs, outputs, plans, container or existing jobs were changed. No replacement seeds or automatic resubmissions are permitted. All 15 primary cells remain in the denominator, and bootstrap is reported separately.

Results remain pending. Final PASS still requires actual executable exit 0, exact three-file baseline bytes and exact COMPLETED/0:0 Slurm task accounting. Resource interpretation stays descriptive. No FDR or quantification claim.
