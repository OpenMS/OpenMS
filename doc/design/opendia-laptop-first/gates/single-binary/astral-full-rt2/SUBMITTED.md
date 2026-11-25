# Frozen rt2 experiment submission, 25 November 2026

Status: submitted, no candidate result yet.

- Source pushed: `bdb0159ff440a46d47990c6a7fb29273cbe2edcc`.
- Frozen protocol pushed: `6533717a83405fc2a5fbed6a7db9dbf2b8ba2e12`.
- Hashed preflight approval pushed: `c409adb65be0c53ccda1478375173cef88c27c20`, before any remote write for this experiment.
- Attempt 1 (08:17 UTC) stopped before any remote write: the EBI VPN dropped at 08:11:38 UTC and `codon-slurm-login.ebi.ac.uk` no longer resolved. Log retained (`submission-evidence/transfer_submit_attempt1_dns.log`). A read-only check after reconnecting confirmed that neither the SIF nor the WORK path existed.
- Attempt 2 (08:38 to 08:40 UTC), the same unchanged script: the new SIF was transferred to its unique path and verified by SHA256 (`775a7bcd493820b513abee1b4f81b2f7e692e5862649009e583fd0be17092e1f`). The ten approved files were placed in the fresh `/nfs/research/juan/singj/opendia/astral-rt2-20261125`. The full preflight (container, both Parquet caches, three raws, all scripts and plans) passed before the first `sbatch`, and so did the two-cache `library.sha256`.

| Stage | Slurm ID | Frozen dependency and resources |
|---|---|---|
| Direct-Parquet primary array | `42535684` | array `0-14%3`, 16 CPUs, 256G, standard, 24h |
| Read-only identity and resource evaluation | `42535685` | `afterany:42535684`, 2 CPUs, 32G, standard, 4h |

A read-only check at 08:40:16 UTC found tasks 0 to 2 RUNNING on `hl-codon-cn-059`, `-063` and `-065`, tasks 3 to 14 pending on the array limit, and evaluation pending on its dependency. This is scheduler state only, not evidence of completion or identity.

No baseline or compact inputs, outputs, plans, containers or jobs were changed. No replacement seeds or automatic resubmissions are permitted. All 15 cells stay in the denominator. Final PASS still requires actual executable exit 0, exact three-file baseline bytes and exact COMPLETED 0:0 task accounting. Resource interpretation stays descriptive. No FDR or quantification claim.
