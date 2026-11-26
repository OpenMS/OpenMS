# SUBMITTED_RTNULL

26 November 2026, 20:18 UTC. Development-only causal control on used D/H, not held-out validation.

- Frozen plan/scripts/manifest were committed and pushed as `83eed8ac1a` before submission. Remote `sha256sum -c PRECHECK_RTNULL.sha256`, container and library checks all passed.
- Existing binary and reference outputs are unchanged. Remote workspace: `/nfs/research/juan/singj/opendia/rtnull-20261126`.
- All ten cells were RUNNING at 20:19:47 UTC. No results are claimed yet.

| Set | Seed | Job |
|---|---:|---:|
| D | 0 | 44052409 |
| D | 1 | 44052410 |
| D | 2 | 44052411 |
| D | 3 | 44052412 |
| D | 4 | 44052413 |
| H | 0 | 44052414 |
| H | 1 | 44052415 |
| H | 2 | 44052416 |
| H | 3 | 44052417 |
| H | 4 | 44052418 |

The first archive upload to `/homes/sing` failed before workspace creation or submission. A direct small-archive transfer as pst_prd then succeeded and verified every frozen checksum. No jobs were submitted by the failed transfer.

Persistent local monitor task `386640n0mt` watches only these ten IDs, then fetches precursor tables/logs/exit files, verifies remote file SHA256 checksums, and runs the frozen local evaluator. Its files are under `/home/singjc/.jcode/scratch/rtnull-20261126-monitor/`. The monitor does not resubmit failed jobs or change parameters. Its source is archived with the calibration diagnostic, outside the immutable freeze.
