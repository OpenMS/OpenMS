# rt2 preflight approval, 25 November 2026

Status: APPROVED FOR THE FROZEN SUBMISSION ONLY AFTER THIS AMENDMENT IS PUSHED.
This is not a result. No Codon workload for this experiment has run at approval preparation.

## Release gates met

- Source `bdb0159ff440a46d47990c6a7fb29273cbe2edcc` pushed to `singjc/OpenMS`, branch `feature/opendia-clean-rewrite`.
- Frozen plan commit `6533717a83405fc2a5fbed6a7db9dbf2b8ba2e12` pushed. The 15 primary cells, unchanged options and acceptance rules remain in force.
- All 19 local OpenDIA CTest gates pass. Raw01 seed0 direct-Parquet runs of mem5, rt1 and rt2 each match the three retained baseline outputs byte-for-byte. See [local result](../runtime-round-20261125/RESULT_LOCAL.md).
- Container built from the clean, pushed source tree. Local, Docker and SIF OpenDIA SHA256 all equal `1587283d3c6fe1d705d145576b44d60ec2d9803a4a8eca1061419dc7515b6794`, the rt2 identity-run binary. Docker and SIF `OpenDIA --help` exit 0. Container-native `OpenDIAWorkflow_test`, `OpenDIALibraryParquet_test` and `OpenDIAScorer_test` pass. The container was built before the plan push. `source_provenance.txt` records this.
- New SIF SHA256 `775a7bcd493820b513abee1b4f81b2f7e692e5862649009e583fd0be17092e1f`. Remote path `/nfs/research/juan/singj/opendia/opendia-rt2-20261125.sif` does not yet exist. No existing SIF is overwritten.
- The two Parquet caches written by the compact experiment's bootstrap were verified read-only on Codon against `bootstrap/library.sha256` on 25 November (both OK). They are used unchanged and only read.
- Sixteen evaluator and shell-contract unit tests pass, and all three submission shell files pass `bash -n`. `compare_resources.py` is byte-identical to the frozen compact evaluator.
- `board_sync.py` ran before the plan commit (4 actions, evidence retained) and runs again before this amendment's commit.

## Exact immutable inputs and scripts

`preflight.sha256` lists this approval, the frozen PLAN, source provenance, the evaluator tests, `runtime.sha256`, the four runtime scripts, the new SIF, both Parquet caches at their compact WORK paths, and all three raw files (hashes from the retained baseline `execution_inputs.sha256`). `submit.sh` verifies it once before the first `sbatch`. `runtime.sha256` lists only `run.sbatch`, `evaluate.sbatch`, `submit.sh` and `compare_resources.py`. Each primary task verifies it and the two-cache `library.sha256` that `submit.sh` writes and checks.

Only the fresh `/nfs/research/juan/singj/opendia/astral-rt2-20261125` work directory and the uniquely named new SIF may be created. Baseline and compact plans, inputs, outputs and jobs remain untouched. A partial submission is inspected, not cleaned or repeated. Actual executable exits, exact three-file baseline bytes and exact COMPLETED 0:0 task accounting remain mandatory. Accounting-lag reevaluation is read-only, never a replacement seed.

Record the pushed amendment commit and job IDs in a later submission record. Do not rewrite this approval or either manifest after submission. No FDR or quantification claim is authorized.
