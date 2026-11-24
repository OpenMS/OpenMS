# Parquet preflight approval, 24 November 2026

Status: APPROVED FOR THE FROZEN SUBMISSION ONLY AFTER THIS AMENDMENT IS PUSHED.
This is not a result. No Codon candidate workload has run at approval preparation.

## Release gates met

- Source `8c0a89f053284c9aebb620d36e351b14a642d15c` pushed to `singjc/OpenMS`, branch `feature/opendia-clean-rewrite`.
- Frozen plan commit `045aae4160719d1a56b127d5ab093bd3efb71705` pushed before container construction. Fixed bootstrap plus 15 primary cells and unchanged options remain in force. No seeds or acceptance rules changed.
- All 19 local OpenDIA CTest gates pass. Astral raw01 seed0 direct-Parquet and exact final-binary automatic-cache runs each match the three retained baseline outputs byte-for-byte. See [local evidence](../compact-library-20261124/RESULT_LOCAL.md), including initial missing-path failure and system-suspend timing contamination. The local source and cache files were unchanged.
- Local symlink regression confirms automatic cache creation beside the work symlink, without modifying the resolved baseline target cache.
- New container built only after the plan push, at 18:15:32-18:16:53 UTC. Existing baseline SIF was not overwritten.
- Local, Docker and SIF executable SHA256 all equal `83c408e5ad42e4865283a0f8586f18b223139f0f9a1b7d5c9f3ea491c9235eb7`. Docker and SIF `OpenDIA --help` actually exit 0. Container-native `OpenDIAWorkflow_test` and `OpenDIALibraryParquet_test` both pass with a read-only PQP fixture.
- New SIF SHA256: `035cd7004896d636dc83800230876b4a504fbc247433ea180c02e7ddaffb9caf`. Remote path: `/nfs/research/juan/singj/opendia/opendia-compact-parquet-20261124.sif`.
- Sixteen evaluator and shell-contract unit tests pass, and all four submission shell files pass `bash -n`.
- Real `board_sync.py` completed before the frozen-plan commit (four actions, retained evidence). It is run again before this amendment's commit, with output retained in container evidence.

## Exact immutable inputs and scripts

`source_provenance.txt` records source and image attribution. `preflight.sha256` supplies exact hashes for this approval, the frozen PLAN, source provenance, all seven original bundle files, the runtime manifest, the new SIF, the two unchanged baseline PQPs and all three raw files. Raw hashes come from the retained baseline `execution_inputs.sha256`, and submit/bootstrap verify their current bytes. The small `runtime.sha256` contains only the five executable workload/evaluation scripts. It excludes container/PQP/raw hashing per primary task, exactly as frozen.

The baseline search6/scoring12 PQP hashes remain respectively `200dea366ea358f9f97674715c070a6bd4ade8a839d2310fbf214fbf845847c7` and `95a19dec2d169d6b19f973b27eda9763df26fcad69d1a3847115cb76186c49a0`.

Only the fresh `/nfs/research/juan/singj/opendia/astral-parquet-20261124` work directory and the uniquely named new SIF may be created. Full preflight must pass before scheduling. Baseline plans, inputs, outputs and all existing jobs remain untouched. A partial submission is inspected, not automatically cleaned or repeated. Bootstrap output identity, actual executable exits and final exact Slurm accounting remain mandatory. Accounting-lag reevaluation is read-only, never a replacement seed.

Record the pushed amendment commit and job IDs in a later submission record. Do not rewrite this approval or either manifest after submission. No FDR or quantification claim is authorized.
