# Running and resuming the frozen benchmark

Use Codon `become pst_prd`, work directory `/nfs/research/juan/singj/opendia/astral-full-20261123`. Every SSH must have an outer timeout60 and BatchMode/ConnectTimeout20/ControlMaster=no/ControlPath=none.

1. Verify `PRECHECK_ASTRAL_FULL.sha256` and the committed precheck revision before `bash submit.sh`. Never resubmit over `job_ids.txt`, `library`, or run directories.
2. `job_ids.txt` records prediction, preparation, search-array and evaluation job IDs in that order. Query `squeue` and `sacct` for those IDs. A dependency failure is not permission to change the recipe. Preserve accounting and logs before retry decisions.
3. Prediction and CPU decoys have separate `library/*.time.txt` and logs. `library/structure.json` authenticates master/search/scoring hashes, counts, pairing and fragment/sibling limitations. No `.complete` file is a substitute for reading exit/accounting state.
4. The search array index maps raw = floor(index/5)+1 and seed=index%5. The exact 15 outputs and their original exit codes remain in the denominator, including failed or absent output. Do not replace a failed seed with the best of additional seeds.
5. After all searches terminate, evaluation creates run status, full calibration grid/curve plots and quantification JSONs for all targets and human-only grading at fixed q=.005,.01,.02. If not all15 succeeded, quantitative completeness is explicitly blocked rather than silently analyzing an easier subset.
6. Fetch summaries and original logs, compare all checksums, publish a result with every failure, resource exceedance and statistical limitation. Run board_sync before every docs commit. Do not claim FDR control or quantitative accuracy from descriptive curves or repeatability.
