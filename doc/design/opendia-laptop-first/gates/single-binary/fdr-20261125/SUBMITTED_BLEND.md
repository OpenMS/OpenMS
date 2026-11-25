# PRECHECK_BLEND submission record (25 Nov 2026)

- **Frozen precheck:** docs `02caf9aa13`, pushed before the container build and before any remote action. The manifest is `PRECHECK_BLEND.sha256`.
- **Source:** OpenMS `b3391c92e8` (pushed to singjc/OpenMS `feature/opendia-clean-rewrite`, clean tree at build).
- **Container:** `opendia-blend-20261125.sif`, sha256 `8ba71d014510183d503727750a3a2bedf5f402585f91fb545c093e898c19e1b4`.
  - The local, Docker and SIF OpenDIA sha256 are all `5b5ca5aecb283e653b5e4cdb57427d0a3a49773f38b96fa5c2f748404f30276f`. This is the binary of the local blend check.
  - The remote sha256 was verified after transfer.
  - `-gbm_blend` is present in the SIF help.
- **Work dir:** `/nfs/research/juan/singj/opendia/blend-20261125/` (new). `inputs.sha256` there covers the four frozen files and the SIF, and `blend_submit.sh` verifies it before each `sbatch`.
- **Gate jobs** (no flag): 42885674 (raw01_s0), 42885675 (D s0), 42885676 (H s0).
- **B5 jobs** (`-gbm_blend 0.5`), 25 in all:
  - 42885680-42885697: Astral raw 01, 02, 03 × seeds 0-4 (15 jobs);
  - 42885698-42885719: HeLa D and H × seeds 0-4 (10 jobs).
  - The full list is in `main_job_ids.txt` on codon.
- **Order:** the gate and main jobs were submitted together. Under the frozen gate rule, the gate only decides which G0 is the reference (existing outputs, or fresh G0 via `blend_submit.sh freshg0`). It does not decide whether B5 runs.
- **Evaluation:** `blend_eval.sh` (written before any output; it applies the frozen gate rule and the unchanged `eval_pdec.py`).
