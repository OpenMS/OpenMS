# PRECHECK_EVSHARE submission record (26 Nov 2026)

- **Frozen precheck:** docs `c8ba28183d`, pushed before any remote action with A and before any job. The manifest is `PRECHECK_EVSHARE.sha256`.
- **Source:** OpenMS `4c6ca367fd` (singjc/OpenMS `feature/opendia-clean-rewrite`, clean tree at build).
- **Container:** `opendia-evshare-20261126.sif`, sha256 `c43320a6e0d350efa29d80db2cafe2470bd026edd5a4fbafe6eba6f5792f191a`.
  - OpenDIA sha256 `733516ac…` locally, in Docker, in the SIF on the laptop and in the SIF on codon. This is the DEV_AB binary.
  - `--helphelp` in the SIF on codon lists `-evidence_share`.
- **New sets:** downloaded by job 43329657 (`evshare/dl.sh`) from the PRIDE FTP. Both zips match the PRIDE SHA-1 (`pride.sha1`: OK, OK).
  - N1 25186.d: 341 MB, 3,176 frames, 337 s.
  - N2 7546.d: 15 GB, 66,775 frames, 7,202 s.
  - Both have 24 diaPASEF windows over 412.5-987.5 m/z, 25 Th wide, as D and H.
  - Only the `analysis.tdf` frame and window tables were read before submission. No search output existed.
- **Work dir:** `/nfs/research/juan/singj/opendia/evshare-20261126/`. `inputs.sha256` there covers the SIF and the three job scripts. `evs_submit.sh` checks it and `pride.sha1` before each `sbatch`.
- **Gate jobs** (no flag): 43339940 (raw03_s0), 43339941 (D s0), 43339942 (H s0).
- **Main jobs:** 35 in all, 43339974-43340008. The full list is in `main_job_ids.txt` on codon.
  - 43339974-43339978: A on raw03, seeds 0-4;
  - 43339979-43339998: A on D, H, N1 and N2, seeds 0-4;
  - 43339999-43340008: G0 on N1 and N2, seeds 0-4.
- **Order:** the gate and main jobs were submitted together. As in PRECHECK_BLEND, the gate only decides which G0 is the reference for raw03, D and H. It does not decide whether A runs.
- **Evaluation:** `evshare/evs_eval.sh`, frozen with the precheck.
