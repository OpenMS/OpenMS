# PRECHECK_LIBSIZE: how much does library size cost OpenDIA against DIA-NN? (21 Nov 2026, diagnostic)

**Question.** An earlier comparison found that moving from the plain library to the r ≈ 7 entrapment library cost DIA-NN about 10 % of its IDs but cost OpenDIA about 49 %. That was measured on older code. If the penalty is still this lopsided on the current default (matched filter plus RT refinement), the gap is mostly about robustness to a larger search space (score separation, competition, FDR), not about signal extraction. If it has shrunk, the gap is already present on the plain library and the levers are elsewhere.

## Setup (fixed before any result)

- **Libraries:** nested subsets of the r ≈ 7 library `entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp`. Every HeLa (non-entrapment) target is kept. An entrapment peptide is kept when `crc32(unmodified sequence) % 10000 < f x 10000`, so each smaller library is a subset of the next. Each kept target keeps its own `DECOY_` partner, so OpenDIA sees exactly the decoys it saw at r ≈ 7. DIA-NN gets the same targets as a targets-only TSV, cut from `entlib-diann5-20261119/targets_r5.tsv`, and makes its own decoys as before.
  - r0: f = 0, no entrapment (the "plain" library). It gives N only, with no FDP estimate.
  - r1: f = 1/7.0539 (r ≈ 1).
  - r3: f = 3/7.0539 (r ≈ 3).
  - r7: the existing runs, OpenDIA `entlib-rtr-20261121/runs/{D,H}_rtr_s{0,1}` and DIA-NN `entlib-diann5-20261119/{D,H}`.
- **OpenDIA:** the current default, container `opendia-pd10-rtrdef.sif` (sha256 `f334b5fc…`, OpenMS `14c3c54113`), `-threads 8 -seed {0,1}`, no extra flags.
  - The r7 reference runs used `opendia-pd10-rtrefine2.sif -rt_refine`. This is the same code path, since the default flip was byte-checked both ways.
- **DIA-NN 2.2:** flags as in `diann_r5.sbatch` (`--qvalue 1 --gen-decoys --unimod4 --min-pr-charge 1 --max-pr-charge 5`).
- **Scripts:** `libsize.py` (build and eval, sha256 `67bc1334…`) and `libsize_submit.sh` (sha256 `4694e60b…`), in `$C/libsize-20261121`.
  - The builder was smoke-tested on the r ≈ 1.56 library: no orphan decoys, and all tables were consistent.
- **Measures, per run (D, H), tool and level:**
  - HeLa IDs at the tool's own q ≤ 0.01 (`hela_q01`).
  - N at combined FDP 1 % for r1, r3 and r7.
  - The lower and combined entrapment bounds.
  - Retention of the seed-0 r0 HeLa ID set at each level.
  - OpenDIA values are the mean of seeds 0 and 1.

## Decision rule (written before any result)

The penalty at r7 for a tool is `P = 1 - hela_q01(r7) / hela_q01(r0)`. The excess penalty is `X = P_OpenDIA - P_DIA-NN`, computed on D and on H.

1. **Library-size robustness is the lead lever** if X ≥ 15 points on D, or X ≥ 10 points on H. The next step is then the scorer ceiling (step 2) with labels taken across levels. The question there is whether a supervised scorer on our features keeps r0 IDs at r7.
2. **The gap is mostly present without entrapment** if OpenDIA / DIA-NN at r0 (`hela_q01`) is ≤ 0.75 on D, or ≤ 0.85 on H. That share of the gap is extraction or scoring at a fixed search space, so steps 2 and 3 are judged at r0 as well as r7.
3. Both can hold. If neither holds, the old 49 % figure no longer applies and this line of work is closed.
4. **Validity check:** OpenDIA's combined bound at r1 and r3 must be ≤ 1.5 %. If it is not, the r1 and r3 numbers are reported but not used for the trend.

No code change follows from this precheck alone.
