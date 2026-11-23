# Astral on the current default, laptop profile (23 Nov 2026)

Diagnostic, no pass rule. Raw 03 is in the set, as in the 16 Nov snapshot (the README holds raw 03 back for design; nothing here was tuned on it).

## Setup
- Inputs as 16 Nov: arm-E search library (`gap-20261114/home/lib/library.pqp`, md5 `e0461beb…`), scoring library `library12_full.pqp`, self-calibrating, raws 01-03, seeds 0-4.
- Laptop profile: `-threads 16`, slurm `-c 16 --mem=32G`, codon AMD EPYC 9555 nodes (shared). Cost from `/usr/bin/time -v`.
- Evaluation: frozen `gates/single-binary/workflow/eval_run.py`.
- Batch 1: OpenMS `46243afae4`, `opendia-pd10-rename.sif` (sha256 `4614e25e…`), `astral_laptop.sbatch` (`EVAL_batch1.txt`, `COST_batch1.txt`).
- Batch 2: OpenMS `88e5da3023` (GBM fix below), `opendia-pd10-gbmfix.sif` (sha256 `2270438a…`), `astral_gbmfix.sbatch` (`EVAL_gbmfix.txt`, `COST_gbmfix.txt`).

## Batch 1: 5 of 15 runs crashed (FAIL, robustness)
raw01_s2, raw02_s0, raw02_s1, raw03_s0 and raw03_s1 stopped with exit 8: `pooledPairedTDC: a fold has zero decoy spread`.

**Cause.** The GBM rescore runs 3 rounds. A round can leave a fold whose model keeps no tree (logged `trees 0/…`) or is otherwise degenerate. All rows of that fold then get one score, so the fold's decoy MAD is 0. For rounds 1 and 2, the next round's start catches this and keeps the previous scores. After round 3 there is no next round, so the throw came from the final q-value step instead. All 5 crashes were this round 3 case.

**Fix (OpenMS `88e5da3023`, local).** After the last round, run the same check and fall back to the round 2 scores. A run that completed before passes the check unchanged.
- HeLa byte-check (codon `gbmfix-20261123`, `BYTECHECK_gbmfix.txt`): D r ≈ 1.56, D r ≈ 7 V7M6 and H r ≈ 7 are **SAME on 9/9 files**.
- Astral: the 10 runs that completed in batch 1 are **byte-identical** in batch 2 (all output TSVs). The 5 crashed runs now complete; each logs `GBM rescore r3: discarded …; keeping round 2 scores`.

## Identifications, batch 2 (15 of 15 runs)

| Raw | N @ FDP_c 0.5 / 1 / 2 % (median of 5) | N @ 1 % range | 16 Nov median @ 1 % | DIA-NN (.raw) @ 1 % | share of DIA-NN | change vs 16 Nov |
|---|---|---|---:|---:|---:|---:|
| 01 | 1770 / **2124** / 2331 | 1950-2275 | 1804 | 2204 | **96 %** | +17.7 % |
| 02 | 1890 / **2144** / 2441 | 1301-2215 | 1909 | 2356 | **91 %** | +12.3 % |
| 03 | 2003 / **2246** / 2498 | 2160-2343 | 1907 | 2259 | **99 %** | +17.8 % |

On 16 Nov the share was 81-84 %. The HeLa levers adopted since then (10 ppm trace tensor, matched filter, RT refinement, sibling and GBM rescoring) carry over to Astral.

**Own q ≤ 0.01 calibration got worse (the main cost):**

| Raw | IDs at own q ≤ 0.01 (median) | FDP_c at own q ≤ 0.01: median (range) | runs above 1 % |
|---|---:|---|---:|
| 01 | 2214 | 1.17 % (0.81-1.77) | 4 of 5 |
| 02 | 2072 | 0.95 % (0.77-1.54) | 2 of 5 |
| 03 | 2275 | 1.32 % (0.79-1.45) | 4 of 5 |

- Overall, 10 of 15 runs are above 1 % FDP_c at their own q ≤ 0.01, with a maximum of 1.77 %. On 16 Nov no run was above 1 %. DIA-NN's own q ≤ 0.01 on the same benchmark is 1.44 % (raw 01) and 1.12 % (raw 02). This matches the HeLa H loosening after the GBM (RESULT_GBT_H), but it is larger here.
- **Seed instability:**
  - raw02_s4 reaches only 1301 at FDP_c 1 %, against at least 2128 for the other raw 02 seeds.
  - raw01_s4 reaches only 497 at 0.5 %, against at least 1740 for the other raw 01 seeds.
  - raw02_s1 reports only 1412 at its own q ≤ 0.01, although it reaches 2215 at FDP_c 1 %.

## Cost at 16 threads / 32 GB

| Raw | batch 2 wall, min (s0-s4) | fastest completed run (both batches) | peak RSS |
|---|---|---:|---:|
| 01 | 7.1, 16.1, 16.9, 28.8, 6.7 | 5.8 min | 16.9-17.2 GB |
| 02 | 27.1, 7.2, 12.0, 7.3, 18.4 | 6.7 min | 16.2-16.6 GB |
| 03 | 12.2, 7.2, 10.7, 11.7, 12.8 | 6.0 min | 16.2-16.5 GB |

- **Memory:** peak RSS is about 17 GB, so the run fits in 32 GB with margin.
- **Wall time:** uncontended, a run takes about 6-7 min. Batch 2 ran 4× slower on some nodes (up to 28.8 min) for the same output. The CPU model is the same everywhere (EPYC 9555). Batch 1 raw01_s1 was slower in every stage (calibration 436 s against 106 s). So this is contention on shared codon nodes, not the code. A clean wall figure needs an exclusive node or the laptop itself.

## Next
- **Own-q looseness on Astral** is now the main problem (R-WP4). It needs a precheck, for example of whether the GBM or the sibling stage moves the q ≤ 0.01 FDP_c, before any change.
- **Seed instability** (raw02_s4, raw01_s4) is worth a look: which stage diverges.
