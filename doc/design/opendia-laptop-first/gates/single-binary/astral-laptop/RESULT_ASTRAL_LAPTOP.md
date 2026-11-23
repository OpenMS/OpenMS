# Astral on the current default, laptop profile (23 Nov 2026), interim

Diagnostic, no pass rule. Raw 03 is in the set, as in the 16 Nov snapshot (the README holds raw 03 back for design; nothing here was tuned on it).

## Setup
- Code: OpenMS `46243afae4` (tool renamed `OpenDIA`, `-scoring_lib` restored), container `opendia-pd10-rename.sif` (sha256 `4614e25e…`).
- Inputs as 16 Nov: arm-E search library (`gap-20261114/home/lib/library.pqp`, md5 `e0461beb…`), scoring library `library12_full.pqp`, self-calibrating, raws 01-03, seeds 0-4.
- Laptop profile: `-threads 16`, slurm `-c 16 --mem=32G`, codon AMD EPYC 9555 nodes. Script `astral_laptop.sbatch`, cost from `/usr/bin/time -v`.
- Evaluation: frozen `gates/single-binary/workflow/eval_run.py` (`EVAL_batch1.txt`, `COST_batch1.txt`).

## Robustness: 5 of 15 runs crashed (FAIL)
raw01_s2, raw02_s0, raw02_s1, raw03_s0 and raw03_s1 stopped with exit 8: `pooledPairedTDC: a fold has zero decoy spread`.

**Cause.** The GBM rescore runs 3 rounds. When a round leaves a fold whose model keeps no tree (logged `trees 0/…`, or a degenerate fit), all rows of that fold get the same score, so the fold's decoy MAD is 0. Rounds 2 and 3 already catch this at the start of the next round and keep the previous scores. After round 3 there is no next round, so the throw happens in the final q-value step instead. All 5 crashes are round 3 cases. The 10 completed runs either had a non-degenerate round 3 or stopped earlier through the existing guard.

**Fix (OpenMS `88e5da3023`, local).** After the last round, apply the same check and fall back to the round 2 scores. This only changes runs that used to crash. A completed run passes the check unchanged, so its output is the same by construction. Byte-check on HeLa D156, D7V and H7 submitted (`gbmfix-20261123`), and all 15 Astral runs re-submitted with the fix (`astral-gbmfix-20261123`, container `opendia-pd10-gbmfix.sif`, sha256 `2270438a…`). **Results pending:** codon became unreachable (DNS) before they finished.

## Identifications, the 10 completed runs

| Raw | Seeds done | N @ FDP_c 0.5 / 1 / 2 % per seed | median N @ 1 % | 16 Nov median | DIA-NN (.raw) @ 1 % | share |
|---|---|---|---:|---:|---:|---:|
| 01 | s0 s1 s3 s4 | 1946/2124/2313, 1810/2275/2474, 1770/1950/2220, 497/1958/2331 | 2041 | 1804 | 2204 | 93 % |
| 02 | s2 s3 s4 | 1987/2144/2431, 1939/2128/2445, 653/1301/2441 | 2128 | 1909 | 2356 | 90 % |
| 03 | s2 s3 s4 | 1869/2160/2503, 2003/2252/2451, 2057/2343/2521 | 2252 | 1907 | 2259 | 100 % |

- **Own q ≤ 0.01 is now anti-conservative on Astral:** FDP_c 0.77-1.77 %, and 6 of 10 runs are above 1 %. Those runs are raw01_s0 1.77 %, raw01_s3 1.25 %, raw01_s4 1.10 %, raw02_s4 1.41 %, raw03_s2 1.32 % and raw03_s3 1.06 %. On 16 Nov all runs were at or below 1 %. This matches the HeLa H loosening seen with the GBM (RESULT_GBT_H), but it is larger here.
- **Seed instability:** raw02_s4 reaches only 1301 at FDP_c 1 %. At 0.5 %, raw01_s4 has 497 and raw02_s4 has 653, against about 1,800-2,000 for the other seeds.
- **Caveat:** these are survivors. The crashed seeds are missing, so the medians may be biased. Do not quote them as the result until the 15-run rerun is in.

## Cost at 16 threads / 32 GB

| Raw | wall (completed runs), min | median wall | peak RSS |
|---|---|---:|---:|
| 01 | 9.4, 18.3, 6.4, 5.8 | 7.9 | 16.9-17.2 GB |
| 02 | 7.8, 7.1, 6.7 | 7.1 | 16.2-16.3 GB |
| 03 | 6.2, 6.0, 9.5 | 6.2 | 16.2 GB |

- Peak RSS is about 17 GB, so the run fits in 32 GB with margin.
- raw01_s1 (18.3 min) ran on a different node (cn-061) with every stage about 1.5-4× slower (calibration 436 s vs 106 s). This looks like node contention, not the code.
