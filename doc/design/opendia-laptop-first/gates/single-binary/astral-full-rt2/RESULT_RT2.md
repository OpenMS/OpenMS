# Result: Astral full-library runtime round (mem5 + rt2) identity and resources, 25 November 2026

**Verdict: PASS, 15/15.** The frozen rule required, per cell: actual executable exit 0, baseline exit 0, exact Slurm task accounting COMPLETED 0:0, and byte-identical `precursors.tsv`, run precursor TSV and calibration JSON against the same raw and seed baseline. All 15 cells meet every condition. No run failed, was replaced or was excluded. The compact-reference evaluation against `astral-parquet-20261124` also gives 15/15, as it must.

Frozen plan `6533717a83`, preflight approval `c409adb65b`, submission record `bfafd873aa`. Source `bdb0159ff4` (contains mem5 `a11a95ec05`). Container SHA256 `775a7bcd493820b513abee1b4f81b2f7e692e5862649009e583fd0be17092e1f`, OpenDIA SHA256 `1587283d3c6fe1d705d145576b44d60ec2d9803a4a8eca1061419dc7515b6794`. Same 17,161,813-precursor library (the compact experiment's Parquet caches, hash-verified), raws, seeds, options, 16 cores and 256G as the baseline and compact arms.

## Jobs

| Stage | Slurm ID | State | Note |
|---|---|---|---|
| Primary array `0-14%3`, direct Parquet | `42535684` | 15/15 COMPLETED 0:0 | 15/15 PASS |
| Evaluation, `afterany` | `42535685` | COMPLETED 0:0 | `primary=0`, `compact_reference=0` |

## Primary cells (GNU time, matched raw and seed)

| Run | Wall s, rt2 | Wall s, compact | Wall s, baseline | Peak RSS GiB, rt2 | Compact | Baseline | RSS reduction vs compact | rt2 node |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| raw01_s0 | 2284 | 2207 | 5903 | 51.7 | 63.2 | 85.7 | 18.2% | cn-059 |
| raw01_s1 | 2288 | 2160 | 5940 | 52.1 | 63.1 | 86.6 | 17.3% | cn-063 |
| raw01_s2 | 2013 | 2182 | 6528 | 51.2 | 63.2 | 86.3 | 18.9% | cn-065 |
| raw01_s3 | 2053 | 2174 | 5343 | 51.9 | 63.4 | 86.0 | 18.1% | cn-087 |
| raw01_s4 | 2332 | 2186 | 4838 | 52.0 | 63.3 | 86.5 | 18.0% | cn-041 |
| raw02_s0 | 2351 | 2257 | 5676 | 53.0 | 66.0 | 87.3 | 19.6% | cn-052 |
| raw02_s1 | 2088 | 2385 | 4038 | 53.3 | 65.6 | 87.9 | 18.7% | cn-075 |
| raw02_s2 | 2115 | 2555 | 7025 | 52.8 | 65.3 | 88.6 | 19.0% | cn-118 |
| raw02_s3 | 2262 | 4195 | 6567 | 53.0 | 65.5 | 88.3 | 19.1% | cn-052 |
| raw02_s4 | 2044 | 2492 | 5607 | 53.7 | 65.1 | 88.1 | 17.5% | cn-071 |
| raw03_s0 | 2014 | 2268 | 5524 | 52.6 | 64.8 | 88.6 | 18.8% | cn-037 |
| raw03_s1 | 2028 | 2274 | 5978 | 52.8 | 64.8 | 88.6 | 18.6% | cn-039 |
| raw03_s2 | 2497 | 2317 | 4826 | 53.3 | 64.5 | 88.6 | 17.4% | cn-011 |
| raw03_s3 | 2106 | 2568 | 7148 | 53.2 | 64.4 | 88.6 | 17.4% | cn-028 |
| raw03_s4 | 2025 | 2227 | 5133 | 53.3 | 64.4 | 86.4 | 17.3% | cn-039 |

**Peak RSS:** median 52.8 GiB against 64.5 GiB compact (median reduction 18.2%, range 17.3% to 19.6%) and 87.9 GiB baseline (−39.9%). The reduction is consistent in every cell.

**Wall:** median 35.1 min against 37.8 min compact (median per-cell ratio 1.10x, range 0.93x to 1.85x) and 94.6 min baseline (2.60x). Total primary wall was 9.0 h against 10.1 h compact. Four cells were slower than their compact match (0.93x to 0.97x). The 1.85x cell is the compact arm's slow outlier (raw02_s3, 4195 s on cn-007), so the range is not a code effect.

## Per-stage seconds (OpenDIA logs, both passes summed, median and range over 15 cells)

| Stage | rt2 | Compact | Median per-cell difference |
|---|---:|---:|---:|
| Calibration | 124 (49-228) | 103 (47-207) | +19 |
| Search | 192 (183-242) | 191 (189-256) | −8 |
| Extraction | 922 (864-1132) | 977 (892-2071) | −54 |
| Features | 261 (251-269) | 261 (245-285) | 0 |
| Scoring | 385 (344-590) | 541 (491-1151) | −159 |
| Run total (OpenDIA `Combined` seconds) | 2094 (2002-2485) | 2256 (2148-4184) | −202 |

## Interpretation

- **Memory.** In the compact arm the peak came during the first scoring pass (about 62 to 65 GB), above the end-of-Features level (about 54 GB). With extraction state released and features moved into the scorer, scoring now stays at about 28 to 35 GB resident. The run peak is now set at the end of the first pass's Features stage: 52.9 / 54.3 / 53.9 GB on raw01/02/03 seed 0, against 53.9 / 55.0 / 54.5 GB at the same point in the compact arm. Further reductions must come from extraction and Features (the open-grid and held-tensor working set), not from scoring.
- **Speed.** Scoring is the only stage with a clear shift: about 159 s per run, roughly 29% of the compact scoring time. This matches the local profile, where the scorer network's scalar forward pass was the largest hot spot. Calibration is 19 s slower in the median. The local calibration change measured faster, so this is attributed to node and filesystem variability (the calibration range overlaps fully), not claimed as a regression or a gain. Extraction and the total wall differences are within the spread between nodes.
- **Limits.** Nodes were assigned by Slurm and differ between arms and cells, and filesystem caches were not controlled. Wall and stage comparisons are descriptive, not an isolated speedup measurement. Byte identity implies unchanged identifications, q-values and calibration for these 15 cells. No new FDR or quantification claim.

## Evidence

`evidence/` holds the official evaluator output (`results/comparison.json`, `results/comparison_compact.json`, `results/slurm_accounting.tsv`, `results/evaluation_exit_codes.txt`), job IDs, the two-cache `library.sha256`, all Slurm logs, and per-cell exit code, GNU time, node and OpenDIA logs. `evidence/compact_reference/` holds the compact arm's OpenDIA logs used for the stage table. Output TSVs are not copied. They remain in `/nfs/research/juan/singj/opendia/astral-rt2-20261125/runs`. `evidence/EVIDENCE.sha256` lists every copied file. No remote file was modified while collecting this record.
