# Result: Astral full-library direct-Parquet identity and resources, 24 November 2026

**Verdict: PASS, 15/15 primary runs plus a separate bootstrap PASS.** The frozen rule required, per cell, actual executable exit 0, baseline exit 0, exact Slurm task accounting COMPLETED 0:0, and byte-identical `precursors.tsv`, run precursor TSV and calibration JSON against the same raw and seed baseline. All 15 cells meet every condition. No run failed, was replaced or was excluded.

Frozen plan `045aae4160`, preflight approval `97a803dc71`, submission record `cc5652e498`. Source `8c0a89f053`, container SHA256 `035cd7004896d636dc83800230876b4a504fbc247433ea180c02e7ddaffb9caf`, OpenDIA SHA256 `83c408e5ad42e4865283a0f8586f18b223139f0f9a1b7d5c9f3ea491c9235eb7`. Baseline is the frozen PLAN_ASTRAL_FULL attempt 3 (`astral-full-20261123`), the same 17,161,813-precursor predicted library, the same raws, seeds, options, 16 cores and 256G.

## Jobs

| Stage | Slurm ID | State | Note |
|---|---|---|---|
| Bootstrap, raw01 seed0, uncached PQP input | `42109936` | COMPLETED 0:0 | 3/3 files identical; wrote both Parquet caches |
| Primary array `0-14%3`, direct Parquet | `42109937` | 15/15 COMPLETED 0:0 | 15/15 PASS |
| Evaluation, `afterany` | `42109938` | COMPLETED 0:0 | `bootstrap=0`, `primary=0` |

## Primary cells (GNU time, matched raw and seed)

| Run | Wall s, compact | Wall s, baseline | Wall ratio | Peak RSS GiB, compact | Peak RSS GiB, baseline | RSS reduction | Compact node |
|---|---:|---:|---:|---:|---:|---:|---|
| raw01_s0 | 2207 | 5903 | 2.67x | 63.2 | 85.7 | 26.3% | hl-codon-cn-043 |
| raw01_s1 | 2160 | 5940 | 2.75x | 63.1 | 86.6 | 27.2% | hl-codon-cn-045 |
| raw01_s2 | 2182 | 6528 | 2.99x | 63.2 | 86.3 | 26.8% | hl-codon-cn-120 |
| raw01_s3 | 2174 | 5343 | 2.46x | 63.4 | 86.0 | 26.2% | hl-codon-cn-045 |
| raw01_s4 | 2186 | 4838 | 2.21x | 63.3 | 86.5 | 26.8% | hl-codon-cn-051 |
| raw02_s0 | 2257 | 5676 | 2.51x | 66.0 | 87.3 | 24.4% | hl-codon-cn-043 |
| raw02_s1 | 2385 | 4038 | 1.69x | 65.6 | 87.9 | 25.4% | hl-codon-cn-064 |
| raw02_s2 | 2555 | 7025 | 2.75x | 65.3 | 88.6 | 26.3% | hl-codon-cn-065 |
| raw02_s3 | 4195 | 6567 | 1.57x | 65.5 | 88.3 | 25.9% | hl-codon-cn-007 |
| raw02_s4 | 2492 | 5607 | 2.25x | 65.1 | 88.1 | 26.1% | hl-codon-cn-064 |
| raw03_s0 | 2268 | 5524 | 2.44x | 64.8 | 88.6 | 26.9% | hl-codon-cn-065 |
| raw03_s1 | 2274 | 5978 | 2.63x | 64.8 | 88.6 | 26.8% | hl-codon-cn-009 |
| raw03_s2 | 2317 | 4826 | 2.08x | 64.5 | 88.6 | 27.2% | hl-codon-cn-065 |
| raw03_s3 | 2568 | 7148 | 2.78x | 64.4 | 88.6 | 27.4% | hl-codon-cn-064 |
| raw03_s4 | 2227 | 5133 | 2.30x | 64.4 | 86.4 | 25.4% | hl-codon-cn-096 |

Median wall 37.8 min against 94.6 min (median per-cell ratio 2.46x, range 1.57x to 2.99x). Median peak RSS 64.5 GiB against 87.9 GiB (median reduction 26.3%, range 24.4% to 27.4%). Total primary wall 10.1 h against 23.9 h. The RSS reduction is consistent in every cell. raw02_s3 is the slow outlier (4195 s on `hl-codon-cn-007`, other compact cells 2160-2568 s) with identical output, so it is attributed to the node, not the code.

## Bootstrap (reported separately, not run 16)

Uncached PQP input, first automatic cache creation: 5874 s wall against the matched baseline's 5903 s, peak RSS 62.8 GiB against 85.7 GiB, 3/3 identical. Conversion wrote `search6.opendia.parquet` (2100 MiB, 94 s) and `scoring12.opendia.parquet` (2971 MiB, 128 s). Cache SHA256 in `evidence/bootstrap/library.sha256`; every primary task verified them before running.

## Interpretation

- **Where the time went (raw01 seed0 stage logs).** Library loading fell from about 3,220 s (baseline PQP read, total 5,885 s minus 2,665 s run time after load) to about 183 s (direct Parquet, 2,196 s minus 2,013 s). That is most of the saving. The analysis after load was also faster, 2,013 s against 2,665 s: calibration 103 against 234 s, extraction 493 and 438 against 678 and 649 s, scoring 266 and 247 against 302 and 274 s. Plausible causes are the compact indexed library (better cache locality) and no calibrated library copies, but the two arms ran on different nodes, so that part is not attributed. The bootstrap still pays the PQP read once and finished in about the baseline's time.
- **Where the memory went.** Resident memory after library load fell from 47,306 MB to 12,391 MB. Peak RSS fell less, from 87,750 MB to 64,684 MB, because it is now set during the first scoring pass (50,315 MB peak at the end of extraction, 64,684 MB after scoring). In the baseline the peak came during extraction, on top of the resident library copies. Further peak reductions need work in scoring and extraction, not in the library.
- **Limits.** Nodes were assigned by Slurm and differ between arms and between cells, and filesystem caches were not controlled. Wall comparisons are descriptive, not a causal or isolated speedup measurement. Byte identity implies unchanged identifications, q-values and calibration for these 15 cells. No new FDR or quantification claim is made.

## Evidence

`evidence/` holds the official evaluator output (`results/comparison.json`, `results/bootstrap.json`, `results/slurm_accounting.tsv`, `results/evaluation_exit_codes.txt`), bootstrap accounting, job IDs, all Slurm logs, and per-cell exit code, GNU time, node and OpenDIA logs. Output TSVs are not copied; they remain in `/nfs/research/juan/singj/opendia/astral-parquet-20261124/runs`. `EVIDENCE.sha256` lists every copied file. No remote file was modified while collecting this record.
