# DIA-NN 2.2 on the r ≈ 7 entrapment library (19 Nov 2026)

This is a descriptive reference, not a frozen plan. The user asked how OpenDIA compares with DIA-NN on the same entrapment
library.

## Setup

**Library.** This is the targets-only library that every PLAN_PDEC and PLAN_ODP run used.
- File: `entlib-r5-20261118/targets.pqp`, sha256 b29d8b43.
- 1,390,163 target precursors, r = 7.05.
- It was converted to TSV with TargetedFileConverter (pd3 container), giving `targets_r5.tsv` (sha256 697f3f6b).
- Script: `diann_r5_conv.sbatch`.

**DIA-NN.**
- DIA-NN 2.2.0 singularity image, 16 threads, codon.
- Flags are the same as `bench/bench_hela.sh`: `--qvalue 1 --unimod4 --min-pr-charge 1 --max-pr-charge 5`, mass
  accuracy automatic. Script: `diann_r5.sbatch`.
- DIA-NN 2.2 does not recognise `--gen-decoys` (its log says so). It generates its own decoys either way, so this is the
  same as in the earlier HeLa benchmark.
- It reads the library's predicted RT, IM and intensities.

**Evaluation.**
- `eval_diann_r5.py` uses the same entrapment labels, r and bounds as `eval_pdec.py`.
- All 58,769 (D) and 84,276 (H) report precursors matched the library.
- DIA-NN is deterministic here, so there is one run per file. OpenDIA numbers are pooled or median over 5 seeds
  (EVAL_ODP_D/H).

## Result (`EVAL_DIANN_R5.txt`)

| | DIA-NN D | OpenDIA ODP D, CUR / NEW | DIA-NN H | OpenDIA ODP H, CUR / NEW |
|---|---|---|---|---|
| IDs at own q ≤ 0.01 | 26,118 | ~8,650 / ~10,420 | 41,649 | ~21,340 / ~25,270 |
| True FDP at own q ≤ 0.01, lower bound | **0.91 %** (0.81-1.01) | 1.06 % / 1.06 % | **1.00 %** (0.92-1.09) | 1.06 % / 1.10 % |
| True FDP at own q ≤ 0.01, combined bound | **1.04 %** (0.93-1.15) | 1.21 % / 1.21 % | **1.15 %** (1.06-1.24) | 1.21 % / 1.26 % |
| N at true FDP 1 % (combined) | **26,006** | 8,233 / 10,100 | **41,272** | 20,745 / 24,549 |
| OpenDIA / DIA-NN at true FDP 1 % | | 32 % / 39 % | | 50 % / 59 % |
| Wall, peak memory (seed 0) | 2:38, 4.6 GB (16 threads) | 2:28 h / 1:56 h, 89 / 92 GB (8 threads) | 4:08, 6.8 GB | 1:01 h / 1:07 h, 110 / 115 GB |

The intervals are 90 % Poisson intervals on the entrapment count. OpenDIA IDs are pooled N divided by 5.

**Overlap (seed 0, non-entrapment IDs at each tool's q ≤ 0.01).**
- 95-97 % of OpenDIA's IDs are also in DIA-NN's list.
- DIA-NN's IDs found by OpenDIA:
  - D: 32 % CUR, 39 % NEW.
  - H: 49 % CUR, 58 % NEW.

## Reading

- **DIA-NN is calibrated on this library and OpenDIA is slightly not.**
  - DIA-NN's lower bound is at or below 1 % on both runs, and its combined bound is 1.04 % (D) and 1.15 % (H). Both
    would pass the PLAN_ODP limits (lower ≤ 1.0 %, combined ≤ 1.2 %), H with no margin.
  - OpenDIA with `-decoys predict` sits about 0.05-0.15 pp higher.
  - On Astral the order was the reverse: DIA-NN's combined FDP there was 1.1-1.4 %.
- **The ID gap is much larger than calibration.**
  - At 1 % true FDP OpenDIA reaches 32-39 % of DIA-NN on D and 50-59 % on H.
  - The plain HeLa library gave 56 % / 63 % (default) and 65 % / 74 % (`-pg_scores`) at each tool's own q. So adding
    seven foreign peptides per real one hurts OpenDIA much more than DIA-NN.
  - DIA-NN loses about 10 % of its plain-library IDs on D (28,950 to 26,006 at true 1 %). OpenDIA loses about 49 %
    (16,116 at own q to 8,233).
- **What OpenDIA reports is almost all in DIA-NN's list** (95-97 %). The deficit is missed IDs, not different IDs.
- **DIA-NN's q is accurate in its tail as well.**
  - At q ≤ 0.005 its combined bound is 0.54 % / 0.67 %. At q ≤ 0.02 it is 1.82 % / 2.10 %.
  - OpenDIA's lower bound is about 1.7-1.9× its own decoy estimate at q ≤ 0.005 (PRECHECK_TDC).
- **Speed and memory.**
  - DIA-NN takes 2.6-4.1 min and 4.6-6.8 GB.
  - OpenDIA on this 1.39 M-precursor library takes 1.0-2.5 h and 89-115 GB, which is 15-55× slower and 15-25× more
    memory. Most of that cost is the tensor extraction per window (13.6 M windows on H).
- **Caveat.** DIA-NN used 16 threads and OpenDIA 8. Wall time is therefore not like-for-like, but a gap of 15× or more
  is not a thread effect.
