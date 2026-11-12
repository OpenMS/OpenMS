# DIA-NN 2.2.0 on our exact library, raw01: IDs, overlap, runtime and memory

12 November 2026. Diagnostic only, no OpenDIA change. Scripts and logs in `diann-ourlib-20261112/` (benchmark tree).

## Setup
- DIA-NN 2.2.0 (docker `diann:2.2.0`; ProteoBench used 2.3.1 with an 8.65M predicted library, 3 runs, MBR).
- **Same library as OpenDIA**: the 300,223 target precursors of our exact TSV (6 fragments, our library RT and intensities). DIA-NN makes its own decoys (`--gen-decoys`).
- Input `raw01.mzML` (FileConverter of the same .raw). Single run, no MBR, no `--reanalyse`, auto mass accuracy, `--unimod4`, `--qvalue 1`, 14 threads.
- Evaluation with our entrapment labels and the same r = 148983/151239.

## IDs

| | DIA-NN, our library | D11FWT |
|---|---:|---:|
| q ≤ 0.01 precursors | **2357** | 1614 |
| of which entrapment | 22 | 11 |
| FDP_c at q ≤ 0.01 | 0.0188 | 0.014 |
| N at FDP_c 1 % | **2161** | 1496 |
| N at FDP_c 2 % | **2425** | 1716 |

- At q ≤ 0.01: 1487 human IDs in both, 116 only ours, **848 only DIA-NN**. 744 of those 848 are in our scoring population.
- ProteoBench reference (2743 human): DIA-NN on our library finds 2121. By our V16 stage: E 1347 / 1413, **D 548 / 866**, C 64 / 113, B 87 / 166, A 75 / 185.
- So most of the gap to ProteoBench is our pipeline, not the library. With our exact library and one run, DIA-NN gets 2161 at FDP_c 1 % against our 1496.

## Where we lose DIA-NN-only IDs (in population, `DIFF_raw01.txt`)
- 91 % of DIA-NN q ≤ 0.01 targets in our population have our winner within 10 s of DIA-NN's RT. 95 % have a candidate there.
- Of the 551 DIA-NN-accepted precursors we rank below decoy p99.9, 359 have our winner at DIA-NN's RT. Median own percentile 0.53, but linear part at 0.95.
- DIA-NN marks them as weaker too: Evidence 2.59 vs 3.75, MS1 profile corr 0.57 vs 0.80, quantity 0.35M vs 2.1M.

## Which model input does the penalising (`attr/attr.py`, exact V16 round-2 checkpoints)
For the 350 DIA-NN-only windows at the DIA-NN RT, the nonlinear part is negative for all of them (median −29 decoy SD, against −20 for decoys). Replacing one input at a time with the decoy median (effects are not additive, so this ranks inputs, it does not split the total):
- **`log_peak_width` has the largest mean effect** (−8.4 SD for DN, about 0 for decoys and entrapment, −3.6 for confident targets). It comes from a subset: median effect −0.3, but a quarter of DN windows lose more than 8 SD and a tenth more than 37 SD. The effect tracks narrow peaks (correlation 0.72 with z). DN widths are not unusual for real peptides (median z −0.02 against −0.24 for confident targets), so the network appears to penalise a narrow peak when it comes with only moderate support.
- `transition_support` (−2.1) and `log_matched_scans` (−2.0) are next. Both are measures the 32-bin grid makes coarse for a 3 s peak.
- This is a diagnostic of the trained model, not a tested fix.

## Runtime and memory, raw01 (same machine: i7-1370P, 20 logical CPUs)

| Tool | Wall | CPU time | Peak memory | Threads | Notes |
|---|---:|---:|---:|---:|---|
| DIA-NN 2.2.0, full search + scoring | **10:13** | 7580 s (12.6 cores avg) | **4.34 GB** (cgroup memory.peak) | 14 (`--cpus 14`) | of which 3:29 reading mzML, 5:42 calibration and mass-accuracy optimisation, 0:48 search + NN |
| OpenDIACalibrate (window300k) | 3:37 | 219 s | 16.5 GB (max RSS) | 1 effective (100 % CPU) | reads `.raw` directly |
| OpenDIA search, 6 fragments | 2:47 | 165 s | 4.33 GB | 1 effective | `astral-design-raw01/pipeline` |
| OpenDIA search + `.xipm` export, 6 fragments | 17:28 | 1032 s | 7.77 GB | 1 effective | voxel population export used by D11FWT |
| OpenDIA search + export, 12 fragments | 18:42 | 1100 s | 8.53 GB | 1 effective | lib12 export |
| D11F training (Python, offline) | 1:15–2:38 | | 2.4 GB | 10 | per arm, plus W and T |

- DIA-NN: identical IDs on a repeated timed run (2357, same set). Memory is the container cgroup peak (includes page cache charged to the container), measured every second. Wall from container start/finish times.
- OpenDIA numbers are `/usr/bin/time -v` from earlier runs, single-threaded in practice. The scoring that gives D11FWT is offline Python and not yet in C++, so there is no end-to-end OpenDIA number comparable to DIA-NN's yet.
- DIA-NN spends 3:29 reading mzML (it would read `.raw` directly on Windows). Excluding that, search plus calibration is about 6:40.
