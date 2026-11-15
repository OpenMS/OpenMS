# PLAN_SINGLE_BINARY steps 2-5: OpenDIAWorkflow

`OpenDIAWorkflow` (OpenMS 794ee1565e, README e612e34d9c) runs the whole analysis in one C++ binary on the `.raw` files: calibration, search, population, one evidence pass, 24 features, the scorer and q-values. A .NET 8 runtime is needed for raw reading.

## Feature parity (raw01, arm-E calibration reused, against the arm-E arrays)

- **Population:** identical, 1,628,266 windows in the same order with equal folds.
- **Search, boundary and b/y (11 features):** equal, with a maximum absolute difference of 8e-6.
- **Interference (4):** 5 windows differ.
- **12-fragment correlation features:** differ, because of the apex-scan bin.
  - The xipm export put the apex scan in the lower central bin in 77 % of windows. The C++ grid always puts it in the upper bin, the rule validated for the interference port.
  - With that bin moved in the reference, 7-16 of 200,000 windows differ.
- **Deficit noise floor:** the median over all windows (1119.7) instead of a 200,000-window sample (1122.0).

## Acceptance (frozen criteria of PLAN_SINGLE_BINARY; arm-E calibration reused, seeds 0-4)

N@FDP_c 1 %:

| Raw | OpenDIAWorkflow | median | Python D11FDC | median | ratio | criterion ≥ 0.98 |
|---|---|---|---|---|---|---|
| raw01 | 1819, 1854, 1817, 1821, 1803 | 1819 | 1837-1857 | 1848 | 0.984 | pass |
| raw02 | 1816, 1902, 1940, 1860, 1870 | 1870 | 1903-1976 | 1925 | 0.971 | **fail** |

- **FDP_c at q ≤ 0.01:** 0.0048-0.0082 on raw01 and 0.0024-0.0037 on raw02, so ≤ 0.015 passes.
  - The tool's own q-values give the same N and FDP_c as the evaluator's TDC on its scores.
- **N at q ≤ 0.01:** 1641-1713 (raw01) and 1610-1676 (raw02).
- **Reading:**
  - The C++ scorer on the arm-E arrays had a raw02 median of 1905 (step 1). The further drop to 1870 comes with the feature differences above. The apex-bin rule is the likely cause, but this is not tested.
  - The raw02 seed spread (1816-1940) is about as wide as the gap.
  - Step 5 therefore fails on raw02 by 0.9 percentage points of the criterion. No retuning was done.

## Self-calibrating two-run analysis (one command, both raws, seed 0)

| Run | rt tolerance | population | N@FDP_c 1 % | q ≤ 0.01 (FDP_c) |
|---|---|---|---|---|
| raw01 | 150.9 s (arm E: 90.8) | 1,858,769 | 1718 | 1723 (0.0104) |
| raw02 | 123.6 s (arm E: 85.6) | 1,710,242 | 1941 | 1777 (0.0045) |

- **Global q:** 1839 target precursors at global q ≤ 0.01, in 3317 run-level rows.
- **Calibration:** the current calibration code needed the windowed second pass on both raws.
  - It gave wider RT tolerances than the arm-E calibration, and so larger populations.
  - raw01 is lower at 1 % (1718, against 1803-1854 with the arm-E calibration), and its q ≤ 0.01 set has FDP_c 0.0104.
  - This is the open "scoring not rerun with the new calibrations" item. It is now measured.

## Cost

- **Single run on 8 threads (raw01, calibration reused):** 7:41 min and 9.6 GB peak.
  - Search 31 s; evidence pass 285 s; features 25 s; scorer 72 s.
  - Acceptance runs took 10-20 min each, because they shared the machine with the two-run job.
- **Two-run self-calibrating job:** 52 min and 22 GB peak, with the concurrent load. Calibration takes 5-6 min per run (2 passes).

## Files

- `EVAL_accept.txt`: raw01 seeds 1-4 and raw02 seeds 0-4. raw01 seed 0 is in `runs_summary.txt`, the w1 run.
- `runs_summary.txt`
- `accept.sh`
- `eval_run.py`
- `eval_workflow.py`
