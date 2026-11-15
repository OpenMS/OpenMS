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

## Diagnosis after the verdict (descriptive; the frozen verdict above is unchanged)

### raw02 shortfall: block swap (`swap.sh`, `SWAP_EVAL.txt`)
The C++ scorer was run offline on the arm-E arrays, with one block replaced by the OpenDIAWorkflow block (seeds 0-4).

| Raw | arm-E arrays (step 1) | all workflow blocks | workflow f12 only | workflow deficit only | in-binary (C++ seed rows) | Python D11FDC |
|---|---|---|---|---|---|---|
| raw01 median N@1 % | 1822 | 1787 | 1808 | 1831 | 1819 | 1848 |
| raw02 median N@1 % | 1905 | 1881 | 1889 | 1927 | 1870 | 1925 |

- **Size of the differences:** all medians lie within 25-55 IDs of each other.
  - The per-seed range of every arm is 50-200 IDs (raw02 WF 1849-1916, step 1 1826-2021).
- **No single cause:** no block explains the raw02 shortfall on its own.
  - The f12 apex-bin change costs about 15 IDs at the median, and the deficit change none.
  - Building the training rows in C++ costs about 10 more.
  - These are smaller than seed noise, so the suspected cause ("apex-bin rule") is not confirmed.
- **What remains:** a scorer-level gap to the Python trainer of about 20-55 IDs at the median, comparable to its seed spread.

### Self-calibration: inflated RT spread (`spread_fix*.py`, `SPREAD_FIX*.txt`)
- **Map accuracy:** the workflow's arm-E raw01 calibration is as accurate as before against DIA-NN 1 % IDs.
  - |err| p50 is 18.3 s (16.9 s before).
  - p95 is 95.2 s (80.8 s before).
- **Inflated spread:** the selection-free spread fit gave 37.7 s on raw01, against a robust spread of 27.2 s for DIA-NN around the same map.
  - On raw02 it gave 30.9 s against 25.4 s.
  - Hence the tolerances of 151 s and 124 s, the larger populations and the raw01 FDP_c of 0.0104.
- **Cause:** the fit takes target residuals relative to the final map, but decoy (background) residuals relative to the stage-3 window centre.
  - On arm E raw01 the two differ by a median 11.4 s, and the half window is 110.7 s.
  - The background edges then do not line up with the targets', and the Gaussian absorbs the mismatch.
  - On the earlier calibration runs the shift is 0.5-4 s and the spread matched DIA-NN (23-28 s against 24-27 s).
- **Two offline fixes tried, neither adopted** (coverage of DIA-NN IDs is the criterion):
  - Decoys relative to the final map: the spread falls to 17.7-23.3 s and coverage to 0.93-0.98, below the current fit on every run.
  - Window frame with a per-target shift: unstable, collapsing to 0-5 s on two runs.
- **Status:** the calibration code is unchanged and needs a new design (for example a background that is shifted per member), with its own plan.

## Further validation (15 November)

- **Inventory checked against the chain scripts.** The validated chain (arm E + PLAN_GAP) runs 4 C++ tools and reads the raw file 6 times:
  - OpenDIA three times;
  - OpenDIACNN once;
  - OpenDIAInterference once;
  - OpenDIACalibrate once, on the PAIRLIB library, whose calibration arm E reused.

  It also runs 7 Python data steps, the Python interference prototype and 3 Python scoring steps.
  OpenDIAWorkflow reads the raw file 3-4 times: calibration takes 1-2 passes, then one for the search and one for the evidence.
- **Determinism:**
  - Four self-calibrated runs give run tables byte-identical to runs that reuse the same calibration file. The calibration itself is byte-identical to a standalone OpenDIACalibrate run.
  - A repeated raw02 run is byte-identical, and so are 1 and 8 threads (raw02, arm-E calibration).
- **Report:**
  - The run table has one row per precursor. RT, area and support equal the search feature table exactly (candidate_row join).
  - The run q-values equal an independent pooled paired TDC recompute, with 0 differences on 2 runs.
  - The global q recomputed from the two run tables gives 1839 targets, equal to the tool's. The 3317 combined rows are exactly the expected set, with no decoys.
  - Calibrated scores share one decoy-tail scale across runs: the decoy median is 0.693 (ln 2) and p99 4.602 (−ln 0.01) in both runs. So pooling best scores over runs is sound.
- **Fail-closed CLI ctests** (OpenMS 9d6bf16f10):
  - a library too small to calibrate;
  - duplicate run file names;
  - `-calibration` without `-anchors`.
