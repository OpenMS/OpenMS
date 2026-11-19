# Why OpenDIA is slower, heavier and less sensitive than DIA-NN on a large library (19 Nov 2026)

This is a diagnosis from existing logs: the PLAN_ODP / BOUNDS_R5 seed-0 runs and RESULT_DIANN_R5. No new runs were made,
and `funnel_r5.py` is the only new analysis. DIA-NN labels are used for evaluation only.

## Stage cost, OpenDIA (8 threads, codon)

| | plain HeLa lib D | r ≈ 7 lib D | r ≈ 7 lib H |
|---|---|---|---|
| library precursors (targets + decoys) | 0.73 M | 2.78 M | 2.78 M |
| search candidates / scored windows | 2.6 M / 2.49 M | 10.4 M / 9.50 M | 13.4 M / 12.4 M |
| calibration used anchors (panel 80,000) | 967 | 140-218 | 441 |
| extraction s | 223 | 738 | 990 |
| **features s** | **170** | **7,491** | 1,926 |
| scorer s | 21 | 39 | 42 |
| wall | 7:54 | 2:28 h | 1:01 h |
| peak RSS | 24 GB | 87 GB | 107 GB |
| 6-fragment trace tensor (6 x 32 x 8 float per window) | 15 GB | 58 GB | 76 GB |

DIA-NN 2.2 on the same r ≈ 7 library takes 2:38 and 4.6 GB on D, and 4:08 and 6.8 GB on H.

## Findings

**1. Memory.** The design keeps every scored window in RAM until scoring.
- Each window holds a dense 6 KB trace tensor plus grids and evidence.
- The tensor alone is 60-70 % of peak RSS, so memory is O(windows).
- There are about 3.8 windows per precursor, and all of them are kept.
- DIA-NN streams per precursor and keeps only the best peak's scores.

**2. Runtime.** The features stage is 84 % of the D wall time, and it grows faster than linearly with crowding.
- Windows grow 3.8x (2.49 M to 9.50 M), but features time grows 44x (170 s to 7,491 s).
- H has 30 % more windows than D, but its gradient is about 4x longer (1,319 against 346 cycles), so windows are about
  3x less crowded per second. Its features stage is 4x faster than D's.
- The per-window feature loop is linear. Inside this stage, only the pairwise interference resolver
  (`InterferenceResolver::compute`) scales with co-eluting density: every overlapping window in the same isolation
  window is tested.
- This is the likely cause but it has not been profiled yet. It needs a timer around the resolver.

**3. Calibration and seeding collapse with library size.**
- The calibration panel is a fixed 80,000 targets. With seven foreign peptides per real one, about 7/8 of the panel cannot
  be present.
- Used anchors therefore drop from 967 to 140-218 on D, and scorer positive seeds from 943 to 138-213.
- DIA-NN trains on its whole population (56,555 targets and 33,881 decoys on D).

**4. Sensitivity: the loss is in scoring, not in finding the peak** (`funnel_r5.py`, seed 0, CUR).

| | D | H |
|---|---|---|
| DIA-NN non-entrapment q ≤ 0.01 | 25,881 | 41,231 |
| with no OpenDIA window at all | 0.0 % | 0.0 % |
| accepted by OpenDIA (q ≤ 0.01) | 31.6 % | 49.3 % |
| scored, not accepted | 68.4 % | 50.7 % |
| ... of which the best window is within 6 s of DIA-NN's RT | 60 % | 51 % |
| ... of which OpenDIA q ≤ 0.05 | 14 % | 13 % |
| accepted by DIA-NN quantity quartile (low to high) | 7 / 16 / 32 / 71 % | 11 / 33 / 63 / 90 % |

- Every DIA-NN ID reached OpenDIA's scorer, so candidate generation is not the limit.
- About 60 % of the missed IDs had the right peak but scored far too low. Most of them are not even near the threshold.
- The other 40 % picked a different peak.
- The loss is concentrated at low abundance.
- This gap also exists on the plain library (56-63 % of DIA-NN), and library size makes it worse.

**5. Why it looked fine before.** The laptop-first design was validated on Astral (302 k targets, 8 threads).
- There, speed matched DIA-NN (7:18 against 7:23) at 14.4 GB.
- A full predicted library (1.4 M targets here, 8.65 M in the ProteoBench DIA-NN reference) is outside what that design
  was tested on. The O(windows) memory and the density-dependent interference cost then dominate.

## What would fix it (each needs its own precheck and frozen plan)

1. **Do not keep the tensor.**
   - Compute the per-window features while extracting, isolation window by isolation window, then release the tensor.
     Keep only the 24-34 floats per window (about 1 GB at 9.5 M windows).
   - Expected peak: about 30 GB, which still includes grids and evidence. With pruning (item 2) and streaming the grid, it
     should come down to a few GB.
   - Outputs should be identical if the features are unchanged. That is a byte-identity check, not a plan.
2. **Prune windows per precursor before the expensive features.** Keep the top-k by search support, for example k = 2.
   - This cuts windows, memory and interference pairs about 2x.
   - It changes IDs, so it needs a plan.
3. **Bound the interference resolver.** Index fragments by m/z so that only windows sharing a fragment are paired,
   instead of all overlapping windows. That is a profile first, then a byte-identity check.
4. **Size the calibration from the file, not a fixed panel.** Grow the panel until a target number of anchors is found,
   for example about 1,000. This is label-free, and it restores seeds on large libraries.
5. **Sensitivity at low abundance.** This is the main gap on every library. It is a scorer and feature question: the
   right peak is found but scored low. It is the largest and least certain item.

## Result of the output-identical fixes (19 Nov, OpenMS 8edbda8426, container pd4 sha256 7e0a574f)

These are the same runs as PLAN_ODP CUR seed 0 on codon (8 threads). All outputs, `precursors.tsv` (both) and
`calibration.json`, are **byte-identical** to the PLAN_ODP runs.

| | before (PLAN_ODP) | after | DIA-NN 2.2 |
|---|---|---|---|
| D wall | 2:28 h | **29:04** | 2:38 |
| D peak RSS | 87 GB | **41.9 GB** | 4.6 GB |
| D features stage | 7,491 s | **268 s** (interference 267 s) | |
| H wall | 1:01 h | **33:13** | 4:08 |
| H peak RSS | 107 GB | **36.7 GB** | 6.8 GB |

- **Tensors are held only while a window's RT range is being read:** at most 1.79 M of 9.50 M on D and 0.66 M of 12.4 M
  on H.
- **D peak RSS comes from extraction (41 GB),** now set by the 1.79 M held tensors (11 GB) plus grids, evidence and
  the library.
- **The interference resolver is still 267 s on D.** The fragment index removed the all-pairs scan, but crowded
  isolation windows still give many shared-fragment partners.
- **Remaining wall time is extraction (807 s / 1,074 s).** That is the single 6-fragment plus 12-fragment pass over the
  raw file.
- **Plain HeLa D (laptop, 8 threads):** byte-identical, 5:25 to 2:38. Peak memory there (4.5 GB) is not tensor-bound.

## OpenSWATH removal check (19 Nov, OpenMS 1f29ce2da8)

**Library reader equivalence (`libreader_equiv.cpp`, validation only).** For each library, `OpenDIA::readLibrary` was
compared with `TransitionPQPFile`, in order. The fields compared were id, RT, 1/K0, charge, sequence, proteins,
transition name, precursor ref, intensity, m/z, fragment charge, decoy flag, annotation, and the TraML id map.

| Library | Precursors | Transitions | Mismatches |
|---|---|---|---|
| plain HeLa | 146,788 | 880,728 | 0 |
| test fixture with compounds | 8 | 44 | 0 |
| r ≈ 1.56 entrapment | 732,674 | 4,310,998 | 0 |
| r ≈ 7 predicted-decoy library | 2,778,189 | 15,049,076 | 0 |

**No OpenSWATH in the core.**
- `ninja -t deps` lists 0 transitive OpenSWATH or OPENSWATHALGO headers in every `OpenDIA_impl` object.
- `libOpenDIA_impl.a` has 0 undefined OpenSwath, TransitionPQPFile, DIAHelpers or PeakMapExtractor symbols.
- The remaining symbols sit in tool objects only.
  - `OpenDIAWorkflow.cpp.o` has 7 (`-decoys predict`).
  - `OpenDIA.cpp.o` has the XIPM writer (`-out_xipm`).
  - `OpenDIACalibrate` and `OpenDIAInterference` have 0.
- **The guard ctest catches a reintroduction.** Adding the DIAHelper include to `OpenDIAStage2.cpp` makes
  `OpenDIA_no_openswath` fail, and it passes again once the include is removed.

## Read vs process split, and parallel per-spectrum work (19 Nov, OpenMS fc24ce9efc / 00cd8bb987)

`fc24ce9efc` logs per raw pass the time spent reading (decoding spectra) and the time spent in OpenDIA's consumer.
Before `00cd8bb987` the consumers ran on one thread. All values below are in seconds, 8 threads on codon, r ≈ 7
library, PLAN_ODP CUR seed 0, with outputs byte-identical to PLAN_ODP throughout.

| Run D | read | process before | process after |
|---|---|---|---|
| calibration pass | 4.7 | 17.9 | **3.7** |
| search pass | 4.5 | 22.0 | 25.0 (already parallel per batch) |
| extraction pass | 7.1 | 715.4 | **205.2** |
| interference resolver | | 215.1 | 266.2 (unchanged code, node noise) |
| run total | | 1,195 (pd4) | **622** |
| wall incl. `-decoys predict` | | 26:12 | **18:34** |
| peak RSS | | 37.7 GB | 38.3 GB |

| Run H | read | process after |
|---|---|---|
| calibration pass | 13.5 | 11.3 |
| search pass | 12.2 | 61.6 |
| extraction pass | 17.3 | 294.1 |
| interference resolver | | 151.8 |
| run total / wall | | 677 s / 19:40 (from 33:13) |
| peak RSS | | 34.0 GB |

**Reading.** Reading the timsTOF file is 4.5-17 s per pass, under 3 % of the time before. The cost was OpenDIA's
serial per-spectrum work.

**Parallelisation (`00cd8bb987`).**
- Each spectrum's independent per-member or per-window work now runs on the workers.
- Order-dependent steps stay serial: tensor slot assignment, grid additions, and isolation windows in order. The
  results therefore do not depend on the thread count.
- Plain HeLa D, default and `-pg_scores`: byte-identical.

**What remains on D (18.5 min wall).**
- About 8 min is `-decoys predict`, which predicts 1.39 M decoys and writes a library before the run. DIA-NN's
  2.6 min does not include an equivalent step.
- The interference resolver takes 266 s.
- The extraction consumer takes 205 s, about 3.5x from 8 threads.
- Search takes 25 s.

## Interference and m/z lookup (19 Nov, OpenMS 571940fcb5 / 36011de90a)

A profile on the r ≈ 1.56 library (local, plain HeLa D, 8 threads) put the interference resolver at 16 % and
`MSSpectrum::MZBegin` at 16 % of samples. Changes:

- **Resolver prefilter.** Index entries carry precursor, RT end and 1/K0, so the cheap rejects of the pair test run
  on the index entry without loading the window or the precursor. `consider()` still repeats every test.
- **Bucketed m/z lower bound.** `SpectrumMzIndex` (0.25 Th buckets) replaces `MZBegin` in traces, the fragment grid
  and calibration. `OpenDIAMzIndex_test` checks that it equals `MZBegin`, including at bucket edges.
- **Bin table (36011de90a).** The resolver index keeps a table of the distinct bins, so each lookup binary-searches
  only inside one bin.
- **Tried and reverted.** Testing the shared fragment during the index scan was output-identical but no faster
  (38.9 s against 36.0 s).

Local r ≈ 1.56 D, interference seconds: 53.0 → 42.4 (prefilter) → 36.0 (bin table). Outputs are byte-identical
after each step, on plain HeLa D (default and `-pg_scores`) and on r ≈ 1.56 D. There are 25/25 ctests.

On codon (r ≈ 7, 8 threads, PLAN_ODP CUR seed 0, container pd7-571940fcb5 without the bin table), outputs are
byte-identical to PLAN_ODP on D and H:

| | D pd6 | D pd7 | H pd6 | H pd7 |
|---|---|---|---|---|
| extraction consume (s) | 205.2 | 181.9 | 294.1 | 270.2 |
| interference resolver (s) | 266.2 | **159.2** | 151.8 | **108.0** |
| run total (s) | 622 | **482** | 677 | **601** |
| wall incl. `-decoys predict` | 18:34 | **16:16** | 19:40 | **18:16** |
| peak RSS (GB) | 38.3 | 37.6 | 34.0 | 33.2 |

On D, roughly 8 of the 16 minutes of wall time are still `-decoys predict`. Without that step the run itself is
8 min, against DIA-NN's 2.6-4.1 min on 16 threads. Memory (33-38 GB against DIA-NN's under 7 GB) is now the larger
gap.

## Predicted decoys once per library (19 Nov, OpenMS 08cda2c3fb)

Within one `OpenDIAWorkflow` call, `-decoys predict` already ran only once, before the loop over the `-in` runs. But
every scale and plan run here was a separate call (one per raw file and seed), and each call wrote the predicted
library into its own `out_dir`. So the same library was predicted again every time, about 8 min on r ≈ 7.

- **Reusing the library is output-identical.** Running the r ≈ 1.56 HeLa D search with `-lib` set to an earlier
  call's predicted library gives byte-identical outputs (16:13 → 5:23).
- **New `-decoy_cache <dir>` (default `out_dir`).** The predicted library is written there with a `.key` file holding
  the input library's SHA-1, the model directory and the OpenMS version, and a call with a matching key reuses it.
  Writes go to a temporary name (tagged with the process ID) that is renamed at the end, and the temporaries are
  removed if prediction fails.
- **Test, r ≈ 1.56 HeLa D, 8 threads, local.**
  - The first call predicts and takes 15:56. The second reuses the library and takes 5:53.
  - Both calls' outputs, and the cached library, are byte-identical to an uncached run.
  - All 25 ctests pass.
- **Caveat.** The key's OpenMS revision is the version string ("exported" in this build, not a git hash). A rebuild
  that changes the prediction code without a version bump would therefore reuse a stale library, so clear the cache
  in that case.

On r ≈ 7 D this should bring a reused call to about the 8 min run time (482 s on pd7), against DIA-NN's 2.6-4.1 min
on 16 threads. This has not been measured on codon yet.
