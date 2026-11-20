# MEM_ROUNDS: output-identical peak-RSS work (20 Nov 2026)

All four rounds are byte-identical in `precursors.tsv`, `<run>.precursors.tsv`, `calibration.json` and `features.tsv`, with 25/25 ctests. None of them needs a plan.

| Round | OpenMS | Change | r≈1.56 D peak | r≈7 D peak |
|---|---|---|---|---|
| 1 | 85c794b211 | FragmentGridCollector: chunked grid pool instead of a growing vector | 11.28 → 10.85 GB | |
| 2 | d0bdf24540 | release the calibrated library copy before extraction | 10.85 → 10.64 GB | |
| 3 | f3f1fc046f | TraceCollector::Sample keeps only apex_rt, and the window id strings are freed before extraction | 10.64 → 10.17 GB | |
| 4 | 0a49d8aa8e | streaming TraceEvidence: slots only while a window is open, keeping the b/y counts after close | about 0 (10.53 → 10.56 GB, local) | 35.57 → 34.11 GB (codon) |

Round 4 r≈7 check: codon jobs 38475115 (base) and 38496478 (new) in `entlib-mem4-20261120`. I re-checked the byte comparison myself (cmp SAME on all four files). Locally, r≈1.56 D against `gap7/ev/newdef` is also byte-identical.

## Peak breakdown at r≈7 D (round 4 base, extraction peak 34.7 GB)

| Item | Size |
|---|---|
| trace tensors, 1.98 M open × 6 KiB | 12.2 GB |
| library and precursors | about 7 GB |
| fragment grids, 1.96 M open × 3 KiB (double) | 6.0 GB |
| InterferenceResolver windows, 10.1 M × 384 B | 3.9 GB |
| TraceEvidence, 10.1 M × 184 B (removed by round 4) | 1.9 GB |
| RunWindow | 1.1 GB |

## Wall time

- Round 4 timings are noisy. On codon the new run took 12:55 against 9:31 for the base. Every stage was slower in proportion, including calibration, which round 4 does not touch (20 → 30 s), so that points to the node.
- The local r≈1.56 run took 6:49 against 5:45. Calibration was the same (27.3 against 28.0 s), but search, extraction and features were 20-25 % slower. That run is a single one on a laptop with a load average of about 7 from other processes.
- Extraction does more per-hit work now (serial slot allocation), so some cost is possible. A same-node A/B on codon is needed before the cost is stated. It is not claimed as zero.

## What is left

- No easy output-identical win remains. Tensors are 97 % dense. The grids and InterferenceResolver windows accumulate in double, so float storage would change outputs.
- Reaching DIA-NN's under 7 GB needs a change that alters outputs (fewer open windows, tensor compression) or a structural one (per-isolation-window extraction passes). Either needs its own precheck and frozen plan.
