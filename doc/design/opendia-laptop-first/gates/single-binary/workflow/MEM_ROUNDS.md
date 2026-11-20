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

## Wall time (same-node A/B, codon 38539596)

This is one sbatch job on one node, running base, new, base, new in turn on r≈7 D s0 (default, `-threads 8`, no `-features_out`). All four outputs are byte-identical.

| Run | Wall | Extraction consume | Features | Peak RSS |
|---|---|---|---|---|
| base 1 | 9:12.6 | 192.2 s | 157.4 s | 35.28 GB |
| new 2 | 9:22.4 | 198.4 s | 160.6 s | 34.03 GB |
| base 3 | 9:15.1 | 191.2 s | 162.1 s | 35.60 GB |
| new 4 | 9:20.1 | 198.7 s | 161.7 s | 33.77 GB |

- **Cost:** wall time is +1.3 % (9:21.2 against 9:13.8 on the means). Extraction consume is +3.6 % (198.6 against 191.7 s), which is the serial slot allocation. Other stages are within noise.
- **Saving:** peak RSS is −1.5 GB (35.44 → 33.90 GB on the means).
- The earlier gaps (codon 12:55 against 9:31, local 6:49 against 5:45) were node and laptop load, not the change.

## What is left

- No easy output-identical win remains. Tensors are 97 % dense. The grids and InterferenceResolver windows accumulate in double, so float storage would change outputs.
- Reaching DIA-NN's under 7 GB needs a change that alters outputs (fewer open windows, tensor compression) or a structural one (per-isolation-window extraction passes). Either needs its own precheck and frozen plan.
