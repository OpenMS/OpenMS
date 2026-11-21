# RESULT_CANDREC: why DIA-NN's peaks get no OpenDIA candidate (21 Nov 2026, diagnostic)

The gate is in `PRECHECK_CANDREC.md` (`1ecf4e88f3`). The part A script is `precheck_candrec_a.py`, and the full output is in `EVAL_CANDREC_A.txt`.

## Part A: causes

C IDs are DIA-NN HeLa q ≤ 0.01 IDs that the tool misses and that have no window within 6 s of DIA-NN's apex.

| Class | D r0 (2,809) | D r7 (2,995) | H r0 (5,906) | H r7 (4,565) |
|---|---:|---:|---:|---:|
| near: the nearest window is 6-15 s away | **71.3 %** | **66.2 %** | 38.8 % | 44.7 % |
| picker: inside the window, < 5 windows, none near | 22.7 % | 32.0 % | 1.9 % | 7.1 % |
| cap: inside the window, 5 windows, none near | 0.0 % | 0.3 % | **57.7 %** | 44.7 % |
| rt_out: outside the refined RT window | 2.9 % | 0.9 % | 0.7 % | 3.3 % |
| no window at all | 3.0 % | 0.6 % | 0.8 % | 0.2 % |

**Supporting numbers**

| | D r0 | D r7 | H r0 | H r7 |
|---|---:|---:|---:|---:|
| Refined RT half-window | ±23.6 s | ±31.7 s | ±115.5 s | ±99.0 s |
| All DIA-NN IDs outside the RT window | 0.6 % | 0.2 % | 0.1 % | 0.5 % |
| Windows per precursor, C-missed (p50) | 2 | 2 | 5 | 5 |
| C-missed with 5 windows | 2.6 % | 6.7 % | **96.3 %** | 87.5 % |
| Nearest window for C-missed (p25 / p50 / p75) | 8.0 / 10.9 / 14.9 s | 8.7 / 11.8 / 17.4 s | 10.5 / 18.1 / 31.5 s | 10.5 / 16.2 / 26.7 s |
| DIA-NN peak width (p50), C-missed vs accepted | 7.9 vs 6.9 s | 7.0 vs 6.9 s | 8.6 vs 8.6 s | 9.5 vs 8.6 s |

**Reading**

- **RT windows are not the cause.** At most 0.6 % of DIA-NN IDs fall outside the refined window, and rt_out is at most 3.3 % of C.
- **D: merged or shifted peaks.** The C IDs have few windows (median 2), and the nearest one sits about 11 s away, roughly 1.5 DIA-NN peak widths. This is the PRECHECK_NOCAND pattern: the valley rule (0.5) merges the right maximum into a neighbour and keeps the neighbour's apex. On D, the C share falls steeply with windows per precursor (25 % at n = 1 down to 1.4 % at n = 5). The picker is simply not splitting.
- **H: the 5-candidate cap.** H's half-window is about 4-5× D's, because H's RT error is larger (p95 56-59 s). Almost every precursor fills all 5 candidate slots (44,172 of 44,866 DIA-NN IDs at r0). The C IDs are crowded out by up to 5 stronger co-elution maxima elsewhere in a ±115 s window. The cap ranks by fragment co-elution, so with 5 slots in a 230 s window a faint true peak loses to brighter interference.

## Part B: arms fixed before any part-B result

PRECHECK_CANDREC names `-valley_fraction 0.8 -max_peaks 8` for picker- or cap-dominant classes. "near" on D is the valley-merge signature, which the same arm targets. Because the cap alone dominates H, a cap-only arm is added to separate the two effects. It is fixed here, before any part-B run:

| Arm | Flags | Targets |
|---|---|---|
| def | none (the PRECHECK_CEIL0 `-features_out` runs) | |
| V8M8 | `-valley_fraction 0.8 -max_peaks 8` | D near/picker and H cap |
| M10 | `-max_peaks 10` | H cap only |

- All arms run with `-features_out`, seed 0, on D and H at r0 and r7: 8 new runs.
- **Recovery** is measured as the share of part-A C IDs that have a window within 6 s in the arm.
- **Gains** are HeLa IDs at own q ≤ 0.01 (r0 and r7) and N at FDP_c 1 % with bounds (r7). Wall time is compared with the def `-features_out` runs:

| Default wall time | D r0 | H r0 | D r7 | H r7 |
|---|---:|---:|---:|---:|
| def | 2:54 | 7:22 | 21:29 | 32:25 |

The rule is as written in `PRECHECK_CANDREC.md`. "The best arm" means the arm with the higher mean r0 gain across D and H.
