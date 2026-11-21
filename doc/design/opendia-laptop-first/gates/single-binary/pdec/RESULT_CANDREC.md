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

## Part B: result

The script is `precheck_candrec_b.py`, and the full output is in `EVAL_CANDREC_B.txt`. All 8 runs exited 0.

| Run | Arm | HeLa q ≤ 0.01 | Δ | N@FDP_c 1 % | Δ | Lower bound Δ | C recovered (near / picker / cap) | C accepted | Wall | RSS |
|---|---|---:|---:|---:|---:|---:|---|---:|---:|---:|
| D r0 | V8M8 | 19,883 | **+4.2 %** | | | | **87.8 %** (91 / 93 / n/a) | 539 | 1.45× | 8.6 GB |
| D r0 | M10 | 19,035 | −0.2 % | | | | 0.5 % | 2 | 2.06× | |
| H r0 | V8M8 | 35,661 | **+2.8 %** | | | | 51.5 % (56 / 50 / 50) | 1,286 | 1.24× | 13.3 GB |
| H r0 | M10 | 35,112 | +1.2 % | | | | 30.1 % (15 / 0 / 42) | 463 | 1.13× | |
| D r7 | V8M8 | 16,081 | +1.4 % | 15,952 | **+3.2 %** | **−0.14 pp** | 91.3 % | 495 | **1.78×** | 59.4 GB |
| D r7 | M10 | 15,269 | −3.7 % | 15,133 | −2.2 % | −0.10 pp | 0.9 % | 5 | 1.12× | |
| H r7 | V8M8 | 32,078 | +3.9 % | 31,678 | **+3.7 %** | **−0.05 pp** | 58.6 % | 998 | 1.15× | 40.1 GB |
| H r7 | M10 | 30,960 | +0.2 % | 30,604 | +0.2 % | −0.09 pp | 21.8 % | 200 | 1.22× | |

The baselines (def) are D r0 19,082, H r0 34,702, D r7 15,855 / 15,465 (lower 1.12 %), and H r7 30,885 / 30,543 (lower 1.20 %).

**Decision (best arm V8M8; mean r0 gain +3.5 %, against M10 +0.5 %)**

| Rule 1 criterion | D | H | |
|---|---|---|---|
| r0 HeLa ≥ +3 % | +4.2 % | **+2.8 %** | fails on H |
| r7 N@FDP_c 1 % ≥ +2 % | +3.2 % | +3.7 % | yes |
| r7 lower bound ≤ +0.2 pp | −0.14 pp | −0.05 pp | yes (it improves) |
| Wall ≤ 1.5× | 1.45× (r0), **1.78× (r7)** | 1.24× / 1.15× | fails on D r7 |

Rule 1 fails, narrowly, on two criteria: H r0 at +2.8 % against the +3 % bar, and D r7 wall time at 1.78×.

Rule 2 (≥ 70 % of C recovered) holds on D (88-91 %) but not on H (52-59 %). Rule 3 (< 50 %) does not hold anywhere. **Verdict: no plan from this precheck as written; record.**

**Reading**

- **D: candidates are recovered but few are accepted.** V8M8 recovers 88-91 % of D's C IDs. Only about 19 % of the recovered ones (539 of about 2,470) then pass q ≤ 0.01. This is the PRECHECK_NOCAND result again: once the right peak is a candidate, it scores below the decoy tail. It is an evidence problem (A_low), not a candidate problem.
- **H: cap and valley interact.** Raising the cap alone (M10) recovers only 42 % of cap-class IDs and nothing for picker-class IDs. Adding the valley change (V8M8) recovers about 50-60 % of every class. So on H many cap-class peaks are also merged, and a wider cap alone mostly adds more interference candidates. About 40-50 % of H's C IDs stay unrecovered even at V8M8. With a ±100-115 s window, the right peak may still not rank in the top 8 by co-elution.
- **V8M8 improves calibration at r7** on both runs (lower bound −0.14 / −0.05 pp). At r7 it gains +3.2 % (D) and +3.7 % (H) at true FDP 1 %. That is a real effect, near the size of the RT refinement. The cost is mostly D r7: 1.78× wall time and 59 GB RSS, from 1.9× the windows.
- **Wall time is the binding constraint, and it is tunable.** D's window count doubles because the valley change splits every precursor's trace, not just the missed ones. A cheaper variant (for example valley 0.7 with max 6, or the valley change only where a precursor has fewer than 3 windows) could keep most of the gain. It is not tested here, because any new arm needs its own gate before results.

**Next:** the Bt scorer port is in progress (RESULT_BT156 passed). Bt and V8M8 act on different classes (decoy-tail sharpening and candidate recovery), so their gains may stack. A combined arm is the natural follow-up once Bt is in C++. It needs its own precheck with the wall-time criterion fixed in advance.
