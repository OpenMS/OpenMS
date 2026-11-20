# PRECHECK_NOCAND: why 17.8 % of DIA-NN IDs have no candidate, and the cheap fixes (19-20 Nov 2026, diagnostic)

These are diagnostics, not a plan. DIA-NN and entrapment labels are used for evaluation only. This follows up on
PRECHECK_SCORER36, which left two open questions: 4,605 DIA-NN IDs at r ≈ 7 D with no window near DIA-NN's RT, and a
+5 % seed-0 lead from hard-negative training.

## 1. The cause: the valley rule merges the right peak into a neighbour

**Data.**
- r ≈ 7 HeLa D, seed 0, default 07131138b8: the features and `calibration.json` from PRECHECK_SCORER36.
- DIA-NN 2.2: 25,881 non-entrapment IDs at q ≤ 0.01.
- Script: `precheck_nocand.py`.

| Cause of the 4,605 missed IDs | n | share |
|---|---|---|
| DIA-NN RT outside the calibrated RT window (±41.4 s) | 19 | 0.4 % |
| DIA-NN 1/K0 outside the mobility window (±0.056) | 1 | 0.0 % |
| Inside both, and the precursor already has the maximum 5 candidates | 521 | 11.3 % |
| Inside both, fewer than 5 candidates, none at DIA-NN's RT | 4,064 | 88.3 % |

The windows are not the problem: only 0.2 % of all DIA-NN IDs fall outside the RT window. Most misses are faint,
covering 25 % of Q1, 25 % of Q2, 16 % of Q3 and 5 % of Q4.

**What the search actually sees.**
- Data: a mini library of sampled DIA-NN IDs (600 picker misses, 150 cap misses, 300 found), searched with the dev tool `OpenDIA` using the run's calibration and `-traces_out`. Scripts: `nocand_mini.py`, `precheck_nocand2.py`.
- The dev tool reproduces the workflow's candidate apexes exactly for all 1,050 precursors.

| Measure | Picker misses | Cap misses | Found |
|---|---|---|---|
| Consensus is a smoothed local maximum at DIA-NN's RT | 0.92 | 0.86 | 0.99 |
| Fragments with signal at that apex ±1 scan (median) | 6 of 6 | 6 | 6 |
| Height at DIA-NN's RT / best picked apex (median) | 0.63 | 0.63 | 1.00 |
| Valley to the nearest picked apex above 0.5 × the lower maximum (so merged by the picker) | **0.875** | 0.65 | |
| DIA-NN's RT inside a picked peak region whose apex is elsewhere | 0.04 | 0.01 | 0.93 |

The right peak is in the trace, with all 6 fragments, as a local maximum at about 0.6 of the neighbouring picked apex,
typically about 15 scans away. The picker's valley rule (`valley_fraction` 0.5) merges the two maxima into one peak and
keeps the taller apex. So the right peak never becomes a candidate.

**Picker settings on the same mini library** (share of precursors with a right candidate):

| Search setting | Picker misses | Cap misses | Found | Candidates per precursor |
|---|---|---|---|---|
| Default (valley 0.5, max 5) | 0.00 | 0.00 | 1.00 | 3.1 |
| valley 0.7, max 5 | 0.70 | 0.63 | 0.99 | 4.8 |
| valley 0.8, max 5 | 0.76 | 0.73 | 0.96 | 5.0 |
| valley 0.8, max 8 | 0.90 | 0.90 | 1.00 | 7.6 |
| valley 1.0, max 8 | 0.89 | 0.90 | 0.99 | 8.0 |

## 2. End to end: more right candidates, but few become IDs

**Data.**
- r ≈ 1.56 HeLa D, seed 0, local, cached predicted library.
- New advanced options `-valley_fraction` and `-max_peaks` (822ef8e362). At their defaults the output is byte-identical (vdef = def2), and 25/25 ctests pass.
- Evaluator `eval_pdec.py`. Script: `precheck_valley.py`, runs in `gap7/ev/{vdef,v08m8,v05m8}`.

| Arm | N at FDP_c 1 % | q ≤ 0.01 lower / combined bound | Pair win | Windows | Wall / RSS |
|---|---|---|---|---|---|
| Default | 12,357 | 0.79 % / 1.49 % | 0.524 | 2.39 M | 5:03 / 11.0 GB |
| valley 0.8, max 8 | **12,578 (+1.8 %)** | 0.80 % / 1.52 % | 0.524 | 5.17 M | 9:37 / 18.3 GB |
| valley 0.5, max 8 | 12,146 (-1.7 %) | 0.74 % / 1.40 % | 0.508 | | |

Outcomes for the DIA-NN IDs:

| | Default | valley 0.8, max 8 |
|---|---|---|
| Right candidate present | 84.1 % | **98.4 %** |
| Q1 / Q2 / Q3 / Q4 | 77.7 / 77.8 / 85.2 / 95.7 % | 97.5 / 98.1 / 98.6 / 99.7 % |
| Accepted with the right window | 47.3 % | 48.7 % |

- 3,744 IDs gain a right candidate, but only 11 % of those are then accepted.
- Accepted-right IDs: 779 gained and 421 lost, since the extra candidates are also extra competition.
- The candidate problem is mostly solved. What remains is the same as in PRECHECK_SCORER36: once recovered, these faint peaks score below the decoy tail.
- This matches PLAN_IM2D, where a different way of recovering candidates also gave about +1-3 %.
- At about twice the wall time and 1.7× the memory for +1.8 %, a plan is not worth freezing on this alone.

## 3. Hard-negative training does not hold up at r ≈ 1.56

The same `precheck_scorer36_tailfit.py` was run on r ≈ 1.56 D s0 with the default features: 2.39 M windows, the
hardest 2 % of decoys (6,904 precursors).

| Scorer | N at own q ≤ 0.01 | lower / combined | N at FDP_c 1 % |
|---|---|---|---|
| Tool | 13,060 | 0.80 / 1.52 % | 12,357 |
| B_tail round 1 / 2 / 3 | 13,762 / 13,973 / 13,786 | 0.86 / 0.99 / 1.02 % lower | 12,505 / 12,251 / 12,311 |
| C_tail (DIA-NN labels) | 12,634 | 0.82 / 1.54 % | 11,890 |

- The best round gains +1.2 %, and later rounds lose IDs.
- Calibration gets worse each round: the lower bound goes 0.80 → 1.02 %. More of the extra N at q ≤ 0.01 is false.
- The +5 % at r ≈ 7 does not reproduce, so the lead is **dropped**. No plan.

## Decision

- **Candidate generation is no longer the main gap.** A wider valley rule gets 98 % of DIA-NN IDs a right candidate, but they still score too low. Both remaining levers point at **evidence for faint peaks**:
  - per-scan fragment agreement across the peak;
  - fragment ratios on the best mobility band;
  - MS1 isotope fit.
- `-valley_fraction 0.8 -max_peaks 8` is a candidate companion for any future faint-peak feature work, since the features need the right candidate in order to score it. It is kept as an advanced option at the old defaults.
- Hard-negative training is dropped.
