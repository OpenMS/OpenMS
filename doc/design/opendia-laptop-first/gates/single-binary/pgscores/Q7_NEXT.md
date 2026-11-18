# Q7: what should the follow-up plan for `-pg_scores` be? (diagnostic, before any new plan)

Written 18 November 2026.

## 1. The raw 02 "failure" is within one run's counting noise (`Q7_POWER.txt`)
- Seeds of one raw share the same data. The unit of evidence is therefore the run, and one Astral run at q ≤ 0.01 has
  about 1,950 IDs.
- At an exact 1 % FDR, that run is expected to hold about 9.7 entrapments.

| run | arm | median N | median entrapments | FDP_c | one-run 90 % interval |
|---|---|---|---|---|---|
| raw 01 | CUR | 1,709 | 7 | 0.0082 | 0.0038-0.0154 |
| raw 01 | NEW | 1,862 | 5 | 0.0054 | 0.0021-0.0113 |
| raw 02 | CUR | 1,811 | 4 | 0.0044 | 0.0015-0.0101 |
| raw 02 | NEW | 1,949 | 10 | 0.0103 | 0.0056-0.0174 |
| raw 03 | CUR | 1,849 | 7 | 0.0076 | 0.0036-0.0142 |
| raw 03 | NEW | 1,930 | 6 | 0.0062 | 0.0027-0.0123 |

Seed 0, pooled over raws 01-03: NEW 0.0073 (0.0049-0.0105), CUR 0.0088 (0.0061-0.0124).

- NEW's raw 02 value of 10 entrapments is what an exactly calibrated 1 % run gives half of the time (P = 0.51).
- Pooled over the three raws, NEW is no worse calibrated than CUR, and is nominally better.
- **The plan's rule "median FDP_c ≤ 0.010 on each raw" cannot tell an exactly calibrated scorer from a bad one.**
  - At a true FDR of 1.0 %, it fails 51 % of the time.
  - At 0.7 % it fails 15 %, and at 1.3 % it passes 19 %.
  - PRECHECK_PGQ's "pooled over seeds" interval treated seeds as independent. They are not, so that interval is too
    narrow.
- **Correction to PRECHECK_PGQ:** the raw 02 excess is not shown to be real. It is consistent with chance.

## 2. diaPASEF has no ground-truth FDR
- The HeLa library has no entrapment proteins: 0 of 14,216.
- The TimSim fixtures (`OpenMS-r1-benchmark-20261106/timsim-control0{1,2}`) hold only 1,000 target precursors, too few
  for a 1 % measurement.
- On D, the false count among NEW's 694 non-DIA-NN additions is bounded at about 270-500 (Q6), with a lower bound from
  DIA-NN's PEP and an estimate from target-decoy competition. The real number is not known.

## 3. What a follow-up plan needs
- A calibration safeguard (for example, a smaller hidden layer at higher dimension) has nothing to fix that the data can
  show. Testing it against the same per-raw rule would be a coin flip, so it is **not** proposed.
- What would settle the question is **a measurement with enough power**:
  1. **Astral:** pool entrapments over raws 01-03 as the frozen calibration criterion, with seed 0 per raw. That gives
     about 5,800 IDs and an interval of about ±0.0028. Add raws if more ProteoBench Astral files become available.
  2. **diaPASEF:** build an entrapment HeLa library, the same recipe as PLAN_PAIRLIB: an FDRBench shuffle of each target
     with its parent kept, predicted fragments and RT/IM. Then run D and H with it. This is the only way to get a true FDR
     on diaPASEF.
     - It needs the library predictor used for the Astral entrapment library, plus a .d-compatible IM prediction.
- **Proposed next plan, to be frozen only after the entrapment HeLa library exists:** `-pg_scores` as the default if
  all of these hold:
  - (a) N@FDP_c 1 % on Astral raws 01-03, pooled, ≥ CUR;
  - (b) pooled FDP_c at q ≤ 0.01 ≤ 0.010 on Astral;
  - (c) on HeLa D with the entrapment library, FDP_c at q ≤ 0.01 ≤ 0.012 and N@FDP_c 1 % ≥ CUR × 1.05;
  - with H and a fresh HeLa run held out.
- The DIA-NN-agreement proxy is dropped as a criterion. Q6 showed it penalises true IDs that DIA-NN misses.

## Status
- `-pg_scores` stays opt-in (a2aa8fc56c), and the default is unchanged.
- No new plan is frozen yet. It waits on the entrapment HeLa library, which is the blocking item.
