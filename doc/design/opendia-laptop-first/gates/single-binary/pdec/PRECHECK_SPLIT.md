# PRECHECK_SPLIT: does peak splitting cut real peaks, and is there a cheaper split? (21 Nov 2026, diagnostic)

**Question.** RESULT_CANDREC part B showed that `-valley_fraction 0.8 -max_peaks 8` (V8M8) recovers merged true peaks for +3.2 / +3.7 % at r7 FDP_c 1 %. It missed its gate on wall time (D r7, 1.78×) and on H r0 (+2.8 %). The valley rule cuts wherever the smoothed fragment-sum dips below `valley_fraction` × the lower maximum. At 0.8, a single real peak with a noisy dip can therefore be cut into two candidates. Only the best window per precursor is scored, so an accidental split shows up as a reported window that is narrower and nested inside the default window, with about the same apex and a smaller area. It would cost quantification, and it could cost IDs when the half-peak features score lower.

## Part A: split harm of V8M8 (no new runs)

- **Data:** the def runs (`ceil0-20261121/{D,H}_{r0,r7}`) and the V8M8 runs (`candrec-20261121/*_V8M8`), seed 0. Each run's `precursors.tsv` gives the best window per precursor (rt_apex, rt_start, rt_end, area).
- **Reference:** DIA-NN HeLa q ≤ 0.01 IDs, using RT, RT.Start, RT.Stop and Precursor.Quantity. Used for evaluation only.
- **Script:** `precheck_split_a.py`.
- **Common set:** HeLa IDs accepted (q ≤ 0.01) by both arms and by DIA-NN, whose def apex is within 6 s of DIA-NN's apex.
- **Split** (V8M8 relative to def) means all three of:
  - the apexes are within 3 s of each other;
  - the V8M8 width is < 0.7 × the def width;
  - the V8M8 window sits inside the def window (±1 s).
- **Measures:**
  - the split share in the common set;
  - for split IDs, the DIA-NN extent covered by each arm's window, and the median V8M8/def area ratio;
  - the Spearman correlation of log area with log DIA-NN quantity, per arm, on the common set;
  - for IDs lost by V8M8 (accepted by def and DIA-NN, not by V8M8), the share whose V8M8 best window is a split of the def window;
  - for IDs gained, the share whose V8M8 apex is within 6 s of DIA-NN and whose def apex is more than 6 s away (true recovery).

**Rule A (written before any result).** Splitting is harmless enough to use without a guard if **all** of these hold on each of the 4 run × level cells:

1. the split share in the common set is ≤ 5 %;
2. the area-quantity Spearman drops by ≤ 0.01;
3. split IDs make up ≤ 25 % of the lost IDs;
4. ≥ 60 % of gained IDs are true recoveries.

If any fails, splitting needs a guard. That would be a C++ change: keep two maxima apart only when the fragment composition differs across the valley, for example when the cosine of the fragment vectors at the two apexes is below a threshold. It gets its own precheck. Part B is then still read, but no arm from part B is planned without the guard.

## Part B: cheaper split arms (fixed before any part-B result)

| Arm | Flags |
|---|---|
| V7M6 | `-valley_fraction 0.7 -max_peaks 6` |
| V8M6 | `-valley_fraction 0.8 -max_peaks 6` |

- Seed 0, `-features_out`, D and H at r0 and r7: 8 runs per arm. They are compared against def and V8M8 with `precheck_candrec_b.py` (extended with the arm names only) and with the part A split measures.
- **Rule B:** the PRECHECK_CANDREC rule 1 bars apply:
  - r0 HeLa ≥ +3 % on D and H;
  - r7 N@FDP_c 1 % ≥ +2 % on D and H;
  - r7 lower bound rises by ≤ +0.2 pp;
  - wall time ≤ 1.5× on every cell.

  Rule A criteria 1-3 must also hold for the arm itself. An arm that passes goes to a frozen 5-seed r7 plan. If both pass, the one with the higher mean r7 gain goes.
