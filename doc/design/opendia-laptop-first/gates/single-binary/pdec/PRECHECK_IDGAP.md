# PRECHECK_IDGAP: where is the ID gap on the r ≈ 7 library? (19 Nov 2026, diagnostic)

This is a diagnostic on existing outputs, not a plan. DIA-NN and entrapment labels are used for evaluation only.

**Data.**
- r ≈ 7 HeLa D, seed 0, CUR, run with `-features_out` (`entlib-odpfeat-20261119`).
- 9.50 M windows, 24 features, the same library (sha256 76f2a6c2…) and the same output as PLAN_ODP CUR s0.
- DIA-NN 2.2 on the same targets gives 25,881 non-entrapment IDs at q ≤ 0.01.
- Scripts are `gap7/precheck_scorer.py`, `precheck_scorer2.py`, `precheck_rtrefit.py` and `calpanel.sh`.

## 1. A better learner on the same features does not help

Each row scores the windows, keeps the best window per precursor, runs paired target-decoy competition and evaluates.

| Scorer | N at own q ≤ 0.01 | lower / combined bound | N at FDP_c 1 % |
|---|---|---|---|
| A: tool score as written | 8,520 | 1.13 % / 1.29 % | 7,968 |
| B: label-free cross-fitted GBM, 3 Percolator rounds | 8,869 | 1.21 % / 1.38 % | 8,052 |
| C: GBM trained on **DIA-NN labels** (22,167 positive windows); a ceiling, not a method | 8,628 | 1.22 % / 1.39 % | 7,798 |

Even the label-cheating ceiling is no better than the tool. These 24 features cannot separate the missed DIA-NN IDs from
decoys, whatever the learner. **The gap is in the features, not in the scorer.**

## 2. Which features fail at low abundance

The table gives the AUC of DIA-NN-confirmed best windows (within 6 s of DIA-NN's RT, 21,276 precursors) against decoy
precursor best windows, by DIA-NN quantity quartile (0.5 means no separation).

| feature | Q1 low | Q2 | Q3 | Q4 high |
|---|---|---|---|---|
| tool score | **0.514** | 0.666 | 0.807 | 0.943 |
| transition_support | 0.766 | 0.800 | 0.816 | 0.823 |
| abs_library_rt_error (lower is better) | 0.243 | 0.234 | 0.233 | 0.239 |
| f12_libcos | 0.758 | 0.795 | 0.831 | 0.887 |
| f12_coelution | 0.671 | 0.766 | 0.844 | 0.946 |
| boundary_coelution | 0.686 | 0.759 | 0.830 | 0.920 |
| i_nclaim (lower is better) | 0.303 | 0.343 | 0.369 | 0.386 |
| byion_b / byion_y | 0.50 | 0.50 | 0.50 | 0.51 |
| f12_extra_support / f12_extra_coelution | 0.500 | 0.500 | 0.500 | 0.500 |

- **The combined score is at chance in the lowest quartile (0.514).** Single features still separate there: RT error
  0.76 and library cosine 0.76. So the combination is tuned to the abundant IDs that it trains on.
- **Four of the 24 features carry nothing on this library.**
  - `f12_extra_*` is constant because the library has 6 fragments, so there are no extra fragments.
  - `byion_b/y` is nonzero in under 1 % of windows.
  - Nothing describes MS1, ion mobility or the core-window mass error, although all three are collected.
- **The right window is not always present.** DIA-NN-confirmed windows exist for 75 % of Q1 IDs and 95 % of Q4 IDs.

## 3. RT is weak, and refitting the map does not fix it

- **The calibration map is built from few anchors.** RT residuals of DIA-NN-confirmed windows are p50 5.8 s and p95 22.6 s.
- **A refit from first-pass IDs barely helps.** Refitting the library-to-run map on the 8,520 first-pass IDs gives
  residuals of p50 4.8 s and p95 17.2 s. The RT-error AUC moves only from 0.76 to 0.78, and the label-free GBM with the
  refit error gets fewer IDs (4,739 at FDP_c 1 %).
- **The library's predicted RT is the limit.** Each target's RT is noisy against this run even with a perfect map.
  DIA-NN's window is about ±29 s (`RT window set to 0.49` min).

## 4. A larger calibration panel restores anchors but not the tolerance

`OpenDIACalibrate` on r ≈ 7 D, local, 8 threads:

| panel targets | used anchors | RT tolerance half width | residual p95 | wall | peak RSS |
|---|---|---|---|---|---|
| 40,000 (as the workflow) | 140 | 41.4 s | 8.2 s | 2:16 | 7.2 GB |
| 160,000 | 703 | 38.3 s | 12.3 s | 4:39 | 16.8 GB |
| 320,000 | 1,521 | 38.4 s | 12.5 s | 6:31 | 29.6 GB |

- More anchors make the map better sampled, but the tolerance stays near 38 s. Candidates and windows would barely
  change.
- A larger panel improves the scorer's seeds (138 positives now), not the search.
- Taken alone this is not the lever for the 3× gap.

## Conclusion

On large libraries OpenDIA is limited by its evidence, not by its learner. Three directions follow:
1. **Add the evidence it already collects but never scores:** MS1 area, ion mobility and core mass error.
2. **Get more out of the same 6 fragments.** DIA-NN was given exactly the same fragments (at most 6 per precursor,
   checked in `targets_r5.tsv`) and finds 3× more. So the fragment count is not the gap; how the evidence is used is.
3. **Weight the scorer towards low abundance.**

The first can be tested at once. The `-extra_scores` experimental block adds 12 such features and is checked as
PRECHECK_EVIDENCE, seed 0 at r ≈ 1.56. A frozen 5-seed plan follows only if it shows a gain.
