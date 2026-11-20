# PRECHECK_IMFEAT: mobility-aware evidence for faint peaks, hand features vs a 2-D trace net (20 Nov 2026, diagnostic)

This is a diagnostic, not a plan. DIA-NN labels are used for evaluation only, and the GBM arms are ceilings. It is
step 1 of the agreed order: test whether "use mobility properly" goes beyond the narrow window (PLAN_IMNARROW), and
whether a net is worth building in C++.

**Data.**
- r ≈ 1.56 HeLa D, seed 0, on the new default (scale 0.6, `ev/ims06` features).
- A mini library of 6,000 DIA-NN non-entrapment q ≤ 0.01 targets (about 1,500 per quantity quartile) and 30,000 random decoy precursors.
- The dev-tool search uses the run calibration with the 1/K0 half width scaled to 0.6. Its peaks are identical to the workflow's for all 36,000 precursors (139,070 windows), and 139,068 join to the 36 workflow features.
- Raw MS2 points of the 6 search fragments (±0.08 1/K0) come from a `.xipm` export.
- Ran on codon (job 38339545). An earlier run (38306769) used peaks searched at the unscaled width, so only half the windows joined. It is superseded and its numbers are not used.
- Scripts: `faint2_mini.py`, `precheck_imfeat_extract.py`, `precheck_imfeat.py`, `imfeat.sbatch`. Output in `EVAL_IMFEAT.txt`.

**Candidate inputs, per candidate window.**
- **faint2:** `p_cos_apex` (library-pattern cosine at the apex) and `p_local_sn` (apex against off-peak background), on the scaled window.
- **band6:** the window's *own* mobility apex m\*, the mode of the intensity-weighted 1/K0 histogram over apex ±2 scans. Six features describe it:
  - `b_offset`: |m\* − calibrated 1/K0|;
  - `b_frag_sd`: spread across fragments of their 1/K0 centroids near m\*;
  - `b_cos_apex`: library cosine within ±0.01 of m\*;
  - `b_support`: fragments present within ±0.01 of m\* at the apex;
  - `b_contrast`: signal in ±0.01 of m\* against 0.02-0.05 away;
  - `b_local_sn`: apex against off-peak background on the ±0.01 band trace.
- **CNN1d / CNN2d:** label-free CNNs on 6 fragments × 33 scans, either mobility-summed (1-D) or × 8 mobility bins over ±0.031 (2-D). Positives are the tool's q ≤ 0.01 target winners and negatives are decoy windows. Cross-fitted on the workflow folds.

## Results

TPR of the 5,350 DIA-NN targets' right windows against the decoy precursors' best windows, at decoy FPR
1e-3 / 3e-3 / 1e-2. Pairs are two GBM seeds.

| Arm | Seed 0 | Seed 1 | Q1 @ 3e-3 |
|---|---|---|---|
| Tool score | 0.586 / 0.646 / 0.682 | | 0.352 |
| CNN1d alone | 0.161 / 0.261 / 0.384 | | 0.096 |
| CNN2d alone | 0.439 / 0.538 / 0.646 | | 0.251 |
| GBM 36 | 0.623 / 0.688 / 0.755 | 0.607 / 0.685 / 0.751 | 0.378 / 0.371 |
| 36 + faint2 | 0.625 / 0.684 / 0.759 | 0.613 / 0.687 / 0.758 | 0.368 / 0.379 |
| **36 + band6** | **0.656 / 0.711 / 0.779** | **0.652 / 0.709 / 0.776** | **0.417 / 0.415** |
| 36 + faint2 + band6 | 0.646 / 0.710 / 0.778 | 0.659 / 0.710 / 0.773 | 0.414 / 0.413 |
| 36 + CNN1d | 0.617 / 0.683 / 0.758 | 0.619 / 0.688 / 0.761 | 0.385 / 0.395 |
| 36 + CNN2d | 0.622 / 0.706 / 0.794 | 0.598 / 0.700 / 0.792 | 0.429 / 0.424 |
| 36 + all (faint2 + band6 + CNN2d) | 0.619 / 0.717 / 0.804 | 0.622 / 0.725 / 0.804 | 0.448 / 0.461 |

Proceed rule, fixed in the script before running: a component goes on if it adds ≥ +0.020 TPR at 3e-3 over 36 in both
seeds.

| Arm | Δ seed 0 / seed 1 | Result |
|---|---|---|
| 36 + faint2 | −0.004 / +0.003 | no |
| **36 + band6** | **+0.023 / +0.024** | **proceed** |
| 36 + faint2 + band6 | +0.022 / +0.025 | proceed |
| 36 + CNN1d | −0.005 / +0.003 | no |
| 36 + CNN2d | +0.018 / +0.015 | no, below the bar |
| 36 + all | +0.029 / +0.040 | proceed |

## Reading

- **Once the window is narrow, the faint2 features carry nothing new.** Their PRECHECK_FAINT gain at the wide window came from the background the narrow window now removes.
- **The window's own mobility apex is new information.** band6 adds +2.3-2.4 points at 3e-3 and +4 points in Q1. This is "use mobility properly" beyond the narrower window, and it is the hand-feature route that passes.
- **A trace net needs mobility.** CNN1d alone is near useless (0.26) and adds nothing. CNN2d alone reaches 0.54, closer to the tool than any earlier net. Stacked, it adds +1.5-1.8 points, mostly at the looser 1e-2 cut (+0.04), and Q1 +0.05.
- **Everything together** gives the largest gain (+2.9-4.0 points, Q1 +0.07-0.09), but only part of it is beyond band6.
- **Caveats.**
  - These are DIA-NN-labelled ceilings on one run. A label-free scorer realises less.
  - The band6 features here are computed from raw points (xipm), which the workflow does not keep. A check of the C++-feasible version, a per-fragment 1/K0 histogram over the core window built from the binned points, is running (codon 38358442) before any code.

## Decision

- **Step 2 candidate:** the band6 evidence, implemented as additional features in the extraction pass (per-fragment 1/K0 histogram in the core window), if the C++-feasible variant keeps the gain. It goes into one frozen plan.
- **Net:** not yet. CNN2d alone is below the bar, and its stacked gain partly overlaps band6. A C++ 2-D net would be a later, separate plan if band6 lands and the remaining gap justifies it.
