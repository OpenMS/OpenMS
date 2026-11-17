# Why HeLa diaPASEF regressed and where the gap to DIA-NN is (2026-11-17, diagnostic)

Run, library and references as in `RESULT_HELA.md`. Scripts (in this folder): `funnel_hela.py`,
`diag2_hela.py` (per-candidate features, `-features_out`) and `coord_hela.py` (RT/1/K0 against DIA-NN).
Counts are each tool's q ≤ 0.01. Funnel groups are DIA-NN's q ≤ 0.01 library precursors (28,568), sorted by our outcome
(right peak = our apex within 3.4 s, half DIA-NN's median peak width, of DIA-NN's RT).

| OpenDIA variant | IDs | accepted | right peak, below cut | wrong peak | not scored |
|---|---|---|---|---|---|
| as benchmarked (eaa2c95cd0) | 11,653 | 10,696 | 9,923 | 7,877 | 72 |
| + calibrated 1/K0 in the evidence pass (0e9e4fc3a5) | **13,290** | 12,189 | 8,971 | 7,336 | 72 |
| 0e9e4fc3a5 with the 7 Nov calibration (`-calibration p40000`) | 13,984 | 13,049 | 10,821 | 4,678 | 20 |
| reference: DIA-NN / OpenSWATH SVM / XGBoost | 28,950 / 19,132 / 19,725 | | | | |

## Cause 1, fixed: uncalibrated 1/K0 in the evidence pass (bug)
- The search used the calibrated precursor 1/K0. The evidence pass (TraceCollector, 6-fragment traces, MS1 and by-ion
  evidence) used the **uncalibrated** library value inside a ±0.03 1/K0 window.
- HeLa's calibrated offset is +0.022 1/K0, so the window sat off-centre: 16-43 % of DIA-NN's IDs fell outside it.
- It entered with the single-binary workflow (794ee1565e). OpenDIACNN's `-apply_calibration_coordinates`, used for the
  7 Nov result, did apply it.
- Fix 0e9e4fc3a5: +1,637 IDs. No effect without library 1/K0 (Astral `has_im` is false).

## Cause 2, open: the RT map at the gradient start (calibration, 8d2771bcec, 12 Nov)
- 6,493 of the wrong-peak DIA-NN IDs never became a candidate. 49 % of them lie outside the RT window, all on the
  early side, and 2,735 sit in the lowest library-RT decile.
- The current isotonic map puts the earliest peptides about 19 s late. Library RT 0-10 maps to 110-122 s, but DIA-NN
  finds them at 79-103 s. The first knot is at 110.7 s, and DIA-NN IDs start at 71 s.
- The 7 Nov 20-bin map is not biased there: DIA-NN-ID |ΔRT| p95 is 5.9 s against 19.3 s now, and the lowest-decile
  median is −1.7 s against −19.3 s.
- With the old calibration, wrong-peak IDs fall from 7,336 to 4,678 and the count rises to 13,984.
- Cause: on HeLa the single pass succeeds, so there is no windowed pass. Stage 3 searches around the EM map with a
  ±30 s window, and no anchors are found before about 110 s. The isotonic fit then extrapolates flat over a 30 s
  stretch with no anchors. The run-edge extension (8d2771bcec) only acts in the windowed pass, and the calibration
  tolerance (11.3 s, from anchors inside the support) does not see this bias.
- This was designed on Astral. It is an RT-calibration issue at the gradient ends, of the same kind as
  PRECHECK_RTPROFILE.

## Cause 3, the main remaining gap: scoring on diaPASEF
Even with both fixes, about 10,800 DIA-NN IDs have our right peak below the cut. Per-feature AUC, right peak against
score-matched decoys:
- The 12-fragment block is empty here. There is no 12-fragment scoring library for HeLa, so only the 6 library fragments
  are used: `f12_extra_*` is always 0, `f12_support` is nonzero in 33 % of windows, and `d_lib` in 0.8 %.
- `byion_b`/`byion_y` are nonzero in < 1 % of windows. Their intensity floor (300, a fixed constant set on Astral) is
  above most diaPASEF fragment intensities (DIA-NN median precursor quantity about 10,000 summed over the peak).
- The grid, interference and deficit features pool all ion mobility (`OpenDIAFragmentGrid.h`: "mobility is ignored").
  On diaPASEF, precursors co-isolated at other 1/K0 add to these sums. This is likely to blur them, but it is not
  measured yet.
- So 6 of the 24 inputs are nearly constant on HeLa (`f12_extra_*` x3, `byion_*` x2, `d_lib`), and the grid and
  interference block is mobility-blind. DIA-NN and OpenSWATH both score in the mobility dimension.

## What it means
- The regression from 13,253 (7 Nov) came from cause 1 plus cause 2. With both fixed the current code is at about
  13,300-14,000, level with or slightly above 7 Nov. It is still at 70-73 % of OpenSWATH and 46-48 % of DIA-NN.
- The rest of the gap is most likely diaPASEF-blind evidence (cause 3). That is a hypothesis from the feature
  coverage, not yet tested.
- Candidate next steps (sizes not yet measured):
  1. Mobility-resolved grid and interference (extract 12-fragment grids within the calibrated 1/K0 window).
  2. Derive the by-ion floor from the run's intensity distribution, not 300.
  3. Fix the early-gradient RT map in the single-pass calibration (run-edge search or an anchor-free extrapolation
     check).
  4. A 12-fragment scoring library for HeLa.
  Each needs a frozen plan with a held-out diaPASEF run (no entrapment exists for HeLa).
