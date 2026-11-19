# PRECHECK_SCORER36: is the scorer the limit now that the default has 36 features? (19 Nov 2026, diagnostic)

This is a diagnostic on one existing default run, not a plan. DIA-NN and entrapment labels are used for evaluation
only. It repeats PRECHECK_IDGAP section 1 (24 features) on the current default, which includes the PLAN_EXTRA evidence
block.

**Data.**
- r ≈ 7 HeLa D, seed 0, default of 07131138b8 (container pd9), run with `-features_out` (`entlib-feat36-20261119`).
- 10:13 wall time and 40 GB on 8 threads, using the cached predicted library (sha256 76f2a6c2…).
- The `precursors.tsv` is byte-identical to PLAN_IM2D CUR s0, which equals PLAN_EXTRA NEW s0.
- 9.50 M windows and 36 features. The features file has sha256 d7f54f50….
- DIA-NN 2.2 on the same targets gives 25,881 non-entrapment IDs at q ≤ 0.01.
- Scripts: `precheck_scorer36.py`, `precheck_scorer36_tail.py`, `precheck_scorer36_tailfit.py`. The method is the same as `precheck_scorer.py`: keep the best window per precursor, run paired target-decoy competition, use 3-fold cross-fitting, and use HistGBM with 300 rounds and 31 leaves.

## 1. A standard learner on the same features does not help

| Scorer | N at own q ≤ 0.01 | lower / combined bound | N at FDP_c 1 % | vs tool |
|---|---|---|---|---|
| A: tool score as written | 11,913 | 1.25 % / 1.43 % | **11,258** | |
| B: label-free GBM, 3 Percolator rounds (rounds 1/2/3 give 11,366 / 11,280 / 11,161) | 11,679 | 1.19 % / 1.36 % | 11,161 | -0.9 % |
| B2: label-free GBM with twice the capacity (600 rounds, 63 leaves) | 11,923 | 1.15 % / 1.31 % | 11,368 | +1.0 % |
| C: GBM trained on **DIA-NN labels** (22,167 positive windows); a ceiling, not a method | 11,648 | 1.10 % / 1.25 % | 11,122 | -1.2 % |

This is the same picture as with 24 features (tool 7,968, ceiling 7,798). The evidence block moved the tool score and
the ceiling together, by about +41 %, and the learner adds nothing. Even DIA-NN's labels do not lift a standard learner
above the tool.

## 2. Control: learners that target the decoy tail

Section 1 trains against all decoy windows, but IDs at 1 % are decided only by the few decoys that score highest. In
this control the negatives are only the windows of the hardest 2 % of decoy precursors, ranked by tool score: 24,589
precursors and 96,134 windows.

| Scorer | N at own q ≤ 0.01 | lower / combined bound | N at FDP_c 1 % | vs tool |
|---|---|---|---|---|
| C_tail: DIA-NN labels against hard decoys (ceiling) | 11,410 | 0.98 % / 1.12 % | 11,314 | +0.5 % |
| B_tail: label-free against hard decoys, round 1 / 2 / 3 | 12,475 / 12,652 / 12,590 | 1.30 / 1.23 / 1.37 % lower | 11,821 / 11,880 / 11,856 | **+5.0 to +5.5 %** |

**Reading.**
- The tail-focused label-free learner gains about 5 % at a true FDP of 1 %.
- That gain is a property of the decoy tail, not new information. With DIA-NN labels the tail-focused ceiling gains only +0.5 %.
- A plausible cause: training on the tool's own winners against near-winner decoys sharpens the boundary where the 1 % cut falls, while DIA-NN's faint positives add only noise there.
- Calibration does not improve. The lower bound is 1.23-1.37 %, against 1.25 % for the tool.
- This is seed 0 on one run, with no seed spread. The +5 % is a candidate for its own precheck on r ≈ 1.56 plus a frozen plan. It is not evidence for anything larger.

## 3. Where the remaining gap is

AUC of DIA-NN-confirmed best windows against decoy best windows, by DIA-NN quantity quartile:

| Score | Q1 low | Q2 | Q3 | Q4 high |
|---|---|---|---|---|
| tool, 24 features (PRECHECK_IDGAP) | 0.514 | 0.666 | 0.807 | 0.943 |
| tool, 36 features | 0.543 | 0.700 | 0.831 | 0.949 |
| B label-free GBM | 0.815 | 0.893 | 0.940 | 0.984 |
| C DIA-NN-labelled (trained on these positives, so optimistic) | 0.958 | 0.982 | 0.990 | 0.998 |

The large AUC gains from B do not become IDs. The 1 % cut sits at a decoy false-positive rate of about 1e-4, and the
bulk AUC says nothing about that region. TPR of the tool score at a fixed decoy FPR:

| Decoy FPR | Q1 | Q2 | Q3 | Q4 | All confirmed |
|---|---|---|---|---|---|
| 1e-4 | 0.220 | 0.408 | 0.618 | 0.867 | 0.528 |
| 1e-3 | 0.320 | 0.522 | 0.713 | 0.914 | 0.617 |
| 1e-2 | 0.413 | 0.609 | 0.779 | 0.935 | 0.684 |
| 5e-2 | 0.454 | 0.644 | 0.800 | 0.939 | 0.709 |

The count at FPR 1e-4 is 11,242 confirmed precursors, which matches the tool's 11,258 at FDP_c 1 %. Two separate
losses explain the 25,881 → 11,258 gap:

1. **Candidate generation: 4,605 DIA-NN IDs (17.8 %).** These have no window within 6 s of DIA-NN's RT, although every one has some window elsewhere. Most are faint: 25 % of Q1, 25 % of Q2, 16 % of Q3 and 5 % of Q4. No scorer can recover them.
2. **Faint-peak evidence: about 10,000 IDs.** The right window exists but scores below the decoy tail (TPR at FPR 1e-4 is 0.22 in Q1 and 0.41 in Q2). DIA-NN's labels cannot teach a learner to separate them with these features (section 1, C).

## Decision

- **Main lever: features and candidates, not the scorer.** The next work is faint-peak evidence. The candidates from the earlier summary are per-scan fragment agreement across the peak, fragment ratios on the best mobility band, and MS1 isotope fit. Also needed is a look at why 17.8 % of DIA-NN IDs have no window near the right RT; this is a different question from PRECHECK_IMPICK's 16 % with no right candidate at r ≈ 1.56.
- **Secondary, cheap:** tail-focused (hard-negative) training of the existing scorer, +5 % at seed 0. It needs an r ≈ 1.56 precheck and then a frozen 5-seed plan before any default change. Calibration must not get worse.
- Nothing in this precheck changes code or defaults.
