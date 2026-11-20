# PLAN_TRPPM: ppm fragment window for the trace tensor, on the r ≈ 7 entrapment (frozen 20 Nov 2026)

**Question.** Does matching the 6-fragment trace tensor within ±10 ppm, instead of the fixed ±0.02 Da (22-58 ppm at this library's fragments), raise IDs at equal true FDR on the large library without breaking own-q calibration?

**Why.** GAPDIAG showed where DIA-NN's IDs go:
- 34 % of DIA-NN's IDs (8,787) have our best window on the right peak but score low.
- On those windows our co-elution features fall to near-decoy values.
- Our trace areas track DIA-NN's quantity much worse than on our IDs (Spearman 0.54 against 0.78).

The trace window is about 4× the calibrated 8.3 ppm half width, so faint peaks collect foreign ions.

**What was seen before this plan.** All runs are r ≈ 1.56 HeLa D, with `-trace_tolerance_ppm` on the current default (OpenMS e19261618f).

| Arm | Seed | N @ q ≤ 1 % | E | Lower | Combined | N @ FDP_c 1 % |
|---|---|---|---|---|---|---|
| default (0.02 Da) | 0 | 13,629 | 78 | 0.57 % | 1.08 % | 13,388 |
| default (0.02 Da) | 1 | 13,962 | 90 | 0.64 % | 1.22 % | 13,545 |
| 8 ppm | 0 | 15,453 | 135 | 0.87 % | 1.65 % | 14,087 |
| **10 ppm** | 0 | 15,899 | 130 | 0.82 % | 1.55 % | **14,948** |
| **10 ppm** | 1 | 15,798 | 131 | 0.83 % | 1.57 % | **14,867** |
| 20 ppm | 0 | 15,290 | 124 | 0.81 % | 1.53 % | 14,253 |

- At 10 ppm, N at FDP_c is +11.7 % (s0) and +9.8 % (s1). The gain holds at FDP_c 0.5 % (11,576 → 13,126) and 2 % (14,536 → 16,603).
- **The own-q bound is worse.** The pooled combined bound goes from 1.15 to 1.56 %, and the tool accepts more at q ≤ 1 % than it can support. This is the same pattern as the earlier feature gains. So criterion 2 below is set as a limit on the bound and is not waived.
- **Choice of 10 ppm.**
  - It is about 1.2 × the run's calibrated half width (8.3 ppm here, 7.9 on r ≈ 1.56). So it is run-derived in spirit, but it is fixed in ppm for this plan.
  - It was chosen after seeing 8, 10 and 20 ppm on seed 0. It had the most IDs, so the choice is not blind. This is a disclosed weakness of the freeze.
- No r ≈ 7 run with `-trace_tolerance_ppm` exists before this plan. The limits below were set knowing the r ≈ 1.56 results above.

## Setup (fixed)

- **Binary:** OpenMS e19261618f, local build. It is deployed to codon `entlib-trppm-20261120/{bin,lib}`, with sha256 in `SHA256SUMS`:
  - OpenDIAWorkflow ef38c049…
  - libOpenMS 2aae4ce5…
  - libOpenMS_CLI a9a63661…
  - libOpenSwathAlgo 2e71282c…

  It runs inside `opendia-pd10-0942ffbdf3.sif` (sha256 fca04be9…) with these libraries first on `LD_LIBRARY_PATH`, as in the round-4 memory A/B. At `-trace_tolerance_ppm 0` the output is byte-identical to `gap7/ev/newdef` on r ≈ 1.56 D. 25/25 ctests pass.
- **Library:** the PLAN_ODP predicted library (sha256 76f2a6c2…), as `-lib`.
- **Arms:**
  - CUR is the current default.
  - NEW adds `-trace_tolerance_ppm 10`.
- **Runs:**
  - Seeds 0-4, D first (2737.d), 8 threads, 128 GB.
  - H (22474.d) is submitted only after the D verdict is committed, and is reported either way.
- **Evaluation:** `pdec/eval_pdec.py` (sha256 bb03387d…). Entrapment labels are read only there.

## Criteria (on D, pooled or median over the 5 seeds)

1. **Sensitivity:** NEW median N at FDP_c 1 % ≥ 1.05 × CUR median.
2. **Calibration not materially worse:**
   - NEW pooled lower bound at own q ≤ 0.01 ≤ CUR pooled lower bound + 0.0030.
   - NEW pooled combined bound ≤ CUR pooled combined bound + 0.0040.
3. **Pair fairness:** NEW P(entrapment beats its own decoy | best q ≤ 0.05) ≤ 0.56.
4. **Cost:** NEW median wall ≤ 1.3 × CUR median wall, and peak RSS ≤ 64 GB.

## Outcome

- **All pass:** report H. Making 10 ppm the default is then a separate commit, with ctests and a byte check.
- **Criterion 2 fails and the rest pass:** the sensitivity gain is real but the tool's q is too loose. Report it as is, keep the option at 0, and route the issue to the open own-q calibration item (b9). Nothing is retuned.
- **Any other fail:** report as is. The option stays at 0 and nothing is retuned.
