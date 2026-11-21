# PRECHECK_RTIM: Sage-style per-run linear correction of library RT and 1/K0 (20 Nov 2026, diagnostic)

**Question.** H's RT error is about 10 % of its gradient, the same share as D's, and both runs' maps sit within about 1 point of the best monotone map (RTCAL D/H diagnostic, `rtcal_dh.py`). So a better map cannot help, but a better per-peptide prediction might. As Sage does, fit per run a ridge regression on a featurised sequence using that run's own confident IDs, and use it to correct the library RT and 1/K0.

## Setup

- **Runs:** r ≈ 7 library, new default (OpenMS 01c5891ce0), seed 0, D and H (`entlib-mt-20261120/runs/{D,H}_mt_s0`).
- **Training pool:** the run's own OpenDIA target IDs at q ≤ 0.01. The fit is 5-fold cross-fitted by stripped sequence, so no precursor is predicted by a model that saw it.
- **RT target:** OpenDIA apex RT minus the run's current map. OpenDIA's apex agrees with DIA-NN on shared IDs to a p50 of 0.0-0.95 s and a p95 of 1.0-2.9 s, so the target is a valid measurement.
- **IM target:** DIA-NN's reported 1/K0 minus the calibrated library 1/K0. This is a caveat: OpenDIA does not write an observed 1/K0 per ID, so the measurement comes from DIA-NN, on IDs OpenDIA itself selected.
- **Features (68):** residue counts (oxidised Met and carbamidomethyl Cys as their own tokens), N- and C-terminal residue one-hots, length, log length, charge one-hot, and precursor m/z. Ridge with CV-chosen alpha on standardised features.
- **Evaluation:** DIA-NN non-entrapment IDs of the same run, all and missed (not in OpenDIA's q ≤ 0.01 set).
- **Script:** `precheck_rtim.py` (sha256 `30b56451…`). The gate was written in the header before the run.
- **Gate:** on H, the corrected p95 on missed IDs must be at least 15 % below the current map's. On D it must be no more than 5 % worse. RT and IM are judged separately.

## Result: RT passes

| p95 RT error, missed DIA-NN IDs | D | H |
|---|---:|---:|
| Current map | 22.6 s | 102.2 s |
| Best monotone map fitted on DIA-NN's own IDs (the earlier "floor") | 20.5 s | 86.5 s |
| **Map + ridge correction** | **16.4 s (−27.4 %)** | **63.1 s (−38.3 %)** |
| Window covering 99 % of all DIA-NN IDs: map, then corrected | 32.8 → 24.2 s | 134.2 → 89.6 s |

Adding the mapped value as a feature ("hyb") gives the same result. The correction gets below the monotone floor because it fixes each peptide's own prediction error, which no single map can.

**Training size** (top-N OpenDIA IDs by score; anchors are 140 on D and 441 on H):

| Training IDs | 140 | 441 | 1,000 | 3,000 | 10,000 | all |
|---|---:|---:|---:|---:|---:|---:|
| D, missed p95 change | −12.8 % | −19.6 % | −22.5 % | −25.7 % | −27.2 % | −27.4 % |
| H, missed p95 change | −18.3 % | −29.8 % | −35.0 % | −37.1 % | −37.5 % | −38.3 % |

Calibration anchors alone give part of the gain, and about 3,000 IDs give nearly all of it. So the fit belongs after a first search, not in the calibration step.

## Result: IM fails

| p95 1/K0 error, missed DIA-NN IDs | D | H |
|---|---:|---:|
| Current calibrated line | 0.0353 | 0.0357 |
| Best correction ("hyb") | 0.0326 (−7.6 %) | 0.0339 (−4.9 %) |

Sequence composition explains little of the mobility residual, so the 1/K0 correction is not taken.

## Decision

The RT correction proceeds. Proposed form, opt-in at first: calibrate, search, fit the ridge on the run's q ≤ 0.01 targets, then re-search with the corrected RT and a window from the corrected residuals. Ridge on 68 features is a closed-form solve, so it needs no libtorch and no Python. It goes through the usual gates: C++ opt-in with a default byte check, a two-seed r ≈ 1.56 workflow check, then a frozen r ≈ 7 plan. These residuals measure RT accuracy, not IDs. Whether the narrower window and better RT centre turn into IDs at equal true FDP is what the workflow check tests.
