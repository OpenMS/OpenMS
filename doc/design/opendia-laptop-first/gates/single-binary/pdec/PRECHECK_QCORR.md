# PRECHECK_QCORR: shared-signal q correction does not move the r ≈ 7 cut

Written 20 Nov 2026. This is a diagnostic only. Labels are used only to evaluate the candidate correction. No C++ change or plan was made.

## Candidate 1: add shared-target excess to the TDC numerator

For winners with the tool's q ≤ 0.05, define a shared-signal event as at least three library fragments within 20 ppm of a stronger co-eluting winner within ±10 s and ±12.5 m/z. At each score rank, estimate the null shared-event rate from decoy winners and add the positive excess among target winners to the TDC numerator:

`q' = (D + 1 + gamma * max(0, S_target - p_decoy * T)) / T`

The correction is run-derived and uses no entrapment labels. The diagnostic swept `gamma = 0.25, 0.5, 1, 2, 5, 10, 20` on one D and one held-out H seed using existing feature tables.

| run | baseline N / lower / combined | all gamma values |
|---|---:|---:|
| D s0 | 8,483 / 0.0111 / 0.0127 | 8,483 / 0.0111 / 0.0127 |
| H s0 | 21,386 / 0.0107 / 0.0122 | 21,386 / 0.0107 / 0.0122 |

The correction did not move the q ≤ 0.01 cut, even with `gamma = 20`. It is rejected as ineffective.

## Candidate 2: subtract a shared-signal score penalty

As a second diagnostic, re-rank winners with

`score' = score - lambda * max(0, shared_fragments - 2)`

and recompute plain TDC. This is still label-free, but it is a scoring change rather than a q-estimator change.

| run | baseline lower / combined | lambda sweep |
|---|---:|---:|
| D s0 | 0.0111 / 0.0127 | λ 0.05-0.80: 0.0111 / 0.0127, N -0.01 to -0.19 % |
| H s0 | 0.0107 / 0.0122 | λ 0.05-0.80: 0.0107 / 0.0122, N -0.00 to +0.21 % |

The penalty changes a few ranks but does not improve calibration. It is also rejected.

## Decision

Do not implement either candidate. A simple label-free correction of the shared-signal residual is not enough to make q = 0.01 nominal on this data. The residual should remain documented as a scoring/interference problem, not a TDC formula bug. The default OpenDIA binary and `-trace_tolerance_ppm 10` are unchanged.

The next viable calibration work would require a stronger model of the claimant graph or acceptance of the measured approximately 1.1× residual. That requires a new design, not a blind C++ patch.
