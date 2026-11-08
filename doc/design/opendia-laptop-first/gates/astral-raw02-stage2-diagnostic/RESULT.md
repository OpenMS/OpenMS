# Raw 02 stage-2 re-score: diagnostic, not a gate

**Why this is not a gate.** Raw 02 was spent by the Option B gate, which FAILED (`../astral-entrapment-raw02-20261108/RESULT.md`). The stage-2 features were chosen after the raw 02 post-mortem looked at its entrapment hits. This re-score is therefore optimistic by construction and does not count as evidence of error control. It only checks that the frozen raw 01 design (code `014cb0dbf1`, `astral-design-raw01-20261107/stage2-design/RESULT.md`) moves the spent run in the expected direction.

## Setup
- **Inputs reused read-only from the raw 02 evaluation:** calibration, seeds, windows and all CNN flags. The only addition is `-stage2 -stage2_features_out` at the frozen defaults.
- **Unchanged columns:** the CNN logit, fold and window columns are byte-identical to the raw 02 evaluation (`CMP.txt`).
- **Labels:** entrapment labels come from the unchanged `evaluate.py`.

## Numbers (pooled paired q; FDP_c is the combined entrapment estimate)

| | accepted q ≤ 0.01 | N_e | FDP_c | 95% CP for FDP_c |
|---|---|---|---|---|
| CNN only (raw 02 gate) | 360 | 5 | 0.028 | 0.009–0.065 |
| stage 2 (this re-score) | 518 | 3 | 0.0117 | 0.002–0.034 |

- **Per fold at q ≤ 0.01:** 264 / 97 / 157, with entrapment hits 2 / 0 / 1.
- **Overlap:** 345 precursors are accepted by both. 15 are accepted only by CNN, 2 of them entrapment. The other 3 entrapment precursors the CNN accepted are still accepted.
- **Other q cuts:**

  | q ≤ | accepted | N_e | FDP_c |
  |---|---|---|---|
  | 0.005 | 372 | 2 | 0.011 |
  | 0.02 | 607 | 3 | 0.010 |
  | 0.03 | 776 | 12 | 0.031 |
  | 0.05 | 967 | 26 | 0.054 |

- **Fold weights:** these were less stable than on raw 01, with `library_correlation` between −0.55 and +0.13. The `cnn_z`, co-elution and MS1 terms were consistent.

## Reading
- **Direction:** count and FDP_c both move the right way. The 95% interval still spans the 0.015 limit, and it rests on 3 entrapment hits.
- **Above q 0.02:** FDP_c rises above the nominal q (0.031 at 0.03, 0.054 at 0.05). The decoy-count estimate is still anti-conservative there.
- **Independent check:** only raw 03, under a hashed gate, can give one.
