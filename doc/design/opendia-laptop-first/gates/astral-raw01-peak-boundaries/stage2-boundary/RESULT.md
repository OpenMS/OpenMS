# Stage-2 boundary group: raw 01 design result (descriptive, not a gate)

`PLAN.md` was written before the run (`PLAN.written_at`, sha256 in `PLAN.md.sha256`). Binary is commit `5e23606c15` (`binary.sha256`). Labels are target/decoy only; entrapment labels were not joined.

## Equivalence
- The CNN path is unchanged: sample, logit, fold and window columns are byte-identical to `../stage2-search-design/search` (`CMP.txt`).
- Features columns 1 to 9 (ids, cnn_z, search) are identical.
- `precursor_best` and the stage-2 columns differ by design.

## Counts (pooled paired q, decoy-count estimate)

| q ≤ | 0.005 | 0.01 | 0.02 | 0.05 |
|---|---|---|---|---|
| search (frozen) | 331 | 717 | 1015 | 1412 |
| search+boundary | 360 | **746** | 1000 | 1441 |

- **At q ≤ 0.01:** +29 (+4.0%), with 700 shared, 46 gained and 17 lost.
- **At q ≤ 0.02:** −15.

## Weights (standardised, all three folds agree in sign)

| feature | weight range |
|---|---|
| boundary co-elution | +1.87 to +1.94 |
| boundary library correlation | +0.82 to +1.21 |
| boundary log FWHM | −1.56 to −1.62 |
| apex agreement | −0.75 to −1.06 |
| log base width | ≈ 0 |

The negative apex-agreement weight is conditional on the other features: it is collinear with co-elution and does not mean agreement is bad.

## Null-band symmetry (medians, target / decoy precursor winners)

| CNN z band | co-elution | library correlation | log FWHM |
|---|---|---|---|
| ≤ 0 (152K / 157K) | 0.066 / 0.053 | 0 / 0 | 1.61 / 1.60 |
| 4 to 6 (7.9K / 6.4K) | 0.179 / 0.166 | 0.223 / 0.150 | 1.86 / 1.87 |
| > 6 (2.6K / 1.1K) | 0.477 / 0.267 | 0.641 / 0.242 | 1.77 / 1.84 |

- The features barely differ in the mostly-null bands and separate where real signal exists.
- The z 4 to 6 band holds a mix of true and false targets, so a small difference there is expected.

## Gained versus lost (`gain.py`)
- **Gained (46):** apex log10 4.82, about as bright as the shared IDs (4.90).
- **Lost (17):** weaker, at 4.09.
- So the gain is not concentrated on weak peaks.

## Decision under the predeclared reading
- **The pre-set threshold (≥ 753, +5%) is not met:** 746 is +4.0%.
- Under the plan, the rule-boundary approach therefore stops here and **no gate is proposed** for it on an unused run.
- The learned segmenter is not ported either, because the gain is not at low SNR.
- The group stays implemented and opt-in, since it is cheap (no measurable time cost: 12:52 wall, 24.7 GB, against about 14:19 for the frozen run).

## Reading
- **Small effect:** a few extra features over correct boundaries help a little but do not close the gap to OpenSWATH+Percolator (1441).
- **The boundaries themselves are not the limit:** the stage-2 width and integration still come from the search's wide window, and the CNN and most features remain unchanged.
- **The larger levers per `MISSES.md`:**
  - the CNN's bias towards bright interference;
  - a linear discriminant that cannot combine the features non-linearly.
