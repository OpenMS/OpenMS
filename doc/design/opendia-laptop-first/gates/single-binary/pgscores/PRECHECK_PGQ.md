# PRECHECK_PGQ: which peak-group score breaks q calibration? (diagnostic, after RESULT_PGSCORES)

Written 18 November 2026. Entrapment and DIA-NN labels are used here for evaluation only.

- **Inputs, Astral:** NEW seed 0 of raws 01 and 02 rerun with `-features_out` (same container). The outputs are
  byte-identical to the design seed 0. We use the precursor winners.
- **Inputs, HeLa D:** NEW seed 0 features.
- **Band:** winners with run q from 0.002 to 0.05 (about 65 entrapment and 118 decoy winners per raw).

## 1. Raw 02's FDP excess is real

**Correction (Q7_NEXT.md):**
- The intervals below pool seeds of the same raw as if they were independent. They are not: the seeds share the run's
  data, so the run is the unit.
- With one run, NEW's raw 02 count of 10 entrapments is what an exactly calibrated 1 % scorer gives half of the time.
- The excess is therefore **not** shown to be real. The heading is kept for the record.
Entrapments at q ≤ 0.01 pooled over 5 seeds, with Poisson 90 % intervals (`Q4_FDP_POOLED.txt`). Entrapments per seed:
raw 01 CUR 7 7 4 8 4, NEW 4 4 7 5 5; raw 02 CUR 8 4 3 4 1, NEW 12 10 7 7 13; raw 03 CUR 9 5 6 8 7, NEW 5 7 5 6 7.

| raw | CUR FDP_c pooled (90 %) | NEW FDP_c pooled (90 %) |
|---|---|---|
| 01 | 0.0070 (0.0050-0.0094) | 0.0054 (0.0037-0.0075) |
| 02 | 0.0045 (0.0030-0.0066) | **0.0100 (0.0078-0.0127)** |
| 03 | 0.0077 (0.0057-0.0102) | 0.0062 (0.0045-0.0084) |
| all three | 0.0064 (0.0053-0.0077) | 0.0072 (0.0061-0.0085) |

- On raw 02, NEW has 49 entrapments across seeds against CUR's 20, and the intervals do not overlap.
- On raws 01 and 03, NEW is better calibrated than CUR.
- Pooled over the three raws, both arms are below 1 %.
- So the excess is specific to one raw. It is not a general bias.

## 2. No single score explains it
- **Per score:** AUC(entrapment > decoy) in the band (`Q4_FEATURE_ED.txt`). A score keeps decoys a valid model of false
  targets only if this is near 0.5.
  - Every pg score is within 0.45-0.55 on both raws.
  - The largest are `pg_massdev` on raw 01 (0.552) and `pg_isotope_overlap` on raw 02 (0.448).
  - Several old features are further from 0.5: `d6_max` 0.594 on raw 01, `boundary_log_fwhm` 0.411 on raw 02.
  - With about 65 entrapments, the noise sd of such an AUC is about 0.045. Nothing stands out.
- **Final score:** AUC(entrapment > decoy) in the band (`Q4_SCORE_ED.txt`), median over 5 seeds:
  - raw 01: CUR 0.414, NEW 0.414;
  - raw 02: CUR 0.492, NEW 0.521;
  - raw 03: CUR 0.452, NEW 0.462.
  - On raw 02, NEW moves entrapments slightly above decoys, but the size is within seed spread (0.502-0.587 against
    0.434-0.519).

## 3. D: the added IDs missing from DIA-NN
`Q4_D_ADDED.txt`, seed 0. NEW adds 3,053 precursors over CUR: 2,359 are in DIA-NN and 694 are not.
- No pg score separates the two groups well. The best AUCs are `pg_massdev` 0.385 (lower in DIA-NN) and `pg_log_sn`
  0.594.
- `transition_support` (0.606) and `f12_libcos` (0.591) do as well as any pg score.
- Both groups sit between the top decoys and the old accepted IDs on almost every feature.
- So they look like a weaker tail, with no feature-specific artefact.

## Reading
- The failure is not one bad score that could be dropped. Two effects add up:
  1. **Raw 02 q calibration.** The excess is real on that raw and absent on the other two. No feature separates
     entrapments from decoys beyond noise, so the cause is more likely the scorer's fit on that raw (a larger model,
     34 inputs) than a biased input.
  2. **D agreement.** NEW goes deeper into the low-evidence tail, and DIA-NN agrees there less often. Part of that may
     be real IDs, which this proxy cannot separate.
- **Not supported:** a follow-up plan that drops specific pg scores. The data give no candidate to drop, and choosing a
  subset now would be fitting to these labels.
- **Candidates for a new frozen plan, each needing its own precheck:**
  - Make the scorer more conservative as the input dimension grows (for example, regularisation or the size of the
    hidden layer scaled with dimension), and check q calibration on raws 01/02.
  - Keep the 10 scores with the current scorer and test the calibration claim on more raws. Raws 01-03 pooled are at
    0.0072, so the raw 02 seeds may be one raw's tail and not a systematic bias. This needs new Astral data, which we
    do not have.
- The block stays removed. `-pg_scores` is reverted (b35655266a).
