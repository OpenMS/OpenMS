# Q6: do the peak-group scores reward "a peak is here" or "it is this peptide"? (diagnostic, D seed 0)

Written 18 November 2026. Labels (DIA-NN, decoy status) are used here for evaluation only.

**Groups** (precursor winners of NEW, `-pg_scores`; `Q6_PEAK_IDENTITY.txt`):
- **T, true:** CUR-accepted targets confirmed by DIA-NN 1 % within ±3 s. n = 14,799.
- **Dpk, peak present but wrong identity:** decoy winners whose apex is within 3 s of their paired target's DIA-NN peak.
  n = 10,200, 15 % of decoy winners.
- **Doff, no real peak:** decoy winners more than 30 s from it. n = 39,329.
- **X:** NEW's 694 additions that are not in DIA-NN at 1 %.

**Measures:**
- "identity" = AUC(T > Dpk).
- "peak" = AUC(Dpk > Doff). This would be far from 0.5 if a score rewarded a real peak whatever the identity.

## Result
- **No score is peak-only.**
  - For all 10 pg scores, the "peak" AUC is 0.48-0.53. Decoys sitting on a real peak look no better than decoys in empty
    regions.
  - Their identity AUCs are strong: 0.85-0.90 for shape, S/N, isotope and intensity fraction, and 0.05-0.18 (inverted)
    for lag and mass deviation.
  - The only feature with a clear peak-only effect is the old `abs_library_rt_error` (peak AUC 0.27). Its cause is the
    shared RT, as expected.
- **So the hypothesis from PLAUS_D is wrong.** The new scores do not let a decoy (or a false target) score well just
  because it sits on another peptide's peak.

## The PLAUS_D estimate was too pessimistic
Two independent cross-checks of the 694 (`Q6_PEP.txt`):

1. **DIA-NN's own posterior.**
   - 64 % of the 694 are seen by DIA-NN at the same peak, with a mean (1 − PEP) of 0.43.
   - Σ(1 − PEP) = 189 expected true, counting unseen IDs as false.
   - CUR's own 1,111 non-DIA-NN IDs give the same 193.
   - That is a lower bound: unseen IDs are not all false.
2. **Our decoys in the same score range.**
   - Between the lowest and highest X score, there are 269 decoy winners.
   - Under target-decoy competition, that estimates about 269 false targets in that range, CUR's included.
   - This is not independent: it is the decoy model this whole question is about. Astral entrapment supports it on raws
     01/03 and shows it under-counting on raw 02.
   - So **if target-decoy competition holds, at most about 270 of the 694 are false**.

- PLAUS_D's decoy-null estimate (~500 false) assumed that a false target gets outside support as often as a well-scored
  decoy. The decoys sit on real peaks only 29 % of the time in this band, so that null was not comparable.
- **Revised estimate:** of NEW's 3,053 additions on D, about 2,360 are DIA-NN-confirmed.
  - Of the other 694, between about 190 (DIA-NN's PEP, a lower bound) and about 420 (the target-decoy estimate) are
    real.
  - So **about 270-500 are false**.
- The extra false share is 1.4-2.7 % of NEW's 18,639. NEW's q ≤ 0.01 on D is probably over its nominal rate, by an
  amount we cannot pin down without entrapment.

## What this leaves
- The new scores are identity evidence, not peak artefacts.
- The agreement drop on D (0.930 → 0.910) is a mix of true IDs that DIA-NN does not report at 1 % and a false excess of
  about 270-500.
- A frozen criterion cannot be reinterpreted afterwards, so the verdict stands.
- The Astral raw 02 excess (FDP_c 0.0100 pooled) is the one clean, entrapment-measured problem. Raws 01 and 03 were
  better calibrated with NEW.
- **Recommended next frozen plan.** Keep the 10 scores and add a calibration safeguard that comes only from the run,
  with entrapment FDP_c on raws 01-03 as the acceptance test. For example:
  - a more conservative scorer as the input dimension grows;
  - or a decoy-count-based q floor.
- Replace the DIA-NN-agreement proxy on diaPASEF with a pre-registered estimate that has both a lower and an upper
  bound, such as the PEP and target-decoy pair used here. That criterion has to be fixed before any run.
