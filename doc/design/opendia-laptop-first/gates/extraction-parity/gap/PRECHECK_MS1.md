# PRECHECK_MS1: is MS1 precursor evidence specific to real peptides? (16 November 2026)

Diagnostic, written before any plan.
- **Inputs:** arm-E population and arrays. `ms1.npy` holds M+0..M+3 over apex ± 15 s in 32 bins, and `xic.npy` the 6
  search fragments; both come from the decoygen arm-E chain, whose population is md5-identical to the gap arm E.
  The C++ scorer at seeds 0-4 is used.
- **Script:** `ms1_precheck.py`. Features are label-free and identical for targets and decoys:
  - M+0 / fragment-sum trace correlation;
  - observed-against-theoretical isotope envelope cosine;
  - M+1 / M+2 correlation with M+0;
  - M+0 apex offset from the fragment apex;
  - log M+0 area.

## 1. Specificity at matched score (seeds 0 and 2; raw01 / raw02)
- **R:** DIA-NN-supported right-peak IDs rejected at today's cut.
- **D:** score-matched decoy winners.
- **E:** entrapment winners in the same band, with their own matched decoys ED.

| feature | AUC R > D | AUC E > ED |
|---|---|---|
| m1_corr_xic | 0.65-0.67 / 0.65 | 0.53-0.54 / 0.57 |
| m1_iso_cos | 0.68-0.69 / 0.66 | 0.53-0.54 / 0.51-0.53 |
| m1_apex_dt (lower = R) | 0.34-0.35 / 0.35 | 0.46-0.50 / 0.46 |
| m1_log_m0 | 0.62-0.63 / 0.61-0.62 | 0.53-0.56 / 0.51-0.54 |
| m1_iso_corr | 0.57-0.58 / 0.54-0.55 | 0.53-0.54 / 0.50-0.54 |

- The features separate rejected real peptides from score-matched decoys at AUC 0.61-0.69.
- Entrapments stay near 0.5: 0.50-0.57, the highest being m1_corr_xic on raw02.
- So MS1 evidence is largely specific to real peptides. Entrapments share their parent's m/z and envelope, and that
  is the residual 0.53-0.57.

## 2. Null-window tell check (windows > 3σ from the RT map, where no peptide is expected)
AUC(target > decoy) is 0.498-0.501 for every feature on both raws, so there is no library or sequence tell.

## 3. Offline scorer probe (`MS1_PROBE.txt`; C++ CrossFitScorer, current 24 inputs + 5 MS1 inputs, seeds 0-4)

| | raw01 CUR | raw01 +MS1 | raw02 CUR | raw02 +MS1 |
|---|---|---|---|---|
| N@FDP_c 1 %, median (range) | 1822 (1804-1864) | **1850** (1840-1860) | 1905 (1826-2021) | **1906** (1868-1955) |
| N@FDP_c 0.5 %, median | 1616 | 1734 | 1666 | 1678 |
| IDs at q ≤ 0.01, median | 1691 | **1738** | 1621 | **1711** |
| FDP_c at q ≤ 0.01, max | 0.0082 | 0.0080 | 0.0062 | 0.0080 |

- **Ranking at equal true FDR improves on raw01 only:** +28 (+1.5 %) against +1 on raw02.
- **The reported IDs at q ≤ 0.01 rise on both raws:** +47 / +90 (+2.8 % / +5.6 %). FDP_c stays ≤ 0.008, so the
  q-values become less conservative while true error remains under 1 %.
- **The seed spread narrows** on raw01 (60 → 20) and raw02 (195 → 87).

## Why the earlier MS1 arm (D11M, v4) failed
- D11M ran on the old library decoys, where every arm had entrapment excess. D11 itself reached FDP_c 0.015 there.
- Its failure (0.016) is therefore not specific to MS1.
- With predicted decoys (arm E), MS1 features show no specificity problem here.

## Reading
- MS1 is the first new evidence block since block I that is specific and label-free.
- Its main effect is on reported IDs at q ≤ 0.01 and on seed stability. Its effect on equal-FDR ranking is small.
- **Worth a frozen plan.** The block must be implemented in C++ in the workflow: the evidence pass already collects the
  M+0 trace and an M+1 core sum, and it needs M+1..M+3 traces and the theoretical envelope added.
