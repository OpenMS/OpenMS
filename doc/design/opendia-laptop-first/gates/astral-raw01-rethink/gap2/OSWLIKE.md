# A teacher-free OSW-like scorer on the shared XICs (raw 01, replicated on raw 02)

This is our own implementation, not OSW code. It runs on the same OSW-extracted traces OSW+Percolator scored (`../gap2/build_xic.py`). The scripts are `oswlike_feats.py`, `oswlike_score.py`, `oswlike_followup.py` and `oswlike_single.py`.

- **Labels:** target/decoy only.
- **Entrapment:** used only for FDP_c ("combined", r = 0.9851).
- **DIA-NN:** used only for the "confirmed right" column.

**Metric:** N target precursors at entrapment FDP_c ≤ 0.5 / 1 / 1.5 / 2%. Precursors are ranked by score, after rank-1 peak-group selection per precursor and paired TDC. This is the same metric as `../matched_fdp.txt`.

## Setup

- **Candidates:** per precursor, local maxima of the Gaussian-smoothed summed fragment trace over the whole XIC (±97 s), top 5 by height. Boundaries walk to a valley or to 5% of the apex, and overlapping candidates are dropped. This gives 2.70M candidates. The right peak is present for 2528 of 2754 DIA-NN-confirmed precursors (92%).
- **20 subscores**, OSW-style but our own code:
  - library-weighted and unweighted lagged-xcorr co-elution and shape;
  - minimum pairwise xcorr;
  - library correlation, dot product, Manhattan, RMSD and spectral angle;
  - intensity share and log S/N;
  - fragments above 3× noise;
  - minimum fragment share;
  - log width and |apex − calibrated library RT|;
  - MS1 xcorr co-elution and shape, MS1 isotope correlation, MS1 log S/N.
- **Learner:**
  - PyProphet-style semi-supervised LDA, cross-fitted on pair-hash folds;
  - initialised on library-weighted xcorr shape;
  - positives are top-ranked targets at q ≤ 0.15, then 0.02 for 4 more rounds;
  - negatives are top-ranked decoys;
  - per-fold decoy-winner median/MAD scaling.
- **MLP variant:** the same loop, then 2 rounds of a 32-16 MLP.

## Results, raw 01 (N at FDP_c ≤ 0.5 / 1 / 1.5 / 2%)

| scorer | N | q ≤ 0.01 (FDP_c) |
|---|---|---|
| **OSW+Percolator (reference)** | **1015 / 1151 / 1251 / 1286** | – |
| OpenDIA tool (search+boundary stage 2) | 641 / 757 / 923 / 950 | 746 (0.008) |
| OSW-like LDA, 20 subscores, XIC maxima | 850 / 1011 / 1096 / 1124 | 993 (0.010) |
| OSW-like MLP, seed 0 | **1116 / 1221 / 1256 / 1326** | 1222 (0.012) |
| OSW-like MLP, seed 1 | 1061 / 1214 / 1271 / 1292 | 1220 (0.012) |
| OSW-like LDA on **our search's** candidate apexes | 962 / 1023 / 1028 / 1080 | 1023 (0.010) |
| OSW-like MLP on our search's candidate apexes | 920 / 1110 / 1196 / 1242 | 1092 (0.009) |

Ablations (LDA, one group dropped at a time), N at FDP_c 1%:

| dropped | N |
|---|---|
| none | 1011 |
| library | 606 |
| signal | 977 |
| shape/RT | 974 |
| MS1 | 1015 |
| xcorr | 7, a bootstrap failure (see below) |

## Reading

1. **The gap closes on raw 01 with a teacher-free scorer on the same traces** (raw 02: mostly; see the replication below). A small MLP on 20 OSW-style subscores, with decoy-only labels, gives 1214–1221 at a true FDP_c of 1%. That is above OSW+Percolator (1151), with decoy-estimated q ≤ 0.01 calibrated to FDP_c 0.012. The LDA alone gives 1011. So most of the 757 → 1151 gap is not missing information and not extraction. It is **which subscores we compute plus a non-linear combiner.**
2. **Our candidate generator is not the main problem either.** The same subscores on our streaming search's apexes give 1023 (LDA) and 1110 (MLP) at FDP_c 1%, against 1011 and about 1215 on the XIC maxima. Coverage of the right peak is 2424 of 2754 (88%) against 2528 (92%). This costs at most about 9% with the MLP and nothing with the LDA.
3. **The subscores that matter, compared with what OpenDIA uses today:**
   - **Library agreement is the largest single group.** Dropping it loses 40%. We have library correlation, but not dot product, Manhattan, RMSD or spectral angle on the integrated fragment areas.
   - **Lagged xcorr is the backbone.** Without it the semi-supervised loop finds no positives at the first step (0 at q ≤ 0.15, even when initialised on the library dot product). No single non-xcorr subscore reaches 15% FDR alone. Every single subscore alone is weak: the best is log S/N at 510, then xcorr co-elution at 304. The strength comes from the combination.
   - **MS1 adds nothing here,** and signal and shape/RT add about 3–4% each.
4. **Peak boundaries and baseline subtraction are not the lever** (see `GAP2.md`). The valley/5% boundaries used here give a median width of about 11 s, close to OSW's 13 s and DIA-NN's 9.8 s. Our search groups' 46 s spans should still not feed features.

## Raw 02 replication (`raw02_score.txt`; same code and settings, XIC cap raised to 144 points)

The first raw 02 pass truncated the traces at 128 of about 134 points. That pass is kept in `trunc128/`; its results are similar.

The right peak is among the candidates for 2504 of 2745 confirmed precursors (91%).

| scorer, raw 02 | N at FDP_c ≤ 0.5 / 1 / 1.5 / 2% | q ≤ 0.01 (FDP_c) |
|---|---|---|
| **OSW+Percolator (reference)** | **1159 / 1287 / 1351 / 1389** | – |
| OpenDIA tool, stage-2 re-score (`astral-raw02-stage2-check`; optimistic, see its RESULT.md) | 169 / 298 / 454 / 494 | 518 (0.012) |
| OSW-like LDA, 20 subscores | 613 / 885 / 954 / 1071 | 1115 (**0.023**) |
| OSW-like MLP, seed 0 | 907 / 1040 / 1093 / 1232 | 1231 (**0.020**) |
| OSW-like MLP, seed 1 | 910 / 1119 / 1158 / 1178 | 1221 (**0.025**) |

Ablations reproduce on raw 02. Measured as N at FDP_c 1%, from 885 with all groups:

| dropped | N |
|---|---|
| library | 481 |
| signal | 841 |
| shape/RT | 852 |
| MS1 | 868 |
| xcorr | 6 (bootstrap failure) |

**Reading of the replication:**

1. **The direction holds.** The OSW-like scorer is far above the current tool on both runs (raw 02: 1040–1119 against 298 at FDP_c 1%). Library agreement and lagged xcorr are again the carrying families.
2. **It does not beat OSW+Percolator on raw 02.** The MLP reaches 81–87% of OSW's 1287, and the LDA 69%. The raw 01 "above OSW" result does not replicate. The honest summary across both runs is **about 85–105% of OSW+Percolator at matched true error.**
3. **Decoy-estimated q is anti-conservative on raw 02** (FDP_c 0.020–0.025 at q ≤ 0.01, against 0.010–0.012 on raw 01). For reference, OSW+Percolator's own m_score ≤ 0.01 export runs at FDP_c 0.032 on raw 01 and 0.034 on raw 02 (`matched_fdp02.txt`). Anti-conservative decoy q is therefore common to both pipelines on this data, and is not a raw 02 peculiarity of OSW. It still fails our gate rule (1.5×). Diagnosis is in `CALIB.md`.
4. Seed variance at 1% is about 8% (1040 against 1119).

## Caveats

- Two runs. Raw 01 is heavily reused and raw 02 is spent, so both are design runs. The OSW reference on raw 02 is 1159 / 1287 / 1351 / 1389.
- XIC extraction is OSW's (about 34 ppm, ±97 s, 1.6 s cycle). Our native extraction must reproduce it. Bright-peak agreement is shown; weak-peak agreement is not.
- The MLP's FDP_c at q ≤ 0.01 is 0.012 on raw 01, within the 1.5× gate rule, but 0.020–0.025 on raw 02, which fails it. Only about 8–15 entrapment hits lie behind each number.
- Two seeds agree at 1%; the 0.5% and 5% points vary more.

## Implication for the design (the user's single-model direction)

The richer learned model should **take these subscore families as inputs alongside the trace encoder**, not replace them:

- lagged xcorr co-elution and shape, unweighted and library-weighted;
- library agreement on integrated areas (correlation, dot product, Manhattan, RMSD, spectral angle);
- intensity share and S/N.

It should be computed per candidate on valley-bounded peaks from a per-precursor XIC.

Section 5 of `GAP2.md` showed that a small trace CNN alone does not recover them. A single cross-fitted MLP head over [CNN embedding + these about 20 subscores], with decoy-only semi-supervised labels, is the minimum design this evidence supports. It replaces stage 2 and its opt-in feature sets. **It is not sufficient on its own:** decoy-estimated q was anti-conservative on raw 02 for both LDA and MLP, so the error-control path needs its own fix (`CALIB.md`).

## Update (2026-11-09): head choice and error control, on the re-extracted XICs (`../decoyregen/`)

The decoy-family test re-extracted both runs with library RT calibrated per run (`CALIB.md` section 6). On those XICs, family O (the library decoys), entrapment-matched N at FDP_c 0.5 / 1 / 1.5 / 2% is:

| Head | raw 01 | raw 02 |
|---|---|---|
| MLP, 5-seed ensemble | 1035 / 1101 / 1136 / 1296 | 1008 / 1083 / 1132 / 1223 |
| Boosted trees as final head, 3 seeds (`head_hgb.py`) | 287 / 419 / 432 / 436 | 257 / 297 / 329 / 405 |
| Rank mean of both | 817 / 908 / 1126 / 1160 | 559 / 591 / 683 / 1027 |

**Boosted trees as the ranking head fail** the predeclared rule (≥ 3% better on both runs) by a wide margin. They overfit the semi-supervised positives, which are noisy and few (about 1000 per fold). The MLP stays the ranking head.

Trees are still the right model for the **bulk-null calibration weights** (`CALIB.md` section 9), where they are trained on about 180k bulk winners per class. With the MLP ensemble plus bulk-null weighted TDC, the q ≤ 0.01 FDP_c is 0.006 (raw 01) and 0.004 (raw 02), and the q ≤ 0.05 FDP_c is 0.050 and 0.055. Plain TDC gave 0.016/0.013 and 0.067/0.082.

The raw 02 count gap to OSW+Percolator (1083 against 1287 at FDP_c 1%) is unchanged by this. It sits in ranking, not calibration.

### Where the remaining count gap sits (`../decoyregen/gap_osw.py`, both runs, re-extracted XICs)

This covers the human targets that OSW+Percolator exports at m_score ≤ 0.01 but our 5-seed MLP does not accept at entrapment FDP_c 1%. OSW's own export is about 3% anti-conservative (section 1 of `CALIB.md`), so part of this set is false.

| Category | raw 02 | raw 01 | DIA-NN also IDs |
|---|---|---|---|
| Accepted by both | 954 | 961 | 0.98 |
| Our winner is OSW's peak, scored below our cut | 398 | 354 | 0.92 |
| OSW's peak is a candidate, but we pick another peak | 58 | 64 | 0.72–0.90 |
| No candidate within 10 s of OSW's RT | 54 | 58 | 0.81–0.86 |
| Not extracted | 2 | 2 | 0 |

**Reading.**
- About 75% of the misses (74–78%) have OSW's peak as our rank-1 candidate. They sit a median 3.5 decoy-MAD below our cut, and only 11–15% are within 1 MAD of it.
- Most of these are real: 92% are DIA-NN IDs.
- So the gap is in **score formation on the right peak**, not in candidate generation (about 11%) or peak choice (about 12%).
- The 20 subscores, or how the model combines them, do not separate these peaks from decoys as well as OSW's scores do.

**Next.** A diagnostic re-extraction of raw 02 that keeps OSW's feature table (`extract_feat02.sh`). It will compare OSW's subscores with ours on these 398 peaks, to locate the missing information. OSW's scores are used only to diagnose, never in our model.

### Which information is missing: b/y ion series in the apex spectrum (raw 02; raw 01 running)

**Setup.** This uses a diagnostic re-extraction of raw 02 that keeps OSW's feature table (`extract_feat02.sh`). Each of our candidates is matched to OSW's nearest peak group of the same precursor, within 5 s; 47% match. OSW's subscores are used here only to find which *information* our features lack. Nothing from OSW's scoring code goes into our model.

**Conditional separation (`osw_info.py`).** Within bins of our score, over the top 5% of decoy winners and above, no single OSW subscore separates targets from decoys by more than AUC 0.54. Library-match variants, MI, mass deviation and b/y series lead.

**Augmentation**, 5-seed MLP ensemble on raw 02. Entrapment-matched N at FDP_c 0.5 / 1 / 1.5 / 2%, against our 20 subscores at 1008 / 1083 / 1132 / 1223:

| Added | N at FDP_c 0.5 / 1 / 1.5 / 2% |
|---|---|
| mass deviation (3 columns) | 1048 / 1119 / 1152 / 1324 |
| MI (4) | 731 / 944 / 1056 / 1125 |
| isotope (4) | 974 / 1119 / 1135 / 1176 |
| **b/y series (2)** | **1053 / 1281 / 1316 / 1377** |
| all 30 | 971 / 1295 / 1329 / 1416 |

**Confound check (`series_check.py`, reading declared before running).** The first run also added a flag for "OSW reported a peak group here". Results:
- The flag alone does not help: 1068.
- b/y series without the flag: **1291** (+19%). That matches OSW+Percolator's 1287 at FDP_c 1% on this run.
- Restricted to the candidates that match an OSW feature: 1028 → 1165 (+13%).

Both arms meet the ≥ 5% reading **on raw 02**. With b/y series plus the flag, plain TDC gives q ≤ 0.01 FDP_c 0.007 and q ≤ 0.05 FDP_c 0.076. Bulk-null TDC gives 0.006 and **0.040**, so calibration holds with the new features.

**Interpretation.** OSW's b/y series scores count how many theoretical b and y ions of the precursor are present in the apex MS2 spectrum. This uses the full theoretical ion ladder, not just the 6 library fragments that our XIC features see. That evidence is **outside the top-6 XICs**, which explains why no XIC-derived feature, CNN or ranker change recovered it.

**For the native pipeline:**
- Own implementation: at the candidate apex (and ±1 scan), match all theoretical b/y ions (charge 1 and 2, inside the scan range) against the centroided MS2 spectrum within the extraction tolerance.
- Features: matched count per series, longest consecutive run, matched intensity share, and decoys scored the same way.
- This is the next native change to test, once raw 01 confirms.

**Raw 01 replication (`series01.txt`).** OSW's b/y series scores also add information on raw 01:
- Without the flag: 1101 → **1242** (+13%).
- Restricted to OSW-matched candidates: 1130 → 1241 (+10%).

The reading passes on both runs. With series plus flag, bulk-null TDC q ≤ 0.05 FDP_c is 0.063 on raw 01 (0.040 on raw 02). Raw 01 thus sits just above the 0.06 bar with these features.

**Our own b/y features, first version: FAIL (`by_feats.py`, `by_test.py`).** The features are fractions of theoretical b and y ions (charge 1 and 2) matched within 17 ppm, longest covered run, and matched-intensity share, averaged over the apex spectrum and its two neighbours.
- They separate on their own: by_ny is 0.374 for confirmed right peaks against 0.109 for decoys.
- Added to the 20 subscores, N at FDP_c 1% moves 1083 → 1069 (raw 02) and 1101 → 1028 (raw 01). The rule (≥ 5% on both runs) fails.
- Bulk-null q ≤ 0.05 FDP_c stays valid (0.042 and 0.043).

**Why our version adds less (`by_vs_osw.py`, raw 02).**
- On matched candidates, our fractions correlate only moderately with OSW's counts (Spearman 0.53–0.54).
- Conditional on our score, OSW's y-series reaches AUC 0.530 against at most 0.523 for ours.
- From OpenSwath's documented behaviour, its score is an **integer count** of charge-1 b or y ions at ≤ 10 ppm above an intensity floor (300) in the apex spectrum. Ours is a length-normalised fraction over charge 1 and 2 at 17 ppm, averaged over 3 spectra.
- A raw count grows with peptide length when the match is real. Charge-2 ions and the wider tolerance add random matches.

Next: our own count-based variant (charge 1, 10 ppm, intensity floor, apex spectrum only). The information is the same as OSW's; the implementation is our own.

**Own b/y count features, v2: PASS (`by_feats.py` v2, `BY_SET=counts by_test.py`, `by_cnt01.txt`, `by_cnt02.txt`).** Two features: integer counts of matched b and y ions, fragment charge 1 only, nearest peak within 10 ppm, intensity > 300, apex spectrum only. They are computed by our own numba code from the mzML. No OSW scores are used.

| run | ours, FDP_c 1% | + own counts | OSW series (diagnostic) | bulk-null q≤0.05 FDP_c | bulk-null q≤0.01 N / FDP_c |
|---|---|---|---|---|---|
| raw 02 | 1083 | **1338 (+23.5%)** | 1291 | 0.054 | 1248 / 0.005 |
| raw 01 | 1101 | **1269 (+15.3%)** | 1242 | 0.040 | 1226 / 0.007 |

- The predeclared rule (≥ 5% on both runs, bulk-null q ≤ 0.05 FDP_c ≤ 0.06) passes.
- On both runs, our counts beat both OSW's series features and OSW itself (1287 on raw 02).
- Honest limit: at FDP_c 0.5%, raw 01 drops slightly (1035 → 1000) while raw 02 rises (1008 → 1286).
- Plain TDC stays anti-conservative at q ≤ 0.05 (0.070 / 0.057), so bulk-null stays the primary q.
- Cost: about 2.5 min per run to stream the mzML.
- Lesson: for the ranking head, the information is in **absolute counts at a tight tolerance on the apex spectrum**. Length-normalised fractions with charge-2 ions and a wider tolerance dilute it.

Native recommendation: add apex-spectrum b/y count subscores (fragment charge 1, 10 ppm, intensity floor) to OpenDIA's scoring as part of the single model.
