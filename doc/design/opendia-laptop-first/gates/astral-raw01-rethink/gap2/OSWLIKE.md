# A teacher-free OSW-like scorer on the shared XICs (raw 01; raw 02 replication pending)

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

1. **The gap closes with a teacher-free scorer on the same traces.** A small MLP on 20 OSW-style subscores, with decoy-only labels, gives 1214–1221 at a true FDP_c of 1%. That is above OSW+Percolator (1151), with decoy-estimated q ≤ 0.01 calibrated to FDP_c 0.012. The LDA alone gives 1011. So the 757 → 1151 gap is not missing information and not extraction. It is **which subscores we compute plus a non-linear combiner.**
2. **Our candidate generator is not the main problem either.** The same subscores on our streaming search's apexes give 1023 (LDA) and 1110 (MLP) at FDP_c 1%, against 1011 and about 1215 on the XIC maxima. Coverage of the right peak is 2424 of 2754 (88%) against 2528 (92%). This costs at most about 9% with the MLP and nothing with the LDA.
3. **The subscores that matter, compared with what OpenDIA uses today:**
   - **Library agreement is the largest single group.** Dropping it loses 40%. We have library correlation, but not dot product, Manhattan, RMSD or spectral angle on the integrated fragment areas.
   - **Lagged xcorr is the backbone.** Without it the semi-supervised loop finds no positives at the first step (0 at q ≤ 0.15, even when initialised on the library dot product). No single non-xcorr subscore reaches 15% FDR alone. Every single subscore alone is weak: the best is log S/N at 510, then xcorr co-elution at 304. The strength comes from the combination.
   - **MS1 adds nothing here,** and signal and shape/RT add about 3–4% each.
4. **Peak boundaries and baseline subtraction are not the lever** (see `GAP2.md`). The valley/5% boundaries used here give a median width of about 11 s, close to OSW's 13 s and DIA-NN's 9.8 s. Our search groups' 46 s spans should still not feed features.

## Caveats

- Raw 01 only so far; a raw 02 replication is running (`raw02_*.txt`). The OSW reference on raw 02 is 1159 / 1287 / 1351 / 1389.
- XIC extraction is OSW's (about 34 ppm, ±97 s, 1.6 s cycle). Our native extraction must reproduce it. Bright-peak agreement is shown; weak-peak agreement is not.
- The MLP's FDP_c at q ≤ 0.01 is 0.012. This is mildly anti-conservative, as with the earlier stage-2 MLP (0.018) but less so, and within the 1.5× gate rule. It is a single run, with only about 8 entrapment hits behind it.
- Two seeds agree at 1%; the 0.5% and 5% points vary more.

## Implication for the design (the user's single-model direction)

The richer learned model should **take these subscore families as inputs alongside the trace encoder**, not replace them:

- lagged xcorr co-elution and shape, unweighted and library-weighted;
- library agreement on integrated areas (correlation, dot product, Manhattan, RMSD, spectral angle);
- intensity share and S/N.

It should be computed per candidate on valley-bounded peaks from a per-precursor XIC.

Section 5 of `GAP2.md` showed that a small trace CNN alone does not recover them. A single cross-fitted MLP head over [CNN embedding + these about 20 subscores], with decoy-only semi-supervised labels, is the minimum design this evidence supports. It replaces stage 2 and its opt-in feature sets.
