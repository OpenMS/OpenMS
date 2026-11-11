# D11: which IDs are missed, and why (descriptive, raws 01 and 02)

This is a descriptive analysis only. Nothing is trained, and no pipeline changes. The scripts are `misses.py` and `misses_b.py` to `misses_f.py`, and the output is `MISSES_D11.txt` plus `misses_D11_raw0{1,2}.tsv`.

**Reference set:** DIA-NN 2.3.1 precursors at 1 % that match a library precursor exactly (`astral-diann-gap-20261108/ref_exact.parquet`), human only.
- This gives 2743 precursors on raw 01 and 2731 on raw 02.
- DIA-NN is a reference only.
- For this library, the reference set is the realistic ceiling (GAP.md).

**D11 setting:** the scores of PLAN_V3_P2b. The ranking is the same as in `evaluate_v2.py`: precursor winners, per-fold decoy median/MAD. The cut is at entrapment FDP_c 1 %.

## Funnel

| Stage | raw 01 | raw 02 | Where the loss is |
|---|---|---|---|
| A: not in the scoring population | 185 | 151 | Search: the pair's best transition_support is ≤ 2 (the population needs ≥ 3) |
| B: no candidate window within 10 s of the DIA-NN RT | 166 | 180 | Search / RT: the DIA-NN RT is far from the calibrated library RT |
| C: a right window exists, but a wrong window wins | **421** | **383** | **Window choice by the D11 score** |
| D: the right window wins, but it is below the 1 % cut | **681** | **744** | **Scoring power at the top** |
| E: accepted | 1290 | 1273 | |
| D11 accepted at FDP_c 1 % (all) | 1348 | 1334 | 96 % / 95 % are DIA-NN-confirmed |

The reference IDs D11 misses at 1 % are 1453 on raw 01 and 1458 on raw 02. Of these, 55 % are group D and 27 % are group C. The search groups A and B together are 24 %.

## C: wrong window wins (421 / 383)

**This is a D11 defect, because BASE picks better here.**
- For the precursors with a right window, D11 picks a right window for 1971 / 2017 and BASE for 2278 / 2306.
- BASE is right where D11 is wrong for 332 / 307 precursors. D11 is right where BASE is wrong for only 25 / 18.

**On every feature the right window is better than the wrong one that D11 picked.**

| | transition_support | coelution | library correlation | byion_y | abs RT error (s) |
|---|---|---|---|---|---|
| right window (median) | 4 | 0.15 | 0.31 | 3 | 20 |
| wrong winner (median) | 3 | 0.00 | −0.09 | 1 | 44 |

The wrong winner has the higher transition_support in only 3 % of cases.

**Cause:**
- **The D11 score carries little order below its top 1 %.**
  - In every fold, the window scores sit in a narrow band up to the 99th percentile (raw 01, fold 0: −876 to −833). The top 0.1 % then jumps to +1733.
  - Inside that band, the Spearman correlation of the D11 score with transition_support, or with a simple sum of the features, is about −0.1 to +0.03. In the top 1 % it is +0.74 to +0.82.
  - A linear fit of the score on the 11 features gives R² = 0.08.
- **Almost all group C windows are in that band:** 410 / 372 of the right windows, and 358 / 332 of the wrong winners.
- **So the choice between two weak windows of one precursor is close to random.**

A probable mechanism, not proven: the tanh features saturate for weak windows. The decoy-whitened readout gives them large weights, so small saturated differences outweigh the linear features. The objective only rewards separation of positives and decoys in the mean. Nothing in it asks the model to order windows within a precursor.

**Possible fix**, in a future plan: train the window order of the same single model on purpose. For example, add a within-precursor ranking term: anchor window above the other windows of the same precursor. Or limit the saturation.
- This keeps one model.
- It could recover up to about 330 precursors per raw at the window step. Not all of them would pass the cut.

## D: right window wins, but scores below the cut (681 / 744)

**This is the largest group. Most of these are far below the cut, not just under it.**
- The cut is about 544 on the normalised score. The group's median is only 96 / 105, so half of them sit below a fifth of the cut.
- 187 / 201 are within half of the cut, and 87 / 80 are within 0.8 of it.

**DIA-NN is very confident about about half of them** (q ≤ 1e-3 for 51 % / 55 %, median q 0.001). So these are real peptides with weaker evidence in our 6-transition extraction.

**Against the accepted IDs (E), they are weaker on co-elution and fragment agreement.**

| | coelution | library correlation | byion_y | byion_b | log intensity |
|---|---|---|---|---|---|
| E accepted (median) | 0.79 | 0.86 | 7 | 4 | 5.1 |
| D (median) | 0.37 / 0.40 | 0.57 / 0.59 | 4 | 2 | 4.7 |
| top 1 % decoy winners (median) | 0.31 / 0.34 | 0.29 / 0.33 | 4 | 2 / 3 | 5.0 |

- **Against the top 1 % of decoy winners, group D is separated mainly by one feature.** transition_support has AUC 0.82. Every other feature is at 0.62–0.75, and the RT error is worse than for the decoys (AUC 0.36).
- **D's median features sit close to the top decoys' medians on several measures.** These are coelution, byion and intensity. The 11 features cannot rank these precursors above the top decoys.
- **The extraction carries less evidence for them.** D has a lower intensity than E (log 4.7 against 5.1) and fewer matching fragments (byion_y 4 against 7).
- **The cut is set by entrapment targets with good features.** The first entrapment target winners have transition_support 5–6, coelution 0.3–0.9 and library correlation 0.5–0.98. They come at ranks 797, 946 and 1145 on raw 01, and at 464, 663 and 981 on raw 02.
- **Beyond the cut:** in the next 1000 ranks, there are 185 / 172 entrapment targets, 499 / 558 DIA-NN-confirmed IDs and 490 / 430 targets with no DIA-NN support.
- **Relaxing the cut does not help much:**

| FDP_c | 1 % | 2 % | 5 % |
|---|---|---|---|
| D11 N | 1348 / 1334 | 1572 / 1572 | 1860 / 1861 |
| DIA-NN-confirmed among them | 1290 / 1273 | 1463 / 1476 | 1632 / 1668 |

**Cause:** this is an evidence limit in the extraction, with 6 transitions, as also found in GAP.md. It is not mainly a modelling gap. A different model on the same 11 features is unlikely to separate them well, because the features themselves overlap with the top decoys. Larger gains need more evidence per precursor:
- more fragments per precursor (the library has 6);
- MS1 isotope features;
- better co-elution features at low intensity;
- predicted-property deviation features (see the on-hold ReDeeM direction).

## A and B: upstream losses (336 / 331)

**A, not in the population (185 / 151):**
- The search finds a candidate for 172 / 140 of them. But the pair's best transition_support is 1–2, which is below the population threshold of 3.
- So DIA-NN sees these peptides, but our search matches only 1–2 of their 6 transitions.

**B, no window near the DIA-NN RT (166 / 180):**
- The DIA-NN RT is far from our calibrated library RT: median 130 / 127 s. The calibration σ is 21 s, and 135 / 151 of them are beyond 3σ.
- For comparison, the median for accepted IDs is 14 s.
- These are library RT outliers. The search window cannot reach them without a wider RT window, at a cost in decoys.

## Summary

| Loss | Size per raw | Cause | Fixable inside one model? |
|---|---|---|---|
| C: wrong window | about 400 | D11 does not order weak windows of one precursor (flat score below its top 1 %); BASE picks better | Yes. Train the window order on purpose (ranking term, less saturation) |
| D: below the cut | about 700 | Weak evidence in the extraction; its features are close to those of the top decoys | Partly. Needs more evidence per precursor (fragments, MS1, prediction deviations), not only a better model |
| A: population | about 170 | Only 1–2 transitions match | Search / extraction |
| B: RT | about 170 | Library RT outliers beyond 3σ | RT calibration or a wider search window |
