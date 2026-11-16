# Percolator rescoring in NASE

`rescoring:enable` rescores NASE's spectrum matches with the in-process
Percolator of OpenMS (`Percolator::rescore`). It is off by default and needs
decoy sequences (`fdr:decoy_pattern`). The search, the candidate methods and
the native scores are unchanged; rescoring only changes the score that the
q-values are calculated on.

```sh
NucleicAcidSearchEngine -in spectra.mzML -database targets_decoys.fasta \
  -out results.mzTab -fdr:decoy_pattern DECOY_ \
  -rescoring:enable -rescoring:out_features features.tsv
```

## What is rescored

Percolator gets one row per spectrum: its best hit by the native score. The
other hits that the search retained (`rescoring:top_hits`, default 2, at least
`report:top_hits`) only describe how far ahead that hit is (`delta_next`,
`delta_best`). This is the usual Percolator input of search engines (one PSM
per spectrum with delta scores), so Percolator's PEPs and NASE's q-values
refer to the same competition: the best hit of each spectrum.

As a consequence:

- The order of the hits of a spectrum stays the native order; rescoring never
  promotes a runner-up.
- Hits tied with the best native score share its Percolator score and PEP.
- With `report:top_hits` above 1, the lower-ranked hits get no Percolator
  score and therefore no q-value; with `fdr:cutoff` below 1 they are removed.

Runner-ups are not rescored as rows of their own. Every RNA target has an
isobaric reversed decoy, so most spectra have a strong decoy runner-up.
Without Percolator's target-decoy competition, its training q-values and PEPs
count these runner-ups as independent chances; with it, the hits that lose the
competition come back as score 0 / PEP 1 and cannot be told from real values.

## Features

Rows come from `RNARescoringFeatures`. Their columns depend only on the
parameters (and on the charges of the rescored hits), never on the label of a
hit, so targets and decoys are described identically.

| Feature | Meaning |
| --- | --- |
| `native_score` | Hyperscore, MVH or route score of the hit |
| `delta_next` | Score minus the best lower or equal score of a hit with another sequence (other charges or isotopes of the same sequence are no competitors) |
| `delta_best` | Score minus the best score of the spectrum (0 for the rescored hit) |
| `abs_ppm_error`, `isotope_offset` | Precursor mass error and isotope offset |
| `charge_<z>` | One-hot precursor charge over the charges of the rescored hits |
| `length`, `modifications`, `adduct` | Oligonucleotide length, modified residues, adduct flag |
| `frac_<series>` | Distinct matched ordinals of each `fragment:ions` series / (length - 1) |
| `matched_ions` | Matched (series, ordinal) pairs / (series × (length - 1)) |
| `matched_intensity` | Intensity of the distinct matched peaks / total intensity of the scored spectrum |
| `longest_series` | Longest run of consecutive ordinals in one series / (length - 1) |
| `coverage` | Backbone positions explained by any series / (length - 1), as `report:require_coverage` |

Optional groups, present when the corresponding mode is on:

| Mode | Features |
| --- | --- |
| `scoring:method=route` | `route_score`, `route_nodes`, `route_intensity_ratio`, `route_prefix_length`, `route_suffix_length` |
| graph candidate methods except `graph_routes` | `graph_tag_length` (0 without a supported tag) |
| `diagnostic:filter` | `diagnostic_supported_mods`, `diagnostic_unsupported_mods` |
| `precursor:correct_with_flash` | `flash_corrected`, `flash_cosine`, `flash_snr` |

Neutral losses other than a-B are not counted as ladder rungs.

## Output

The PSM section gets `percolator:score` (MS:1001492) and `percolator:PEP`
(MS:1001493) next to the native score; the PSM-level q-value is calculated on
the Percolator score by NASE's usual FDR path. `rescoring:out_features`
writes one line per rescored spectrum: native ID, scan index, sequence,
charge, decoy label, the features and Percolator's score, q-value and PEP
(Percolator's own q-values are for inspection only).

## Fallback

NASE uses the native score, with a warning, if

- fewer than `rescoring:min_decoys` spectra (default 100) have a decoy as
  their best hit, or no target does; Percolator's model and PEPs are
  unreliable on fewer decoys, and it does not reliably refuse such input, or
- Percolator's sanity check rejects the learned model (no target at q < 1% on a
  held-out fold, or fewer than its best single feature), or
- Percolator throws.

The output is then the same as without rescoring. Percolator itself would score
a rejected model with its best single feature, which can rank far worse than
the native score: on a HAP1 file searched with an entrapment database, a
high-scoring decoy left no target at q < 1% on a fold, and the single-feature
score kept 16 / 102 / 329 spectra at 1% / 5% / 10% FDR, against 0 / 2,095 /
2,474 for the native score.

`rescoring:train_fdr` (default 0.01) sets the FDR at which Percolator selects
positive training examples and checks the learned model (its `train_fdr` and
`test_fdr`). A higher value can let hard searches train: with 0.05, Percolator
trained on all three HAP1 clone1 entrapment searches, accepting 9.9% / 5.7%
more spectra than the native score at 5% / 10% q-value, with estimated
entrapment FDPs of 4.4% / 8.8% (native: 4.0% / 7.6%), pooled over the three
files.

## Reproducibility and limits

Results depend on `rescoring:seed` (Percolator's cross-validation folds), the
input order and `threads`. They are reproducible for the same three.

The test data are synthetic (target and reversed-decoy spectra plus poor
spectra with a few true fragments) and check the wiring, not the model: on one
set Percolator trains a model (which keeps fewer targets at 1% FDR than the
native score there), on the other it rejects it. Reversed RNA
decoys also share many fragment masses with their targets (on the synthetic
data, a decoy matched a median 32% of its possible ions in its target's
spectrum, against 40% for the target); shuffled decoys
(`DecoyDatabase -type RNA -method shuffle`) may separate better. Neither a
nominal q-value yield nor this test shows that the FDR is controlled on real
data; an entrapment experiment is the appropriate check.
