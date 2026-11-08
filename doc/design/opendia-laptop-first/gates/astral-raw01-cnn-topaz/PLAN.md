# Raw 01 design: TOPAZ-style modality-flexible CNN (`-cnn_arch topaz`)

Written 2026-11-08, before any code or run. Raw 01 is design data with target/decoy labels only. No entrapment labels and no DIA-NN or OpenSWATH output are used. Raw 02 is spent and raw 03 stays untouched.

## Why

The default CNN (`grid2d`, 6141 parameters) takes a 6 × 32 RT × 8 IM grid. On Astral 7/8 of that grid is empty, MS1 is never seen, and it has no explicit co-elution. Stage 2 gets 717 at q ≤ 0.01 with `cnn_z` plus the search features. Same-library OpenSWATH+Percolator reaches about 1441.

## What is borrowed from TOPAZ (redeem `add/topaz`)

- 1D conv branch over RT on a dual input view: per-channel max-normalised plus log1p (`trace_input.rs`, `conv_encoder.rs`). Layers 1×1→32, k3→64, k3→64, max+mean pooling, then a projection to 64.
- An optional MS1 branch on M+0, M+1, M+2 traces, plus the 4 MS1-vs-MS2 agreement statistics (mean/max cosine, apex difference, apex std).
- Co-elution head (`coelution.rs`): mean/min cosine at lag 0 and with up to ±2 bins of lag, soft-apex std/range/MAD, entropy mean/min, peak-like fraction, and a learned 8-d embedding of the pairwise similarities.
- Candidate MLP 128→64→1 with dropout 0.1, and masked-max MIL over a precursor's candidates (`topaz.rs`, `losses.rs`).

## Changes from TOPAZ

- **IM is a separate optional branch with an explicit mask.** The XIM is the fragment mobilogram summed over the core RT bins 12..19 of 32 (±3.75 s; written before implementation). Its embedding and co-elution features are multiplied by `has_im`, so a run without IM adds exactly zero, not learned biases on zero-filled input. `has_im` and `has_ms1` are also scorer inputs. One architecture covers Astral (XIC + MS1) and timsTOF (XIC + MS1 + XIM).
- **No OSW heuristic features in the CNN.** Search features enter in stage 2 as before, so the CNN stays trace-only plus library intensities. Two library-agreement statistics are added: cosine of sqrt core-summed fragment intensities against sqrt library intensities, and the same without the sqrt.
- The search's 32 RT bins (±15 s) are kept. No re-extraction.

## Options (opt-in, default byte-identical)

- `-cnn_arch grid2d|topaz`: `grid2d` is the default and its code path is untouched.
- `-cnn_training seeds|mil`:
  - `seeds` (default) trains on the `-seeds` windows as now.
  - `mil` trains on every in-RT-window precursor bag: target = 1, decoy = 0, masked-max over its candidates. Each epoch draws `-mil_bags_per_epoch` (default 40000) training bags, and validation uses a fixed subsample of the same size.
  - Requires `-cnn_arch topaz`. Excludes `-ss_iterations`.
  - In `mil`, seed labels are not used for training. They are still used for the stage-2 mass offsets.

## Runs (raw 01, the frozen stage-2 search command otherwise unchanged)

Base command: `../stage2-search-design/run.sh` with `-stage2_feature_set search`.
- **A:** + `-cnn_arch topaz` (seeds training). Tests the architecture alone.
- **B:** + `-cnn_arch topaz -cnn_training mil`. Tests MIL training on top of the architecture.

## Decision rule (fixed now)

- **Primary metric:** target precursors at stage-2 pooled paired q ≤ 0.01. The reference is 717 (`stage2-search-design/search`).
- A variant **passes** if it reaches ≥ 732 (+2%, the rule used for the isotope group).
- **Null-symmetry guard (required to pass):** among in-window CNN precursor winners scoring below the decoy-winner median in their fold, count targets T_low and decoys D_low. R = (T_low / D_low) / (N_T / N_D) must be ≥ 0.95.
  - Under a decoy-like null for absent targets, R ≈ π0 ≈ 1. R < 0.95 means the whole target distribution is shifted, i.e. the network learned a target-vs-decoy artifact rather than peptide evidence.
  - This matters most for MIL, where every target bag is labelled 1.
  - R is also reported for the `grid2d` reference.
  - **Amendment (2026-11-08 15:58 UTC, after computing R on the frozen `grid2d` reference only; no variant output existed yet):** the 0.95 threshold was miscalibrated. Present targets also leave the low region, so R estimates π0, not 1.
    - Reference R by decoy-quantile cut: 0.861 (q 0.05), 0.861 (0.10), 0.889 (0.25), 0.926 (0.50), 0.961 (0.75).
    - The 0.05-quantile value 0.86 is a Storey-style π0 estimate for this run.
    - **Replacement guard:** R at the decoy median must be ≥ 0.83, the reference π0 estimate minus 0.03.
    - A model that pushes more targets out of the below-median region than the run's estimated present fraction allows is treated as having learned a target/decoy artifact. It fails regardless of count.
    - R at the 0.05 cut is also reported.
- If both pass, freeze the higher count. On a tie within 1%, freeze A (simpler).
- If neither passes, `grid2d` stays frozen and the topaz code stays opt-in and experimental.
- **Reported but not part of the decision:**
  - CNN-only pooled paired count at q ≤ 0.01 (reference 321);
  - counts at q ≤ 0.02 and q ≤ 0.05;
  - per-fold counts;
  - stage-2 weights;
  - decoy z 0.999 quantile per fold.
- **Default equivalence** is required before any run: the HeLa frozen CNN command (`../hela-cnn-default-cmp-stage2/run.sh`) must give a byte-identical `crossfit.tsv`.
- **HeLa `-cnn_arch topaz`** (seeds, `-ss_iterations 1`, fold-wise paired) is an IM-path diagnostic only, reported against 13253. It is not part of the decision.

## Limits

- Design-run counts are decoy-count estimates, not calibrated FDR.
- With two variants on one design run, the winner's count is optimistically biased. Only raw 03 under a hashed gate can evaluate it.
