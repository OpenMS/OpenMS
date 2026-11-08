# Raw 01 TOPAZ design: result

Rule in `PLAN.md`: at least 732 targets at stage-2 pooled paired q ≤ 0.01 (the reference is 717), with R ≥ 0.83 at the decoy median.

## Arm A: `-cnn_arch topaz`, trained on seeds. FAIL

| | grid2d (reference) | topaz A |
|---|---|---|
| Stage 2, q ≤ 0.01 / 0.02 / 0.05 | 717 / 1015 / 1412 | **0** / 152 / 853 |
| CNN only, q ≤ 0.01 / 0.05 | 321 / 634 | 0 / 207 |
| Decoy z 0.999 quantile per fold | 7.6 / 7.1 / 7.1 | 9.5 / 9.2 / 12.8 |
| Validation loss per fold | 0.051 / 0.034 / 0.078 | 0.019 / 0.0005 / 0.041 |
| R at decoy median / at 0.05 cut | 0.926 / 0.861 | 0.934 / 0.882 (guard passes) |
| CNN time (train + score, 3 folds) | 243 s | 228 s |

The seed separation is better: validation loss is lower and held-out seed AUC is 1.0. The decoy tail, though, is much heavier, so few targets survive.
- The model has about 100K parameters and is fitted to about 600 training seeds per fold. It overfits the seeds and scores some decoy windows with extreme confidence.
- The architecture is not slower than grid2d: the CNN stage took 228 s against 243 s. The 16 min wall time is dominated by MS1+MS2 extraction (about 10 min) and stage 2.

## Arm B: `-cnn_training mil`. Stopped, not evaluated

I stopped it (by PID) 19 s after it started, at extraction. The design draws 40000 bags × about 5.2 windows per epoch, with forward and backward passes, for up to 40 epochs over 3 folds. That is about 8M window updates per fold. Going by arm A's inference throughput, it would take on the order of an hour or more. That breaks the requirement that OpenDIA stay fast (user, 2026-11-08 16:16 UTC).

## Decision

`grid2d` stays frozen (717). `-cnn_arch topaz` stays opt-in and experimental. Any next variant must be predeclared with an explicit compute budget.

## HeLa timsTOF diagnostic (IM path; not a gate, not part of the decision)

This is the frozen HeLa teacher-free command (`-ss_iterations 1`, fold-wise paired) with `-cnn_arch topaz` at width 32, in `hela-topaz/`.

| | grid2d (frozen) | topaz |
|---|---|---|
| Iteration 0 targets at q ≤ 0.01 | 12906 | 13882 |
| Iteration 1 targets at q ≤ 0.01 | **13253** | **15119** (+14%) |
| CNN seconds (2 rounds × 3 folds) | about 145 | 342 |

On a strong IM run with 6000 seeds, the same architecture that collapsed on the 600-seed Astral design run gains 14%. This fits the diagnosis of over-capacity on few seeds, rather than an architectural defect.

These are decoy-count estimates only. HeLa was earlier used for default-equivalence checks, not as a held-out evaluation.

## Round 2 (`PLAN2.md`): compact width 8 with weight decay 0.01, under a CNN budget of ≤ 365 s. Both arms FAIL

The HeLa default `crossfit.tsv` stayed byte-identical after the rebuild (`hela-default-cmp2/`).

| | grid2d (reference) | C: topaz w8 (seeds) | D: topaz w8 MIL (4000 bags × 10 epochs) |
|---|---|---|---|
| Stage 2, q ≤ 0.01 (rule: ≥ 732) | 717 | **473** | **653** |
| Stage 2, q ≤ 0.02 / 0.05 | 1015 / 1412 | 733 / 1177 | 1028 / **1539** |
| CNN only, q ≤ 0.01 / 0.05 | 321 / 634 | 0 / 623 | 204 / 582 |
| Decoy z 0.999 quantile per fold | 7.6 / 7.1 / 7.1 | 7.5 / 8.7 / 8.9 | 4.1 / 4.1 / 4.7 |
| R at decoy median (guard ≥ 0.83) / at 0.05 cut | 0.926 / 0.861 | 0.922 / 0.871 | 0.900 / **0.710** |
| CNN seconds (budget ≤ 365) | 243 | 146 | 253 |
| Wall time, peak RSS | 14:20, 23.7 GB | about 16 min, 22.7 GB | 14:14, 25.7 GB |

**C (compact, seed training).** It is fast, but it still overfits the 600 seeds: validation loss is 0.003 to 0.02 and the decoy tail is heavier than grid2d's.

**D (MIL).** Training on every precursor bag removed the heavy decoy tail: the decoy z 0.999 quantile fell from about 7 to about 4. D is the best arm at q ≤ 0.05 (1539 against 1412), but not at the predeclared q ≤ 0.01.
- Its bag validation loss stayed at 0.69 (about ln 2). That is expected, since most target bags hold absent peptides, so a bag label is a weak signal.
- R at the 0.05 cut fell to 0.71, against a π0 estimate of about 0.86. Targets are pushed out of the extreme low tail more than for grid2d. Under MIL this can be a mild target-vs-decoy distribution shift rather than peptide evidence.
- The predeclared guard (R at the median) passes, but this is a warning to keep in mind.

## Decision (rounds 1 and 2)

No TOPAZ arm passed. `grid2d` with stage-2 `search` (717) stays frozen.
- `-cnn_arch topaz`, `-topaz_width`, `-weight_decay` and `-cnn_training mil` stay opt-in and experimental.
- All four TOPAZ arms (A, B, C, D) have now been used on raw 01. Further raw 01 iteration on this architecture risks fitting design noise.
