# Raw 01 design, round 2: compact TOPAZ under a compute budget

Written 2026-11-08 16:22 UTC, before any round-2 code is built or run, after the round-1 result (`RESULT.md`: arm A gave 0, arm B was stopped for cost). Raw 01 only, with target/decoy labels only. Raw 03 stays untouched.

## User requirement
OpenDIA must stay fast. Training and scoring must be quick and still learn.

## Diagnosis of arm A
- About 100K parameters on about 600 seeds per fold overfit. Validation loss reached 0.0005, but the decoy tail was heavy (decoy z 0.999 quantile 9 to 13, against grid2d's 7).

## Variable: model capacity plus regularisation, then MIL
- `-topaz_width 8`: conv 8-16-16, embeddings 16, MLP 32-16. About 13K parameters, about 2× grid2d's 6141 and about 8× fewer than arm A.
- `-weight_decay 0.01` (AdamW).

## Arms (the frozen stage-2 search command, otherwise unchanged)
- **C:** `-cnn_arch topaz -topaz_width 8 -weight_decay 0.01` (seeds training).
- **D:** C + `-cnn_training mil -mil_bags_per_epoch 4000 -epochs 10`. About 21K windows per epoch, so about 0.2M window updates per fold.

## Compute budget (part of the pass rule)
- The CNN stage (the tool's "Native CNN crossfit seconds", covering training and scoring of all 3 folds) must be ≤ 365 s, i.e. 1.5× grid2d's 243 s on raw 01.
- If an arm exceeds 2× the budget while running (730 s of CNN time), it is stopped and counted as a fail.

## Pass rule (unchanged from PLAN.md)
- ≥ 732 targets at stage-2 pooled paired q ≤ 0.01;
- R ≥ 0.83 at the decoy median;
- within the budget.
- If both arms pass, freeze the higher count. On a tie within 1%, freeze C.
- If neither passes, `grid2d` stays frozen.

Reported but not part of the decision: the same diagnostics as `PLAN.md`.

## Default equivalence
After the rebuild and before raw 01, HeLa frozen `crossfit.tsv` must be byte-identical.
