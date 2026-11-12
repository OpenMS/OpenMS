# Original-thread reconstruction amendment

Approved by coordinator rat at 2026-11-12T09:12:08Z: restore original 10 Torch threads without CPU affinity restriction, retaining the failed four-thread attempt. This changes execution numerics only, not the model, seeds, optimizer, selection thresholds, fold assignment, validation split, or reference arrays.

New outputs are exclusively under `score/original_threads/`. The source `train_v4.py` is imported, never edited, and its `run()` is never called. Original OMP settings are restored (6 raw01, 10 raw02), with `TRAIN_THREADS=10` and no MKL thread override. As before, checkpoints and round memberships are saved without evaluation labels. Additional r1/r2 snapshots record actual full-score training selectors and the linear selector of the identical prior-round model. Evaluation joins are gated on both raws passing Spearman >=0.99999, maximum absolute score difference <=0.001, and W-selected row agreement 1.

The prior four-thread gate failure and raw01 artifacts remain untouched. Its reporting-only W agreement correction remains documented in the parent `AMENDMENTS.md`.
