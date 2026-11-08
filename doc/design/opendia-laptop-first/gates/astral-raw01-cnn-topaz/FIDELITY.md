# Port fidelity: OpenDIACNN `-cnn_arch topaz` against redeem-topaz (`add/topaz`)

Checked 2026-11-08 17:04 UTC, after the four raw 01 arms. This audits whether the port is correct, so that the failed arms are not blamed on the wrong cause. It changes nothing about the decisions above.

## Numeric parity (scratch harness `$JCODE_SCRATCH_DIR/topaz-parity/`)

- **Setup.** The Rust side is a `git archive` of `add/topaz` (redeem worktree untouched) with a small example (`parity_coelution.rs`). It is built offline against a workspace limited to `redeem-io` and `redeem-topaz`, with features `io-sqlite,io-parquet`. The C++ side compiles the exact `OpenDIACNN.cpp` source blocks, extracted by marker, against LibTorch 2.12.1.
- **Inputs.** The same deterministic pseudo-random traces on both sides: MS2 [3,6,32] and MS1 [3,3,32]. Learned similarity weights are zeroed on both sides, matching candle's `VarBuilder::zeros`.
- **Result.** The CoelutionHead output (10 statistics plus the 8-d embedding) for C = 6 and C = 3, and `ms12_features`, agree to within **1e-6** in every column. That is float32 rounding.
- **Parameter counts per block are identical:**
  - `ConvBranch` (dual input): 27232 for C = 6 and 27040 for C = 3;
  - `CoelutionHead`: 1256 and 488.

## Differences from TOPAZ, beyond those declared in `PLAN.md`

These were found by reading `train/trainer.rs`, `config.rs`, `infer/pipeline.rs` and `trace_window.rs`.

| | redeem TOPAZ | OpenDIACNN topaz | Effect |
|---|---|---|---|
| Input scaling | `normalize_max: true` in every CLI and run path: each [C,L] block is divided by its global max before the dual view | raw intensities; the dual view's log1p channel sees absolute log intensity | the log1p view carries absolute intensity, which can separate runs and peptides by abundance instead of shape |
| Optimiser | AdamW, weight decay 0.01 by default | Adam (arm A); AdamW 0.01 only in arms C and D | arm A had no weight decay |
| LR schedule | cosine warm-up (5%), min LR ratio 0.05 | constant 1e-3 | none |
| Gradient clipping | max norm 5 | none | none |
| Losses | bag BCE + 0.1 × pairwise pos-neg softplus + 0.1 × in-bag ranking (margin 1) | bag BCE only (D), window BCE (A, C) | in-bag ranking pushes the winner above the other candidates in a target bag; missing in D |
| Patience | 3 | 6 (40-epoch cap; D 10 epochs) | none |
| Heuristic features | OSW feature-table scores enter the MLP | none (by design; they enter stage 2) | declared in `PLAN.md` |
| Window length | `l` = 64 (XIC), 258 (XIM) | 32 RT bins, 8 IM bins | declared |

## Reading

- The architecture itself was ported faithfully, and its features are numerically identical.
- The training regime was not. The largest undeclared difference is **input scaling**: TOPAZ max-normalises each trace block, OpenDIACNN does not. A log1p channel on raw Astral intensities, with a free bias, makes overfitting to abundance easier, and that fits the heavy decoy tail seen in arms A and C.
- TOPAZ's auxiliary ranking losses are also missing from the MIL arm.
- Any further TOPAZ work should first adopt `normalize_max` and the TOPAZ loss and optimiser defaults, as one declared variable. It must be designed on a run that is not yet spent; raw 01 has now been used for four arms.
