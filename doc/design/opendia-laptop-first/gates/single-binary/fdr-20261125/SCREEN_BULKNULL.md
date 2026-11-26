# Dev screen: opt-in `-bulk_null_q` (label-free bulk-null weighted TDC q)

Status: dev screen on used data. This is not a gate. The rule below is fixed and pushed before any job runs.

## Change under test

`bulkNullQ` (OpenDIAGBM) is CALIB §9 ported to the native pipeline. The final ranking (pooled paired TDC of the final score)
is not changed. Only the decoy count in q is reweighted:
- The bulk is the precursor winners below their fold's median decoy-winner score.
- Per fold f, a HistGBM with default options learns target versus decoy on the bulk winners of the other folds. It uses the 49
  feature columns, not the score.
- Each decoy winner of f gets the weight w = exp(margin) * n_d / n_t, clipped to [0.2, 5] and normalised to mean 1 over f's bulk decoys.
- The q-value is q = min(1, (sum of w over decoys + 1) / targets), with a reverse cumulative minimum.

No species, entrapment or label information is used. The option is advanced and off by default. When it is off nothing changes, because the
code path is not entered. Unit test (e) checks three things on synthetic data with known truth:
- The ranking is kept.
- A shifted null goes from FDP 12.7 % at q <= 5 % to 4.3 %.
- An exchangeable null stays neutral (5153 versus 5155 IDs).

## Inputs (read-only)

- `calsearch-20261126/runs/{D,H}_G0f_s{0,1}/*.features.tsv` (default OpenDIA, the evshare container, with `-features_out`).
- A key table exported read-only from `entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp`. It carries the
  DIAG_STEPA label correction: 2912 exact-human "entrapment" precursors are relabelled human, giving r = 6.9176. This is identical to the RTNULL gate.
- The harness `bulknull/bulknull_replay.cpp` is linked statically against fork `306a9ed402` plus `bulknull/bulknull.diff`. It
  rebuilds the input as `gbm_replay` does and calls `bulkNullQ` on the exported final score. Hashes are in `BULKNULL.sha256`.

## Metric (gate scale, `bulknull/eval_bulknull.py`)

At q <= 0.01: N, E, combined bound E(1+1/r)/N, with a 90 % Poisson interval on pooled E. N@FDP_c 1 % is taken on the exported score
order. It is the same for both q columns because the ranking is identical.

## Rule (fixed)

The screen passes only if all of these hold:
1. **Fidelity**: for every kept row, the replayed plain q equals the exported q (log line). This shows the ranking and
   competition were rebuilt exactly.
2. **Calibration**: the pooled combined bound under `bn_q` is <= 1.00 % on D and on H (2 seeds each).
3. **No ID loss at true FDR**: N@FDP_c 1 % under `bn_q` is >= 1.00x of default. This holds by construction once 1 holds.

N at reported q <= 1 % will drop. It is reported, not gated, as agreed with the user on 26 Nov. If the screen passes, the next step is a
frozen precheck on held-out PXD017703 runs. All four runs are reported. There are no reruns and no retuning after results are seen.
