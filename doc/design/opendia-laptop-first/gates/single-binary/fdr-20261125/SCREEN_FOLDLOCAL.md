# Dev screen: opt-in `-gbm_fold_local` (fold-local GBM training)

Status: dev screen. This is not a gate. These rules are written down before any replay runs.

## Change under test

In `OpenDIAGBM.cpp`, the default GBM picks its training rows (TDC positives and the top-2 % hard negatives) from pooled
scores over all folds. The held-out fold is dropped only after that choice is made. DIAG_CALIBRATION.md shows that held-out scores can then affect their own model.
With `-gbm_fold_local`, each outer fold f selects positives and negatives only from rows where fold != f, with the same
rule as the default. Between rounds those rows are rescored by inner cross-fit, using only rows where fold != f. The final model
for f is trained on the last selection and predicts f. Weights, the 1.5 M negative cap, the validation split, calibration and
the kept-window logic do not change. The option is advanced and off by default. Off output must be byte-identical.
The rule does not use species or entrapment information.

## Inputs (read-only)

- Features: `calsearch-20261126/runs/{D,H}_G0f_s{0,1}/*.features.tsv`, read in place on codon.
- Library: `entlib-odp5-20261118/runs/D_cur_s0/hela_r5.predicted_decoys.pqp`. This is the same library the runs used.

## Arms (same harness, same seed 0 for s0 and 1 for s1, 8 threads, on codon)

- P: `gbm_replay` pooled (the default code path).
- L: `gbm_replay ... fold_local`.
- Each run's own exported `score` column (written by the replay as `*.self.tsv`) is evaluated as a check on P. It is not an arm.

The harness is `foldlocal/gbm_replay.cpp`, linked statically against `libOpenDIA_impl.a` built from fork `4c6ca367fd` plus
`foldlocal/foldlocal.diff`. It writes kept rows only. The local download was too slow
(about 2.4 MB/s), so the gz copies are not used. The key table `hela_keys_e.tsv` (id, traml, decoy, entrapment) was
exported read-only from the library above. Hashes are in `FOLDLOCAL.sha256`.

## Metrics

Same as `devAB/eval_replay.py`: precursor-level paired TDC, q <= 1 %. N is the number of accepted targets. E is the number of accepted
entrapments (`instr(accession, '_p_target')`), used for evaluation only. The combined bound is 2E/N.

## Promotion rule (fixed)

L moves to a frozen remote precheck only if both of these hold:
1. The pooled combined bound for L is lower than for P on both D and H (2 seeds each).
2. Pooled N for L is >= 0.995x of P on both D and H.

A bound <= 1 % is the end goal. It is not required at this screen. All four runs are reported. No reruns or retuning are done after
results are seen.
