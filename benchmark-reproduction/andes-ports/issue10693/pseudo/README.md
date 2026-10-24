# #10693: pseudo-transfers as ground truth for PIP-ECHO's transfer FDR

`pseudo_transfers.patch` applies on top of `../extraction/proto_extraction_v2.patch` (develop `5edf8394`); it is for local
investigation only. With it, `PIPECHO_PSEUDO=<percent>` hides, for that share of the precursors identified in at least
2 runs, the identification in one of those runs (chosen by a hash of the precursor key). The hidden pairs and the
accepted transfers go to `PIPECHO_PSEUDO_OUT` (`.accepted.tsv`).

There are two modes:
- **present (default):** the hidden feature stays as an unidentified acceptor. A transfer that lands on it is right;
  any other accepted transfer into the pair is a wrong peak.
- **`PIPECHO_PSEUDO_MODE=absent`:** the feature is removed, as if the precursor's peak had not been detected. Every
  accepted transfer into the pair is then a wrong peak.

`run_pseudo.sh` and `run_pseudo_absent.sh` combine the develop multiplex checkpoints of `../seeding` on HYE with 15 %
hidden, i.e. 9,684 pairs. They run once with the random-donor decoys (default) and once with probe decoys
(`PIPECHO_PROBE_DECOYS=120`, v2). Each run uses PIP-ECHO fdr 0.05; the evaluation applies the lower q thresholds to its
transfers. `pseudo_eval.py` evaluates a run; `results/pseudo_eval.txt` has all four. `fdp_estimate.py` extrapolates to
the real transfers (`results/fdp_estimate.txt`).

Caveats:
- The hidden pairs come from precursors identified in at least 2 runs, so they are biased towards abundant precursors.
- The wrong-peak counts are small (14–126).
- The extrapolation assumes that the real pairs without a detected feature of the precursor (status `no_feature` /
  `outside_rt` in `../candidates.py`, 115,519) behave like the absent pseudo pairs.
