# #10693, step 2: extraction evidence as transfer-FDR predictors (prototype, local only)

Develop `5edf8394` + `predictors_v3.patch`. The patch contains `../extraction` v2 (extra seeds, per-seed apex tolerance,
probe decoys), `../pseudo` (pseudo-transfers) and the following. It is not for a PR. The build stamp is the same as
develop's, so the multiplex checkpoints of `../seeding` and the extraction checkpoints of `../extraction` are reused.

## What the patch adds

- **`PIPECHO_EXTRA_PREDICTORS=1`** adds three SVM predictors (`Score`, `RunStatistics::score`, `predictors_t::encode`):
  - `ffid_quality`: the acceptor's overall quality from FeatureFinderIdentification;
  - `width_ratio`: |log2(acceptor / donor elution width of the monoisotopic trace)|;
  - `n_traces`: the acceptor's number of mass traces.

  The apex offset to the predicted RT is not new: `rt_diff_error` / `rt_score` already measure it.
- **`PIPECHO_SCORE_DUMP=<file>`**: every scored acceptor before the FDR cut, with all predictors, PEP and q.
- **`PIPECHO_PROBE_FIXED=1`**: probe decoys with the fixed v1 shift. The extraction checkpoints carry decoy seeds at
  ±120 s, so the extraction arms need v1. The pseudo-transfer runs use v2, as in `../pseudo`.

`width_ratio` and `n_traces` are constant 0. ProteomicsLFQ drops the convex hulls after detection
(`f.setConvexHulls({})` in `detectRun_`), so the checkpoints carry neither widths nor traces, and SimpleSVM ignores
both as uninformative. **Only `ffid_quality` is effective.** Widths or shapes need re-detection, as with the #10637
snapshots.

## Arms

`master_pred.sh` runs all arms on HYE with the plain database and the #10585 settings. The baselines are the arms of
`../extraction` and `../pseudo` without the extra predictors.

| arm | checkpoints | decoys | `max_training_points` | baseline |
|---|---|---|---|---|
| ctrl_x_q05 / q01 | develop multiplex | probes v1 (fixed 120 s) | 15000 | ctrl15 |
| proto_x_q05 / q01 | multiplex + extra seeds v1 | probes v1 | 15000 | proto15 |
| rand_x / rand_x_absent | develop multiplex, 15 % pseudo | random donor (default) | 50000 | rand / rand_absent |
| probe_x / probe_x_absent | develop multiplex, 15 % pseudo | probes v2 | 50000 | probe / probe_absent |

Each arm takes 0.75–2.2 h on one core and about 2 GB.

## Results

- `results/proteobench_hye.tsv`: raw (`out_msstats.csv`) and rt_local (`norm.csv`); min 3 runs, proteobench 0.18.6.
- `results/pseudo_eval.txt`: the `../pseudo/pseudo_eval.py` tables.
- `results/transfer_counts.txt`: accepted transfers per arm.
- `results/score_dump_summary.txt`: predictor medians of targets and decoys, and the lowest target q.
