# #10693, lever 3: extraction at the predicted position (prototype, local only)

Develop `5edf8394` + `proto_extraction_v1.patch` (v2 adds the window-scaled decoy shift). This is not for a PR. The
build stamp is the same as develop's, so the develop multiplex checkpoints of `../seeding` are reused. Recipe as in
`../seeding` (HYE plain database, #10585 settings).

## What the prototype does

- **`PLFQ_EXTRA_SEEDS_DIR=<dir>`** (ProteomicsLFQ `detectRun_`): adds the seeds in `<dir>/<run>.featureXML` to the
  multiplex seeds, so FFID extracts features there.
- **Per-seed apex tolerance:** a seed may carry `seed_apex_rt_tolerance` (FFID `addSeeds_` / `filterFeatures_`).
  The prototype uses 20 s, because the predicted RT is within 5 s of the real apex for only 53 % of the found acceptors
  (median 4.6 s, 90 % within 17 s).
- **`PIPECHO_PROBE_DECOYS=<s>`** (PIP-ECHO): the decoy looks at the donor's m/z at a fixed RT shift instead of at a
  random other donor's RT. The sign follows from the precursor key.
  - v1: shift = `<s>` (120 s);
  - v2: shift = max(`<s>`, 4 × local half-window), the separation PIP-ECHO uses for random donors.
- **`make_probe_seeds.py`**: for every reported precursor in every run where it is not identified, one target seed at
  the predicted RT (the PIP-ECHO local-model centre from the `PIPECHO_DIAG` dump, mapped to the run's own RT) and one
  decoy seed at the probe position. Both are at the theoretical m/z. This gives 21k–25k targets per run.

## Arms

All are on HYE. `master_proto.sh` runs the first pass, `master_proto2.sh` and `run_ctrl_v2.sh` the second.

| arm | checkpoints | decoys | `max_training_points` |
|---|---|---|---|
| ctrl | develop multiplex | probes v1 | 50000 (default) |
| ctrl15 | develop multiplex | probes v1 | 15000 |
| ctrl2 | develop multiplex | probes v2 | 15000 |
| protorand | multiplex + extra seeds (v1) | random donor (default) | 50000 |
| proto15 | multiplex + extra seeds (v1) | probes v1 | 15000 |

With probe decoys and the default limit of 50000, the prototype's SVMs trained on 50k points each and were still in round
3 of 10 after 8 h. The fair arms therefore use 15000 for both control and prototype.

## Results

- `results/proteobench_hye.tsv`: all arms, raw (`out_msstats.csv`) and rt_local (`norm.csv`).
- `results/transfer_counts.txt`: candidates, decoy candidates and accepted transfers per arm.
- `results/spots_extraction.tsv`: the `../candidates.py` spot table with the extra seeds. Pairs with a candidate rise from
  18.5 % to 48.7 % (acceptance from the protorand arm).
- `results/seed_features_A02.txt`: what FFID extracted at target vs decoy seeds in A_02.
