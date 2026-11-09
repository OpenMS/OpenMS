# k11: commonness-matched decoys, native OpenDIA on raws 01 and 02 (declared before any run)

## Why
CALIB sections 8 and 11: no pipeline-side input carries the 5% excess, but fragment commonness is lower in decoys than in
targets (library means 0.668 vs 0.792). Re-weighting the existing decoys to the target commonness marginal
(`../astral-raw02-byion-design-20261109/com_gate.txt`) lowers plain q<=0.05 FDP_c from 0.065/0.074 to 0.050/0.062. A library
whose decoys match target commonness by construction is the deployable form (no weights, no labels, no run data).

## Library (label-free; uses the library only)
- Targets: the 300,222 target precursors of `astral-entrapment-20261107/lib/library.pqp`, rows unchanged and in the same order.
- Five candidate decoys per target: the current pseudo-reverse decoy (switchKR) and four `OpenSwathDecoyGenerator -method shuffle
  -shuffle_seed 1..4` (switchKR on) from `targets.pqp`.
- Per candidate: commonness = mean over its top-6 fragments of log1p(number of top-6 fragments of other target precursors in its
  2 Th window within 17 ppm), excluding its own target; own-target overlap = number of its top-6 fragments within 34 ppm of its
  target's top-6.
- Eligible: own-target overlap <= 1. Choose the eligible candidate with the smallest |commonness(decoy) - commonness(target)|;
  ties go to the pseudo-reverse. No eligible candidate: keep the pseudo-reverse.
- Name the chosen decoy `DECOY_<target TRAML_ID>`; convert with TargetedFileConverter.

## Runs (one variable: the decoys)
Per raw, the same binaries and options as the byion arm (`stage2-byion-design/run.sh`, `astral-raw02-byion-design-20261109/run.sh`):
reuse cal.json and anchors.tsv (target IDs checked identical by TRAML_ID), re-run OpenDIA search, prepare_trace_windows.py and
OpenDIACNN crossfit with `-stage2_feature_set search+boundary+byion`. Sequential, one job at a time (memory).

## Reading (tool stage2_score, the base evaluate.py metric; entrapment labels joined only after the runs)
Base: raw 01 N@FDP_c1% 1029, q<=0.01 FDP_c 0.013, q<=0.05 FDP_c 0.065; raw 02 1298, 0.007, 0.074.
- PASS if on BOTH raws plain pooled q<=0.05 FDP_c <= 0.06 AND q<=0.01 FDP_c <= 0.015 AND N@FDP_c1% >= 95% of base.
- PARTIAL if q<=0.05 FDP_c drops by >= 0.01 on both raws but another condition fails.
- FAIL otherwise. Also reported: library commonness means and AUC (target vs chosen decoy), candidate choice shares.
