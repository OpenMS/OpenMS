# Stage-2 byion group: raw 01 native design run (written before the run)

Binary: OpenDIACNN at commit 8c49a5bbb2. Same command as ../stage2-boundary-design/run.sh, with
-stage2_feature_set search+boundary+byion (defaults -stage2_byion_ppm 10, -stage2_byion_intensity_min 300).

Equivalence check: the sample, logit, fold and window columns of crossfit.tsv are byte-identical to
../stage2-boundary-design/boundary, and the first 15 feature columns are identical.

Reading (entrapment-matched N at FDP_c <= 1%, the matched_fdp.py metric; reference search+boundary 757):
- PASS if N >= 795 (+5%) and pooled q <= 0.05 has FDP_c <= 0.06.
- Otherwise FAIL. That would mean the stage-2 LDA cannot use the counts the way the offline MLP does (+15%),
  which points to the head rather than the features.
Raw 02 replication follows only on a raw 01 PASS.
