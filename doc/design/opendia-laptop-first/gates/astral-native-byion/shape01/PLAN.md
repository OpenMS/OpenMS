# Stage-2 shape group: raw 01 native design run (written before the run)

Binary: OpenDIACNN at f98c0dca2c. Same command as ../stage2-byion-design/run.sh, with
-stage2_feature_set search+boundary+byion+shape.

Check: crossfit columns 1-6,9 and feature columns 1-16 (through byion_y) are byte-identical to ../stage2-byion-design/byion.
This shows the picker refactor left the boundary features unchanged.

Reading (entrapment-matched N at FDP_c 1%; reference byion 1029):
- PASS if N >= 1060 (+3%) and the bulk-null (bulknull_native.py method) q<=0.05 FDP_c <= 0.06.
- On PASS, replicate on raw 02 against 1298 with the same rule.
