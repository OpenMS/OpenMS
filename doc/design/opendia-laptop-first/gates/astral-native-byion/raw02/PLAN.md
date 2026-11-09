# Raw 02 native byion replication (written before the run; raw 02 is a design run, not a gate)

Binary: OpenDIACNN at 8c49a5bbb2. Raw 02 inputs (calibration, seeds, windows, search features) read-only from
../astral-entrapment-raw02-20261108. Two arms with identical CNN flags (as the raw 02 stage-2 check):
A = -stage2_feature_set search+boundary, B = search+boundary+byion.
Check: crossfit columns 1-6,9 of A and B are byte-identical to each other.
Reading (as raw 01): PASS if B's entrapment-matched N at FDP_c <= 1% is >= 1.05 x A's AND bulk-null weighted
(bulknull_native.py method) q<=0.05 FDP_c <= 0.06 and q<=0.01 FDP_c <= 0.015.
