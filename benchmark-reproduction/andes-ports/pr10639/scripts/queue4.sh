cd /home/user/bench
for e in v2_frozen v2_multi3 v2_phos2_E v2_mhc v2_whole; do python3 -I scripts/bench.py experiments.json $e 3; done
echo QUEUE4_DONE
