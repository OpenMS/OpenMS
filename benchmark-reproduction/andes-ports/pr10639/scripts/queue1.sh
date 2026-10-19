cd /home/user/bench
for e in frozen_hfx_A2 frozen_hfx_B1 frozen_astral_B1 frozen_velos_125_R1 multi3 chunk calib frozen_hfx_A2_t1; do python3 -I scripts/bench.py experiments.json $e 3; done
echo QUEUE1_DONE
