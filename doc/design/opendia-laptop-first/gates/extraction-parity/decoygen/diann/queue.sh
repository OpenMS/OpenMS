#!/bin/bash
# Queue (rerun without --report-decoys): five DIA-NN runs, two at a time (9 threads each).
cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/diann-decoys-20261114
run() { rm -rf "$1/raw$3"; bash run.sh "$@"; echo "JCODE_CHECKPOINT {\"message\":\"done $1 raw$3 exit $(cat $1/raw$3/EXIT)\"}"; }
( run pairlib pairlibdec 01 9; run armD armDdec 02 9; run own owndec 02 9 ) &
( run armD armDdec 01 9; run pairlib pairlibdec 02 9 ) &
wait
