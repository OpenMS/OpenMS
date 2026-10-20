#!/bin/bash
# Per-run feature detection (-detect_only) into lfq/<db>/feat_pe, one run at a time; logs peak memory.
# usage: detect_all.sh <build> <db>
set -u
source /home/user/bench/scripts/env.sh
for r in $RUNS; do
  s=$(date +%s)
  DETECT_RUN=$r /home/user/bench/scripts/lfq_arm.sh $1 $2 detect || { echo "FAILED $r $?"; exit 1; }
  echo "$2 $r detected $(( $(date +%s) - s )) s" >> $BENCH/lfq/timing.txt
done
echo "DETECT_DONE $2"
