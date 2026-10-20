#!/bin/bash
# after a container restart: finish P_base (reuses its written checkpoints), then the Agilent set
B=/home/user/bench; L=/home/user/work/lfq
echo "$(date +%T) start P_base (resume)"
BIN=/opt/openms-base/bin $B/plfq.sh $L/P_base $L/feat_P_base -targeted_only false -pip_echo true > $L/P_base.log 2>&1; echo "$(date +%T) P_base exit=$?"
$B/run_agilent.sh > /home/user/logs/agilent.log 2>&1
echo RESUME_DONE
