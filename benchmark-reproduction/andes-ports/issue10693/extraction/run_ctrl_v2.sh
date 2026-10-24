#!/bin/bash
# control with window-scaled decoy probes (prototype v2), from the develop multiplex checkpoints
S=/home/user/bench/scripts; D=/home/user/data; R=/home/user/runs; P=/opt/openms-proto2
export OPENMS_DATA_PATH=/home/user/OpenMS/share/OpenMS PIPECHO_PROBE_DECOYS=120
MX="-targeted_only false -pip_echo true -PipEcho:max_training_points 15000"
H="$D/db/HYE_td.fasta $D/ids_hye"
( $S/plfq.sh $P $H $R/seed/hye/feat_mx $R/proto/hye/ctrl2_q05 $MX && echo "done ctrl2_q05" ) &
( $S/plfq.sh $P $H $R/seed/hye/feat_mx $R/proto/hye/ctrl2_q01 $MX -PipEcho:fdr 0.01 && echo "done ctrl2_q01" ) &
wait; echo CTRL2_DONE
