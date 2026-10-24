#!/bin/bash
# #10693: pseudo-transfer ground truth for the transfer FDR, random-donor decoys (default) vs probe decoys (v2).
S=/home/user/bench/scripts; D=/home/user/data; R=/home/user/runs; P=/opt/openms-pseudo
export OPENMS_DATA_PATH=/home/user/OpenMS/share/OpenMS PIPECHO_PSEUDO=15
MX="-targeted_only false -pip_echo true"
H="$D/db/HYE_td.fasta $D/ids_hye"; mkdir -p $R/pseudo
( PIPECHO_PSEUDO_OUT=$R/pseudo/rand.tsv $S/plfq.sh $P $H $R/seed/hye/feat_mx $R/pseudo/rand_q05 $MX && echo "done rand" ) &
( PIPECHO_PSEUDO_OUT=$R/pseudo/probe.tsv PIPECHO_PROBE_DECOYS=120 $S/plfq.sh $P $H $R/seed/hye/feat_mx $R/pseudo/probe_q05 $MX && echo "done probe" ) &
wait; echo PSEUDO_DONE
