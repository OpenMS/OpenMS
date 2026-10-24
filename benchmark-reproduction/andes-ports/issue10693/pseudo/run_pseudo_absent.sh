#!/bin/bash
# #10693: pseudo-transfers with the hidden precursor's feature removed ("not detected"): every accepted transfer is wrong.
S=/home/user/bench/scripts; D=/home/user/data; R=/home/user/runs; P=/opt/openms-pseudo2
export OPENMS_DATA_PATH=/home/user/OpenMS/share/OpenMS PIPECHO_PSEUDO=15 PIPECHO_PSEUDO_MODE=absent
MX="-targeted_only false -pip_echo true"
H="$D/db/HYE_td.fasta $D/ids_hye"
( PIPECHO_PSEUDO_OUT=$R/pseudo/rand_absent.tsv $S/plfq.sh $P $H $R/seed/hye/feat_mx $R/pseudo/rand_absent_q05 $MX && echo "done rand_absent" ) &
( PIPECHO_PSEUDO_OUT=$R/pseudo/probe_absent.tsv PIPECHO_PROBE_DECOYS=120 $S/plfq.sh $P $H $R/seed/hye/feat_mx $R/pseudo/probe_absent_q05 $MX && echo "done probe_absent" ) &
wait; echo PSEUDO_ABSENT_DONE
