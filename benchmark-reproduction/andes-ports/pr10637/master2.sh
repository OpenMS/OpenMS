#!/bin/bash
# Operating curve: extra transfer-q points on the HYE checkpoints, after master.sh has finished.
S=/home/user/bench/scripts; D=/home/user/data; R=/home/user/runs
export OPENMS_DATA_PATH=/home/user/OpenMS/share/OpenMS
log() { echo "$(date -u +%FT%TZ) $*"; }
until grep -q MASTER_DONE /home/user/logs/master.log; do sleep 60; done
arm() { local out=$5; [ -s $out/norm.csv ] && return 0; log "start $out"; $S/plfq.sh "$@" && log "done $out" || log "FAILED $out"; }
P="-targeted_only false -pip_echo true"
for q in 0.02 0.03; do
  t=${q#0.}
  arm /opt/openms-fix $D/db/HYE_td.fasta $D/ids_hye $R/hye/feat_pe $R/hye/off_q$t $P -PipEcho:shape_features false -PipEcho:fdr $q &
  arm /opt/openms-fix $D/db/HYE_td.fasta $D/ids_hye $R/hye/feat_pe $R/hye/on_q$t $P -PipEcho:shape_features true -PipEcho:fdr $q &
done
wait
log "MASTER2_DONE"
