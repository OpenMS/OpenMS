#!/bin/bash
# Full benchmark chain (runs unattended; every step skips what is already done).
S=/home/user/bench/scripts; D=/home/user/data; R=/home/user/runs
export OPENMS_DATA_PATH=/home/user/OpenMS/share/OpenMS
log() { echo "$(date -u +%FT%TZ) $*"; }
until grep -q SEARCH_ALL_DONE /home/user/logs/search_hye.log; do sleep 30; done
log "HYE search done"
mkdir -p $D/ids_entrap
$S/search_all.sh /opt/openms-fix/bin $D/db/hye_entrap_td.fasta $D/ids_entrap
log "entrapment search done"
arm() { # arm <prefix> <db> <ids> <featdir> <outdir> args...
  local out=$5; [ -s $out/norm.csv ] && return 0
  log "start $out"; $S/plfq.sh "$@" && log "done $out" || log "FAILED $out"
}
# HYE: detection alone (11-13 GB peak), then the combine arms from the same checkpoints
arm /opt/openms-fix $D/db/HYE_td.fasta $D/ids_hye $R/hye/feat_pe $R/hye/detect -targeted_only false -pip_echo true -detect_only
log "HYE detection done"
P="-targeted_only false -pip_echo true"
arm /opt/openms-base $D/db/HYE_td.fasta $D/ids_hye $R/hye/feat_pe $R/hye/base_q05 $P &
arm /opt/openms-fix  $D/db/HYE_td.fasta $D/ids_hye $R/hye/feat_pe $R/hye/off_q05 $P -PipEcho:shape_features false &
arm /opt/openms-fix  $D/db/HYE_td.fasta $D/ids_hye $R/hye/feat_pe $R/hye/on_q05 $P -PipEcho:shape_features true &
arm /opt/openms-fix  $D/db/HYE_td.fasta $D/ids_hye $R/hye/feat_pe $R/hye/off_q01 $P -PipEcho:shape_features false -PipEcho:fdr 0.01 &
wait
arm /opt/openms-fix  $D/db/HYE_td.fasta $D/ids_hye $R/hye/feat_pe $R/hye/on_q01 $P -PipEcho:shape_features true -PipEcho:fdr 0.01 &
arm /opt/openms-base $D/db/HYE_td.fasta $D/ids_hye $R/hye/feat_pe $R/hye/base_q01 $P -PipEcho:fdr 0.01 &
wait
log "HYE arms done"
# Entrapment database
arm /opt/openms-fix $D/db/hye_entrap_td.fasta $D/ids_entrap $R/entrap/feat_pe $R/entrap/detect -targeted_only false -pip_echo true -detect_only
log "entrapment detection done"
arm /opt/openms-fix $D/db/hye_entrap_td.fasta $D/ids_entrap $R/entrap/feat_pe $R/entrap/off_q05 $P -PipEcho:shape_features false &
arm /opt/openms-fix $D/db/hye_entrap_td.fasta $D/ids_entrap $R/entrap/feat_pe $R/entrap/on_q05 $P -PipEcho:shape_features true &
arm /opt/openms-fix $D/db/hye_entrap_td.fasta $D/ids_entrap $R/entrap/feat_pe $R/entrap/off_q01 $P -PipEcho:shape_features false -PipEcho:fdr 0.01 &
arm /opt/openms-fix $D/db/hye_entrap_td.fasta $D/ids_entrap $R/entrap/feat_pe $R/entrap/on_q01 $P -PipEcho:shape_features true -PipEcho:fdr 0.01 &
wait
log "MASTER_DONE"
