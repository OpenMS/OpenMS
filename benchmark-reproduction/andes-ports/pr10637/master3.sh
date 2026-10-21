#!/bin/bash
# shape_features x donor gate (#10633), at -proteinFDR 1 (where ProteomicsLFQ turns the gate on).
# /opt/openms-fix = PR (no gate), /opt/openms-gate = PR + #10633 (combined_10633.patch, same build stamp so the checkpoints are reused).
S=/home/user/bench/scripts; D=/home/user/data; R=/home/user/runs
export OPENMS_DATA_PATH=/home/user/OpenMS/share/OpenMS PROTEIN_FDR=1
log() { echo "$(date -u +%FT%TZ) $*"; }
arm() { local out=$5; [ -s $out/norm.csv ] && return 0; log "start $out"; $S/plfq.sh "$@" && log "done $out" || log "FAILED $out"; }
P="-targeted_only false -pip_echo true"
E="$D/db/hye_entrap_td.fasta $D/ids_entrap $R/entrap/feat_pe"; H="$D/db/HYE_td.fasta $D/ids_hye $R/hye/feat_pe"
arm /opt/openms-fix  $E $R/entrap/p1_nogate_off_q05 $P -PipEcho:shape_features false &
arm /opt/openms-fix  $E $R/entrap/p1_nogate_on_q05  $P -PipEcho:shape_features true &
arm /opt/openms-gate $E $R/entrap/p1_gate_off_q05   $P -PipEcho:shape_features false &
arm /opt/openms-gate $E $R/entrap/p1_gate_on_q05    $P -PipEcho:shape_features true &
wait
for q in 0.01 0.05; do t=${q#0.}
  arm /opt/openms-fix  $H $R/hye/p1_nogate_off_q$t $P -PipEcho:shape_features false -PipEcho:fdr $q &
  arm /opt/openms-fix  $H $R/hye/p1_nogate_on_q$t  $P -PipEcho:shape_features true -PipEcho:fdr $q &
  arm /opt/openms-gate $H $R/hye/p1_gate_off_q$t   $P -PipEcho:shape_features false -PipEcho:fdr $q &
  arm /opt/openms-gate $H $R/hye/p1_gate_on_q$t    $P -PipEcho:shape_features true -PipEcho:fdr $q &
  wait
done
log "MASTER3_DONE"
