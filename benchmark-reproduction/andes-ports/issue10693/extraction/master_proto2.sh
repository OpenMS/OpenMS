#!/bin/bash
# #10693 lever 3 prototype, second pass: fair arms with PipEcho:max_training_points 15000 for control and prototype
# (with the default 50000 the prototype's SVMs train on 50k points each and need about a day).
S=/home/user/bench/scripts; D=/home/user/data; R=/home/user/runs; P=/opt/openms-proto; SHIFT=120
export OPENMS_DATA_PATH=/home/user/OpenMS/share/OpenMS
log() { echo "$(date -u +%FT%TZ) $*"; }
arm() { local out=$5; [ -s $out/norm.csv ] && return 0; log "start $out"; $S/plfq.sh "$@" && log "done $out" || log "FAILED $out"; }
det() { local out=$5; grep -q "ProteomicsLFQ took" $out/plfq.log 2>/dev/null && return 0; log "start $out"; $S/plfq.sh "$@" && log "done $out" || log "FAILED $out"; }
MX="-targeted_only false -pip_echo true"; C15="-PipEcho:max_training_points 15000"
H="$D/db/HYE_td.fasta $D/ids_hye"; E="$D/db/hye_entrap_td.fasta $D/ids_entrap"
PIPECHO_PROBE_DECOYS=$SHIFT arm $P $H $R/seed/hye/feat_mx $R/proto/hye/ctrl15_q05 $MX $C15 &
PIPECHO_PROBE_DECOYS=$SHIFT arm $P $H $R/seed/hye/feat_mx $R/proto/hye/ctrl15_q01 $MX $C15 -PipEcho:fdr 0.01 &
PIPECHO_PROBE_DECOYS=$SHIFT arm $P $H $R/proto/hye/feat_proto $R/proto/hye/proto15_q05 $MX $C15 &
PIPECHO_PROBE_DECOYS=$SHIFT arm $P $H $R/proto/hye/feat_proto $R/proto/hye/proto15_q01 $MX $C15 -PipEcho:fdr 0.01 &
wait
log "HYE done"
grep -q "ProteomicsLFQ took" $R/proto/entrap/det_proto/plfq.log 2>/dev/null || { PLFQ_EXTRA_SEEDS_DIR=$R/proto/seeds_entrap THREADS=1 det $P $E $R/proto/entrap/feat_proto $R/proto/entrap/det_proto $MX -detect_only; }
grep -q "ProteomicsLFQ took" $R/proto/entrap/det_proto/plfq.log || { log "STOP: entrapment detection incomplete"; exit 1; }
PIPECHO_PROBE_DECOYS=$SHIFT arm $P $E $R/seed/entrap/feat_mx $R/proto/entrap/ctrl15_q05 $MX $C15 &
PIPECHO_PROBE_DECOYS=$SHIFT arm $P $E $R/proto/entrap/feat_proto $R/proto/entrap/proto15_q05 $MX $C15 &
arm $P $E $R/proto/entrap/feat_proto $R/proto/entrap/protorand15_q05 $MX $C15 &
wait
log "MASTER_PROTO2_DONE"
