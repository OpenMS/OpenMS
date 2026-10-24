#!/bin/bash
# #10693 lever 3 prototype: extraction at predicted positions with RT-shift decoy probes (local build /opt/openms-proto
# = develop 5edf8394 + proto_extraction.patch, same build stamp as /opt/openms-dev, so the multiplex checkpoints are reused).
S=/home/user/bench/scripts; D=/home/user/data; R=/home/user/runs; P=/opt/openms-proto; SHIFT=120; TOL=20
export OPENMS_DATA_PATH=/home/user/OpenMS/share/OpenMS
log() { echo "$(date -u +%FT%TZ) $*"; }
arm() { local out=$5; [ -s $out/norm.csv ] && return 0; log "start $out"; $S/plfq.sh "$@" && log "done $out" || log "FAILED $out"; }
det() { local out=$5; grep -q "ProteomicsLFQ took" $out/plfq.log 2>/dev/null && return 0; log "start $out"; $S/plfq.sh "$@" && log "done $out" || log "FAILED $out"; }
MX="-targeted_only false -pip_echo true"
H="$D/db/HYE_td.fasta $D/ids_hye"; E="$D/db/hye_entrap_td.fasta $D/ids_entrap"
# a) controls: develop multiplex checkpoints, probe decoys
PIPECHO_PROBE_DECOYS=$SHIFT PIPECHO_DIAG=$R/proto/hye/link_ctrl_q05.tsv arm $P $H $R/seed/hye/feat_mx $R/proto/hye/ctrl_q05 $MX &
PIPECHO_PROBE_DECOYS=$SHIFT arm $P $H $R/seed/hye/feat_mx $R/proto/hye/ctrl_q01 $MX -PipEcho:fdr 0.01 &
PIPECHO_PROBE_DECOYS=$SHIFT arm $P $E $R/seed/entrap/feat_mx $R/proto/entrap/ctrl_q05 $MX &
until grep -q -E "DIAG_DONE|DIAG_FAILED" $R/seed/entrap/run_diag.log; do sleep 60; done
wait
# b) entrapment seeds from its diag dump
if [ ! -d $R/proto/seeds_entrap ]; then
  /home/user/pbenv/bin/python -I $S/candidates.py $R/seed/entrap/link_mx_q05.tsv $R/seed/entrap/mx_q05_diag/out.consensusXML.gz $R/proto/entrap/cand_mx > $R/proto/entrap/cand_mx.log 2>&1
  /home/user/pbenv/bin/python -I $S/make_probe_seeds.py $R/seed/entrap/link_mx_q05.tsv $R/proto/entrap/cand_mx_pairs.tsv.gz $R/seed/entrap/feat_mx $R/proto/seeds_entrap $SHIFT $TOL > $R/proto/entrap/seeds.log 2>&1
fi
# c) HYE: detection with the extra seeds (2 threads: memory), then combines
PLFQ_EXTRA_SEEDS_DIR=$R/proto/seeds_hye THREADS=2 det $P $H $R/proto/hye/feat_proto $R/proto/hye/det_proto $MX -detect_only
PIPECHO_PROBE_DECOYS=$SHIFT PIPECHO_DIAG=$R/proto/hye/link_proto_q05.tsv arm $P $H $R/proto/hye/feat_proto $R/proto/hye/proto_q05 $MX &
PIPECHO_PROBE_DECOYS=$SHIFT arm $P $H $R/proto/hye/feat_proto $R/proto/hye/proto_q01 $MX -PipEcho:fdr 0.01 &
arm $P $H $R/proto/hye/feat_proto $R/proto/hye/protorand_q05 $MX &
arm $P $H $R/proto/hye/feat_proto $R/proto/hye/protorand_q01 $MX -PipEcho:fdr 0.01 &
wait
log "HYE done"
# d) entrapment
PLFQ_EXTRA_SEEDS_DIR=$R/proto/seeds_entrap THREADS=1 det $P $E $R/proto/entrap/feat_proto $R/proto/entrap/det_proto $MX -detect_only
PIPECHO_PROBE_DECOYS=$SHIFT arm $P $E $R/proto/entrap/feat_proto $R/proto/entrap/proto_q05 $MX &
arm $P $E $R/proto/entrap/feat_proto $R/proto/entrap/protorand_q05 $MX &
wait
log "MASTER_PROTO_DONE"
