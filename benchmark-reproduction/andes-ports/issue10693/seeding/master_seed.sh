#!/bin/bash
# #10693 lever 2: Seeding multiplex vs biosaur2 on develop (5edf8394), #10585 recipe, -proteinFDR 0.01.
# Detection runs alone (memory; biosaur2 with -threads 1: with 4 and 2 threads it was OOM-killed at ~12 GB); combines from the checkpoints run 4 at a time. PIP-ECHO arms write the
# PIPECHO_DIAG dump (/opt/openms-dev-diag = develop + diag_link patch, same build stamp).
S=/home/user/bench/scripts; D=/home/user/data; R=/home/user/runs/seed
export OPENMS_DATA_PATH=/home/user/OpenMS/share/OpenMS
log() { echo "$(date -u +%FT%TZ) $*"; }
arm() { local out=$5; [ -s $out/norm.csv ] && return 0; log "start $out"; $S/plfq.sh "$@" && log "done $out" || log "FAILED $out"; }
det() { local out=$5; grep -q "ProteomicsLFQ took" $out/plfq.log 2>/dev/null && return 0; log "start $out"; $S/plfq.sh "$@" && log "done $out" || log "FAILED $out"; }
P=/opt/openms-dev; PD=/opt/openms-dev-diag
MX="-targeted_only false -pip_echo true"; B2="$MX -Seeding:algorithm biosaur2"
H="$D/db/HYE_td.fasta $D/ids_hye"; E="$D/db/hye_entrap_td.fasta $D/ids_entrap"
det $P $H $R/hye/feat_mx $R/hye/det_mx $MX -detect_only
THREADS=1 det $P $H $R/hye/feat_b2 $R/hye/det_b2 $B2 -detect_only
arm $P $H $R/hye/feat_id $R/hye/id_only -targeted_only true
PIPECHO_DIAG=$R/hye/link_mx_q05.tsv arm $PD $H $R/hye/feat_mx $R/hye/mx_q05 $MX &
PIPECHO_DIAG=$R/hye/link_b2_q05.tsv arm $PD $H $R/hye/feat_b2 $R/hye/b2_q05 $B2 &
arm $P $H $R/hye/feat_mx $R/hye/mx_q01 $MX -PipEcho:fdr 0.01 &
arm $P $H $R/hye/feat_b2 $R/hye/b2_q01 $B2 -PipEcho:fdr 0.01 &
wait
log "HYE done"
det $P $E $R/entrap/feat_mx $R/entrap/det_mx $MX -detect_only
THREADS=1 det $P $E $R/entrap/feat_b2 $R/entrap/det_b2 $B2 -detect_only
arm $P $E $R/entrap/feat_mx $R/entrap/mx_q05 $MX &
arm $P $E $R/entrap/feat_b2 $R/entrap/b2_q05 $B2 &
wait
log "MASTER_SEED_DONE"
