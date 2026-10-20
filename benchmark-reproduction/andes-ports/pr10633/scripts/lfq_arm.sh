#!/bin/bash
# One ProteomicsLFQ arm (MBR with PIP-ECHO) on the six HF-X runs as in #10585, then rt_local normalization and
# MSstats export. All arms of one database share the per-run feature checkpoints in lfq/<db>/feat_pe (-feat_dir),
# so only alignment, linking (PIP-ECHO), inference and quantification run per arm.
# usage: lfq_arm.sh <build: pr|base> <db: entrap|hye> <arm name> [extra ProteomicsLFQ options ...]
# env: MZTAB=1 also writes the mzTab; DETECT_RUN=<run> only detects that run (-detect_only).
set -euo pipefail
source /home/user/bench/scripts/env.sh
build=$1; db=$2; arm=$3; shift 3
use_openms /opt/openms-$build
case $db in entrap) fasta=$BENCH/db/HYE_entrap_td.fasta ;; hye) fasta=$BENCH/db/HYE_td.fasta ;; esac
fd=$BENCH/lfq/$db/feat_pe
in=(); ids=()
for r in ${DETECT_RUN:-$RUNS}; do in+=($DATA_PB/$r.mzML); ids+=($BENCH/ids/$db/$r.idXML); done
common=(-design $BENCH/run/design.tsv -fasta $fasta -feat_dir $fd \
  -psmFDR 0.01 -proteinFDR 0.01 -picked_proteinFDR true -protein_inference aggregation \
  -targeted_only false -pip_echo true -threads 4 -no_progress)
if [ -n "${DETECT_RUN:-}" ]; then
  mkdir -p $BENCH/lfq/$db/detect
  ProteomicsLFQ -in "${in[@]}" -ids "${ids[@]}" "${common[@]}" -detect_only "$@" > $BENCH/lfq/$db/detect/$DETECT_RUN.log 2>&1
  exit 0
fi
out=$BENCH/lfq/$db/$arm; mkdir -p $out
extra_out=()
[ "${MZTAB:-0}" = 1 ] && extra_out=(-out $out/lfq.mzTab)
start=$(date +%s)
ProteomicsLFQ -in "${in[@]}" -ids "${ids[@]}" "${common[@]}" \
  -out_cxml $out/lfq.consensusXML -out_msstats $out/msstats.csv "${extra_out[@]}" "$@" > $out/lfq.log 2>&1
ConsensusMapNormalizer -in $out/lfq.consensusXML -out $out/rtlocal.consensusXML -algorithm_type rt_local -no_progress > $out/norm.log 2>&1
MSstatsConverter -in $out/rtlocal.consensusXML -in_design $BENCH/run/design.tsv -out $out/rtlocal.msstats.csv -no_progress > $out/msstats.log 2>&1
echo "$build $db $arm $* : $(( $(date +%s) - start )) s" >> $BENCH/lfq/timing.txt
