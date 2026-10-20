#!/bin/bash
# ProteomicsLFQ arm (#10585 recipe): plfq.sh <out_dir> <feat_dir> [extra args]; BIN selects the build
set -e; source $(dirname $0)/env.sh
out=$1; fd=$2; shift 2; mkdir -p $out $fd
IN=""; IDS=""
for r in $RUNS; do IN="$IN $DATA/$r.mzML"; IDS="$IDS $WORK/search/$r.idXML"; done
s=$(date +%s)
$BIN/ProteomicsLFQ -in $IN -ids $IDS -design $(dirname $0)/design.tsv -fasta $WORK/db/HYE_td.fasta -feat_dir $fd \
  -psmFDR 0.01 -proteinFDR 0.01 -picked_proteinFDR true -protein_inference aggregation \
  -out_msstats $out/out_msstats.csv -out_cxml $out/out.consensusXML -out $out/out.mzTab -threads 4 -no_progress "$@"
echo "ProteomicsLFQ wall time: $(( $(date +%s)-s )) s"
# rt_local normalization route documented by #10585
$BIN/ConsensusMapNormalizer -in $out/out.consensusXML -out $out/norm.consensusXML -algorithm_type rt_local -no_progress
$BIN/MSstatsConverter -in $out/norm.consensusXML -in_design $(dirname $0)/design.tsv -out $out/norm_msstats.csv -no_progress
rm $out/norm.consensusXML
echo ARM_DONE
