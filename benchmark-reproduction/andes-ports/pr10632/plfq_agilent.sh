#!/bin/bash
# ProteomicsLFQ arm on the Agilent set: plfq_agilent.sh <out_dir> <feat_dir> [extra args]
set -e; source $(dirname $0)/env.sh
out=$1; fd=$2; shift 2; mkdir -p $out $fd
R2="202505616 202505617 202505618 202505622 202505623 202505624"; IN=""; IDS=""
for r in $R2; do IN="$IN $DATA2/$r.mzML"; IDS="$IDS $WORK2/search/$r.idXML"; done
$BIN/ProteomicsLFQ -in $IN -ids $IDS -design $(dirname $0)/design_agilent.tsv -fasta $DATA2/2026-07-17-decoys-contam-GCF_000695525.1_BOL_protein.fas -feat_dir $fd \
  -psmFDR 0.01 -proteinFDR 0.01 -picked_proteinFDR true -protein_inference aggregation \
  -out_msstats $out/out_msstats.csv -out_cxml $out/out.consensusXML -threads 4 -no_progress "$@"
$BIN/ConsensusMapNormalizer -in $out/out.consensusXML -out $out/norm.consensusXML -algorithm_type rt_local -no_progress
$BIN/MSstatsConverter -in $out/norm.consensusXML -in_design $(dirname $0)/design_agilent.tsv -out $out/norm_msstats.csv -no_progress
rm $out/norm.consensusXML
echo ARM_DONE
