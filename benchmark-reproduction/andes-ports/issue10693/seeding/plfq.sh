#!/bin/bash
# One ProteomicsLFQ arm as in #10585 (scripts/plfq.sh there). Usage:
#   plfq.sh <openms_prefix> <db> <ids_dir> <feat_dir> <out_dir> [extra ProteomicsLFQ args]
# Writes <out_dir>/out.consensusXML, out_msstats.csv (no normalization), then rt_local-normalized
# norm.consensusXML and norm.csv (ConsensusMapNormalizer -algorithm_type rt_local -> MSstatsConverter).
set -euo pipefail
P=$1; DB=$2; IDS=$3; FD=$4; OUT=$5; shift 5
export LD_LIBRARY_PATH=$P/lib OPENMS_DATA_PATH=${OPENMS_DATA_PATH:-/home/user/OpenMS/share/OpenMS}
MZ=/home/user/data/mzml_pb; DESIGN=/home/user/data/design_pb.tsv
RUNS="A_Sample_Alpha_01 A_Sample_Alpha_02 A_Sample_Alpha_03 B_Sample_Alpha_01 B_Sample_Alpha_02 B_Sample_Alpha_03"
IN=(); ID=()
for r in $RUNS; do IN+=($MZ/LFQ_Orbitrap_DDA_Condition_$r.mzML); ID+=($IDS/LFQ_Orbitrap_DDA_Condition_$r.idXML); done
mkdir -p $OUT
$P/bin/ProteomicsLFQ -in "${IN[@]}" -ids "${ID[@]}" -design $DESIGN -fasta $DB -feat_dir $FD \
  -psmFDR 0.01 -proteinFDR ${PROTEIN_FDR:-0.01} -picked_proteinFDR true -protein_inference aggregation \
  -out_msstats $OUT/out_msstats.csv -out_cxml $OUT/out.consensusXML -threads ${THREADS:-4} -no_progress "$@" > $OUT/plfq.log 2>&1
for a in "$@"; do [ "$a" = "-detect_only" ] && exit 0; done
$P/bin/ConsensusMapNormalizer -in $OUT/out.consensusXML -out $OUT/norm.consensusXML -algorithm_type rt_local -no_progress > $OUT/norm.log 2>&1
$P/bin/MSstatsConverter -in $OUT/norm.consensusXML -in_design $DESIGN -out $OUT/norm.csv -no_progress >> $OUT/norm.log 2>&1
# Keep disk lean: hashes (for the option-off vs baseline identity check), transfer counts, then compress.
( cd $OUT && sha256sum out.consensusXML out_msstats.csv norm.csv > sha256.txt )
python3 -I $(dirname "$0")/transfers.py $OUT/out.consensusXML 0.01 0.05 > $OUT/transfers.tsv
rm -f $OUT/norm.consensusXML
gzip -1 -f $OUT/out.consensusXML
