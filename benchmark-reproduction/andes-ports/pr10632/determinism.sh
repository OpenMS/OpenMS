#!/bin/bash
# Is PIP-ECHO-mode detection of run A_01 deterministic? Detect it again with each build (-detect_only).
source $(dirname $0)/env.sh
L=$WORK/lfq; r=LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01
for b in base pr; do
  BIN=/opt/openms-$b/bin; export LD_LIBRARY_PATH=/opt/openms-$b/lib
  $BIN/ProteomicsLFQ -in $DATA/$r.mzML -ids $WORK/search/$r.idXML -design $(dirname $0)/design.tsv -fasta $WORK/db/HYE_td.fasta \
    -feat_dir $L/det_A01_$b -detect_only -psmFDR 0.01 -proteinFDR 0.01 -picked_proteinFDR true -protein_inference aggregation \
    -targeted_only false -pip_echo true -threads 4 -no_progress > $L/det_A01_$b.log 2>&1
  echo "$(date +%T) det $b exit=$?"
done
echo DET_DONE
