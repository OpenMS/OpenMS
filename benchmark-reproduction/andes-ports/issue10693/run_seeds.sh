#!/bin/bash
export OPENMS_DATA_PATH=/home/user/OpenMS/share/OpenMS LD_LIBRARY_PATH=/opt/openms-fix/lib
D=/home/user/data; r=LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01
cd /home/user/runs/diag/seeds
/opt/openms-fix/bin/ProteomicsLFQ -in $D/mzml_pb/$r.mzML -ids $D/ids_hye/$r.idXML -design $D/design_pb.tsv -fasta $D/db/HYE_td.fasta \
  -feat_dir /home/user/runs/diag/seeds/ckpt -detect_only -targeted_only false -pip_echo true \
  -psmFDR 0.01 -proteinFDR 0.01 -picked_proteinFDR true -protein_inference aggregation -threads 4 -debug 667 -no_progress > detect.log 2>&1 && echo SEED_DONE || echo SEED_FAILED
