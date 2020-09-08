#!/bin/bash
./bin/ProteomicsLFQ -in ../share/OpenMS/examples/FRACTIONS/BSA1_F1.mzML ../share/OpenMS/examples/FRACTIONS/BSA1_F2.mzML ../share/OpenMS/examples/FRACTIONS/BSA2_F1.mzML \
../share/OpenMS/examples/FRACTIONS/BSA2_F2.mzML ../share/OpenMS/examples/FRACTIONS/BSA3_F1.mzML ../share/OpenMS/examples/FRACTIONS/BSA3_F2.mzML \
-ids ../share/OpenMS/examples/FRACTIONS/BSA1_F1.idXML ../share/OpenMS/examples/FRACTIONS/BSA1_F2.idXML ../share/OpenMS/examples/FRACTIONS/BSA2_F1.idXML ../share/OpenMS/examples/FRACTIONS/BSA2_F2.idXML ../share/OpenMS/examples/FRACTIONS/BSA3_F1.idXML ../share/OpenMS/examples/FRACTIONS/BSA3_F2.idXML \
-design ../share/OpenMS/examples/FRACTIONS/BSA_design.tsv \
-Alignment:max_rt_shift 0 \
-fasta ../share/OpenMS/examples/TOPPAS/data/BSA_Identification/18Protein_SoCe_Tr_detergents_trace_target_decoy.fasta \
-targeted_only true \
-transfer_ids "false" \
-mass_recalibration "false" \
-out_cxml BSA.consensusXML \
-out_msstats BSA.csv \
-out BSA.mzTab -threads 4 

./bin/FileInfo -in BSA.consensusXML -out BSA.fileinfo

