#!/bin/bash
export PATH="$PATH:/Users/pfeuffer/git/OpenMS-inference-src/cmake-build-debug/bin"
export PATH="$PATH:/Users/pfeuffer/git/THIRDPARTY/MacOS/64bit/percolator"

export MSGF_PATH="/Users/pfeuffer/git/THIRDPARTY/All/MSGFPlus/MSGFPlus.jar"
#input data
export DATA_PATH="/Users/pfeuffer/git/OpenMS-inference-src/share/OpenMS/examples/FRACTIONS/"

# perform search, score calibration, and PSM-level q-value/PEP estimation

# percolator doesn't work on the small file so we use IDPEP instead
for f in ${DATA_PATH}/*.mzML; do
  echo $f
  fn=${f%.mzML} # filename and path without mzML extension
  # search ( we need multiple matches so we can build a proper decoy model on this small dataset)
  MSGFPlusAdapter -in ${f} -out ${fn}.idXML -database ${DATA_PATH}../TOPPAS/data/BSA_Identification/18Protein_SoCe_Tr_detergents_trace_target_decoy.fasta -executable ${MSGF_PATH} -max_precursor_charge 5 -matches_per_spec 10
  # annotate target/decoy and protein links
  PeptideIndexer -fasta ${DATA_PATH}../TOPPAS/data/BSA_Identification/18Protein_SoCe_Tr_detergents_trace_target_decoy.fasta -in ${fn}.idXML -out ${fn}.idXML -enzyme:specificity none
  IDPosteriorErrorProbability -in ${fn}.idXML -out ${fn}.idXML -out_plot ${fn}_fit.txt
done


# changeed protein FDR filtering threshold...

ProteomicsLFQ -in \
${DATA_PATH}/BSA1.mzML ${DATA_PATH}/BSA2.mzML ${DATA_PATH}/BSA3.mzML \
-ids \
${DATA_PATH}/BSA1.idXML ${DATA_PATH}/BSA2.idXML ${DATA_PATH}/BSA3.idXML \
-design ${DATA_PATH}/experimental_design.tsv \
-Alignment:max_rt_shift 0.1 \
-fasta ${DATA_PATH}/../TOPPAS/data/BSA_Identification/18Protein_SoCe_Tr_detergents_trace_target_decoy.fasta -targeted_only "true" \
-transfer_ids "false" \
-proteinFDR 0.3 \
-mass_recalibration "false" \
-out_msstats "BSA_targeted_only.csv" \
-out_cxml "BSA_targeted_only.consensusXML" \
-out BSA_targeted_only.mzTab -debug 0  \

