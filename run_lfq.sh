#!bin/bash
export PATH="$PATH:/home/sachsenb/OpenMS/openms-build/bin"
export PATH="$PATH:/home/sachsenb/OpenMS/THIRDPARTY/Linux/64bit/Percolator"

export MSGF_PATH="/home/sachsenb/OpenMS/THIRDPARTY/All/MSGFPlus/MSGFPlus.jar"
#input data
export DATA_PATH="/home/sachsenb/OpenMS/openms-build/yasset_iPRG2015"

# perform search, score calibration, and PSM-level q-value/PEP estimation
for f in ${DATA_PATH}/*.mzML; do
  echo $f
  # search
  MSGFPlusAdapter -in ${f} -out ${f}.idXML -database ${DATA_PATH}/iPRG2015_decoy.fasta -executable ${MSGF_PATH} -max_precursor_charge 5 -threads 10
  # annotate target/decoy and protein links
  PeptideIndexer -fasta ${DATA_PATH}/iPRG2015_decoy.fasta  -in ${f}.idXML -out ${f}.idXML
  # run percolator so we get well calibrated PEPs and q-values
  PercolatorAdapter -in ${f}.idXML -out ${f}.idXML -percolator_executable percolator -post-processing-tdc -train-best-positive -subset-max-train 100000 
  FalseDiscoveryRate -in ${f}.idXML -out ${f}.idXML -algorithm:add_decoy_peptides -algorithm:add_decoy_proteins 
  # pre-filter to 5% PSM-level FDR to reduce data
  IDFilter -in ${f}.idXML -out ${f}.idXML -score:pep 0.05 
done

###########################
ProteomicsLFQ -in \
${DATA_PATH}/JD_06232014_sample1_A.mzML  ${DATA_PATH}/JD_06232014_sample1_B.mzML  ${DATA_PATH}/JD_06232014_sample1_C.mzML \
${DATA_PATH}/JD_06232014_sample2_A.mzML  ${DATA_PATH}/JD_06232014_sample2_B.mzML  ${DATA_PATH}/JD_06232014_sample2_C.mzML \
${DATA_PATH}/JD_06232014_sample3_A.mzML  ${DATA_PATH}/JD_06232014_sample3_B.mzML  ${DATA_PATH}/JD_06232014_sample3_C.mzML \
${DATA_PATH}/JD_06232014_sample4_A.mzML  ${DATA_PATH}/JD_06232014_sample4_B.mzML  ${DATA_PATH}/JD_06232014_sample4_C.mzML \
-ids \
${DATA_PATH}/JD_06232014_sample1_A.idXML  ${DATA_PATH}/JD_06232014_sample1_B.idXML  ${DATA_PATH}/JD_06232014_sample1_C.idXML \
${DATA_PATH}/JD_06232014_sample2_A.idXML  ${DATA_PATH}/JD_06232014_sample2_B.idXML  ${DATA_PATH}/JD_06232014_sample2_C.idXML \
${DATA_PATH}/JD_06232014_sample3_A.idXML  ${DATA_PATH}/JD_06232014_sample3_B.idXML  ${DATA_PATH}/JD_06232014_sample3_C.idXML \
${DATA_PATH}/JD_06232014_sample4_A.idXML  ${DATA_PATH}/JD_06232014_sample4_B.idXML  ${DATA_PATH}/JD_06232014_sample4_C.idXML \
-design ${DATA_PATH}/experimental_design.tsv \
-Alignment:max_rt_shift 0.1 \
-fasta ${DATA_PATH}/iPRG2015_decoy.fasta -targeted_only "true" \
-transfer_ids "false" \
-mass_recalibration "true" \
-out_msstats "iPRG2015_targeted_only.csv" \
-out_cxml "iPRG2015_targeted_only.consensusXML" \
-out iPRG2015_targeted_only.mzTab  \
> iPRG2015_targeted_only.log 

FileInfo -in iPRG2015_targeted_only.consensusXML > iPRG2015_targeted_only.fileinfo

Rscript ./home/sachsenb/OpenMS/share/OpenMS/SCRIPTS/ProteomicsLFQ.R iPRG2015_targeted_only.csv iPRG2015_targeted_only.mzTab iPRG2015_targeted_only_final.mzTab
