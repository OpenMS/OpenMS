#!bin/bash
export PATH="$PATH:/home/sachsenb/OpenMS/openms-build/bin"
export PATH="$PATH:/home/sachsenb/OpenMS/THIRDPARTY/All/MSGFPlus"
export PATH="$PATH:/home/sachsenb/OpenMS/THIRDPARTY/Linux/64bit/Percolator"

#input data
export DATA_PATH="/home/sachsenb/OpenMS/openms-build/yasset_iPRG2015"

# perform search, score calibration, and PSM-level q-value/PEP estimation
for f in ${DATA_PATH}/*.mzML; do
  echo $f
  # search
  MSGFPlusAdapter -in ${f} -out ${f}.idXML -database ${DATA_PATH}/iPRG2015_decoy.fasta -max_precursor_charge 5
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

exit

./bin/ProteomicsLFQ -in \
./iPRG2015/JD_06232014_sample1_A.mzML  ./iPRG2015/JD_06232014_sample1_B.mzML  ./iPRG2015/JD_06232014_sample1_C.mzML \
./iPRG2015/JD_06232014_sample2_A.mzML  ./iPRG2015/JD_06232014_sample2_B.mzML  ./iPRG2015/JD_06232014_sample2_C.mzML \
./iPRG2015/JD_06232014_sample3_A.mzML  ./iPRG2015/JD_06232014_sample3_B.mzML  ./iPRG2015/JD_06232014_sample3_C.mzML \
./iPRG2015/JD_06232014_sample4_A.mzML  ./iPRG2015/JD_06232014_sample4_B.mzML  ./iPRG2015/JD_06232014_sample4_C.mzML \
-ids \
./iPRG2015/JD_06232014_sample1_A.idXML  ./iPRG2015/JD_06232014_sample1_B.idXML  ./iPRG2015/JD_06232014_sample1_C.idXML \
./iPRG2015/JD_06232014_sample2_A.idXML  ./iPRG2015/JD_06232014_sample2_B.idXML  ./iPRG2015/JD_06232014_sample2_C.idXML \
./iPRG2015/JD_06232014_sample3_A.idXML  ./iPRG2015/JD_06232014_sample3_B.idXML  ./iPRG2015/JD_06232014_sample3_C.idXML \
./iPRG2015/JD_06232014_sample4_A.idXML  ./iPRG2015/JD_06232014_sample4_B.idXML  ./iPRG2015/JD_06232014_sample4_C.idXML \
-design ./iPRG2015/experimental_design.tsv \
-Alignment:max_rt_shift 0.1 \
-fasta ./iPRG2015/iPRG2015_decoy.fasta -targeted_only "true" \
-transfer_ids "false" \
-mass_recalibration "false" \
-out_msstats "iPRG2015_targeted_only.csv" \
-out_cxml "iPRG2015_targeted_only.consensusXML" \
-out iPRG2015_targeted_only.mzTab > iPRG2015_targeted_only.log 

./bin/FileInfo -in iPRG2015_targeted_only.consensusXML > iPRG2015_targeted_only.fileinfo

exit

./bin/ProteomicsLFQ -in \
./iPRG2015/JD_06232014_sample1_A.mzML  ./iPRG2015/JD_06232014_sample1_B.mzML  ./iPRG2015/JD_06232014_sample1_C.mzML \
./iPRG2015/JD_06232014_sample2_A.mzML  ./iPRG2015/JD_06232014_sample2_B.mzML  ./iPRG2015/JD_06232014_sample2_C.mzML \
./iPRG2015/JD_06232014_sample3_A.mzML  ./iPRG2015/JD_06232014_sample3_B.mzML  ./iPRG2015/JD_06232014_sample3_C.mzML \
./iPRG2015/JD_06232014_sample4_A.mzML  ./iPRG2015/JD_06232014_sample4_B.mzML  ./iPRG2015/JD_06232014_sample4_C.mzML \
-ids \
./iPRG2015/JD_06232014_sample1_A.idXML  ./iPRG2015/JD_06232014_sample1_B.idXML  ./iPRG2015/JD_06232014_sample1_C.idXML \
./iPRG2015/JD_06232014_sample2_A.idXML  ./iPRG2015/JD_06232014_sample2_B.idXML  ./iPRG2015/JD_06232014_sample2_C.idXML \
./iPRG2015/JD_06232014_sample3_A.idXML  ./iPRG2015/JD_06232014_sample3_B.idXML  ./iPRG2015/JD_06232014_sample3_C.idXML \
./iPRG2015/JD_06232014_sample4_A.idXML  ./iPRG2015/JD_06232014_sample4_B.idXML  ./iPRG2015/JD_06232014_sample4_C.idXML \
-design ./iPRG2015/experimental_design.tsv \
-Alignment:max_rt_shift 0.1 \
-fasta ./iPRG2015/iPRG2015_decoy.fasta -targeted_only "true" \
-transfer_ids "merged" \
-out_msstats "iPRG2015_targeted_only_merged.csv" \
-out_cxml "iPRG2015_targeted_only_merged.consensusXML" \
-out iPRG2015_targeted_only_merged.mzTab -threads 10 > iPRG2015_targeted_only_merged.log 

./bin/FileInfo -in iPRG2015_targeted_only_merged.consensusXML > iPRG2015_targeted_only_merged.fileinfo

./bin/ProteomicsLFQ -in \
./iPRG2015/JD_06232014_sample1_A.mzML  ./iPRG2015/JD_06232014_sample1_B.mzML  ./iPRG2015/JD_06232014_sample1_C.mzML \
./iPRG2015/JD_06232014_sample2_A.mzML  ./iPRG2015/JD_06232014_sample2_B.mzML  ./iPRG2015/JD_06232014_sample2_C.mzML \
./iPRG2015/JD_06232014_sample3_A.mzML  ./iPRG2015/JD_06232014_sample3_B.mzML  ./iPRG2015/JD_06232014_sample3_C.mzML \
./iPRG2015/JD_06232014_sample4_A.mzML  ./iPRG2015/JD_06232014_sample4_B.mzML  ./iPRG2015/JD_06232014_sample4_C.mzML \
-ids \
./iPRG2015/JD_06232014_sample1_A.idXML  ./iPRG2015/JD_06232014_sample1_B.idXML  ./iPRG2015/JD_06232014_sample1_C.idXML \
./iPRG2015/JD_06232014_sample2_A.idXML  ./iPRG2015/JD_06232014_sample2_B.idXML  ./iPRG2015/JD_06232014_sample2_C.idXML \
./iPRG2015/JD_06232014_sample3_A.idXML  ./iPRG2015/JD_06232014_sample3_B.idXML  ./iPRG2015/JD_06232014_sample3_C.idXML \
./iPRG2015/JD_06232014_sample4_A.idXML  ./iPRG2015/JD_06232014_sample4_B.idXML  ./iPRG2015/JD_06232014_sample4_C.idXML \
-design ./iPRG2015/experimental_design.tsv \
-Alignment:max_rt_shift 0.1 \
-fasta ./iPRG2015/iPRG2015_decoy.fasta -targeted_only "true" \
-transfer_ids "SVM" \
-out_msstats "iPRG2015_targeted_only_SVM.csv" \
-out_cxml "iPRG2015_targeted_only_SVM.consensusXML" \
-out iPRG2015_targeted_only_SVM.mzTab -threads 10 > iPRG2015_targeted_only_SVM.log 

./bin/FileInfo -in iPRG2015_targeted_only_SVM.consensusXML > iPRG2015_targeted_only_SVM.fileinfo
exit 1

./bin/ProteomicsLFQ -in \
./iPRG2015/JD_06232014_sample1_A.mzML  ./iPRG2015/JD_06232014_sample1_B.mzML  ./iPRG2015/JD_06232014_sample1_C.mzML \
./iPRG2015/JD_06232014_sample2_A.mzML  ./iPRG2015/JD_06232014_sample2_B.mzML  ./iPRG2015/JD_06232014_sample2_C.mzML \
./iPRG2015/JD_06232014_sample3_A.mzML  ./iPRG2015/JD_06232014_sample3_B.mzML  ./iPRG2015/JD_06232014_sample3_C.mzML \
./iPRG2015/JD_06232014_sample4_A.mzML  ./iPRG2015/JD_06232014_sample4_B.mzML  ./iPRG2015/JD_06232014_sample4_C.mzML \
-ids \
./iPRG2015/JD_06232014_sample1_A.idXML  ./iPRG2015/JD_06232014_sample1_B.idXML  ./iPRG2015/JD_06232014_sample1_C.idXML \
./iPRG2015/JD_06232014_sample2_A.idXML  ./iPRG2015/JD_06232014_sample2_B.idXML  ./iPRG2015/JD_06232014_sample2_C.idXML \
./iPRG2015/JD_06232014_sample3_A.idXML  ./iPRG2015/JD_06232014_sample3_B.idXML  ./iPRG2015/JD_06232014_sample3_C.idXML \
./iPRG2015/JD_06232014_sample4_A.idXML  ./iPRG2015/JD_06232014_sample4_B.idXML  ./iPRG2015/JD_06232014_sample4_C.idXML \
-design ./iPRG2015/experimental_design.tsv \
-Alignment:max_rt_shift 0 \
-fasta ./iPRG2015/iPRG2015_decoy.fasta \
-out iPRG2015_2.mzTab -debug 667 -threads 10 > iPRG2015.log 

./bin/FileInfo -in debug_consensus.consensusXML > iPRG2015.fileinfo

./bin/ProteomicsLFQ -in \
./iPRG2015/JD_06232014_sample1_A.mzML  ./iPRG2015/JD_06232014_sample1_B.mzML  ./iPRG2015/JD_06232014_sample1_C.mzML \
./iPRG2015/JD_06232014_sample2_A.mzML  ./iPRG2015/JD_06232014_sample2_B.mzML  ./iPRG2015/JD_06232014_sample2_C.mzML \
./iPRG2015/JD_06232014_sample3_A.mzML  ./iPRG2015/JD_06232014_sample3_B.mzML  ./iPRG2015/JD_06232014_sample3_C.mzML \
./iPRG2015/JD_06232014_sample4_A.mzML  ./iPRG2015/JD_06232014_sample4_B.mzML  ./iPRG2015/JD_06232014_sample4_C.mzML \
-ids \
./iPRG2015/JD_06232014_sample1_A.idXML  ./iPRG2015/JD_06232014_sample1_B.idXML  ./iPRG2015/JD_06232014_sample1_C.idXML \
./iPRG2015/JD_06232014_sample2_A.idXML  ./iPRG2015/JD_06232014_sample2_B.idXML  ./iPRG2015/JD_06232014_sample2_C.idXML \
./iPRG2015/JD_06232014_sample3_A.idXML  ./iPRG2015/JD_06232014_sample3_B.idXML  ./iPRG2015/JD_06232014_sample3_C.idXML \
./iPRG2015/JD_06232014_sample4_A.idXML  ./iPRG2015/JD_06232014_sample4_B.idXML  ./iPRG2015/JD_06232014_sample4_C.idXML \
-design ./iPRG2015/experimental_design.tsv \
-Alignment:min_run_occur 9 \
-Alignment:max_rt_shift 0 \
-fasta ./iPRG2015/iPRG2015_decoy.fasta \
-out iPRG2015_mrc9.mzTab -debug 667 -threads 10 > iPRG2015mrc9.log 

./bin/FileInfo -in debug_consensus.consensusXML > iPRG2015mrc9.fileinfo



#IDFilter -score:pep 0.85 -in ./iPRG2015/PeptideAtlas/idXML/1A.idXML -out ./iPRG2015/PeptideAtlas/idXML/1A_IDF.idXML
#IDFilter -score:pep 0.85 -in ./iPRG2015/PeptideAtlas/idXML/2A.idXML -out ./iPRG2015/PeptideAtlas/idXML/2A_IDF.idXML
#IDFilter -score:pep 0.85 -in ./iPRG2015/PeptideAtlas/idXML/3A.idXML -out ./iPRG2015/PeptideAtlas/idXML/3A_IDF.idXML
#IDFilter -score:pep 0.85 -in ./iPRG2015/PeptideAtlas/idXML/4A.idXML -out ./iPRG2015/PeptideAtlas/idXML/4A_IDF.idXML
#IDFilter -score:pep 0.85 -in ./iPRG2015/PeptideAtlas/idXML/1B.idXML -out ./iPRG2015/PeptideAtlas/idXML/1B_IDF.idXML
#IDFilter -score:pep 0.85 -in ./iPRG2015/PeptideAtlas/idXML/2B.idXML -out ./iPRG2015/PeptideAtlas/idXML/2B_IDF.idXML
#IDFilter -score:pep 0.85 -in ./iPRG2015/PeptideAtlas/idXML/3B.idXML -out ./iPRG2015/PeptideAtlas/idXML/3B_IDF.idXML
#IDFilter -score:pep 0.85 -in ./iPRG2015/PeptideAtlas/idXML/4B.idXML -out ./iPRG2015/PeptideAtlas/idXML/4B_IDF.idXML
#IDFilter -score:pep 0.85 -in ./iPRG2015/PeptideAtlas/idXML/1C.idXML -out ./iPRG2015/PeptideAtlas/idXML/1C_IDF.idXML
#IDFilter -score:pep 0.85 -in ./iPRG2015/PeptideAtlas/idXML/2C.idXML -out ./iPRG2015/PeptideAtlas/idXML/2C_IDF.idXML
#IDFilter -score:pep 0.85 -in ./iPRG2015/PeptideAtlas/idXML/3C.idXML -out ./iPRG2015/PeptideAtlas/idXML/3C_IDF.idXML
#IDFilter -score:pep 0.85 -in ./iPRG2015/PeptideAtlas/idXML/4C.idXML -out ./iPRG2015/PeptideAtlas/idXML/4C_IDF.idXML
#
./bin/ProteomicsLFQ -in \
./iPRG2015/JD_06232014_sample1_A.mzML  ./iPRG2015/JD_06232014_sample1_B.mzML  ./iPRG2015/JD_06232014_sample1_C.mzML \
./iPRG2015/JD_06232014_sample2_A.mzML  ./iPRG2015/JD_06232014_sample2_B.mzML  ./iPRG2015/JD_06232014_sample2_C.mzML \
./iPRG2015/JD_06232014_sample3_A.mzML  ./iPRG2015/JD_06232014_sample3_B.mzML  ./iPRG2015/JD_06232014_sample3_C.mzML \
./iPRG2015/JD_06232014_sample4_A.mzML  ./iPRG2015/JD_06232014_sample4_B.mzML  ./iPRG2015/JD_06232014_sample4_C.mzML \
-ids \
./iPRG2015/PeptideAtlas/idXML/1A_IDF.idXML  ./iPRG2015/PeptideAtlas/idXML/1B_IDF.idXML  ./iPRG2015/PeptideAtlas/idXML/1C_IDF.idXML \
./iPRG2015/PeptideAtlas/idXML/2A_IDF.idXML  ./iPRG2015/PeptideAtlas/idXML/2B_IDF.idXML  ./iPRG2015/PeptideAtlas/idXML/2C_IDF.idXML \
./iPRG2015/PeptideAtlas/idXML/3A_IDF.idXML  ./iPRG2015/PeptideAtlas/idXML/3B_IDF.idXML  ./iPRG2015/PeptideAtlas/idXML/3C_IDF.idXML \
./iPRG2015/PeptideAtlas/idXML/4A_IDF.idXML  ./iPRG2015/PeptideAtlas/idXML/4B_IDF.idXML  ./iPRG2015/PeptideAtlas/idXML/4C_IDF.idXML \
-design ./iPRG2015/experimental_design.tsv \
-Alignment:max_rt_shift 0 \
-fasta ./iPRG2015/iPRG2015_decoy.fasta \
-out iPRG2015_PeptideAtlas.mzTab -debug 667 -threads 10 > iPRG2015_PeptideAtlas.log 

./bin/FileInfo -in debug_consensus.consensusXML > iPRG2015PeptideAtlas.fileinfo


