#
# README:
#
# Test names (the first argument of 'add_test') must be distinct.
# If not, the first defintion of the test is overwritten!
#
# Therefore we use the following naming scheme for tests:
# TOPP_<tool>_<num>[_out<num2>]
#
# 'TOPP_' This prefix is used to distinguish TOPP and OpenMS tests in the nightly builds.
# <tool>  is the TOPP tool name used to identify failed TOPP tool tests during the nightly builds.
# <num>   Index that is incresed for each call of a TOPP tool during a test.
#         For each new TOPP tool test this is reset to 1.
# '_out'  When the result file of a test is checked with FuzzyDiff, this suffix is appended.
# <num2>  If several files need to be checked in the same test, the checks are numbered from 1.
#
#
# Warning: When passing parameters avoid usage of "<" and ">" signs, as this will confuse CDASH when 
#          interpreting the test results file (which happens to be XML)


### General settings
set(TOPP_BIN_PATH ${OPENMS_BINARY_DIR})
set(DATA_DIR_TOPP ${PROJECT_SOURCE_DIR}/TOPP)
set(DATA_DIR_TOPP_BIN ${PROJECT_BINARY_DIR}/TOPP)
set(DIFF ${TOPP_BIN_PATH}/FuzzyDiff -ini ${DATA_DIR_TOPP}/FuzzyDiff.ini)

### AdditiveSeries tests
add_test("TOPP_AdditiveSeries_1" ${TOPP_BIN_PATH}/AdditiveSeries -in ${DATA_DIR_TOPP}/AdditiveSeries_1_feat11.featureXML ${DATA_DIR_TOPP}/AdditiveSeries_1_feat20.featureXML ${DATA_DIR_TOPP}/AdditiveSeries_1_feat22.featureXML ${DATA_DIR_TOPP}/AdditiveSeries_1_feat34.featureXML -ini ${DATA_DIR_TOPP}/AdditiveSeries_1_parameters.ini)
add_test("TOPP_AdditiveSeries_1_out1" ${DIFF} -in1 AdditiveSeries_1_output.tmp -in2 ${DATA_DIR_TOPP}/AdditiveSeries_1_output.txt)
add_test("TOPP_AdditiveSeries_1_out2" ${DIFF} -in1 AdditiveSeries_1_gnuplot_tmp.cmd -in2 ${DATA_DIR_TOPP}/AdditiveSeries_1_gnuplot.cmd)
add_test("TOPP_AdditiveSeries_1_out3" ${DIFF} -in1 AdditiveSeries_1_gnuplot_tmp.dat -in2 ${DATA_DIR_TOPP}/AdditiveSeries_1_gnuplot.dat)
add_test("TOPP_AdditiveSeries_1_out4" ${DIFF} -in1 AdditiveSeries_1_gnuplot_tmp.err -in2 ${DATA_DIR_TOPP}/AdditiveSeries_1_gnuplot.err)

### IDMerger tests
add_test("TOPP_IDMerger_1" ${TOPP_BIN_PATH}/IDMerger -ini ${DATA_DIR_TOPP}/IDMerger_1_parameters.ini -out IDMerger_1_output.tmp -in ${DATA_DIR_TOPP}/IDMerger_1_input.idXML ${DATA_DIR_TOPP}/IDMerger_2_input.idXML)
add_test("TOPP_IDMerger_1_out1" ${DIFF} -in1 IDMerger_1_output.tmp -in2 ${DATA_DIR_TOPP}/IDMerger_1_output.idXML)

### BaselineFilter tests
add_test("TOPP_BaselineFilter_1" ${TOPP_BIN_PATH}/BaselineFilter -test -in ${DATA_DIR_TOPP}/BaselineFilter_input.mzML -out BaselineFilter.tmp -struc_elem_length 1.5)
add_test("TOPP_BaselineFilter_1_out1" ${DIFF} -in1 BaselineFilter.tmp -in2 ${DATA_DIR_TOPP}/BaselineFilter_output.mzML)

### MapNormalizer tests
add_test("TOPP_MapNormalizer_1" ${TOPP_BIN_PATH}/MapNormalizer -test -in ${DATA_DIR_TOPP}/MapNormalizer_input.mzML -out MapNormalizer.tmp)
add_test("TOPP_MapNormalizer_1_out1" ${DIFF} -in1 MapNormalizer.tmp -in2 ${DATA_DIR_TOPP}/MapNormalizer_output.mzML)

### DBImporter tests
if (DB_TEST)
	add_test("TOPP_DBImporter_1" ${TOPP_BIN_PATH}/DBImporter -init -ini ${DATA_DIR_TOPP_BIN}/DBImporter.ini)
	add_test("TOPP_DBImporter_2" ${TOPP_BIN_PATH}/DBImporter -test -in ${DATA_DIR_TOPP}/DBImporter_input.mzML -ini ${DATA_DIR_TOPP_BIN}/DBImporter.ini)
endif()

### DBExporter tests
if (DB_TEST)
	add_test("TOPP_DBExporter_1" ${TOPP_BIN_PATH}/DBExporter -test -out DBExporter.tmp -ini ${DATA_DIR_TOPP_BIN}/DBExporter.ini -id 1)
	add_test("TOPP_DBExporter_1_out1" ${DIFF} -in1 DBExporter.tmp -in2 ${DATA_DIR_TOPP}/DBExporter_output.mzML)
	add_test("TOPP_DBExporter_2" ${TOPP_BIN_PATH}/DBExporter -test -out DBExporter.tmp -ini ${DATA_DIR_TOPP_BIN}/DBExporter.ini -query "SELECT 1")
	add_test("TOPP_DBExporter_2_out1" ${DIFF} -in1 DBExporter.tmp -in2 ${DATA_DIR_TOPP}/DBExporter_output.mzML)
endif()

### Decharger tests
add_test("TOPP_Decharger_1" ${TOPP_BIN_PATH}/Decharger -test -in ${DATA_DIR_TOPP}/Decharger_input.featureXML -out_fm Decharger_output_fm.tmp -out_cm Decharger_output.tmp -outpairs Decharger_p_output.tmp)
##tmp disabled:
#add_test("TOPP_Decharger_1_out1" ${DIFF} -in1 Decharger_output.tmp -in2 ${DATA_DIR_TOPP}/Decharger_output.consensusXML)

### DTAExtractor tests
add_test("TOPP_DTAExtractor_1" ${TOPP_BIN_PATH}/DTAExtractor -in ${DATA_DIR_TOPP}/DTAExtractor_1_input.mzML -out DTAExtractor -rt :61)
add_test("TOPP_DTAExtractor_1_out1" ${DIFF} -in1 DTAExtractor_RT60.dta -in2 ${DATA_DIR_TOPP}/DTAExtractor_1_output.dta)
add_test("TOPP_DTAExtractor_2" ${TOPP_BIN_PATH}/DTAExtractor -in ${DATA_DIR_TOPP}/DTAExtractor_1_input.mzML -out DTAExtractor -level 1)
add_test("TOPP_DTAExtractor_2_out1" ${DIFF} -in1 DTAExtractor_RT60.dta -in2 ${DATA_DIR_TOPP}/DTAExtractor_2_output.dta)
add_test("TOPP_DTAExtractor_3" ${TOPP_BIN_PATH}/DTAExtractor -in ${DATA_DIR_TOPP}/DTAExtractor_1_input.mzML -out DTAExtractor -level 2 -mz :1000)
add_test("TOPP_DTAExtractor_3_out1" ${DIFF} -in1 DTAExtractor_RT140_MZ5.dta -in2 ${DATA_DIR_TOPP}/DTAExtractor_3_output.dta)

### FeatureFinder tests
add_test("TOPP_FeatureFinder_1" ${TOPP_BIN_PATH}/FeatureFinder -test -ini ${DATA_DIR_TOPP}/FeatureFinder_1_parameters.ini -in ${DATA_DIR_TOPP}/FeatureFinder_1_input.mzML -out FeatureFinder.tmp)
add_test("TOPP_FeatureFinder_1_out1" ${DIFF} -in1 FeatureFinder.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinder_1_output.featureXML)
add_test("TOPP_FeatureFinder_2" ${TOPP_BIN_PATH}/FeatureFinder -test -ini ${DATA_DIR_TOPP}/FeatureFinder_2_parameters.ini -in ${DATA_DIR_TOPP}/FeatureFinder_1_input.mzML -out FeatureFinder.tmp)
add_test("TOPP_FeatureFinder_2_out1" ${DIFF} -in1 FeatureFinder.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinder_2_output.featureXML)
add_test("TOPP_FeatureFinder_3" ${TOPP_BIN_PATH}/FeatureFinder -test -ini ${DATA_DIR_TOPP}/FeatureFinder_3_parameters.ini -in ${DATA_DIR_TOPP}/FeatureFinder_1_input.mzML -out FeatureFinder.tmp)
add_test("TOPP_FeatureFinder_3_out1" ${DIFF} -in1 FeatureFinder.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinder_3_output.featureXML)
add_test("TOPP_FeatureFinder_4" ${TOPP_BIN_PATH}/FeatureFinder -test -ini ${DATA_DIR_TOPP}/FeatureFinder_4_parameters.ini -in ${DATA_DIR_TOPP}/FeatureFinder_4_input.mzML -out FeatureFinder.tmp)
add_test("TOPP_FeatureFinder_4_out1" ${DIFF} -in1 FeatureFinder.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinder_4_output.featureXML)
add_test("TOPP_FeatureFinder_5" ${TOPP_BIN_PATH}/FeatureFinder -test -in ${DATA_DIR_TOPP}/FeatureFinder_5_input.mzML -ini ${DATA_DIR_TOPP}/FeatureFinder_5.ini -out FeatureFinder.tmp)
add_test("TOPP_FeatureFinder_5_out1" ${DIFF} -in1 FeatureFinder.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinder_5_output.featureXML)

### FileConverter tests
add_test("TOPP_FileConverter_1" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_1_input.mzData -out FileConverter.tmp -out_type mzML)
add_test("TOPP_FileConverter_1_out1" ${DIFF} -in1 FileConverter.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_1_output.mzML)
add_test("TOPP_FileConverter_2" ${TOPP_BIN_PATH}/FileConverter -test  -in ${DATA_DIR_TOPP}/FileConverter_2_input.txt -in_type DTA2D -out FileConverter.tmp -out_type mzML)
add_test("TOPP_FileConverter_2_out1" ${DIFF} -in1 FileConverter.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_2_output.mzML)
add_test("TOPP_FileConverter_3" ${TOPP_BIN_PATH}/FileConverter -test  -in ${DATA_DIR_TOPP}/FileConverter_3_input.featureXML -out FileConverter.tmp -out_type mzML)
add_test("TOPP_FileConverter_3_out1" ${DIFF} -in1 FileConverter.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_3_output.mzML)
add_test("TOPP_FileConverter_4" ${TOPP_BIN_PATH}/FileConverter -test  -in ${DATA_DIR_TOPP}/FileConverter_4_input.mzXML -out FileConverter.tmp -out_type mzML)
add_test("TOPP_FileConverter_4_out1" ${DIFF} -in1 FileConverter.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_4_output.mzML)
if (USE_ANDIMS)
	add_test("TOPP_FileConverter_5" ${TOPP_BIN_PATH}/FileConverter -test  -in ${DATA_DIR_TOPP}/FileConverter_5_input.cdf -out FileConverter.tmp -out_type mzML)
	add_test("TOPP_FileConverter_5_out1" ${DIFF} -in1 FileConverter.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_5_output.mzML)
endif()
add_test("TOPP_FileConverter_6" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_6_input.mzML -out FileConverter.tmp -out_type mzXML)
add_test("TOPP_FileConverter_6_out1" ${DIFF} -in1 FileConverter.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_6_output.mzXML)
add_test("TOPP_FileConverter_7" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_7_input.consensusXML -out FileConverter.tmp -out_type featureXML)
add_test("TOPP_FileConverter_7_out1" ${DIFF} -in1 FileConverter.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_7_output.featureXML)
add_test("TOPP_FileConverter_8" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_8_input.mzML -out FileConverter.tmp -out_type mzData)
add_test("TOPP_FileConverter_8_out1" ${DIFF} -in1 FileConverter.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_8_output.mzData)
add_test("TOPP_FileConverter_9" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_9_input.consensusXML -out FileConverter.tmp -out_type featureXML)
add_test("TOPP_FileConverter_9_out1" ${DIFF} -in1 FileConverter.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_9_output.featureXML)

### FileFilter tests
add_test("TOPP_FileFilter_1" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out FileFilter.tmp -rt :30 -mz :1000 -int :20000)
add_test("TOPP_FileFilter_1_out1" ${DIFF} -in1 FileFilter.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_1_output.mzML)
add_test("TOPP_FileFilter_2" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out FileFilter.tmp -rt 30: -mz 1000: -int 100:)
add_test("TOPP_FileFilter_2_out1" ${DIFF} -in1 FileFilter.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_2_output.mzML)
add_test("TOPP_FileFilter_3" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out FileFilter.tmp -level 2)
add_test("TOPP_FileFilter_3_out1" ${DIFF} -in1 FileFilter.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_3_output.mzML)
add_test("TOPP_FileFilter_4" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_4_input.mzML -out FileFilter.tmp -remove_zoom)
add_test("TOPP_FileFilter_4_out1" ${DIFF} -in1 FileFilter.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_4_output.mzML)
add_test("TOPP_FileFilter_5" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_5_input.featureXML -out FileFilter.tmp -rt :1000 -mz :480 -int :79000 -charge :3 -q :0.6)
add_test("TOPP_FileFilter_5_out1" ${DIFF} -in1 FileFilter.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_5_out.featureXML)
add_test("TOPP_FileFilter_6" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_5_input.featureXML -out FileFilter.tmp -rt 1000: -mz 440: -int 70000: -charge 3: -q 0.51:)
add_test("TOPP_FileFilter_6_out1" ${DIFF} -in1 FileFilter.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_6_out.featureXML)
add_test("TOPP_FileFilter_7" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_7_input.mzML -out FileFilter.tmp -int 7000: -level 1 2 3)
add_test("TOPP_FileFilter_7_out1" ${DIFF} -in1 FileFilter.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_7_output.mzML)
add_test("TOPP_FileFilter_8" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_8_input.consensusXML -out FileFilter.tmp -rt 600:1400 -mz 700:2300 -int 1100:6000)
add_test("TOPP_FileFilter_8_out1" ${DIFF} -in1 FileFilter.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_8_out.consensusXML)
add_test("TOPP_FileFilter_9" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_9_input.mzML -out FileFilter.tmp -remove_mode SelectedIonMonitoring)
add_test("TOPP_FileFilter_9_out1" ${DIFF} -in1 FileFilter.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_9_output.mzML)
add_test("TOPP_FileFilter_10" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_10_input.mzML -out FileFilter.tmp -remove_activation "Collision-induced dissociation")
add_test("TOPP_FileFilter_10_out1" ${DIFF} -in1 FileFilter.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_10_output.mzML)
add_test("TOPP_FileFilter_11" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_11_input.mzML -out FileFilter.tmp -remove_activation "Plasma desorption")
add_test("TOPP_FileFilter_11_out1" ${DIFF} -in1 FileFilter.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_11_output.mzML)
add_test("TOPP_FileFilter_12" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_12_input.mzML -out FileFilter.tmp -remove_chromatograms)
add_test("TOPP_FileFilter_12_out1" ${DIFF} -in1 FileFilter.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_12_output.mzML)

### FileInfo tests
add_test("TOPP_FileInfo_1" ${TOPP_BIN_PATH}/FileInfo -in ${DATA_DIR_TOPP}/FileInfo_1_input.txt -in_type DTA -no_progress -out FileInfo.tmp)
add_test("TOPP_FileInfo_1_out1" ${DIFF} -whitelist "file name" -in1 FileInfo.tmp -in2 ${DATA_DIR_TOPP}/FileInfo_1_output.txt)
add_test("TOPP_FileInfo_2" ${TOPP_BIN_PATH}/FileInfo -in ${DATA_DIR_TOPP}/FileInfo_2_input.dta2d -no_progress -out FileInfo.tmp)
add_test("TOPP_FileInfo_2_out1" ${DIFF} -whitelist "file name" -in1 FileInfo.tmp -in2 ${DATA_DIR_TOPP}/FileInfo_2_output.txt)
add_test("TOPP_FileInfo_3" ${TOPP_BIN_PATH}/FileInfo -in ${DATA_DIR_TOPP}/FileInfo_3_input.featureXML -m -s -p -no_progress -out FileInfo.tmp)
add_test("TOPP_FileInfo_3_out1" ${DIFF} -whitelist "file name" -in1 FileInfo.tmp -in2 ${DATA_DIR_TOPP}/FileInfo_3_output.txt)
add_test("TOPP_FileInfo_4" ${TOPP_BIN_PATH}/FileInfo -in ${DATA_DIR_TOPP}/FileInfo_4_input.mzXML -m -no_progress -out FileInfo.tmp)
add_test("TOPP_FileInfo_4_out1" ${DIFF} -whitelist "file name" -in1 FileInfo.tmp -in2 ${DATA_DIR_TOPP}/FileInfo_4_output.txt)
add_test("TOPP_FileInfo_5" ${TOPP_BIN_PATH}/FileInfo -in ${DATA_DIR_TOPP}/FileInfo_5_input.txt -in_type mzData -m -s -no_progress -out FileInfo.tmp)
add_test("TOPP_FileInfo_5_out1" ${DIFF} -whitelist "file name" -in1 FileInfo.tmp -in2 ${DATA_DIR_TOPP}/FileInfo_5_output.txt)
add_test("TOPP_FileInfo_6" ${TOPP_BIN_PATH}/FileInfo -in ${DATA_DIR_TOPP}/FileInfo_6_input.mzData -d -s -no_progress -out FileInfo.tmp)
add_test("TOPP_FileInfo_6_out1" ${DIFF} -whitelist "file name" -in1 FileInfo.tmp -in2 ${DATA_DIR_TOPP}/FileInfo_6_output.txt)
add_test("TOPP_FileInfo_7" ${TOPP_BIN_PATH}/FileInfo -in ${DATA_DIR_TOPP}/FileInfo_7_input.consensusXML -s -m -p -no_progress -out FileInfo.tmp)
add_test("TOPP_FileInfo_7_out1" ${DIFF} -whitelist "file name" -in1 FileInfo.tmp -in2 ${DATA_DIR_TOPP}/FileInfo_7_output.txt)
if (USE_ANDIMS)
	add_test("TOPP_FileInfo_8" ${TOPP_BIN_PATH}/FileInfo -in ${DATA_DIR_TOPP}/FileInfo_8_input.notype -in_type cdf -m -no_progress -out FileInfo.tmp)
	add_test("TOPP_FileInfo_8_out1" ${DIFF} -whitelist "file name" -in1 FileInfo.tmp -in2 ${DATA_DIR_TOPP}/FileInfo_8_output.txt)
endif()
add_test("TOPP_FileInfo_9" ${TOPP_BIN_PATH}/FileInfo -in ${DATA_DIR_TOPP}/FileInfo_9_input.mzML -m -p -s -no_progress -out FileInfo.tmp)
add_test("TOPP_FileInfo_9_out1" ${DIFF} -whitelist "file name" -in1 FileInfo.tmp -in2 ${DATA_DIR_TOPP}/FileInfo_9_output.txt)

### FileMerger tests
add_test("TOPP_FileMerger_1" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_1_input1.dta2d ${DATA_DIR_TOPP}/FileMerger_1_input2.dta ${DATA_DIR_TOPP}/FileMerger_1_input3.dta2d -out FileMerger.tmp -rt_custom 1 2 3)
add_test("TOPP_FileMerger_1_out1" ${DIFF} -in1 FileMerger.tmp -in2 ${DATA_DIR_TOPP}/FileMerger_1_output.mzML)
add_test("TOPP_FileMerger_2" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_2_input1.txt ${DATA_DIR_TOPP}/FileMerger_2_input2.txt -in_type DTA -out FileMerger.tmp -rt_custom 5 10)
add_test("TOPP_FileMerger_2_out1" ${DIFF} -in1 FileMerger.tmp -in2 ${DATA_DIR_TOPP}/FileMerger_2_output.mzML)
add_test("TOPP_FileMerger_3" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_3_input1.dta ${DATA_DIR_TOPP}/FileMerger_3_input2.dta -out FileMerger.tmp -rt_auto)
add_test("TOPP_FileMerger_3_out1" ${DIFF} -in1 FileMerger.tmp -in2 ${DATA_DIR_TOPP}/FileMerger_3_output.mzML)
add_test("TOPP_FileMerger_4" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_4_input1.dta2d ${DATA_DIR_TOPP}/FileMerger_4_input2.dta2d -out FileMerger.tmp)
add_test("TOPP_FileMerger_4_out1" ${DIFF} -in1 FileMerger.tmp -in2 ${DATA_DIR_TOPP}/FileMerger_4_output.mzML)
add_test("TOPP_FileMerger_5" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_5_input_rt1023.331.dta ${DATA_DIR_TOPP}/FileMerger_5_input_rt2044.334.dta ${DATA_DIR_TOPP}/FileMerger_5_input_rt889.32.dta -rt_filename -user_ms_level -ms_level 2 -out FileMerger.tmp)
add_test("TOPP_FileMerger_5_out1" ${DIFF} -in1 FileMerger.tmp -in2 ${DATA_DIR_TOPP}/FileMerger_5_output.mzML)
add_test("TOPP_FileMerger_6" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_6_input1.mzML ${DATA_DIR_TOPP}/FileMerger_6_input2.mzML -out FileMerger.tmp)
add_test("TOPP_FileMerger_6_out1" ${DIFF} -in1 FileMerger.tmp -in2 ${DATA_DIR_TOPP}/FileMerger_6_output.mzML)

### IDRTCalibration tests
add_test("TOPP_IDRTCalibration_1" ${TOPP_BIN_PATH}/IDRTCalibration -in ${DATA_DIR_TOPP}/IDRTCalibration_1_input.IdXML -out IDRTCalibration_1_output.tmp -calibrant_1_input 10 -calibrant_2_input 90)
add_test("TOPP_IDRTCalibration_1_out1" ${DIFF} -in1 IDRTCalibration_1_output.tmp -in2 ${DATA_DIR_TOPP}/IDRTCalibration_1_output.IdXML)

### IDFileConverter tests
add_test("TOPP_IDFileConverter_1" ${TOPP_BIN_PATH}/IDFileConverter -in ${DATA_DIR_TOPP}/IDFileConverter_1_input_sequestoutdir -mz_file ${DATA_DIR_TOPP}/IDFileConverter_1_input_mzfile.mzML -out IDFileConverter_1_output.tmp -out_type idXML -ignore_proteins_per_peptide )
add_test("TOPP_IDFileConverter_1_out1" ${DIFF} -in1 IDFileConverter_1_output.tmp -in2 ${DATA_DIR_TOPP}/IDFileConverter_1_output.idXML)

### IDFilter tests
add_test("TOPP_IDFilter_1" ${TOPP_BIN_PATH}/IDFilter -in ${DATA_DIR_TOPP}/IDFilter_1_input.idXML -out IDFilter_1_output.tmp -sequences_file ${DATA_DIR_TOPP}/IDFilter_1_input.fas)
add_test("TOPP_IDFilter_1_out1" ${DIFF} -in1 IDFilter_1_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_1_output.idXML)
add_test("TOPP_IDFilter_2" ${TOPP_BIN_PATH}/IDFilter -in ${DATA_DIR_TOPP}/IDFilter_2_input.idXML -out IDFilter_2_output.tmp -pep_fraction 1 -prot_fraction 1)
add_test("TOPP_IDFilter_2_out1" ${DIFF} -in1 IDFilter_2_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_2_output.idXML)
add_test("TOPP_IDFilter_3" ${TOPP_BIN_PATH}/IDFilter -in ${DATA_DIR_TOPP}/IDFilter_3_input.idXML -out IDFilter_3_output.tmp -exclusion_peptides_file ${DATA_DIR_TOPP}/IDFilter_3_2_input.idXML)
add_test("TOPP_IDFilter_3_out1" ${DIFF} -in1 IDFilter_3_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_3_output.idXML)
add_test("TOPP_IDFilter_4" ${TOPP_BIN_PATH}/IDFilter -in ${DATA_DIR_TOPP}/IDFilter_4_input.idXML -out IDFilter_4_output.tmp -rt_filtering -p_value 0.08)
add_test("TOPP_IDFilter_4_out1" ${DIFF} -in1 IDFilter_4_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_4_output.idXML)
add_test("TOPP_IDFilter_5" ${TOPP_BIN_PATH}/IDFilter -in ${DATA_DIR_TOPP}/IDFilter_5_input.idXML -out IDFilter_5_output.tmp -pep_score 32 -prot_score 25)
add_test("TOPP_IDFilter_5_out1" ${DIFF} -in1 IDFilter_5_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_5_output.idXML)
add_test("TOPP_IDFilter_6" ${TOPP_BIN_PATH}/IDFilter -in ${DATA_DIR_TOPP}/IDFilter_6_input.idXML -out IDFilter_6_output.tmp -best_n_peptide_hits 2 -best_n_protein_hits 10)
add_test("TOPP_IDFilter_6_out1" ${DIFF} -in1 IDFilter_6_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_6_output.idXML)
add_test("TOPP_IDFilter_7" ${TOPP_BIN_PATH}/IDFilter -in ${DATA_DIR_TOPP}/IDFilter_7_input.idXML -out IDFilter_7_output.tmp -unique)
add_test("TOPP_IDFilter_7_out1" ${DIFF} -in1 IDFilter_7_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_7_output.idXML)

### MapAligner tests
add_test("TOPP_MapAligner_1" ${TOPP_BIN_PATH}/MapAligner -test -ini ${DATA_DIR_TOPP}/MapAligner_1_parameters.ini -in ${DATA_DIR_TOPP}/MapAligner_1_input1.featureXML ${DATA_DIR_TOPP}/MapAligner_1_input2.featureXML ${DATA_DIR_TOPP}/MapAligner_1_input3.featureXML)
add_test("TOPP_MapAligner_1_out1" ${DIFF} -in1 MapAligner_1_output1.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_1_output1.featureXML)
add_test("TOPP_MapAligner_1_out2" ${DIFF} -in1 MapAligner_1_output2.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_1_output2.featureXML)
add_test("TOPP_MapAligner_1_out3" ${DIFF} -in1 MapAligner_1_output3.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_1_output3.featureXML)
add_test("TOPP_MapAligner_1_out4" ${DIFF} -in1 MapAligner_1_trafo1.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_1_trafo1.trafoXML)
add_test("TOPP_MapAligner_1_out5" ${DIFF} -in1 MapAligner_1_trafo2.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_1_trafo2.trafoXML)
add_test("TOPP_MapAligner_1_out6" ${DIFF} -in1 MapAligner_1_trafo3.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_1_trafo3.trafoXML)
add_test("TOPP_MapAligner_2" ${TOPP_BIN_PATH}/MapAligner -test -ini ${DATA_DIR_TOPP}/MapAligner_2_parameters.ini -in ${DATA_DIR_TOPP}/MapAligner_2_input1.mzML ${DATA_DIR_TOPP}/MapAligner_2_input2.mzML ${DATA_DIR_TOPP}/MapAligner_2_input3.mzML)
add_test("TOPP_MapAligner_2_out1" ${DIFF} -in1 MapAligner_2_output1.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_2_output1.mzML)
add_test("TOPP_MapAligner_2_out2" ${DIFF} -in1 MapAligner_2_output2.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_2_output2.mzML)
add_test("TOPP_MapAligner_2_out3" ${DIFF} -in1 MapAligner_2_output3.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_2_output3.mzML)
add_test("TOPP_MapAligner_3" ${TOPP_BIN_PATH}/MapAligner -test -ini ${DATA_DIR_TOPP}/MapAligner_3_parameters.ini -in ${DATA_DIR_TOPP}/MapAligner_3_input1.mzML ${DATA_DIR_TOPP}/MapAligner_3_input2.mzML ${DATA_DIR_TOPP}/MapAligner_3_input3.mzML)
add_test("TOPP_MapAligner_3_out1" ${DIFF} -in1 MapAligner_3_output1.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_3_output1.mzML)
add_test("TOPP_MapAligner_3_out2" ${DIFF} -in1 MapAligner_3_output2.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_3_output2.mzML)
add_test("TOPP_MapAligner_3_out3" ${DIFF} -in1 MapAligner_3_output3.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_3_output3.mzML)
add_test("TOPP_MapAligner_4" ${TOPP_BIN_PATH}/MapAligner -test -type apply_given_trafo -in ${DATA_DIR_TOPP}/MapAligner_4_input.featureXML ${DATA_DIR_TOPP}/MapAligner_4_input.featureXML -given_transformations ${DATA_DIR_TOPP}/MapAligner_4_trafo1.trafoXML ${DATA_DIR_TOPP}/MapAligner_4_trafo2.trafoXML -out MapAligner_4_output1.tmp MapAligner_4_output2.tmp)
add_test("TOPP_MapAligner_4_out1" ${DIFF} -in1 MapAligner_4_output1.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_4_output1.featureXML)
add_test("TOPP_MapAligner_4_out2" ${DIFF} -in1 MapAligner_4_output2.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_4_output2.featureXML)
add_test("TOPP_MapAligner_5" ${TOPP_BIN_PATH}/MapAligner -test -type apply_given_trafo -in ${DATA_DIR_TOPP}/MapAligner_5_input.mzML -given_transformations ${DATA_DIR_TOPP}/MapAligner_5_trafo.trafoXML -out MapAligner_5_output.tmp)
add_test("TOPP_MapAligner_5_out1" ${DIFF} -in1 MapAligner_5_output.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_5_output.mzML)
add_test("TOPP_MapAligner_6" ${TOPP_BIN_PATH}/MapAligner -test -ini ${DATA_DIR_TOPP}/MapAligner_6_parameters.ini -in ${DATA_DIR_TOPP}/MapAligner_6_input1.featureXML ${DATA_DIR_TOPP}/MapAligner_6_input2.featureXML -out MapAligner_6_output1.tmp MapAligner_6_output2.tmp -transformations MapAligner_6_trafo1.tmp MapAligner_6_trafo2.tmp)
add_test("TOPP_MapAligner_6_out1" ${DIFF} -in1 MapAligner_6_trafo1.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_6_trafo1.trafoXML)
add_test("TOPP_MapAligner_6_out2" ${DIFF} -in1 MapAligner_6_trafo2.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_6_trafo2.trafoXML)


### MascotAdapter tests
add_test("TOPP_MascotAdapter_1" ${TOPP_BIN_PATH}/MascotAdapter -ini ${DATA_DIR_TOPP}/MascotAdapter_1_parameters.ini -mascot_in -in ${DATA_DIR_TOPP}/MascotAdapter_1_input.mzData)
add_test("TOPP_MascotAdapter_1_out1" ${DIFF} -in1 MascotAdapter_1_output.tmp -in2 ${DATA_DIR_TOPP}/MascotAdapter_1_output.mascot_in)
add_test("TOPP_MascotAdapter_2" ${TOPP_BIN_PATH}/MascotAdapter -in ${DATA_DIR_TOPP}/MascotAdapter_2_input.mascotXML -out MascotAdapter_2_output.tmp -mascot_out)
add_test("TOPP_MascotAdapter_2_out1" ${DIFF} -in1 MascotAdapter_2_output.tmp -in2 ${DATA_DIR_TOPP}/MascotAdapter_2_output.idXML)

### NoiseFilter tests
add_test("TOPP_NoiseFilter_1" ${TOPP_BIN_PATH}/NoiseFilter -test -ini ${DATA_DIR_TOPP}/NoiseFilter_1_parameters.ini -in ${DATA_DIR_TOPP}/NoiseFilter_1_input.mzML -out NoiseFilter.tmp)
add_test("TOPP_NoiseFilter_1_out1" ${DIFF} -in1 NoiseFilter.tmp -in2 ${DATA_DIR_TOPP}/NoiseFilter_1_output.mzML)
add_test("TOPP_NoiseFilter_2" ${TOPP_BIN_PATH}/NoiseFilter -test -ini ${DATA_DIR_TOPP}/NoiseFilter_2_parameters.ini -in ${DATA_DIR_TOPP}/NoiseFilter_2_input.mzML -out NoiseFilter.tmp)
add_test("TOPP_NoiseFilter_2_out1" ${DIFF} -in1 NoiseFilter.tmp -in2 ${DATA_DIR_TOPP}/NoiseFilter_2_output.mzML)

### PeakPicker tests
add_test("TOPP_PeakPicker_1" ${TOPP_BIN_PATH}/PeakPicker  -test -ini ${DATA_DIR_TOPP}/PeakPicker_parameters.ini -in ${DATA_DIR_TOPP}/PeakPicker_input.mzML -out PeakPicker.tmp)
add_test("TOPP_PeakPicker_1_out1" ${DIFF} -in1 PeakPicker.tmp -in2 ${DATA_DIR_TOPP}/PeakPicker_output.mzML)
add_test("TOPP_PeakPicker_2" ${TOPP_BIN_PATH}/PeakPicker  -test -ini ${DATA_DIR_TOPP}/PeakPicker_deconv_parameters.ini -in ${DATA_DIR_TOPP}/PeakPicker_deconv_input.mzML -out PeakPicker.tmp)
add_test("TOPP_PeakPicker_2_out1" ${DIFF} -in1 PeakPicker.tmp -in2 ${DATA_DIR_TOPP}/PeakPicker_deconv_output.mzML)
add_test("TOPP_PeakPicker_3" ${TOPP_BIN_PATH}/PeakPicker  -test -ini ${DATA_DIR_TOPP}/PeakPicker_2Dopt_parameters.ini -in ${DATA_DIR_TOPP}/PeakPicker_2Dopt_input.mzML -out PeakPicker.tmp)
add_test("TOPP_PeakPicker_3_out1" ${DIFF} -in1 PeakPicker.tmp -in2 ${DATA_DIR_TOPP}/PeakPicker_2Dopt_output.mzML)
add_test("TOPP_PeakPicker_4" ${TOPP_BIN_PATH}/PeakPicker  -test -ini ${DATA_DIR_TOPP}/PeakPicker_parameters.ini -in ${DATA_DIR_TOPP}/PeakPicker_input.mzML -out PeakPicker.tmp -threads 2)
add_test("TOPP_PeakPicker_4_out1" ${DIFF} -in1 PeakPicker.tmp -in2 ${DATA_DIR_TOPP}/PeakPicker_output.mzML)

### Resampler tests
add_test("TOPP_Resampler_1" ${TOPP_BIN_PATH}/Resampler -test -in ${DATA_DIR_TOPP}/Resampler_1_input.mzML -out Resampler.mzML -sampling_rate 0.3)
add_test("TOPP_Resampler_1_out1" ${DIFF} -in1 Resampler.mzML -in2 ${DATA_DIR_TOPP}/Resampler_1_output.mzML)
add_test("TOPP_Resampler_2" ${TOPP_BIN_PATH}/Resampler -test -in ${DATA_DIR_TOPP}/Resampler_2_input.mzML -out Resampler2.png -width 20 -height 15)
add_test("TOPP_Resampler_2_out1" ${DIFF} -in1 Resampler2.png -in2 ${DATA_DIR_TOPP}/Resampler_2_output.png)
add_test("TOPP_Resampler_3" ${TOPP_BIN_PATH}/Resampler -test -in ${DATA_DIR_TOPP}/Resampler_1_input.mzML -out Resampler3.png -width 20 -height 15 -precursor_size 3 -precursor_color green -log_intensity)
add_test("TOPP_Resampler_3_out1" ${DIFF} -in1 Resampler3.png -in2 ${DATA_DIR_TOPP}/Resampler_3_output.png)

### RTModel tests
add_test("TOPP_RTModel_1" ${TOPP_BIN_PATH}/RTModel -in ${DATA_DIR_TOPP}/RTModel_1_input.idXML -out RTModel_1_output.tmp -ini ${DATA_DIR_TOPP}/RTModel_1_parameters.ini)
add_test("TOPP_RTModel_1_out1" ${DIFF} -in1 RTModel_1_output.tmp -in2 ${DATA_DIR_TOPP}/RTModel_1_output.model)
add_test("TOPP_RTModel_2" ${TOPP_BIN_PATH}/RTModel -in_positive ${DATA_DIR_TOPP}/RTModel_2_input_positive.idXML -in_negative ${DATA_DIR_TOPP}/RTModel_2_input_negative.idXML -out RTModel_2_output.tmp -ini ${DATA_DIR_TOPP}/RTModel_2_parameters.ini)
add_test("TOPP_RTModel_2_out1" ${DIFF} -in1 RTModel_2_output.tmp -in2 ${DATA_DIR_TOPP}/RTModel_2_output.model)
add_test("TOPP_RTModel_2_out2" ${DIFF} -in1 RTModel_2_output.tmp_samples -in2 ${DATA_DIR_TOPP}/RTModel_2_output.samples)
add_test("TOPP_RTModel_2_out3" ${DIFF} -in1 RTModel_2_output.tmp_additional_parameters -in2 ${DATA_DIR_TOPP}/RTModel_2_output.additional_parameters)
add_test("TOPP_RTModel_3" ${TOPP_BIN_PATH}/RTModel -in ${DATA_DIR_TOPP}/RTModel_3_input.idXML -out RTModel_3_output.tmp -ini ${DATA_DIR_TOPP}/RTModel_3_parameters.ini)
add_test("TOPP_RTModel_3_out1" ${DIFF} -in1 RTModel_3_output.tmp -in2 ${DATA_DIR_TOPP}/RTModel_3_output.model)
add_test("TOPP_RTModel_3_out2" ${DIFF} -in1 RTModel_3_output.tmp_samples -in2 ${DATA_DIR_TOPP}/RTModel_3_output.samples)
add_test("TOPP_RTModel_4" ${TOPP_BIN_PATH}/RTModel -in ${DATA_DIR_TOPP}/RTModel_4_input.txt -out RTModel_4_output.tmp -ini ${DATA_DIR_TOPP}/RTModel_4_parameters.ini -textfile_input)
add_test("TOPP_RTModel_4_out1" ${DIFF} -in1 RTModel_4_output.tmp -in2 ${DATA_DIR_TOPP}/RTModel_4_output.model)
add_test("TOPP_RTModel_4_out2" ${DIFF} -in1 RTModel_4_output.tmp_samples -in2 ${DATA_DIR_TOPP}/RTModel_4_output.samples)

### PTModel tests
add_test("TOPP_PTModel_1" ${TOPP_BIN_PATH}/PTModel -in_positive ${DATA_DIR_TOPP}/PTModel_1_input_positive.idXML -in_negative ${DATA_DIR_TOPP}/PTModel_1_input_negative.idXML -out PTModel_1_output.tmp -ini ${DATA_DIR_TOPP}/PTModel_1_parameters.ini)
add_test("TOPP_PTModel_1_out1" ${DIFF} -whitelist "probA" "probB" -in1 PTModel_1_output.tmp -in2 ${DATA_DIR_TOPP}/PTModel_1_output.model)
add_test("TOPP_PTModel_1_out2" ${DIFF} -in1 PTModel_1_output.tmp_samples -in2 ${DATA_DIR_TOPP}/PTModel_1_output.samples)
add_test("TOPP_PTModel_1_out3" ${DIFF} -in1 PTModel_1_output.tmp_additional_parameters -in2 ${DATA_DIR_TOPP}/PTModel_1_output.additional_parameters)

### PTPredict tests
add_test("TOPP_PTPredict_1" ${TOPP_BIN_PATH}/PTPredict -in ${DATA_DIR_TOPP}/PTPredict_1_input.idXML -out PTPredict_1_output.tmp -svm_model ${DATA_DIR_TOPP}/PTPredict_1_input.model)
add_test("TOPP_PTPredict_1_out1" ${DIFF} -in1 PTPredict_1_output.tmp -in2 ${DATA_DIR_TOPP}/PTPredict_1_output.idXML)

### RTPredict tests
add_test("TOPP_RTPredict_1" ${TOPP_BIN_PATH}/RTPredict -in ${DATA_DIR_TOPP}/RTPredict_1_input.idXML -out RTPredict_1_output.tmp -total_gradient_time 3000 -svm_model ${DATA_DIR_TOPP}/RTPredict_1_input.model)
add_test("TOPP_RTPredict_1_out1" ${DIFF} -in1 RTPredict_1_output.tmp -in2 ${DATA_DIR_TOPP}/RTPredict_1_output.idXML)
add_test("TOPP_RTPredict_2" ${TOPP_BIN_PATH}/RTPredict -in ${DATA_DIR_TOPP}/RTPredict_2_input.idXML -out_positive RTPredict_2_output_positive.tmp -out_negative RTPredict_2_output_negative.tmp -svm_model ${DATA_DIR_TOPP}/RTPredict_2_input.model)
add_test("TOPP_RTPredict_2_out1" ${DIFF} -in1 RTPredict_2_output_positive.tmp -in2 ${DATA_DIR_TOPP}/RTPredict_2_output_positive.idXML)
add_test("TOPP_RTPredict_2_out2" ${DIFF} -in1 RTPredict_2_output_negative.tmp -in2 ${DATA_DIR_TOPP}/RTPredict_2_output_negative.idXML)
add_test("TOPP_RTPredict_3" ${TOPP_BIN_PATH}/RTPredict -in ${DATA_DIR_TOPP}/RTPredict_3_input.idXML -out RTPredict_3_output.tmp -total_gradient_time 1 -svm_model ${DATA_DIR_TOPP}/RTPredict_3_input.model)
add_test("TOPP_RTPredict_3_out1" ${DIFF} -whitelist "?xml-stylesheet" "predicted_RT_p_value" -in1 RTPredict_3_output.tmp -in2 ${DATA_DIR_TOPP}/RTPredict_3_output.idXML)
add_test("TOPP_RTPredict_4" ${TOPP_BIN_PATH}/RTPredict -in ${DATA_DIR_TOPP}/RTPredict_4_input.txt -out RTPredict_4_output.tmp -total_gradient_time 1 -svm_model ${DATA_DIR_TOPP}/RTPredict_4_input.model -textfile_input)
add_test("TOPP_RTPredict_4_out1" ${DIFF} -in1 RTPredict_4_output.tmp -in2 ${DATA_DIR_TOPP}/RTPredict_4_output.txt)
add_test("TOPP_RTPredict_5" ${TOPP_BIN_PATH}/RTPredict -in ${DATA_DIR_TOPP}/RTPredict_5_input.idXML -out RTPredict_5_output.tmp -total_gradient_time 1 -svm_model ${DATA_DIR_TOPP}/RTPredict_5_input.model -textfile_output)
add_test("TOPP_RTPredict_5_out1" ${DIFF} -in1 RTPredict_5_output.tmp -in2 ${DATA_DIR_TOPP}/RTPredict_5_output.txt)

### InspectAdapter tests
add_test("TOPP_InspectAdapter_1" ${TOPP_BIN_PATH}/InspectAdapter -ini ${DATA_DIR_TOPP}/InspectAdapter_1_parameters.ini -trie_dbs ${DATA_DIR_TOPP}/Inspect_FASTAFile_test2.trie -in ${DATA_DIR_TOPP}/Inspect.mzXML -dbs ${DATA_DIR_TOPP}/Inspect_FASTAFile_test.fasta -inspect_in -out InspectAdapter_2_output.tmp)

add_test("TOPP_InspectAdapter_1_out1" ${DIFF} -whitelist "spectra," "db," -in1 InspectAdapter_2_output.tmp -in2 ${DATA_DIR_TOPP}/InspectAdapter_2_output.inspect_in)

add_test("TOPP_InspectAdapter_2" ${TOPP_BIN_PATH}/InspectAdapter -ini ${DATA_DIR_TOPP}/InspectAdapter_1_parameters.ini -trie_dbs ${DATA_DIR_TOPP}/Inspect_FASTAFile_test2.trie -in ${DATA_DIR_TOPP}/Inspect.mzData -dbs ${DATA_DIR_TOPP}/Inspect_FASTAFile_test.fasta -inspect_in -out InspectAdapter_3_output.tmp)

add_test("TOPP_InspectAdapter_2_out1" ${DIFF} -whitelist "spectra," "db," -in1 InspectAdapter_3_output.tmp -in2 ${DATA_DIR_TOPP}/InspectAdapter_2_output.inspect_in)

# TODO the following tests are waiting for better implementations of InspectAdapter and
# associated classes
#add_test("TOPP_InspectAdapter_3" ${TOPP_BIN_PATH}/InspectAdapter -ini ${DATA_DIR_TOPP}/InspectAdapter_1_parameters.ini -trie_dbs ${DATA_DIR_TOPP}/Inspect_FASTAFile_test2.trie -in ${DATA_DIR_TOPP}/InspectAdapter.out -dbs ${DATA_DIR_TOPP}/Inspect_FASTAFile_test.fasta -out InspectAdapter_4_output.tmp -inspect_out)
#add_test("TOPP_InspectAdapter_3_out1" ${DIFF} -whitelist "?xml-stylesheet" "IdentificationRun date" -in1 InspectAdapter_4_output.tmp -in2 ${DATA_DIR_TOPP}/InspectAdapter_4_output.idXML)

### SequestAdapter tests
add_test("TOPP_SequestAdapter_1" ${TOPP_BIN_PATH}/SequestAdapter -ini ${DATA_DIR_TOPP}/SequestAdapter_1_parameters.ini -in ${DATA_DIR_TOPP}/Sequest.mzXML -mz_files ${DATA_DIR_TOPP}/Sequest.mzXML -modifications_xml_file ${DATA_DIR_TOPP}/Sequest_PTMs.xml -out SequestAdapter_2_output.tmp -sequest_in -temp_data_directory ${DATA_DIR_TOPP} -db ${DATA_DIR_TOPP}/Sequest_FASTAFile_test.fasta)
add_test("TOPP_SequestAdapter_1_out1" ${DIFF} -in1 SequestAdapter_2_output.tmp -in2 ${DATA_DIR_TOPP}/SequestAdapter_2_output.sequest_in)
add_test("TOPP_SequestAdapter_2" ${TOPP_BIN_PATH}/SequestAdapter -ini ${DATA_DIR_TOPP}/SequestAdapter_1_parameters.ini -in ${DATA_DIR_TOPP}/Sequest.mzData -mz_files ${DATA_DIR_TOPP}/Sequest.mzXML -modifications_xml_file ${DATA_DIR_TOPP}/Sequest_PTMs.xml -out SequestAdapter_3_output.tmp -sequest_in -temp_data_directory ${DATA_DIR_TOPP} -db ${DATA_DIR_TOPP}/Sequest_FASTAFile_test.fasta)
add_test("TOPP_SequestAdapter_2_out1" ${DIFF} -in1 SequestAdapter_3_output.tmp -in2 ${DATA_DIR_TOPP}/SequestAdapter_2_output.sequest_in)

# TODO the following tests are waiting for better implementations of InspectAdapter and
# associated classes
#add_test("TOPP_SequestAdapter_3" ${TOPP_BIN_PATH}/SequestAdapter -ini ${DATA_DIR_TOPP}/SequestAdapter_2_parameters.ini -mz_files ${DATA_DIR_TOPP}/Sequest.mzXML -modifications_xml_file ${DATA_DIR_TOPP}/Sequest_PTMs.xml -in ${DATA_DIR_TOPP}/Sequest.mzXML -out SequestAdapter_4_output.tmp -temp_data_directory ${DATA_DIR_TOPP} -db ${DATA_DIR_TOPP}/Sequest_FASTAFile_test.fasta)
#add_test("TOPP_SequestAdapter_3_out1" ${DIFF} -in1 SequestAdapter_4_output.tmp -in2 ${DATA_DIR_TOPP}/SequestAdapter_4_output.idXML)

### SpecLibSearcher tests
#add_test("TOPP_SpecLibSearcher_1" ${TOPP_BIN_PATH}/SpecLibSearcher -test -ini ${DATA_DIR_TOPP}/SpecLibSearcher_1_parameters.ini -in ${DATA_DIR_TOPP}/SpecLibSearcher_1.MzData -lib $(DATA_DIR_TOPP)/SpecLibSearcher_1.MSP -out SpecLibSearcher_1.tmp)
#add_test("TOPP_SpecLibSearcher_1_out1" ${DIFF} -in1 SpecLibSearcher_1.tmp  -in2 $(DATA_DIR_TOPP)/SpecLibSearcher_1.idXML)
### PepNovoAdapter tests
add_test("TOPP_PepNovoAdapter_1" ${TOPP_BIN_PATH}/PepNovoAdapter -ini ${DATA_DIR_TOPP}/PepNovoAdapter_1_parameters.ini -in ${DATA_DIR_TOPP}/PepNovo.mzXML -pepnovo_in -out PepNovoAdapter_3_output.tmp -dta_list ${DATA_DIR_TOPP}/tmp/dta_list.txt -model_directory ${DATA_DIR_TOPP}/tmp/ -temp_data_directory ${DATA_DIR_TOPP}/tmp/)
add_test("TOPP_PepNovoAdapter_1_out1" ${DIFF} -in1 ${DATA_DIR_TOPP}/tmp/PepNovo_PTMs_.txt -in2 ${DATA_DIR_TOPP}/tmp/PepNovo_PTMs.txt)
#TODO ANDREAS - We have to clean up the /tmp/ directory to run this test multiple times
#add_test("TOPP_PepNovoAdapter_2" ${TOPP_BIN_PATH}/PepNovoAdapter -ini ${DATA_DIR_TOPP}/PepNovoAdapter_1_parameters.ini -in ${DATA_DIR_TOPP}/PepNovo.mzData -pepnovo_in -out PepNovoAdapter_4_output.tmp -temp_data_directory ${DATA_DIR_TOPP}) 
#add_test("TOPP_PepNovoAdapter_2_out1" ${DIFF} -in1 ${DATA_DIR_TOPP}/PepNovo_PTMs_.txt -in2 ${DATA_DIR_TOPP}/PepNovo_PTMs.txt) 
add_test("TOPP_PepNovoAdapter_3" ${TOPP_BIN_PATH}/PepNovoAdapter -ini ${DATA_DIR_TOPP}/PepNovoAdapter_5_parameters.ini -in ${DATA_DIR_TOPP}/PepNovoAdapter_5_output.pepnovo_out -out PepNovoAdapter_5_output.tmp -pepnovo_out -dta_list ${DATA_DIR_TOPP}/tmp/dta_list.txt -model_directory ${DATA_DIR_TOPP}/tmp/ -temp_data_directory ${DATA_DIR_TOPP}/tmp/ -modifications_xml_file ${DATA_DIR_TOPP}/PepNovo_PTMs.xml -mz_files ${DATA_DIR_TOPP}/PepNovo.mzXML)
add_test("TOPP_PepNovoAdapter_3_out1" ${DIFF} -whitelist "?xml-stylesheet" "date_group_1" -in1 PepNovoAdapter_5_output.tmp -in2 ${DATA_DIR_TOPP}/PepNovoAdapter_5_output.idXML)

### SpectraFilter tests
add_test("TOPP_SpectraFilter_1" ${TOPP_BIN_PATH}/SpectraFilter -test -in ${DATA_DIR_TOPP}/SpectraFilter_1_input.mzML -out SpectraFilter.tmp -type SqrtMower)
add_test("TOPP_SpectraFilter_1_out1" ${DIFF} -in1 SpectraFilter.tmp -in2 ${DATA_DIR_TOPP}/SpectraFilter_1_output.mzML)
add_test("TOPP_SpectraFilter_2" ${TOPP_BIN_PATH}/SpectraFilter -test -in ${DATA_DIR_TOPP}/SpectraFilter_2_input.mzML -out SpectraFilter.tmp -type WindowMower)
add_test("TOPP_SpectraFilter_2_out1" ${DIFF} -in1 SpectraFilter.tmp -in2 ${DATA_DIR_TOPP}/SpectraFilter_2_output.mzML)
add_test("TOPP_SpectraFilter_3" ${TOPP_BIN_PATH}/SpectraFilter -test -in ${DATA_DIR_TOPP}/SpectraFilter_3_input.mzML -out SpectraFilter.tmp -ini ${DATA_DIR_TOPP}/SpectraFilter_3_parameters.ini -type WindowMower)
add_test("TOPP_SpectraFilter_3_out1" ${DIFF} -in1 SpectraFilter.tmp -in2 ${DATA_DIR_TOPP}/SpectraFilter_3_output.mzML)

### InternalCalibration tests
add_test("TOPP_InternalCalibration_1" ${TOPP_BIN_PATH}/InternalCalibration -test -in ${DATA_DIR_TOPP}/InternalCalibration_2_input.mzML -out InternalCalibration.tmp -ref_peaks ${DATA_DIR_TOPP}/InternalCalibration_ref_masses -ini ${DATA_DIR_TOPP}/InternalCalibration_parameters.ini -type spectrumwise)
#add_test("TOPP_InternalCalibration_1_out1" ${DIFF} -in1 InternalCalibration.tmp -in2 ${DATA_DIR_TOPP}/InternalCalibration_2_output.mzML)

### TOFCalibration tests
add_test("TOPP_TOFCalibration_1" ${TOPP_BIN_PATH}/TOFCalibration -test -in ${DATA_DIR_TOPP}/TOFCalibration_1_input.mzML -out TOFCalibration.tmp -ref_masses ${DATA_DIR_TOPP}/TOFCalibration_ref_masses -ini ${DATA_DIR_TOPP}/TOFCalibration_parameters.ini -tof_const ${DATA_DIR_TOPP}/TOFCalibration_const -ext_calibrants ${DATA_DIR_TOPP}/TOFCalibration_1_calibrants.mzML)
#add_test("TOPP_TOFCalibration_1_out1" ${DIFF} -in1 TOFCalibration.tmp -in2 ${DATA_DIR_TOPP}/TOFCalibration_1_output.mzML)
add_test("TOPP_TOFCalibration_2" ${TOPP_BIN_PATH}/TOFCalibration -test -in ${DATA_DIR_TOPP}/TOFCalibration_2_input.mzML -out TOFCalibration.tmp -ref_masses ${DATA_DIR_TOPP}/TOFCalibration_ref_masses -ini ${DATA_DIR_TOPP}/TOFCalibration_parameters.ini -tof_const ${DATA_DIR_TOPP}/TOFCalibration_const -ext_calibrants ${DATA_DIR_TOPP}/TOFCalibration_2_calibrants.mzML -peak_data)
#add_test("TOPP_TOFCalibration_2_out1" ${DIFF} -in1 TOFCalibration.tmp -in2 ${DATA_DIR_TOPP}/TOFCalibration_2_output.mzML)

### TextExporter tests
add_test("TOPP_TextExporter_1" ${TOPP_BIN_PATH}/TextExporter -in ${DATA_DIR_TOPP}/TextExporter_1_input.featureXML -no_progress -out TextExporter.tmp)
add_test("TOPP_TextExporter_1_out1" ${DIFF} -in1 TextExporter.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_1_output.txt)

add_test("TOPP_TextExporter_2" ${TOPP_BIN_PATH}/TextExporter -in ${DATA_DIR_TOPP}/TextExporter_2_input.consensusXML -ini ${DATA_DIR_TOPP}/TextExporter_2_parameters.ini)
add_test("TOPP_TextExporter_2_out1" ${DIFF} -whitelist "extracted from" -in1 TextExporter_2_consensus_centroids.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_2_consensus_centroids.wsv)
add_test("TOPP_TextExporter_2_out2" ${DIFF} -whitelist "extracted from" -in1 TextExporter_2_consensus_elements.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_2_consensus_elements.wsv)
add_test("TOPP_TextExporter_2_out3" ${DIFF} -whitelist "extracted from" -in1 TextExporter_2_consensus_features.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_2_consensus_features.wsv)
add_test("TOPP_TextExporter_2_out4" ${DIFF} -whitelist "extracted from" -in1 TextExporter_2_consensus_tsv.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_2_consensus_tsv.tsv)

add_test("TOPP_TextExporter_3" ${TOPP_BIN_PATH}/TextExporter -in ${DATA_DIR_TOPP}/TextExporter_3_input.idXML -no_progress -out TextExporter.tmp)
add_test("TOPP_TextExporter_3_out1" ${DIFF} -in1 TextExporter.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_3_output.txt)

add_test("TOPP_TextExporter_4" ${TOPP_BIN_PATH}/TextExporter -in ${DATA_DIR_TOPP}/TextExporter_3_input.idXML -no_progress -out TextExporter.tmp -proteins_only)
add_test("TOPP_TextExporter_4_out1" ${DIFF} -in1 TextExporter.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_4_output_proteins.txt)

add_test("TOPP_TextExporter_5" ${TOPP_BIN_PATH}/TextExporter -in ${DATA_DIR_TOPP}/TextExporter_3_input.idXML -no_progress -out TextExporter.tmp -peptides_only -first_dim_rt)
add_test("TOPP_TextExporter_5_out1" ${DIFF} -in1 TextExporter.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_5_output_peptides.txt)

add_test("TOPP_TextExporter_6" ${TOPP_BIN_PATH}/TextExporter -in ${DATA_DIR_TOPP}/TextExporter_6_input.featureXML -no_progress -out TextExporter.tmp -no_ids)
add_test("TOPP_TextExporter_6_out1" ${DIFF} -in1 TextExporter.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_6_output.txt)

add_test("TOPP_TextExporter_7" ${TOPP_BIN_PATH}/TextExporter -in ${DATA_DIR_TOPP}/TextExporter_7_input.consensusXML -ini ${DATA_DIR_TOPP}/TextExporter_7_parameters.ini)
add_test("TOPP_TextExporter_7_out1" ${DIFF} -whitelist "extracted from" -in1 TextExporter_7_consensus_tsv.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_7_consensus_tsv.tsv)

### TextImporter tests
add_test("TOPP_TextImporter_1" ${TOPP_BIN_PATH}/TextImporter -test -in ${DATA_DIR_TOPP}/TextImporter_1_input.csv -no_progress -out TextImporter.tmp)
add_test("TOPP_TextImporter_1_out1" ${DIFF} -in1 TextImporter.tmp -in2 ${DATA_DIR_TOPP}/TextImporter_1_output.featureXML)

add_test("TOPP_TextImporter_2" ${TOPP_BIN_PATH}/TextImporter -test -in ${DATA_DIR_TOPP}/TextImporter_2_input.pepList -no_progress -out TextImporter.tmp -mode SpecArray)
add_test("TOPP_TextImporter_2_out1" ${DIFF} -in1 TextImporter.tmp -in2 ${DATA_DIR_TOPP}/TextImporter_2_output.featureXML)

add_test("TOPP_TextImporter_3" ${TOPP_BIN_PATH}/TextImporter -test -in ${DATA_DIR_TOPP}/TextImporter_3_input.peptides.tsv -no_progress -out TextImporter.tmp -mode msInspect)
add_test("TOPP_TextImporter_3_out1" ${DIFF} -in1 TextImporter.tmp -in2 ${DATA_DIR_TOPP}/TextImporter_3_output.featureXML)

add_test("TOPP_TextImporter_4" ${TOPP_BIN_PATH}/TextImporter -test -in ${DATA_DIR_TOPP}/TextImporter_4_input.peptides.dat -no_progress -out TextImporter.tmp -mode Kroenik)
add_test("TOPP_TextImporter_4_out1" ${DIFF} -in1 TextImporter.tmp -in2 ${DATA_DIR_TOPP}/TextImporter_4_output.featureXML)

### SILACAnalyzer tests
add_test("TOPP_SILACAnalyzer_1" ${TOPP_BIN_PATH}/SILACAnalyzer -test -ini ${DATA_DIR_TOPP}/SILACAnalyzer_double.ini -in ${DATA_DIR_TOPP}/SILACAnalyzer_double.mzML -out SILACAnalyzer.tmp -out_visual SILACAnalyzer_vis.tmp)
add_test("TOPP_SILACAnalyzer_1_out1" ${DIFF} -in1 SILACAnalyzer.tmp -in2 ${DATA_DIR_TOPP}/SILACAnalyzer_double.consensusXML)
add_test("TOPP_SILACAnalyzer_1_out2" ${DIFF} -in1 SILACAnalyzer_vis.tmp -in2 ${DATA_DIR_TOPP}/SILACAnalyzer_double.featureXML)

add_test("TOPP_SILACAnalyzer_2" ${TOPP_BIN_PATH}/SILACAnalyzer -test -ini ${DATA_DIR_TOPP}/SILACAnalyzer_triple.ini -in ${DATA_DIR_TOPP}/SILACAnalyzer_triple.mzML -out SILACAnalyzer.tmp -out_visual SILACAnalyzer_vis.tmp)
add_test("TOPP_SILACAnalyzer_2_out1" ${DIFF} -in1 SILACAnalyzer.tmp -in2 ${DATA_DIR_TOPP}/SILACAnalyzer_triple.consensusXML)
add_test("TOPP_SILACAnalyzer_2_out2" ${DIFF} -in1 SILACAnalyzer_vis.tmp -in2 ${DATA_DIR_TOPP}/SILACAnalyzer_triple.featureXML)

### FeatureLinker tests
add_test("TOPP_FeatureLinker_1" ${TOPP_BIN_PATH}/FeatureLinker -test -ini ${DATA_DIR_TOPP}/FeatureLinker_1_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinker_1_input.featureXML)
add_test("TOPP_FeatureLinker_1_out1" ${DIFF} -in1 FeatureLinker_1_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinker_1_output.consensusXML)
add_test("TOPP_FeatureLinker_2" ${TOPP_BIN_PATH}/FeatureLinker -test -ini ${DATA_DIR_TOPP}/FeatureLinker_2_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinker_2_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinker_2_input2.featureXML ${DATA_DIR_TOPP}/FeatureLinker_2_input3.featureXML)
add_test("TOPP_FeatureLinker_2_out1" ${DIFF} -in1 FeatureLinker_2_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinker_2_output.consensusXML)
add_test("TOPP_FeatureLinker_3" ${TOPP_BIN_PATH}/FeatureLinker -test -ini ${DATA_DIR_TOPP}/FeatureLinker_3_parameters.ini  -in ${DATA_DIR_TOPP}/FeatureLinker_3_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinker_3_input2.featureXML ${DATA_DIR_TOPP}/FeatureLinker_3_input3.featureXML)
add_test("TOPP_FeatureLinker_3_out1" ${DIFF} -in1 FeatureLinker_3_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinker_3_output.consensusXML)
add_test("TOPP_FeatureLinker_4" ${TOPP_BIN_PATH}/FeatureLinker -test -ini ${DATA_DIR_TOPP}/FeatureLinker_4_parameters.ini  -in ${DATA_DIR_TOPP}/FeatureLinker_4_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinker_4_input2.featureXML)
add_test("TOPP_FeatureLinker_4_out1" ${DIFF} -in1 FeatureLinker_4_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinker_4_output.consensusXML)
add_test("TOPP_FeatureLinker_5" ${TOPP_BIN_PATH}/FeatureLinker -test -ini ${DATA_DIR_TOPP}/FeatureLinker_5_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinker_5_input.featureXML)
add_test("TOPP_FeatureLinker_5_out1" ${DIFF} -in1 FeatureLinker_5_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinker_5_output.consensusXML)

### IDMapper tests
add_test("TOPP_IDMapper_1" ${TOPP_BIN_PATH}/IDMapper -test -in ${DATA_DIR_TOPP}/IDMapper_1_input.featureXML -out IDMapper_1_output.tmp -id ${DATA_DIR_TOPP}/IDMapper_1_input.idXML)
add_test("TOPP_IDMapper_1_out1" ${DIFF} -in1 IDMapper_1_output.tmp -in2 ${DATA_DIR_TOPP}/IDMapper_1_output.featureXML)
add_test("TOPP_IDMapper_2" ${TOPP_BIN_PATH}/IDMapper -test -in ${DATA_DIR_TOPP}/IDMapper_2_input.consensusXML -out IDMapper_2_output.tmp -id ${DATA_DIR_TOPP}/IDMapper_2_input.idXML)
add_test("TOPP_IDMapper_2_out1" ${DIFF} -in1 IDMapper_2_output.tmp -in2 ${DATA_DIR_TOPP}/IDMapper_2_output.consensusXML)

### ConsensusID tests
add_test("TOPP_ConsensusID_1" ${TOPP_BIN_PATH}/ConsensusID -ini ${DATA_DIR_TOPP}/ConsensusID_1_parameters.ini -in ${DATA_DIR_TOPP}/ConsensusID_1_input.idXML -out ConsensusID_1_output.tmp)
add_test("TOPP_ConsensusID_1_out1" ${DIFF} -whitelist "?xml-stylesheet" "IdentificationRun date" -in1 ConsensusID_1_output.tmp -in2 ${DATA_DIR_TOPP}/ConsensusID_1_output.idXML)
add_test("TOPP_ConsensusID_2" ${TOPP_BIN_PATH}/ConsensusID -ini ${DATA_DIR_TOPP}/ConsensusID_2_parameters.ini -out ConsensusID_2_output.tmp -in ${DATA_DIR_TOPP}/ConsensusID_2_input.featureXML)
add_test("TOPP_ConsensusID_2_out1" ${DIFF} -whitelist "IdentificationRun id" -in1 ConsensusID_2_output.tmp -in2 ${DATA_DIR_TOPP}/ConsensusID_2_output.featureXML)
add_test("TOPP_ConsensusID_3" ${TOPP_BIN_PATH}/ConsensusID -ini ${DATA_DIR_TOPP}/ConsensusID_3_parameters.ini -out ConsensusID_3_output.tmp -in ${DATA_DIR_TOPP}/ConsensusID_3_input.consensusXML)
add_test("TOPP_ConsensusID_3_out1" ${DIFF} -whitelist "?xml-stylesheet" "IdentificationRun id" -in1 ConsensusID_3_output.tmp -in2 ${DATA_DIR_TOPP}/ConsensusID_3_output.consensusXML)


### PrecursorIonSelector tests
add_test("TOPP_PrecursorIonSelector_1" ${TOPP_BIN_PATH}/PrecursorIonSelector -in ${DATA_DIR_TOPP}/PrecursorIonSelector_features.featureXML -sim_results PrecursorIonSelector_1_output.tmp -ini ${DATA_DIR_TOPP}/PrecursorIonSelector_1_parameters.ini -db_path ${DATA_DIR_TOPP}/PrecursorIonSelector_db.fasta -ids ${DATA_DIR_TOPP}/PrecursorIonSelector_ids.IdXML)
add_test("TOPP_PrecursorIonSelector_1_out1" ${DIFF} -in1 PrecursorIonSelector_1_output.tmp -in2 ${DATA_DIR_TOPP}/PrecursorIonSelector_1_output.txt)
add_test("TOPP_PrecursorIonSelector_2" ${TOPP_BIN_PATH}/PrecursorIonSelector -in ${DATA_DIR_TOPP}/PrecursorIonSelector_features.featureXML -out PrecursorIonSelector_2_output.tmp -ini ${DATA_DIR_TOPP}/PrecursorIonSelector_2_parameters.ini -db_path ${DATA_DIR_TOPP}/PrecursorIonSelector_db.fasta -ids ${DATA_DIR_TOPP}/PrecursorIonSelector_ids.IdXML)
add_test("TOPP_PrecursorIonSelector_2_out1" ${DIFF} -in1 PrecursorIonSelector_2_output.tmp -in2 ${DATA_DIR_TOPP}/PrecursorIonSelector_2_output.featureXML)

### CompNovo tests
add_test("TOPP_CompNovo_1" ${TOPP_BIN_PATH}/CompNovo -in ${DATA_DIR_TOPP}/CompNovo_1_input.mzML -ini ${DATA_DIR_TOPP}/CompNovo_1.ini -out CompNovo_1_output.tmp)
add_test("TOPP_CompNovo_1_out1" ${DIFF} -whitelist "?xml-stylesheet" "IdentificationRun date" -in1 CompNovo_1_output.tmp -in2 ${DATA_DIR_TOPP}/CompNovo_1_output.idXML)
add_test("TOPP_CompNovo_2" ${TOPP_BIN_PATH}/CompNovo -in ${DATA_DIR_TOPP}/CompNovo_2_input.mzML -ini ${DATA_DIR_TOPP}/CompNovo_2.ini -out CompNovo_2_output.tmp)
add_test("TOPP_CompNovo_2_out1" ${DIFF} -whitelist "?xml-stylesheet" "IdentificationRun date" -in1 CompNovo_2_output.tmp -in2 ${DATA_DIR_TOPP}/CompNovo_2_output.idXML)


