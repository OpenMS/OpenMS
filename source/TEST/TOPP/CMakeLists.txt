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

### General settings
set(TOPP_BIN_PATH ${PROJECT_SOURCE_DIR}/bin/)
set(DIFF ${TOPP_BIN_PATH}/FuzzyDiff -ini FuzzyDiff.ini)

### AdditiveSeries tests
add_test("TOPP_AdditiveSeries_1" ${TOPP_BIN_PATH}/AdditiveSeries -ini AdditiveSeries_1_parameters.ini)
add_test("TOPP_AdditiveSeries_1_out1" ${DIFF} -in1 AdditiveSeries_1_output.tmp -in2 AdditiveSeries_1_output.txt)
add_test("TOPP_AdditiveSeries_1_out2" ${DIFF} -in1 AdditiveSeries_1_gnuplot_tmp.cmd -in2 AdditiveSeries_1_gnuplot.cmd)
add_test("TOPP_AdditiveSeries_1_out3" ${DIFF} -in1 AdditiveSeries_1_gnuplot_tmp.dat -in2 AdditiveSeries_1_gnuplot.dat)
add_test("TOPP_AdditiveSeries_1_out4" ${DIFF} -in1 AdditiveSeries_1_gnuplot_tmp.err -in2 AdditiveSeries_1_gnuplot.err)

### IDMerger tests
add_test("TOPP_IDMerger_1" ${TOPP_BIN_PATH}/IDMerger -ini IDMerger_1_parameters.ini -out IDMerger_1_output.tmp)
add_test("TOPP_IDMerger_1_out1" ${DIFF} -in1 IDMerger_1_output.tmp -in2 IDMerger_1_output.idXML)

### BaselineFilter tests
add_test("TOPP_BaselineFilter_1" ${TOPP_BIN_PATH}/BaselineFilter  -in RawData_input.mzData -out BaselineFilter_output.tmp -struc_elem_length 1.5)
add_test("TOPP_BaselineFilter_1_out1" ${DIFF} -in1 BaselineFilter_output.tmp -in2 BaselineFilter_output.mzData)

### MapNormalizer tests
add_test("TOPP_MapNormalizer_1" ${TOPP_BIN_PATH}/MapNormalizer  -in MapNormalizer_input.mzData -out MapNormalizer_output.tmp)
add_test("TOPP_MapNormalizer_1_out1" ${DIFF} -in1 MapNormalizer_output.tmp -in2 MapNormalizer_output.mzData)

### DBImporter tests
if (EXISTS source/TEST/TOPP/DBImporter_input.ini)
	add_test("TOPP_DBImporter_1" ${TOPP_BIN_PATH}/DBImporter  -init -ini DBImporter_input.ini)
	add_test("TOPP_DBImporter_2" ${TOPP_BIN_PATH}/DBImporter  -in DBImporter_input.mzData -ini DBImporter_input.ini)
	add_test("TOPP_DBImporter_3" ${TOPP_BIN_PATH}/DBImporter  -in DBImporter_input2.mzData -ini DBImporter_input.ini)
endif()

### DBExporter tests
if (EXISTS source/TEST/TOPP/DBExporter_input.ini)
	add_test("TOPP_DBExporter_1" ${TOPP_BIN_PATH}/DBExporter  -out DBExporter_output.tmp -ini DBExporter_input.ini -id 1)
	add_test("TOPP_DBExporter_1_out1" ${DIFF} -in1 DBExporter_output.tmp -in2 DBImporter_input.mzData)
#TODO	add_test("TOPP_DBExporter_2" ${TOPP_BIN_PATH}/DBExporter  -out DBExporter_output2.tmp -ini DBExporter_input.ini -query 'SELECT 1')
#add_test("TOPP_DBExporter_2_out1" ${DIFF} -in1 DBExporter_output2.tmp -in2 DBImporter_input.mzData)
endif()

### Decharger tests
add_test("TOPP_Decharger_1" ${TOPP_BIN_PATH}/Decharger  -ini Decharger.ini -in Decharger_input.featureXML -out Decharger_output.tmp)
add_test("TOPP_Decharger_1_out1" ${DIFF} -in1 Decharger_output.tmp -in2 Decharger_output.featureXML)

### DTAExtractor tests
add_test("TOPP_DTAExtractor_1" ${TOPP_BIN_PATH}/DTAExtractor -in DTAExtractor_1_input.mzData -out DTAExtractor -rt :61)
add_test("TOPP_DTAExtractor_1_out1" ${DIFF} -in1 DTAExtractor_RT60.dta -in2 DTAExtractor_1_output.dta)
add_test("TOPP_DTAExtractor_2" ${TOPP_BIN_PATH}/DTAExtractor -in DTAExtractor_1_input.mzData -out DTAExtractor -level 1)
add_test("TOPP_DTAExtractor_2_out1" ${DIFF} -in1 DTAExtractor_RT60.dta -in2 DTAExtractor_2_output.dta)
add_test("TOPP_DTAExtractor_3" ${TOPP_BIN_PATH}/DTAExtractor -in DTAExtractor_1_input.mzData -out DTAExtractor -level 2 -mz :1000)
add_test("TOPP_DTAExtractor_3_out1" ${DIFF} -in1 DTAExtractor_RT140_MZ5.dta -in2 DTAExtractor_3_output.dta)

### FeatureFinder tests
add_test("TOPP_FeatureFinder_1" ${TOPP_BIN_PATH}/FeatureFinder  -ini FeatureFinder_1_parameters.ini -in FeatureFinder_1_input.mzData -out FeatureFinder_1_output.tmp)
add_test("TOPP_FeatureFinder_1_out1" ${DIFF} -in1 FeatureFinder_1_output.tmp -in2 FeatureFinder_1_output.featureXML)
add_test("TOPP_FeatureFinder_2" ${TOPP_BIN_PATH}/FeatureFinder  -ini FeatureFinder_2_parameters.ini -in FeatureFinder_1_input.mzData -out FeatureFinder_2_output.tmp)
add_test("TOPP_FeatureFinder_2_out1" ${DIFF} -in1 FeatureFinder_2_output.tmp -in2 FeatureFinder_2_output.featureXML)
add_test("TOPP_FeatureFinder_3" ${TOPP_BIN_PATH}/FeatureFinder  -ini FeatureFinder_3_parameters.ini -in FeatureFinder_1_input.mzData -out FeatureFinder_3_output.tmp)
add_test("TOPP_FeatureFinder_3_out1" ${DIFF} -in1 FeatureFinder_3_output.tmp -in2 FeatureFinder_3_output.featureXML)
add_test("TOPP_FeatureFinder_4" ${TOPP_BIN_PATH}/FeatureFinder  -ini FeatureFinder_4_parameters.ini -in FeatureFinder_4_input.mzData -out FeatureFinder_4_output.tmp)
add_test("TOPP_FeatureFinder_4_out1" ${DIFF} -in1 FeatureFinder_4_output.tmp -in2 FeatureFinder_4_output.featureXML)

### FileConverter tests
add_test("TOPP_FileConverter_1" ${TOPP_BIN_PATH}/FileConverter -in FileConverter_1_input.txt -in_type mzData -out FileConverter.tmp -out_type DTA2D)
add_test("TOPP_FileConverter_1_out1" ${DIFF} -in1 FileConverter.tmp -in2 FileConverter_1_output.dta2d)
add_test("TOPP_FileConverter_2" ${TOPP_BIN_PATH}/FileConverter -in FileConverter_2_input.txt -in_type DTA2D -out FileConverter.tmp -out_type DTA2D)
add_test("TOPP_FileConverter_2_out1" ${DIFF} -in1 FileConverter.tmp -in2 FileConverter_2_output.dta2d)
add_test("TOPP_FileConverter_3" ${TOPP_BIN_PATH}/FileConverter -in FileConverter_3_input.txt -in_type featureXML -out FileConverter.tmp -out_type DTA2D)
add_test("TOPP_FileConverter_3_out1" ${DIFF} -in1 FileConverter.tmp -in2 FileConverter_3_output.dta2d)
add_test("TOPP_FileConverter_4" ${TOPP_BIN_PATH}/FileConverter -in FileConverter_4_input.txt -in_type mzXML -out FileConverter.tmp -out_type DTA2D)
add_test("TOPP_FileConverter_4_out1" ${DIFF} -in1 FileConverter.tmp -in2 FileConverter_4_output.dta2d)
if (USE_ANDIMS)
	add_test("TOPP_FileConverter_5" ${TOPP_BIN_PATH}/FileConverter -in FileConverter_5_input.txt -in_type cdf -out FileConverter.tmp -out_type DTA2D)
	add_test("TOPP_FileConverter_5_out1" ${DIFF} -in1 FileConverter.tmp -in2 FileConverter_5_output.dta2d)
endif()
add_test("TOPP_FileConverter_6" ${TOPP_BIN_PATH}/FileConverter -in FileConverter_6_input.dta2d -out FileConverter.tmp -out_type mzXML)
add_test("TOPP_FileConverter_7" ${TOPP_BIN_PATH}/FileConverter -in FileConverter.tmp -in_type mzXML -out FileConverter2.tmp -out_type DTA2D)
add_test("TOPP_FileConverter_7_out1" ${DIFF} -in1 FileConverter2.tmp -in2 FileConverter_6_input.dta2d)
add_test("TOPP_FileConverter_8" ${TOPP_BIN_PATH}/FileConverter -in FileConverter_6_input.dta2d -out FileConverter.tmp -out_type mzData)
add_test("TOPP_FileConverter_9" ${TOPP_BIN_PATH}/FileConverter -in FileConverter.tmp -in_type mzData -out FileConverter2.tmp -out_type DTA2D)
add_test("TOPP_FileConverter_9_out1" ${DIFF} -in1 FileConverter2.tmp -in2 FileConverter_6_input.dta2d)
add_test("TOPP_FileConverter_10" ${TOPP_BIN_PATH}/FileConverter -in FileConverter_7_input.consensusXML -out FileConverter_7_output.tmp -out_type featureXML)
add_test("TOPP_FileConverter_10_out1" ${DIFF} -in1 FileConverter_7_output.tmp -in2 FileConverter_7_output.featureXML)
add_test("TOPP_FileConverter_11" ${TOPP_BIN_PATH}/FileConverter -in FileConverter_8_input.mzML -out FileConverter_8_output.tmp -out_type mzData)
add_test("TOPP_FileConverter_11_out1" ${DIFF} -in1 FileConverter_8_output.tmp -in2 FileConverter_8_output.mzData)

### FileFilter tests
add_test("TOPP_FileFilter_1" ${TOPP_BIN_PATH}/FileFilter -in FileFilter_1_input.mzData -out FileFilter_1_out.tmp -rt :30 -mz :1000 -int :20000)
add_test("TOPP_FileFilter_1_out1" ${DIFF} -in1 FileFilter_1_out.tmp -in2 FileFilter_1_output.mzData)
add_test("TOPP_FileFilter_2" ${TOPP_BIN_PATH}/FileFilter -in FileFilter_1_input.mzData -out FileFilter_2_out.tmp -rt 30: -mz 1000: -int 100:)
add_test("TOPP_FileFilter_2_out1" ${DIFF} -in1 FileFilter_2_out.tmp -in2 FileFilter_2_output.mzData)
add_test("TOPP_FileFilter_3" ${TOPP_BIN_PATH}/FileFilter -in FileFilter_1_input.mzData -out FileFilter_3_out.tmp -level 2)
add_test("TOPP_FileFilter_3_out1" ${DIFF} -in1 FileFilter_3_out.tmp -in2 FileFilter_3_output.mzData)
add_test("TOPP_FileFilter_4" ${TOPP_BIN_PATH}/FileFilter -in FileFilter_4_input.mzData -out FileFilter_4_out.tmp -remove_mode Zoom)
add_test("TOPP_FileFilter_4_out1" ${DIFF} -in1 FileFilter_4_out.tmp -in2 FileFilter_4_output.mzData)
add_test("TOPP_FileFilter_5" ${TOPP_BIN_PATH}/FileFilter -in FileFilter_5_input.featureXML -out FileFilter_5_out.tmp -rt :1000 -mz :480 -int :79000 -charge :3 -q :0.6)
add_test("TOPP_FileFilter_5_out1" ${DIFF} -in1 FileFilter_5_out.tmp -in2 FileFilter_5_out.featureXML)
add_test("TOPP_FileFilter_6" ${TOPP_BIN_PATH}/FileFilter -in FileFilter_5_input.featureXML -out FileFilter_6_out.tmp -rt 1000: -mz 440: -int 70000: -charge 3: -q 0.51:)
add_test("TOPP_FileFilter_6_out1" ${DIFF} -in1 FileFilter_6_out.tmp -in2 FileFilter_6_out.featureXML)
add_test("TOPP_FileFilter_7" ${TOPP_BIN_PATH}/FileFilter -in FileFilter_7_input.mzData -out FileFilter_7_out.tmp -int 7000: -level 1 2 3)
add_test("TOPP_FileFilter_7_out1" ${DIFF} -in1 FileFilter_7_out.tmp -in2 FileFilter_7_out.mzData)
add_test("TOPP_FileFilter_8" ${TOPP_BIN_PATH}/FileFilter -in FileFilter_8_input.consensusXML -out FileFilter_8_out.tmp -rt 600:1400 -mz 700:2300 -int 1100:6000)
add_test("TOPP_FileFilter_8_out1" ${DIFF} -in1 FileFilter_8_out.tmp -in2 FileFilter_8_out.consensusXML)

### FileInfo tests
#TODO add_test("TOPP_FileInfo_1" ${TOPP_BIN_PATH}/FileInfo -in FileInfo_1_input.txt -in_type DTA -no_progress > FileInfo.tmp)
#add_test("TOPP_FileInfo_1_out1" ${DIFF} -in1 FileInfo.tmp -in2 FileInfo_1_output.txt)
#add_test("TOPP_FileInfo_2" ${TOPP_BIN_PATH}/FileInfo -in FileInfo_2_input.dta2d -no_progress > FileInfo.tmp)
#add_test("TOPP_FileInfo_2_out1" ${DIFF} -in1 FileInfo.tmp -in2 FileInfo_2_output.txt)
#add_test("TOPP_FileInfo_3" ${TOPP_BIN_PATH}/FileInfo -in FileInfo_3_input.featureXML -m -s -no_progress > FileInfo.tmp)
#add_test("TOPP_FileInfo_3_out1" ${DIFF} -in1 FileInfo.tmp -in2 FileInfo_3_output.txt)
#add_test("TOPP_FileInfo_4" ${TOPP_BIN_PATH}/FileInfo -in FileInfo_4_input.mzXML -m -no_progress > FileInfo.tmp)
#add_test("TOPP_FileInfo_4_out1" ${DIFF} -in1 FileInfo.tmp -in2 FileInfo_4_output.txt)
#add_test("TOPP_FileInfo_5" ${TOPP_BIN_PATH}/FileInfo -in FileInfo_5_input.txt -in_type mzData -m -s -no_progress > FileInfo.tmp)
#add_test("TOPP_FileInfo_5_out1" ${DIFF} -in1 FileInfo.tmp -in2 FileInfo_5_output.txt)
#add_test("TOPP_FileInfo_6" ${TOPP_BIN_PATH}/FileInfo -in FileInfo_6_input.mzData -d -s -no_progress > FileInfo.tmp)
#add_test("TOPP_FileInfo_6_out1" ${DIFF} -in1 FileInfo.tmp -in2 FileInfo_6_output.txt)
#add_test("TOPP_FileInfo_7" ${TOPP_BIN_PATH}/FileInfo -in FileInfo_7_input.consensusXML -s -m -no_progress > FileInfo.tmp)
#add_test("TOPP_FileInfo_7_out1" ${DIFF} -in1 FileInfo.tmp -in2 FileInfo_7_output.txt)
#if (USE_ANDIMS)
#	add_test("TOPP_FileInfo_8" ${TOPP_BIN_PATH}/FileInfo -in FileInfo_8_input.notype -in_type cdf -m -no_progress > FileInfo.tmp)
#	add_test("TOPP_FileInfo_8_out1" ${DIFF} -in1 FileInfo.tmp -in2 FileInfo_8_output.txt)
#endif()
#add_test("TOPP_FileInfo_9" ${TOPP_BIN_PATH}/FileInfo -in FileInfo_9_input.mzML -m -no_progress > FileInfo.tmp)
#add_test("TOPP_FileInfo_9_out1" ${DIFF} -in1 FileInfo.tmp -in2 FileInfo_9_output.txt)

### FileMerger tests
add_test("TOPP_FileMerger_1" ${TOPP_BIN_PATH}/FileMerger -in FileMerger_1_input1.dta2d FileMerger_1_input2.dta FileMerger_1_input3.dta2d -out FileMerger.tmp -rt_custom 1 2 3)
add_test("TOPP_FileMerger_2" ${TOPP_BIN_PATH}/FileConverter -in FileMerger.tmp -in_type mzData -out FileMerger2.tmp -out_type DTA2D)
add_test("TOPP_FileMerger_2_out1" ${DIFF} -in1 FileMerger2.tmp -in2 FileMerger_1_output.txt)
add_test("TOPP_FileMerger_3" ${TOPP_BIN_PATH}/FileMerger -in FileMerger_2_input1.txt FileMerger_2_input2.txt -in_type DTA -out FileMerger.tmp -rt_custom 5 10)
add_test("TOPP_FileMerger_4" ${TOPP_BIN_PATH}/FileConverter -in FileMerger.tmp -in_type mzData -out FileMerger2.tmp -out_type DTA2D)
add_test("TOPP_FileMerger_4_out1" ${DIFF} -in1 FileMerger2.tmp -in2 FileMerger_2_output.txt)
add_test("TOPP_FileMerger_5" ${TOPP_BIN_PATH}/FileMerger -in FileMerger_3_input1.dta FileMerger_3_input2.dta -out FileMerger.tmp -rt_auto)
add_test("TOPP_FileMerger_6" ${TOPP_BIN_PATH}/FileConverter -in FileMerger.tmp -in_type mzData -out FileMerger2.tmp -out_type DTA2D)
add_test("TOPP_FileMerger_6_out1" ${DIFF} -in1 FileMerger2.tmp -in2 FileMerger_3_output.txt)
add_test("TOPP_FileMerger_7" ${TOPP_BIN_PATH}/FileMerger -in FileMerger_4_input1.dta2d FileMerger_4_input2.dta2d -out FileMerger.tmp)
add_test("TOPP_FileMerger_8" ${TOPP_BIN_PATH}/FileConverter -in FileMerger.tmp -in_type mzData -out FileMerger2.tmp -out_type DTA2D)
add_test("TOPP_FileMerger_8_out1" ${DIFF} -in1 FileMerger2.tmp -in2 FileMerger_4_output.txt)
add_test("TOPP_FileMerger_9" ${TOPP_BIN_PATH}/FileMerger -in FileMerger_5_input_rt1023.331.dta FileMerger_5_input_rt2044.334.dta FileMerger_5_input_rt889.32.dta -rt_filename -user_ms_level -ms_level 2 -out FileMerger5.tmp)
add_test("TOPP_FileMerger_9_out1" ${DIFF} -in1 FileMerger5.tmp -in2 FileMerger_5_output.txt)
add_test("TOPP_FileMerger_10" ${TOPP_BIN_PATH}/FileMerger -in FileMerger_6_input1.mzML FileMerger_6_input2.mzML -out FileMerger6.tmp)
add_test("TOPP_FileMerger_10_out1" ${DIFF} -in1 FileMerger6.tmp -in2 FileMerger_6_output.mzData)

### IDRTCalibration tests
add_test("TOPP_IDRTCalibration_1" ${TOPP_BIN_PATH}/IDRTCalibration -in IDRTCalibration_1_input.IdXML -out IDRTCalibration_1_output.tmp -calibrant_1_input 10 -calibrant_2_input 90)
add_test("TOPP_IDRTCalibration_1_out1" ${DIFF} -in1 IDRTCalibration_1_output.tmp -in2 IDRTCalibration_1_output.IdXML)

### IDFilter tests
add_test("TOPP_IDFilter_1" ${TOPP_BIN_PATH}/IDFilter -in IDFilter_1_input.idXML -out IDFilter_1_output.tmp -sequences_file IDFilter_1_input.fas)
add_test("TOPP_IDFilter_1_out1" ${DIFF} -in1 IDFilter_1_output.tmp -in2 IDFilter_1_output.idXML)
add_test("TOPP_IDFilter_2" ${TOPP_BIN_PATH}/IDFilter -in IDFilter_2_input.idXML -out IDFilter_2_output.tmp -pep_fraction 1 -prot_fraction 1)
add_test("TOPP_IDFilter_2_out1" ${DIFF} -in1 IDFilter_2_output.tmp -in2 IDFilter_2_output.idXML)
add_test("TOPP_IDFilter_3" ${TOPP_BIN_PATH}/IDFilter -in IDFilter_3_input.idXML -out IDFilter_3_output.tmp -exclusion_peptides_file IDFilter_3_2_input.idXML)
add_test("TOPP_IDFilter_3_out1" ${DIFF} -in1 IDFilter_3_output.tmp -in2 IDFilter_3_output.idXML)
add_test("TOPP_IDFilter_4" ${TOPP_BIN_PATH}/IDFilter -in IDFilter_4_input.idXML -out IDFilter_4_output.tmp -rt_filtering -p_value 0.08)
add_test("TOPP_IDFilter_4_out1" ${DIFF} -in1 IDFilter_4_output.tmp -in2 IDFilter_4_output.idXML)
add_test("TOPP_IDFilter_5" ${TOPP_BIN_PATH}/IDFilter -in IDFilter_5_input.idXML -out IDFilter_5_output.tmp -pep_score 32 -prot_score 25)
add_test("TOPP_IDFilter_5_out1" ${DIFF} -in1 IDFilter_5_output.tmp -in2 IDFilter_5_output.idXML)
add_test("TOPP_IDFilter_6" ${TOPP_BIN_PATH}/IDFilter -in IDFilter_6_input.idXML -out IDFilter_6_output.tmp -best_n_peptide_hits 2 -best_n_protein_hits 10)
add_test("TOPP_IDFilter_6_out1" ${DIFF} -in1 IDFilter_6_output.tmp -in2 IDFilter_6_output.idXML)

### MapAligner tests
add_test("TOPP_MapAligner_1" ${TOPP_BIN_PATH}/MapAligner -ini MapAligner_1_parameters.ini)
add_test("TOPP_MapAligner_1_out1" ${DIFF} -in1 MapAligner_1_output1.tmp -in2 MapAligner_1_output1.featureXML)
add_test("TOPP_MapAligner_1_out2" ${DIFF} -in1 MapAligner_1_output2.tmp -in2 MapAligner_1_output2.featureXML)
add_test("TOPP_MapAligner_1_out3" ${DIFF} -in1 MapAligner_1_output3.tmp -in2 MapAligner_1_output3.featureXML)
add_test("TOPP_MapAligner_1_out4" ${DIFF} -in1 MapAligner_1_trafo1.tmp -in2 MapAligner_1_trafo1.trafoXML)
add_test("TOPP_MapAligner_1_out5" ${DIFF} -in1 MapAligner_1_trafo2.tmp -in2 MapAligner_1_trafo2.trafoXML)
add_test("TOPP_MapAligner_1_out6" ${DIFF} -in1 MapAligner_1_trafo3.tmp -in2 MapAligner_1_trafo3.trafoXML)
add_test("TOPP_MapAligner_2" ${TOPP_BIN_PATH}/MapAligner -ini MapAligner_2_parameters.ini)
add_test("TOPP_MapAligner_2_out1" ${DIFF} -in1 MapAligner_2_output1.tmp -in2 MapAligner_2_output1.mzData)
add_test("TOPP_MapAligner_2_out2" ${DIFF} -in1 MapAligner_2_output2.tmp -in2 MapAligner_2_output2.mzData)
add_test("TOPP_MapAligner_2_out3" ${DIFF} -in1 MapAligner_2_output3.tmp -in2 MapAligner_2_output3.mzData)
add_test("TOPP_MapAligner_3" ${TOPP_BIN_PATH}/MapAligner -ini MapAligner_3_parameters.ini)
add_test("TOPP_MapAligner_3_out1" ${DIFF} -in1 MapAligner_3_output1.tmp -in2 MapAligner_3_output1.mzData)
add_test("TOPP_MapAligner_3_out2" ${DIFF} -in1 MapAligner_3_output2.tmp -in2 MapAligner_3_output2.mzData)
add_test("TOPP_MapAligner_3_out3" ${DIFF} -in1 MapAligner_3_output3.tmp -in2 MapAligner_3_output3.mzData)
add_test("TOPP_MapAligner_4" ${TOPP_BIN_PATH}/MapAligner -type apply_given_trafo -in MapAligner_4_input.featureXML MapAligner_4_input.featureXML -given_transformations MapAligner_4_trafo1.trafoXML MapAligner_4_trafo2.trafoXML -out MapAligner_4_output1.tmp MapAligner_4_output2.tmp)
add_test("TOPP_MapAligner_4_out1" ${DIFF} -in1 MapAligner_4_output1.tmp -in2 MapAligner_4_output1.featureXML)
add_test("TOPP_MapAligner_4_out2" ${DIFF} -in1 MapAligner_4_output2.tmp -in2 MapAligner_4_output2.featureXML)
add_test("TOPP_MapAligner_5" ${TOPP_BIN_PATH}/MapAligner -type apply_given_trafo -in MapAligner_5_input.mzData -given_transformations MapAligner_5_trafo.trafoXML -out MapAligner_5_output.tmp)
add_test("TOPP_MapAligner_5_out1" ${DIFF} -in1 MapAligner_5_output.tmp -in2 MapAligner_5_output.mzData)

### MascotAdapter tests
add_test("TOPP_MascotAdapter_1" ${TOPP_BIN_PATH}/MascotAdapter -ini MascotAdapter_1_parameters.ini -mascot_in)
add_test("TOPP_MascotAdapter_1_out1" ${DIFF} -in1 MascotAdapter_1_output.tmp -in2 MascotAdapter_1_output.mascot_in)
add_test("TOPP_MascotAdapter_2" ${TOPP_BIN_PATH}/MascotAdapter -in MascotAdapter_2_input.mascotXML -out MascotAdapter_2_output.tmp -mascot_out)
add_test("TOPP_MascotAdapter_2_out1" ${DIFF} -in1 MascotAdapter_2_output.tmp -in2 MascotAdapter_2_output.idXML)

### NoiseFilter tests
add_test("TOPP_NoiseFilter_1" ${TOPP_BIN_PATH}/NoiseFilter  -ini NoiseFilter_1_parameters.ini -in NoiseFilter_1_input.mzData -out NoiseFilter_1_output.tmp)
add_test("TOPP_NoiseFilter_1_out1" ${DIFF} -in1 NoiseFilter_1_output.tmp -in2 NoiseFilter_1_output.mzData)
add_test("TOPP_NoiseFilter_2" ${TOPP_BIN_PATH}/NoiseFilter  -ini NoiseFilter_2_parameters.ini -in RawData_input.mzData -out NoiseFilter_2_output.tmp)
add_test("TOPP_NoiseFilter_2_out1" ${DIFF} -in1 NoiseFilter_2_output.tmp -in2 NoiseFilter_2_output.mzData)

### PeakPicker tests
add_test("TOPP_PeakPicker_1" ${TOPP_BIN_PATH}/PeakPicker  -ini PeakPicker_parameters.ini -in RawData_input.mzData -out PeakPicker_output.tmp)
add_test("TOPP_PeakPicker_1_out1" ${DIFF} -in1 PeakPicker_output.tmp -in2 PeakPicker_output.mzData)
add_test("TOPP_PeakPicker_2" ${TOPP_BIN_PATH}/PeakPicker  -ini PeakPicker_deconv_parameters.ini -in PeakPicker_deconv_input.mzData -out PeakPicker_deconv_output.tmp)
add_test("TOPP_PeakPicker_2_out1" ${DIFF} -in1 PeakPicker_deconv_output.tmp -in2 PeakPicker_deconv_output.mzData)
add_test("TOPP_PeakPicker_3" ${TOPP_BIN_PATH}/PeakPicker  -ini PeakPicker_2Dopt_parameters.ini -in RawData_input.mzData -out PeakPicker_2Dopt_output.tmp)
add_test("TOPP_PeakPicker_3_out1" ${DIFF} -in1 PeakPicker_2Dopt_output.tmp -in2 PeakPicker_2Dopt_output.mzData)

### Resampler tests
add_test("TOPP_Resampler_1" ${TOPP_BIN_PATH}/Resampler -in Resampler_1_input.mzData -out Resampler_1_output.tmp -sampling_rate 0.3)
add_test("TOPP_Resampler_1_out1" ${DIFF} -in1 Resampler_1_output.tmp -in2 Resampler_1_output.mzData)
add_test("TOPP_Resampler_2" ${TOPP_BIN_PATH}/Resampler -in Resampler_2_input.mzData -out Resampler_2_output.tmp -image -width 20 -height 15)
add_test("TOPP_Resampler_2_out1" ${DIFF} -in1 Resampler_2_output.tmp -in2 Resampler_2_output.png)

### RTModel tests
add_test("TOPP_RTModel_1" ${TOPP_BIN_PATH}/RTModel -in RTModel_1_input.idXML -out RTModel_1_output.tmp -ini RTModel_1_parameters.ini)
add_test("TOPP_RTModel_1_out1" ${DIFF} -in1 RTModel_1_output.tmp -in2 RTModel_1_output.model)
add_test("TOPP_RTModel_2" ${TOPP_BIN_PATH}/RTModel -in_positive RTModel_2_input_positive.idXML -in_negative RTModel_2_input_negative.idXML -out RTModel_2_output.tmp -ini RTModel_2_parameters.ini)
add_test("TOPP_RTModel_2_out1" ${DIFF} -in1 RTModel_2_output.tmp -in2 RTModel_2_output.model)
add_test("TOPP_RTModel_2_out2" ${DIFF} -in1 RTModel_2_output.tmp_samples -in2 RTModel_2_output.samples)
add_test("TOPP_RTModel_2_out3" ${DIFF} -in1 RTModel_2_output.tmp_additional_parameters -in2 RTModel_2_output.additional_parameters)
#REMOVE(-rf RTModel_2_output.tmp_additional_parameters RTModel_2_output.tmp_samples
add_test("TOPP_RTModel_3" ${TOPP_BIN_PATH}/RTModel -in RTModel_3_input.idXML -out RTModel_3_output.tmp -ini RTModel_3_parameters.ini)
add_test("TOPP_RTModel_3_out1" ${DIFF} -in1 RTModel_3_output.tmp -in2 RTModel_3_output.model)
add_test("TOPP_RTModel_3_out2" ${DIFF} -in1 RTModel_3_output.tmp_samples -in2 RTModel_3_output.samples)
#REMOVE(-rf RTModel_3_output.tmp_additional_parameters RTModel_3_output.tmp_samples points.txt
add_test("TOPP_RTModel_4" ${TOPP_BIN_PATH}/RTModel -in RTModel_4_input.txt -out RTModel_4_output.tmp -ini RTModel_4_parameters.ini -textfile_input)
add_test("TOPP_RTModel_4_out1" ${DIFF} -in1 RTModel_4_output.tmp -in2 RTModel_4_output.model)
add_test("TOPP_RTModel_4_out2" ${DIFF} -in1 RTModel_4_output.tmp_samples -in2 RTModel_4_output.samples)
#REMOVE(-rf RTModel_4_output.tmp_additional_parameters RTModel_4_output.tmp_samples points.txt

### PTModel tests
add_test("TOPP_PTModel_1" ${TOPP_BIN_PATH}/PTModel -in_positive PTModel_1_input_positive.idXML -in_negative PTModel_1_input_negative.idXML -out PTModel_1_output.tmp -ini PTModel_1_parameters.ini)
#SED('s/^.*prob.*$$//g' PTModel_1_output.tmp > PTModel_1_output.tmp_
#TODO add_test("TOPP_PTModel_1_out1" ${DIFF} -in1 PTModel_1_output.tmp_ -in2 PTModel_1_output.model)
add_test("TOPP_PTModel_1_out2" ${DIFF} -in1 PTModel_1_output.tmp_samples -in2 PTModel_1_output.samples)
add_test("TOPP_PTModel_1_out3" ${DIFF} -in1 PTModel_1_output.tmp_additional_parameters -in2 PTModel_1_output.additional_parameters)
#REMOVE(-rf PTModel_1_output.tmp_additional_parameters PTModel_1_output.tmp_samples

### PTPredict tests
add_test("TOPP_PTPredict_1" ${TOPP_BIN_PATH}/PTPredict -in PTPredict_1_input.idXML -out PTPredict_1_output.tmp -svm_model PTPredict_1_input.model)
add_test("TOPP_PTPredict_1_out1" ${DIFF} -in1 PTPredict_1_output.tmp -in2 PTPredict_1_output.idXML)

### RTPredict tests
add_test("TOPP_RTPredict_1" ${TOPP_BIN_PATH}/RTPredict -in RTPredict_1_input.idXML -out RTPredict_1_output.tmp -total_gradient_time 3000 -svm_model RTPredict_1_input.model)
add_test("TOPP_RTPredict_1_out1" ${DIFF} -in1 RTPredict_1_output.tmp -in2 RTPredict_1_output.idXML)
add_test("TOPP_RTPredict_2" ${TOPP_BIN_PATH}/RTPredict -in RTPredict_2_input.idXML -out_positive RTPredict_2_output_positive.tmp -out_negative RTPredict_2_output_negative.tmp -svm_model RTPredict_2_input.model)
add_test("TOPP_RTPredict_2_out1" ${DIFF} -in1 RTPredict_2_output_positive.tmp -in2 RTPredict_2_output_positive.idXML)
add_test("TOPP_RTPredict_2_out2" ${DIFF} -in1 RTPredict_2_output_negative.tmp -in2 RTPredict_2_output_negative.idXML)
add_test("TOPP_RTPredict_3" ${TOPP_BIN_PATH}/RTPredict -in RTPredict_3_input.idXML -out RTPredict_3_output.tmp_ -total_gradient_time 1 -svm_model RTPredict_3_input.model)
#SED('s/^.*predicted_RT_p_value.*$$//g' RTPredict_3_output.tmp_ > RTPredict_3_output.tmp
#TODO add_test("TOPP_RTPredict_3_out1" ${DIFF} -in1 RTPredict_3_output.tmp -in2 RTPredict_3_output.idXML)
add_test("TOPP_RTPredict_4" ${TOPP_BIN_PATH}/RTPredict -in RTPredict_4_input.txt -out RTPredict_4_output.tmp -total_gradient_time 1 -svm_model RTPredict_4_input.model -textfile_input)
add_test("TOPP_RTPredict_4_out1" ${DIFF} -in1 RTPredict_4_output.tmp -in2 RTPredict_4_output.txt)
add_test("TOPP_RTPredict_5" ${TOPP_BIN_PATH}/RTPredict -in RTPredict_5_input.idXML -out RTPredict_5_output.tmp -total_gradient_time 1 -svm_model RTPredict_5_input.model -textfile_output)
add_test("TOPP_RTPredict_5_out1" ${DIFF} -in1 RTPredict_5_output.tmp -in2 RTPredict_5_output.txt)

### InspectAdapter tests
#TODO add_test("TOPP_InspectAdapter_1" ${TOPP_BIN_PATH}/InspectAdapter -ini InspectAdapter_1_parameters.ini -list_modifications > InspectAdapter_1_output.tmp)
#add_test("TOPP_InspectAdapter_1_out1" ${DIFF} -in1 InspectAdapter_1_output.tmp -in2 InspectAdapter_1_output.modification_list.txt)
#add_test("TOPP_InspectAdapter_2" ${TOPP_BIN_PATH}/InspectAdapter -ini InspectAdapter_1_parameters.ini -inspect_in -out InspectAdapter_2_output.tmp)
##SED('s/^spectra,.*$$/spectra,dummy.mzXML/g' InspectAdapter_2_output.tmp > InspectAdapter_2_output.tmp_3
##SED('s/^db,.*$$/db,dummy.trie/g' InspectAdapter_2_output.tmp_3 > InspectAdapter_2_output.tmp
#add_test("TOPP_InspectAdapter_2_out1" ${DIFF} -in1 InspectAdapter_2_output.tmp -in2 InspectAdapter_2_output.inspect_in)
#add_test("TOPP_InspectAdapter_3" ${TOPP_BIN_PATH}/InspectAdapter -ini InspectAdapter_1_parameters.ini -inspect_in -in Inspect.mzData -out InspectAdapter_3_output.tmp)
##SED('s/^spectra,.*$$/spectra,dummy.mzXML/g' InspectAdapter_3_output.tmp > InspectAdapter_3_output.tmp_2
##SED('s/^db,.*$$/db,dummy.trie/g' InspectAdapter_3_output.tmp_2 > InspectAdapter_3_output.tmp
#add_test("TOPP_InspectAdapter_3_out1" ${DIFF} -in1 InspectAdapter_3_output.tmp -in2 InspectAdapter_2_output.inspect_in)
#add_test("TOPP_InspectAdapter_4" ${TOPP_BIN_PATH}/InspectAdapter -ini InspectAdapter_1_parameters.ini -in InspectAdapter.out -out InspectAdapter_4_output.tmp_ -inspect_out)
##SED('s/^.*date.*$$/# date removed/g' InspectAdapter_4_output.tmp_ > InspectAdapter_4_output.tmp
#add_test("TOPP_InspectAdapter_4_out1" ${DIFF} -in1 InspectAdapter_4_output.tmp -in2 InspectAdapter_4_output.idXML)
##REMOVE(-rf Inspect_FASTAFile_test4.tmp.trie Inspect_FASTAFile_test4.tmp.index InspectAdapter_4_output.tmp_

### SequestAdapter tests
#TODO add_test("TOPP_SequestAdapter_1" ${TOPP_BIN_PATH}/SequestAdapter -ini SequestAdapter_1_parameters.ini -list_modifications > SequestAdapter_1_output.tmp)
#add_test("TOPP_SequestAdapter_1_out1" ${DIFF} -in1 SequestAdapter_1_output.tmp -in2 SequestAdapter_1_output.modification_list.txt)
##REMOVE(-rf ./tmp/Sequest*.dta_* ./tmp/sequest_run.bat ./tmp/temp.sequest.in
#add_test("TOPP_SequestAdapter_2" ${TOPP_BIN_PATH}/SequestAdapter -ini SequestAdapter_1_parameters.ini -in Sequest.mzXML -out SequestAdapter_2_output.tmp -sequest_in)
#add_test("TOPP_SequestAdapter_2_out1" ${DIFF} -in1 SequestAdapter_2_output.tmp -in2 SequestAdapter_2_output.sequest_in)
##REMOVE(-rf ./tmp/Sequest*.dta_*
##REMOVE(-rf ./tmp/Sequest*.dta_* ./tmp/sequest_run.bat ./tmp/temp.sequest.in
#add_test("TOPP_SequestAdapter_3" ${TOPP_BIN_PATH}/SequestAdapter -ini SequestAdapter_1_parameters.ini -in Sequest.mzData -out SequestAdapter_3_output.tmp -sequest_in -in Sequest.mzData)
#add_test("TOPP_SequestAdapter_3_out1" ${DIFF} -in1 SequestAdapter_3_output.tmp -in2 SequestAdapter_2_output.sequest_in)
##REMOVE(-rf ./tmp/Sequest*.dta_*
#add_test("TOPP_SequestAdapter_4" ${TOPP_BIN_PATH}/SequestAdapter -ini SequestAdapter_2_parameters.ini -in tmp -out SequestAdapter_4_output.tmp -sequest_out)
#add_test("TOPP_SequestAdapter_4_out1" ${DIFF} -in1 SequestAdapter_4_output.tmp -in2 SequestAdapter_4_output.idXML)
##REMOVE(-rf ./tmp/sequest_run.bat ./tmp/temp.sequest.in

### PepNovoAdapter tests
#TODO add_test("TOPP_PepNovoAdapter_1" ${TOPP_BIN_PATH}/PepNovoAdapter -ini PepNovoAdapter_1_parameters.ini -list_modifications > PepNovoAdapter_1_output.tmp)
#add_test("TOPP_PepNovoAdapter_1_out1" ${DIFF} -in1 PepNovoAdapter_1_output.tmp -in2 PepNovoAdapter_1_output.modification_list.txt)
#add_test("TOPP_PepNovoAdapter_2" ${TOPP_BIN_PATH}/PepNovoAdapter -ini PepNovoAdapter_1_parameters.ini -list_models > PepNovoAdapter_2_output.tmp)
#add_test("TOPP_PepNovoAdapter_2_out1" ${DIFF} -in1 PepNovoAdapter_2_output.tmp -in2 PepNovoAdapter_2_output.model_list.txt)
#add_test("TOPP_PepNovoAdapter_3" ${TOPP_BIN_PATH}/PepNovoAdapter -ini PepNovoAdapter_1_parameters.ini -in PepNovo.mzXML -pepnovo_in -out PepNovoAdapter_3_output.tmp)
#add_test("TOPP_PepNovoAdapter_3_out1" ${DIFF} -in1 ./tmp/PepNovo_PTMs_.txt -in2 ./tmp/PepNovo_PTMs.txt)
##REMOVE(-rf ./tmp/PepNovo*.dta
#add_test("TOPP_PepNovoAdapter_4" ${TOPP_BIN_PATH}/PepNovoAdapter -ini PepNovoAdapter_1_parameters.ini -in PepNovo.mzData -pepnovo_in -out PepNovoAdapter_4_output.tmp)
#add_test("TOPP_PepNovoAdapter_4_out1" ${DIFF} -in1 ./tmp/PepNovo_PTMs_.txt -in2 ./tmp/PepNovo_PTMs.txt)
##REMOVE(-rf ./tmp/PepNovo*.dta ./tmp/PepNovo_PTMs.txt ./tmp/dta_list.txt
#add_test("TOPP_PepNovoAdapter_5" ${TOPP_BIN_PATH}/PepNovoAdapter -ini PepNovoAdapter_5_parameters.ini -in PepNovoAdapter_5_output.pepnovo_out -out PepNovoAdapter_5_output.tmp_ -pepnovo_out)
##SED('s/^.*date_group_1.*$$/# date group removed/g' PepNovoAdapter_5_output.tmp_ > PepNovoAdapter_5_output.tmp
#add_test("TOPP_PepNovoAdapter_5_out1" ${DIFF} -in1 PepNovoAdapter_5_output.tmp -in2 PepNovoAdapter_5_output.idXML)

### SpectraFilter tests
add_test("TOPP_SpectraFilter_1" ${TOPP_BIN_PATH}/SpectraFilter -in SpectraFilter_1_input.mzData -out SpectraFilter_1_output.mzData.tmp -type SqrtMower)
add_test("TOPP_SpectraFilter_1_out1" ${DIFF} -in1 SpectraFilter_1_output.mzData.tmp -in2 SpectraFilter_1_output.mzData)
add_test("TOPP_SpectraFilter_2" ${TOPP_BIN_PATH}/SpectraFilter -in SpectraFilter_2_input.mzData -out SpectraFilter_2_output.mzData.tmp -type WindowMower)
add_test("TOPP_SpectraFilter_2_out1" ${DIFF} -in1 SpectraFilter_2_output.mzData.tmp -in2 SpectraFilter_2_output.mzData)
add_test("TOPP_SpectraFilter_3" ${TOPP_BIN_PATH}/SpectraFilter -in SpectraFilter_3_input.mzData -out SpectraFilter_3_output.mzData.tmp -ini SpectraFilter_3_parameters.ini -type WindowMower)
add_test("TOPP_SpectraFilter_3_out1" ${DIFF} -in1 SpectraFilter_3_output.mzData.tmp -in2 SpectraFilter_3_output.mzData)

### InternalCalibration tests
add_test("TOPP_InternalCalibration_1" ${TOPP_BIN_PATH}/InternalCalibration -in InternalCalibration_1_input.mzData -out InternalCalibration_1_output.tmp -ref_masses InternalCalibration_ref_masses -ini InternalCalibration_parameters.ini)
add_test("TOPP_InternalCalibration_1_out1" ${DIFF} -in1 InternalCalibration_1_output.tmp -in2 InternalCalibration_1_output.mzData)
add_test("TOPP_InternalCalibration_2" ${TOPP_BIN_PATH}/InternalCalibration -in InternalCalibration_2_input.mzData -out InternalCalibration_2_output.tmp -ref_masses InternalCalibration_ref_masses -peak_data)
add_test("TOPP_InternalCalibration_2_out1" ${DIFF} -in1 InternalCalibration_2_output.tmp -in2 InternalCalibration_2_output.mzData)

### TOFCalibration tests
add_test("TOPP_TOFCalibration_1" ${TOPP_BIN_PATH}/TOFCalibration -in TOFCalibration_1_input.mzData -out TOFCalibration_1_output.tmp -ref_masses TOFCalibration_ref_masses -ini TOFCalibration_parameters.ini -tof_const TOFCalibration_const -ext_calibrants TOFCalibration_1_calibrants.mzData)
add_test("TOPP_TOFCalibration_1_out1" ${DIFF} -in1 TOFCalibration_1_output.tmp -in2 TOFCalibration_1_output.mzData)
add_test("TOPP_TOFCalibration_2" ${TOPP_BIN_PATH}/TOFCalibration -in TOFCalibration_2_input.mzData -out TOFCalibration_2_output.tmp -ref_masses TOFCalibration_ref_masses -ini TOFCalibration_parameters.ini -tof_const TOFCalibration_const -ext_calibrants TOFCalibration_2_calibrants.mzData -peak_data)
add_test("TOPP_TOFCalibration_2_out1" ${DIFF} -in1 TOFCalibration_2_output.tmp -in2 TOFCalibration_2_output.mzData)

### TextExporter tests
add_test("TOPP_TextExporter_1" ${TOPP_BIN_PATH}/TextExporter -in TextExporter_1_input.featureXML -no_progress -out TextExporter.tmp)
add_test("TOPP_TextExporter_1_out1" ${DIFF} -in1 TextExporter.tmp -in2 TextExporter_1_output.txt)
#TODO add_test("TOPP_TextExporter_2" ${TOPP_BIN_PATH}/TextExporter -ini TextExporter_2_parameters.ini)
##SED('s:\(.*extracted\ from\ \).*$$:\1\ <omitted>:g' TextExporter_2_consensus_centroids.tmp > TextExporter_2_consensus_centroids.tmp_
#add_test("TOPP_TextExporter_2_out1" ${DIFF} -in1 TextExporter_2_consensus_centroids.tmp_ -in2 TextExporter_2_consensus_centroids.wsv)
##SED('s:\(.*extracted\ from\ \).*$$:\1\ <omitted>:g' TextExporter_2_consensus_elements.tmp > TextExporter_2_consensus_elements.tmp_
#add_test("TOPP_TextExporter_2_out2" ${DIFF} -in1 TextExporter_2_consensus_elements.tmp_ -in2 TextExporter_2_consensus_elements.wsv)
##SED('s:\(.*extracted\ from\ \).*$$:\1\ <omitted>:g' TextExporter_2_consensus_features.tmp > TextExporter_2_consensus_features.tmp_
#add_test("TOPP_TextExporter_2_out3" ${DIFF} -in1 TextExporter_2_consensus_features.tmp_ -in2 TextExporter_2_consensus_features.wsv)
add_test("TOPP_TextExporter_3" ${TOPP_BIN_PATH}/TextExporter -in TextExporter_3_input.idXML -no_progress -out TextExporter3.tmp)
add_test("TOPP_TextExporter_3_out1" ${DIFF} -in1 TextExporter3.tmp -in2 TextExporter_3_output.txt)
add_test("TOPP_TextExporter_4" ${TOPP_BIN_PATH}/TextExporter -in TextExporter_3_input.idXML -no_progress -out TextExporter3.tmp -proteins_only)
add_test("TOPP_TextExporter_4_out1" ${DIFF} -in1 TextExporter3.tmp -in2 TextExporter_3_output_proteins.txt)
add_test("TOPP_TextExporter_5" ${TOPP_BIN_PATH}/TextExporter -in TextExporter_3_input.idXML -no_progress -out TextExporter3.tmp -peptides_only)
add_test("TOPP_TextExporter_5_out1" ${DIFF} -in1 TextExporter3.tmp -in2 TextExporter_3_output_peptides.txt)
add_test("TOPP_TextExporter_6" ${TOPP_BIN_PATH}/TextExporter -in TextExporter_4_input.idXML -no_progress -out TextExporter4.tmp -peptides_only_csv)
add_test("TOPP_TextExporter_6_out1" ${DIFF} -in1 TextExporter4.tmp -in2 TextExporter_4_output.txt)
add_test("TOPP_TextExporter_7" ${TOPP_BIN_PATH}/TextExporter -in TextExporter_5_input.idXML -no_progress -out TextExporter5.tmp -peptides_only_csv -first_dim_rt)
add_test("TOPP_TextExporter_7_out1" ${DIFF} -in1 TextExporter5.tmp -in2 TextExporter_5_output.txt)

### SILACAnalyzer tests
add_test("TOPP_SILACAnalyzer_1" ${TOPP_BIN_PATH}/SILACAnalyzer -ini SILACAnalyzer_1_parameters.ini -in SILACAnalyzer_1_input.mzData -out SILACAnalyzer_1_output.tmp)
add_test("TOPP_SILACAnalyzer_1_out1" ${DIFF} -in1 SILACAnalyzer_1_output.tmp -in2 SILACAnalyzer_1_output.consensusXML)

### FeatureFinderMRM tests
add_test("TOPP_FeatureFinderMRM_1" ${TOPP_BIN_PATH}/FeatureFinderMRM -in FeatureFinderMRM_1_input.mzData -ini FeatureFinderMRM_1.ini -out FeatureFinderMRM.tmp)
add_test("TOPP_FeatureFinderMRM_1_out1" ${DIFF} -in1 FeatureFinderMRM.tmp -in2 FeatureFinderMRM_1_output.featureXML)
add_test("TOPP_FeatureFinderMRM_2" ${TOPP_BIN_PATH}/FeatureFinderMRM -in FeatureFinderMRM_1_input.mzData -ini FeatureFinderMRM_1.ini -d -out FeatureFinderMRM.tmp)
add_test("TOPP_FeatureFinderMRM_2_out1" ${DIFF} -in1 FeatureFinderMRM.tmp -in2 FeatureFinderMRM_1_output.featureXML)
add_test("TOPP_FeatureFinderMRM_2_out2" ${DIFF} -in1 precursor_461.743 -in2 FeatureFinderMRM_2_output.txt)
#REMOVE(-rf precursor_461.743

### FeatureLinker tests
add_test("TOPP_FeatureLinker_1" ${TOPP_BIN_PATH}/FeatureLinker -ini FeatureLinker_1_parameters.ini -in FeatureLinker_1_input.featureXML)
add_test("TOPP_FeatureLinker_1_out1" ${DIFF} -in1 FeatureLinker_1_output.tmp -in2 FeatureLinker_1_output.consensusXML)
add_test("TOPP_FeatureLinker_2" ${TOPP_BIN_PATH}/FeatureLinker -ini FeatureLinker_2_parameters.ini)
add_test("TOPP_FeatureLinker_2_out1" ${DIFF} -in1 FeatureLinker_2_output.tmp -in2 FeatureLinker_2_output.consensusXML)
add_test("TOPP_FeatureLinker_3" ${TOPP_BIN_PATH}/FeatureLinker -ini FeatureLinker_3_parameters.ini)
add_test("TOPP_FeatureLinker_3_out1" ${DIFF} -in1 FeatureLinker_3_output.tmp -in2 FeatureLinker_3_output.consensusXML)
add_test("TOPP_FeatureLinker_4" ${TOPP_BIN_PATH}/FeatureLinker -ini FeatureLinker_4_parameters.ini)
add_test("TOPP_FeatureLinker_4_out1" ${DIFF} -in1 FeatureLinker_4_output.tmp -in2 FeatureLinker_4_output.consensusXML)
add_test("TOPP_FeatureLinker_5" ${TOPP_BIN_PATH}/FeatureLinker -ini FeatureLinker_5_parameters.ini -in FeatureLinker_5_input.featureXML)
add_test("TOPP_FeatureLinker_5_out1" ${DIFF} -in1 FeatureLinker_5_output.tmp -in2 FeatureLinker_5_output.consensusXML)

### IDMapper tests
add_test("TOPP_IDMapper_1" ${TOPP_BIN_PATH}/IDMapper -in IDMapper_1_input.featureXML -out IDMapper_1_output.tmp -id IDMapper_1_input.idXML)
add_test("TOPP_IDMapper_1_out1" ${DIFF} -in1 IDMapper_1_output.tmp -in2 IDMapper_1_output.featureXML)
add_test("TOPP_IDMapper_2" ${TOPP_BIN_PATH}/IDMapper -in IDMapper_2_input.consensusXML -out IDMapper_2_output.tmp -id IDMapper_2_input.idXML)
add_test("TOPP_IDMapper_2_out1" ${DIFF} -in1 IDMapper_2_output.tmp -in2 IDMapper_2_output.consensusXML)

### ConsensusID tests
#TODO add_test("TOPP_ConsensusID_1" ${TOPP_BIN_PATH}/ConsensusID -ini ConsensusID_1_parameters.ini -out ConsensusID_1_output.tmp_)
##SED('s/^.*IdentificationRun date.*$$/#Run removed because of date/g' ConsensusID_1_output.tmp_ > ConsensusID_1_output.tmp
#add_test("TOPP_ConsensusID_1_out1" ${DIFF} -in1 ConsensusID_1_output.tmp -in2 ConsensusID_1_output.idXML)
#add_test("TOPP_ConsensusID_2" ${TOPP_BIN_PATH}/ConsensusID -ini ConsensusID_2_parameters.ini -out ConsensusID_2_output.tmp_ -in ConsensusID_2_input.featureXML)
##SED('s/^.*IdentificationRun id.*$$/#Run removed because of date/g' ConsensusID_2_output.tmp_ > ConsensusID_2_output.tmp
#add_test("TOPP_ConsensusID_2_out1" ${DIFF} -in1 ConsensusID_2_output.tmp -in2 ConsensusID_2_output.featureXML)
#add_test("TOPP_ConsensusID_3" ${TOPP_BIN_PATH}/ConsensusID -ini ConsensusID_3_parameters.ini -out ConsensusID_3_output.tmp_ -in ConsensusID_3_input.consensusXML)
##SED('s/^.*IdentificationRun id.*$$/#Run removed because of date/g' ConsensusID_3_output.tmp_ > ConsensusID_3_output.tmp
#add_test("TOPP_ConsensusID_3_out1" ${DIFF} -in1 ConsensusID_3_output.tmp -in2 ConsensusID_3_output.consensusXML)
