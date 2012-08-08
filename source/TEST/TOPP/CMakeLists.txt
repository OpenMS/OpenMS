#
# README:
#
# This file contains the declarations of all TOPP/UTILS tests that are executed by CTest
#
# In order to add a test do the following:
# - add the name of the tool to the ALL_WRITEINI_TOOLS list (see below), in order to check -write_ini functionality
# - add specific tests for your tool (this usually involves FuzzyDiff) - see below for examples
# - for UTILS we do not require any tests, but you are encouraged to provide them!
#
# Test names (the first argument of 'add_test') must be distinct.
# If not, the first definition of the test is overwritten!
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
# When the result file of a test is checked with FuzzyDiff, we need to set the dependency to the test
# that generates the output file by
#
# set_tests_properties("TOPP_<tool>_<num>" PROPERTIES DEPENDS "TOPP_<tool>_<num>_out<num2>")
#
# Warning: When passing parameters avoid usage of "<" and ">" signs, as this will confuse CDASH when
#          interpreting the test results file (which happens to be XML)
#
# Hint: Most TOPP tools add DataProcessing information to their output.
#       The -test flag ensures that this output is independent of local path settings etc, so use it when *testing* TOPP tools!
# Hint: Provide inidividual names for each test output file e.g., the first test for the FileConverter should have FileConverter_1.tmp
#       as output file.

### General settings
set(TOPP_BIN_PATH ${OPENMS_BINARY_DIR})
set(DATA_DIR_TOPP ${PROJECT_SOURCE_DIR}/TOPP)
set(DATA_DIR_TOPP_BIN ${PROJECT_BINARY_DIR}/TOPP)
set(DIFF ${TOPP_BIN_PATH}/FuzzyDiff -ini ${DATA_DIR_TOPP}/FuzzyDiff.ini)

## configure resource files (for TOPPAS/ExecutePipeline) - there is no other way to specify the input
set (CONFIGURED_CONFIG_H ${PROJECT_BINARY_DIR}/include/OpenMS/config.h)
configure_file(${DATA_DIR_TOPP}/ExecutePipeline_1.trf.in ${DATA_DIR_TOPP_BIN}/ExecutePipeline_1.trf)


### test if -write_wdsl works
add_test("TOPP_WSDL_1" ${TOPP_BIN_PATH}/AdditiveSeries -write_wsdl WSDL_out.tmp)
add_test("TOPP_WSDL_1_out1" ${DIFF} -in1 WSDL_out.tmp -in2 ${DATA_DIR_TOPP}/WSDL_out.wsdl)

#####
## one test for each TOPP tool, checking -write_ini functionality:
#####

set(TOPP_executables_INI ${TOPP_executables})
## GenericWrapper (has a -type; special treatment below)
list(REMOVE_ITEM TOPP_executables_INI GenericWrapper)

foreach (i ${TOPP_executables_INI})
  add_test("TOPPWRITEINI_${i}" ${TOPP_BIN_PATH}/${i} -write_ini "${DATA_DIR_TOPP_BIN}/${i}.ini.tmp")
	add_test("TOPPWRITEINI_${i}_SectionName" ${CMAKE_COMMAND} -DNAME=${i} -DFILE=${DATA_DIR_TOPP_BIN}/${i}.ini.tmp -P ${DATA_DIR_TOPP}/check_ini.cmake)
  set_tests_properties("TOPPWRITEINI_${i}_SectionName" PROPERTIES DEPENDS "TOPPWRITEINI_${i}")
endforeach(i)
add_test("TOPPWRITEINI_GenericWrapper" ${TOPP_BIN_PATH}/GenericWrapper -type RAWFileConvert -write_ini "${DATA_DIR_TOPP_BIN}/GenericWrapper.ini.tmp")


##############################################################################################################################################################################
###### TOPP tests  ############
##############################################################################################################################################################################

### AdditiveSeries tests
add_test("TOPP_AdditiveSeries_1" ${TOPP_BIN_PATH}/AdditiveSeries -in ${DATA_DIR_TOPP}/AdditiveSeries_1_feat11.featureXML ${DATA_DIR_TOPP}/AdditiveSeries_1_feat20.featureXML ${DATA_DIR_TOPP}/AdditiveSeries_1_feat22.featureXML ${DATA_DIR_TOPP}/AdditiveSeries_1_feat34.featureXML -ini ${DATA_DIR_TOPP}/AdditiveSeries_1_parameters.ini)
add_test("TOPP_AdditiveSeries_1_out1" ${DIFF} -in1 AdditiveSeries_1_output.tmp -in2 ${DATA_DIR_TOPP}/AdditiveSeries_1_output.txt )
set_tests_properties("TOPP_AdditiveSeries_1_out1" PROPERTIES DEPENDS "TOPP_AdditiveSeries_1")
add_test("TOPP_AdditiveSeries_1_out2" ${DIFF} -in1 AdditiveSeries_1_gnuplot_tmp.cmd -in2 ${DATA_DIR_TOPP}/AdditiveSeries_1_gnuplot.cmd )
set_tests_properties("TOPP_AdditiveSeries_1_out2" PROPERTIES DEPENDS "TOPP_AdditiveSeries_1")
add_test("TOPP_AdditiveSeries_1_out3" ${DIFF} -in1 AdditiveSeries_1_gnuplot_tmp.dat -in2 ${DATA_DIR_TOPP}/AdditiveSeries_1_gnuplot.dat )
set_tests_properties("TOPP_AdditiveSeries_1_out3" PROPERTIES DEPENDS "TOPP_AdditiveSeries_1")
add_test("TOPP_AdditiveSeries_1_out4" ${DIFF} -in1 AdditiveSeries_1_gnuplot_tmp.err -in2 ${DATA_DIR_TOPP}/AdditiveSeries_1_gnuplot.err )
set_tests_properties("TOPP_AdditiveSeries_1_out4" PROPERTIES DEPENDS "TOPP_AdditiveSeries_1")

### IDMerger tests
add_test("TOPP_IDMerger_1" ${TOPP_BIN_PATH}/IDMerger -test -in ${DATA_DIR_TOPP}/IDMerger_1_input1.idXML ${DATA_DIR_TOPP}/IDMerger_1_input2.idXML -out IDMerger_1_output.tmp)
add_test("TOPP_IDMerger_1_out1" ${DIFF} -in1 IDMerger_1_output.tmp -in2 ${DATA_DIR_TOPP}/IDMerger_1_output.idXML )
set_tests_properties("TOPP_IDMerger_1_out1" PROPERTIES DEPENDS "TOPP_IDMerger_1")
add_test("TOPP_IDMerger_2" ${TOPP_BIN_PATH}/IDMerger -test -pepxml_protxml -in ${DATA_DIR_TOPP}/IDMerger_2_input1.idXML ${DATA_DIR_TOPP}/IDMerger_2_input2.idXML -out IDMerger_2_output.tmp)
add_test("TOPP_IDMerger_2_out1" ${DIFF} -in1 IDMerger_2_output.tmp -in2 ${DATA_DIR_TOPP}/IDMerger_2_output.idXML )
set_tests_properties("TOPP_IDMerger_2_out1" PROPERTIES DEPENDS "TOPP_IDMerger_2")
# conflicting IDs:
add_test("TOPP_IDMerger_3" ${TOPP_BIN_PATH}/IDMerger -test -in ${DATA_DIR_TOPP}/IDMerger_1_input1.idXML ${DATA_DIR_TOPP}/IDMerger_1_input1.idXML -out IDMerger_3_output.tmp)
add_test("TOPP_IDMerger_3_out1" ${DIFF} -in1 IDMerger_3_output.tmp -in2 ${DATA_DIR_TOPP}/IDMerger_3_output.idXML )
set_tests_properties("TOPP_IDMerger_3_out1" PROPERTIES DEPENDS "TOPP_IDMerger_3")
# empty IDs:
add_test("TOPP_IDMerger_4" ${TOPP_BIN_PATH}/IDMerger -test -in ${DATA_DIR_TOPP}/../data/degenerate_cases/empty.idXML ${DATA_DIR_TOPP}/../data/degenerate_cases/empty.idXML -out IDMerger_4_output.tmp)
add_test("TOPP_IDMerger_4_out1" ${DIFF} -in1 IDMerger_4_output.tmp -in2 ${DATA_DIR_TOPP}/IDMerger_4_output.idXML )
set_tests_properties("TOPP_IDMerger_4_out1" PROPERTIES DEPENDS "TOPP_IDMerger_4")

### BaselineFilter tests
add_test("TOPP_BaselineFilter_1" ${TOPP_BIN_PATH}/BaselineFilter -test -in ${DATA_DIR_TOPP}/BaselineFilter_input.mzML -out BaselineFilter.tmp -struc_elem_length 1.5)
add_test("TOPP_BaselineFilter_1_out1" ${DIFF} -in1 BaselineFilter.tmp -in2 ${DATA_DIR_TOPP}/BaselineFilter_output.mzML )
set_tests_properties("TOPP_BaselineFilter_1_out1" PROPERTIES DEPENDS "TOPP_BaselineFilter_1")

### ConsensusMapNormalizer tests
add_test("TOPP_ConsensusMapNormalizer_1" ${TOPP_BIN_PATH}/ConsensusMapNormalizer -test -in ${DATA_DIR_TOPP}/ConsensusMapNormalizer_input.consensusXML -out ConsensusMapNormalizer.tmp)
add_test("TOPP_ConsensusMapNormalizer_1_out1" ${DIFF} -in1 ConsensusMapNormalizer.tmp -in2 ${DATA_DIR_TOPP}/ConsensusMapNormalizer_output.consensusXML )
set_tests_properties("TOPP_ConsensusMapNormalizer_1_out1" PROPERTIES DEPENDS "TOPP_ConsensusMapNormalizer_1")

### MapNormalizer tests
add_test("TOPP_MapNormalizer_1" ${TOPP_BIN_PATH}/MapNormalizer -test -in ${DATA_DIR_TOPP}/MapNormalizer_input.mzML -out MapNormalizer.tmp)
add_test("TOPP_MapNormalizer_1_out1" ${DIFF} -in1 MapNormalizer.tmp -in2 ${DATA_DIR_TOPP}/MapNormalizer_output.mzML )
set_tests_properties("TOPP_MapNormalizer_1_out1" PROPERTIES DEPENDS "TOPP_MapNormalizer_1")

### DBImporter tests
if (DB_TEST)
	add_test("TOPP_DBImporter_1" ${TOPP_BIN_PATH}/DBImporter -init -ini ${DATA_DIR_TOPP_BIN}/DBImporter.ini)
  add_test("TOPP_DBImporter_2" ${TOPP_BIN_PATH}/DBImporter -test -in ${DATA_DIR_TOPP}/DBImporter_input.mzML -ini ${DATA_DIR_TOPP_BIN}/DBImporter.ini)
endif()

### DBExporter tests
if (DB_TEST)
	add_test("TOPP_DBExporter_1" ${TOPP_BIN_PATH}/DBExporter -test -out DBExporter.tmp -ini ${DATA_DIR_TOPP_BIN}/DBExporter.ini -id 1)
  add_test("TOPP_DBExporter_1_out1" ${DIFF} -in1 DBExporter.tmp -in2 ${DATA_DIR_TOPP}/DBExporter_output.mzML )
set_tests_properties("TOPP_DBExporter_1_out1" PROPERTIES DEPENDS "TOPP_DBExporter_1")
	add_test("TOPP_DBExporter_2" ${TOPP_BIN_PATH}/DBExporter -test -out DBExporter.tmp -ini ${DATA_DIR_TOPP_BIN}/DBExporter.ini -query "SELECT 1")
  add_test("TOPP_DBExporter_2_out1" ${DIFF} -in1 DBExporter.tmp -in2 ${DATA_DIR_TOPP}/DBExporter_output.mzML )
set_tests_properties("TOPP_DBExporter_2_out1" PROPERTIES DEPENDS "TOPP_DBExporter_2")
endif()

### Decharger tests
add_test("TOPP_Decharger_1" ${TOPP_BIN_PATH}/Decharger -test -in ${DATA_DIR_TOPP}/Decharger_input.featureXML -ini ${DATA_DIR_TOPP}/Decharger_input.ini -out_fm Decharger_output_fm.tmp -out_cm Decharger_output.tmp -outpairs Decharger_p_output.tmp)
add_test("TOPP_Decharger_1_out1" ${DIFF} -in1 Decharger_output.tmp -in2 ${DATA_DIR_TOPP}/Decharger_output.consensusXML )
set_tests_properties("TOPP_Decharger_1_out1" PROPERTIES DEPENDS "TOPP_Decharger_1")

### DTAExtractor tests
add_test("TOPP_DTAExtractor_1" ${TOPP_BIN_PATH}/DTAExtractor -in ${DATA_DIR_TOPP}/DTAExtractor_1_input.mzML -out DTAExtractor -rt :61)
add_test("TOPP_DTAExtractor_1_out1" ${DIFF} -in1 DTAExtractor_RT60.dta -in2 ${DATA_DIR_TOPP}/DTAExtractor_1_output.dta )
set_tests_properties("TOPP_DTAExtractor_1_out1" PROPERTIES DEPENDS "TOPP_DTAExtractor_1")
add_test("TOPP_DTAExtractor_2" ${TOPP_BIN_PATH}/DTAExtractor -in ${DATA_DIR_TOPP}/DTAExtractor_1_input.mzML -out DTAExtractor -level 1)
add_test("TOPP_DTAExtractor_2_out1" ${DIFF} -in1 DTAExtractor_RT60.dta -in2 ${DATA_DIR_TOPP}/DTAExtractor_2_output.dta )
set_tests_properties("TOPP_DTAExtractor_2_out1" PROPERTIES DEPENDS "TOPP_DTAExtractor_2")
add_test("TOPP_DTAExtractor_3" ${TOPP_BIN_PATH}/DTAExtractor -in ${DATA_DIR_TOPP}/DTAExtractor_1_input.mzML -out DTAExtractor -level 2 -mz :1000)
add_test("TOPP_DTAExtractor_3_out1" ${DIFF} -in1 DTAExtractor_RT140_MZ5.dta -in2 ${DATA_DIR_TOPP}/DTAExtractor_3_output.dta )
set_tests_properties("TOPP_DTAExtractor_3_out1" PROPERTIES DEPENDS "TOPP_DTAExtractor_3")

### MassTraceExtractor test
add_test("TOPP_MassTraceExtractor_1" ${TOPP_BIN_PATH}/MassTraceExtractor -test -ini ${DATA_DIR_TOPP}/MassTraceExtractor.ini -in ${DATA_DIR_TOPP}/MassTraceExtractor_1_input.mzML -out MassTraceExtractor_1.tmp)
add_test("TOPP_MassTraceExtractor_1_out" ${DIFF} -in1 MassTraceExtractor_1.tmp -in2 ${DATA_DIR_TOPP}/MassTraceExtractor_1_output.featureXML )
set_tests_properties("TOPP_MassTraceExtractor_1_out" PROPERTIES DEPENDS "TOPP_MassTraceExtractor_1")


### FeatureFinderCentroided test
add_test("TOPP_FeatureFinderCentroided_1" ${TOPP_BIN_PATH}/FeatureFinderCentroided -test -ini ${DATA_DIR_TOPP}/FeatureFinderCentroided_1_parameters.ini -in ${DATA_DIR_TOPP}/FeatureFinderCentroided_1_input.mzML -out FeatureFinderCentroided_1.tmp)
add_test("TOPP_FeatureFinderCentroided_1_out1" ${DIFF} -in1 FeatureFinderCentroided_1.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinderCentroided_1_output.featureXML )
set_tests_properties("TOPP_FeatureFinderCentroided_1_out1" PROPERTIES DEPENDS "TOPP_FeatureFinderCentroided_1")
### FeatureFinderMRM test
add_test("TOPP_FeatureFinderMRM_1" ${TOPP_BIN_PATH}/FeatureFinderMRM -test -in ${DATA_DIR_TOPP}/FeatureFinderMRM_1_input.mzML -ini ${DATA_DIR_TOPP}/FeatureFinderMRM_1_parameters.ini -out FeatureFinderMRM_1.tmp)
add_test("TOPP_FeatureFinderMRM_1_out1" ${DIFF} -in1 FeatureFinderMRM_1.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinderMRM_1_output.featureXML )
set_tests_properties("TOPP_FeatureFinderMRM_1_out1" PROPERTIES DEPENDS "TOPP_FeatureFinderMRM_1")

### FileConverter tests
add_test("TOPP_FileConverter_1" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_1_input.mzData -out FileConverter_1.tmp -out_type mzML)
add_test("TOPP_FileConverter_1_out1" ${DIFF} -in1 FileConverter_1.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_1_output.mzML )
set_tests_properties("TOPP_FileConverter_1_out1" PROPERTIES DEPENDS "TOPP_FileConverter_1")
add_test("TOPP_FileConverter_2" ${TOPP_BIN_PATH}/FileConverter -test  -in ${DATA_DIR_TOPP}/FileConverter_2_input.txt -in_type DTA2D -out FileConverter_2.tmp -out_type mzML)
add_test("TOPP_FileConverter_2_out1" ${DIFF} -in1 FileConverter_2.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_2_output.mzML )
set_tests_properties("TOPP_FileConverter_2_out1" PROPERTIES DEPENDS "TOPP_FileConverter_2")
add_test("TOPP_FileConverter_3" ${TOPP_BIN_PATH}/FileConverter -test  -in ${DATA_DIR_TOPP}/FileConverter_3_input.featureXML -out FileConverter_3.tmp -out_type mzML)
add_test("TOPP_FileConverter_3_out1" ${DIFF} -in1 FileConverter_3.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_3_output.mzML )
set_tests_properties("TOPP_FileConverter_3_out1" PROPERTIES DEPENDS "TOPP_FileConverter_3")
add_test("TOPP_FileConverter_4" ${TOPP_BIN_PATH}/FileConverter -test  -in ${DATA_DIR_TOPP}/FileConverter_4_input.mzXML -out FileConverter_4.tmp -out_type mzML)
add_test("TOPP_FileConverter_4_out1" ${DIFF} -in1 FileConverter_4.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_4_output.mzML )
set_tests_properties("TOPP_FileConverter_4_out1" PROPERTIES DEPENDS "TOPP_FileConverter_4")
add_test("TOPP_FileConverter_6" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_6_input.mzML -out FileConverter_6.tmp -out_type mzXML)
add_test("TOPP_FileConverter_6_out1" ${DIFF} -in1 FileConverter_6.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_6_output.mzXML )
set_tests_properties("TOPP_FileConverter_6_out1" PROPERTIES DEPENDS "TOPP_FileConverter_6")
add_test("TOPP_FileConverter_7" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_7_input.consensusXML -out FileConverter_7.tmp -out_type featureXML)
add_test("TOPP_FileConverter_7_out1" ${DIFF} -in1 FileConverter_7.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_7_output.featureXML )
set_tests_properties("TOPP_FileConverter_7_out1" PROPERTIES DEPENDS "TOPP_FileConverter_7")
add_test("TOPP_FileConverter_8" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_8_input.mzML -out FileConverter_8.tmp -out_type mzData)
add_test("TOPP_FileConverter_8_out1" ${DIFF} -in1 FileConverter_8.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_8_output.mzData )
set_tests_properties("TOPP_FileConverter_8_out1" PROPERTIES DEPENDS "TOPP_FileConverter_8")
add_test("TOPP_FileConverter_9" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_9_input.consensusXML -out FileConverter_9.tmp -out_type featureXML)
add_test("TOPP_FileConverter_9_out1" ${DIFF} -in1 FileConverter_9.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_9_output.featureXML )
set_tests_properties("TOPP_FileConverter_9_out1" PROPERTIES DEPENDS "TOPP_FileConverter_9")
add_test("TOPP_FileConverter_10" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_10_input.edta -no_progress -out FileConverter_10.tmp -out_type featureXML)
add_test("TOPP_FileConverter_10_out1" ${DIFF} -in1 FileConverter_10.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_10_output.featureXML )
set_tests_properties("TOPP_FileConverter_10_out1" PROPERTIES DEPENDS "TOPP_FileConverter_10")
#SpecArray
add_test("TOPP_FileConverter_11" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_11_input.pepList -no_progress -out FileConverter_11.tmp -out_type featureXML)
add_test("TOPP_FileConverter_11_out1" ${DIFF} -in1 FileConverter_11.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_11_output.featureXML )
set_tests_properties("TOPP_FileConverter_11_out1" PROPERTIES DEPENDS "TOPP_FileConverter_11")
#msInspect
add_test("TOPP_FileConverter_12" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_12_input.peptides.tsv -no_progress -out FileConverter_12.tmp -out_type featureXML)
add_test("TOPP_FileConverter_12_out1" ${DIFF} -in1 FileConverter_12.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_12_output.featureXML )
set_tests_properties("TOPP_FileConverter_12_out1" PROPERTIES DEPENDS "TOPP_FileConverter_12")
#kroenik
add_test("TOPP_FileConverter_13" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_13_input.peptides.kroenik -no_progress -out FileConverter_13.tmp -out_type featureXML)
add_test("TOPP_FileConverter_13_out1" ${DIFF} -in1 FileConverter_13.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_13_output.featureXML )
set_tests_properties("TOPP_FileConverter_13_out1" PROPERTIES DEPENDS "TOPP_FileConverter_13")
# featureXML to consensusXML
add_test("TOPP_FileConverter_14" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_9_output.featureXML -no_progress -out FileConverter_14.tmp -out_type consensusXML)
add_test("TOPP_FileConverter_14_out1" ${DIFF} -in1 FileConverter_14.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_14_output.consensusXML )
set_tests_properties("TOPP_FileConverter_14_out1" PROPERTIES DEPENDS "TOPP_FileConverter_14")
# EDTA to consensusXML
add_test("TOPP_FileConverter_15" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_10_input.edta -no_progress -out FileConverter_15.tmp -out_type consensusXML)
add_test("TOPP_FileConverter_15_out1" ${DIFF} -in1 FileConverter_15.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_15_output.consensusXML )
set_tests_properties("TOPP_FileConverter_15_out1" PROPERTIES DEPENDS "TOPP_FileConverter_15")
# consensus EDTA to consensusXML
add_test("TOPP_FileConverter_16" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_16_input.edta -no_progress -out FileConverter_16.tmp -out_type consensusXML)
add_test("TOPP_FileConverter_16_out1" ${DIFF} -in1 FileConverter_16.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_16_output.consensusXML )
set_tests_properties("TOPP_FileConverter_16_out1" PROPERTIES DEPENDS "TOPP_FileConverter_16")

### FileFilter tests
add_test("TOPP_FileFilter_1" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out FileFilter_1.tmp -rt :30 -mz :1000 -int :20000 -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_1_out1" ${DIFF} -in1 FileFilter_1.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_1_output.mzML )
set_tests_properties("TOPP_FileFilter_1_out1" PROPERTIES DEPENDS "TOPP_FileFilter_1")
add_test("TOPP_FileFilter_2" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out FileFilter_2.tmp -rt 30: -mz 1000: -int 100: -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_2_out1" ${DIFF} -in1 FileFilter_2.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_2_output.mzML )
set_tests_properties("TOPP_FileFilter_2_out1" PROPERTIES DEPENDS "TOPP_FileFilter_2")
add_test("TOPP_FileFilter_3" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out FileFilter_3.tmp -level 2 -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_3_out1" ${DIFF} -in1 FileFilter_3.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_3_output.mzML )
set_tests_properties("TOPP_FileFilter_3_out1" PROPERTIES DEPENDS "TOPP_FileFilter_3")
add_test("TOPP_FileFilter_4" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_4_input.mzML -out FileFilter_4.tmp -remove_zoom -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_4_out1" ${DIFF} -in1 FileFilter_4.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_4_output.mzML )
set_tests_properties("TOPP_FileFilter_4_out1" PROPERTIES DEPENDS "TOPP_FileFilter_4")
add_test("TOPP_FileFilter_5" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_5_input.featureXML -out FileFilter_5.tmp -rt :1000 -mz :480 -int :79000 -charge :3 -q :0.6 -in_type featureXML -out_type featureXML)
add_test("TOPP_FileFilter_5_out1" ${DIFF} -in1 FileFilter_5.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_5_out.featureXML )
set_tests_properties("TOPP_FileFilter_5_out1" PROPERTIES DEPENDS "TOPP_FileFilter_5")
add_test("TOPP_FileFilter_6" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_5_input.featureXML -out FileFilter_6.tmp -rt 1000: -mz 440: -int 70000: -charge 3: -q 0.51: -in_type featureXML -out_type featureXML)
add_test("TOPP_FileFilter_6_out1" ${DIFF} -in1 FileFilter_6.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_6_out.featureXML )
set_tests_properties("TOPP_FileFilter_6_out1" PROPERTIES DEPENDS "TOPP_FileFilter_6")
add_test("TOPP_FileFilter_7" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_7_input.mzML -out FileFilter_7.tmp -int 7000: -level 1 2 3 -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_7_out1" ${DIFF} -in1 FileFilter_7.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_7_output.mzML )
set_tests_properties("TOPP_FileFilter_7_out1" PROPERTIES DEPENDS "TOPP_FileFilter_7")
add_test("TOPP_FileFilter_8" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_8_input.consensusXML -out FileFilter_8.tmp -rt 600:1400 -mz 700:2300 -int 1100:6000 -in_type consensusXML -out_type consensusXML)
add_test("TOPP_FileFilter_8_out1" ${DIFF} -in1 FileFilter_8.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_8_output.consensusXML )
set_tests_properties("TOPP_FileFilter_8_out1" PROPERTIES DEPENDS "TOPP_FileFilter_8")
add_test("TOPP_FileFilter_9" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_9_input.mzML -out FileFilter_9.tmp -remove_mode SelectedIonMonitoring -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_9_out1" ${DIFF} -in1 FileFilter_9.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_9_output.mzML)
set_tests_properties("TOPP_FileFilter_9_out1" PROPERTIES DEPENDS "TOPP_FileFilter_9")
add_test("TOPP_FileFilter_10" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_10_input.mzML -out FileFilter_10.tmp -remove_activation "Collision-induced dissociation" -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_10_out1" ${DIFF} -in1 FileFilter_10.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_10_output.mzML )
set_tests_properties("TOPP_FileFilter_10_out1" PROPERTIES DEPENDS "TOPP_FileFilter_10")
add_test("TOPP_FileFilter_11" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_11_input.mzML -out FileFilter_11.tmp -remove_activation "Plasma desorption" -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_11_out1" ${DIFF} -in1 FileFilter_11.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_11_output.mzML )
set_tests_properties("TOPP_FileFilter_11_out1" PROPERTIES DEPENDS "TOPP_FileFilter_11")
add_test("TOPP_FileFilter_12" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_12_input.mzML -out FileFilter_12.tmp -remove_chromatograms -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_12_out1" ${DIFF} -in1 FileFilter_12.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_12_output.mzML )
set_tests_properties("TOPP_FileFilter_12_out1" PROPERTIES DEPENDS "TOPP_FileFilter_12")
add_test("TOPP_FileFilter_13" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_13_input.consensusXML -out FileFilter_13.tmp -in_type consensusXML -out_type featureXML -map 2)
add_test("TOPP_FileFilter_13_out1" ${DIFF} -in1 FileFilter_13.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_13_output.featureXML )
set_tests_properties("TOPP_FileFilter_13_out1" PROPERTIES DEPENDS "TOPP_FileFilter_13")
add_test("TOPP_FileFilter_14" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_14_input.consensusXML -out FileFilter_14.tmp -in_type consensusXML -out_type consensusXML -map 0 2)
add_test("TOPP_FileFilter_14_out1" ${DIFF} -in1 FileFilter_14.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_14_output.consensusXML )
set_tests_properties("TOPP_FileFilter_14_out1" PROPERTIES DEPENDS "TOPP_FileFilter_14")
add_test("TOPP_FileFilter_15" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_15_input.featureXML -out FileFilter_15.tmp -id:sequences_whitelist Oxidation -id:remove_unassigned_ids)
add_test("TOPP_FileFilter_15_out1" ${DIFF} -in1 FileFilter_15.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_15_output.featureXML )
set_tests_properties("TOPP_FileFilter_15_out1" PROPERTIES DEPENDS "TOPP_FileFilter_15")
add_test("TOPP_FileFilter_16" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_15_input.featureXML -out FileFilter_16.tmp -id:sequences_whitelist Oxidation -id:remove_unassigned_ids -mz 400:600 -rt 3000:4000)
add_test("TOPP_FileFilter_16_out1" ${DIFF} -in1 FileFilter_16.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_16_output.featureXML )
set_tests_properties("TOPP_FileFilter_16_out1" PROPERTIES DEPENDS "TOPP_FileFilter_16")
add_test("TOPP_FileFilter_17" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_15_input.featureXML -out FileFilter_17.tmp -id:remove_annotated_features -mz 400:600 -rt 3000:4000)
add_test("TOPP_FileFilter_17_out1" ${DIFF} -in1 FileFilter_17.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_17_output.featureXML )
set_tests_properties("TOPP_FileFilter_17_out1" PROPERTIES DEPENDS "TOPP_FileFilter_17")
add_test("TOPP_FileFilter_18" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_18_input.consensusXML -out FileFilter_18.tmp -id:sequences_whitelist Oxidation -id:remove_unassigned_ids)
add_test("TOPP_FileFilter_18_out1" ${DIFF} -in1 FileFilter_18.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_18_output.consensusXML )
set_tests_properties("TOPP_FileFilter_18_out1" PROPERTIES DEPENDS "TOPP_FileFilter_18")
add_test("TOPP_FileFilter_19" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_18_input.consensusXML -out FileFilter_19.tmp -id:remove_unannotated_features -mz 400:600 -rt 3000:4000)
add_test("TOPP_FileFilter_19_out1" ${DIFF} -in1 FileFilter_19.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_19_output.consensusXML )
set_tests_properties("TOPP_FileFilter_19_out1" PROPERTIES DEPENDS "TOPP_FileFilter_19")
add_test("TOPP_FileFilter_20" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_15_input.featureXML -out FileFilter_20.tmp -id:accessions_whitelist YDL217C -id:remove_unassigned_ids)
add_test("TOPP_FileFilter_20_out1" ${DIFF} -in1 FileFilter_20.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_20_output.featureXML )
set_tests_properties("TOPP_FileFilter_20_out1" PROPERTIES DEPENDS "TOPP_FileFilter_20")
add_test("TOPP_FileFilter_21" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_15_input.featureXML -out FileFilter_21.tmp -id:remove_unassigned_ids -id:remove_unannotated_features -id:keep_best_score_id)
add_test("TOPP_FileFilter_21_out1" ${DIFF} -in1 FileFilter_21.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_21_output.featureXML )
set_tests_properties("TOPP_FileFilter_21_out1" PROPERTIES DEPENDS "TOPP_FileFilter_21")
add_test("TOPP_FileFilter_22" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_22_input.consensusXML -out FileFilter_22.tmp -remove_meta distinct_charges gt "1,2")
add_test("TOPP_FileFilter_22_out1" ${DIFF} -in1 FileFilter_22.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_22_output.consensusXML )
set_tests_properties("TOPP_FileFilter_22_out1" PROPERTIES DEPENDS "TOPP_FileFilter_22")
add_test("TOPP_FileFilter_23" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_22_input.consensusXML -out FileFilter_23.tmp -remove_meta distinct_charges_size gt 2)
add_test("TOPP_FileFilter_23_out1" ${DIFF} -in1 FileFilter_23.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_22_output.consensusXML )
set_tests_properties("TOPP_FileFilter_23_out1" PROPERTIES DEPENDS "TOPP_FileFilter_23")
add_test("TOPP_FileFilter_24" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_22_input.consensusXML -out FileFilter_24.tmp -remove_meta DOESNOTEXIST lt "whatever")
add_test("TOPP_FileFilter_24_out1" ${DIFF} -in1 FileFilter_24.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_24_output.consensusXML )
set_tests_properties("TOPP_FileFilter_24_out1" PROPERTIES DEPENDS "TOPP_FileFilter_24")


### FileInfo tests
add_test("TOPP_FileInfo_1" ${TOPP_BIN_PATH}/FileInfo -in ${DATA_DIR_TOPP}/FileInfo_1_input.txt -in_type DTA -no_progress -out FileInfo_1.tmp)
add_test("TOPP_FileInfo_1_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_1.tmp -in2 ${DATA_DIR_TOPP}/FileInfo_1_output.txt )
set_tests_properties("TOPP_FileInfo_1_out1" PROPERTIES DEPENDS "TOPP_FileInfo_1")
add_test("TOPP_FileInfo_2" ${TOPP_BIN_PATH}/FileInfo -in ${DATA_DIR_TOPP}/FileInfo_2_input.dta2d -no_progress -out FileInfo_2.tmp)
add_test("TOPP_FileInfo_2_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_2.tmp -in2 ${DATA_DIR_TOPP}/FileInfo_2_output.txt )
set_tests_properties("TOPP_FileInfo_2_out1" PROPERTIES DEPENDS "TOPP_FileInfo_2")
add_test("TOPP_FileInfo_3" ${TOPP_BIN_PATH}/FileInfo -in ${DATA_DIR_TOPP}/FileInfo_3_input.featureXML -m -s -p -no_progress -out FileInfo_3.tmp)
add_test("TOPP_FileInfo_3_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_3.tmp -in2 ${DATA_DIR_TOPP}/FileInfo_3_output.txt )
set_tests_properties("TOPP_FileInfo_3_out1" PROPERTIES DEPENDS "TOPP_FileInfo_3")
add_test("TOPP_FileInfo_4" ${TOPP_BIN_PATH}/FileInfo -in ${DATA_DIR_TOPP}/FileInfo_4_input.mzXML -m -no_progress -out FileInfo_4.tmp)
add_test("TOPP_FileInfo_4_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_4.tmp -in2 ${DATA_DIR_TOPP}/FileInfo_4_output.txt )
set_tests_properties("TOPP_FileInfo_4_out1" PROPERTIES DEPENDS "TOPP_FileInfo_4")
add_test("TOPP_FileInfo_5" ${TOPP_BIN_PATH}/FileInfo -in ${DATA_DIR_TOPP}/FileInfo_5_input.txt -in_type mzData -m -s -no_progress -out FileInfo_5.tmp)
add_test("TOPP_FileInfo_5_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_5.tmp -in2 ${DATA_DIR_TOPP}/FileInfo_5_output.txt )
set_tests_properties("TOPP_FileInfo_5_out1" PROPERTIES DEPENDS "TOPP_FileInfo_5")
add_test("TOPP_FileInfo_6" ${TOPP_BIN_PATH}/FileInfo -in ${DATA_DIR_TOPP}/FileInfo_6_input.mzData -d -s -no_progress -out FileInfo_6.tmp)
add_test("TOPP_FileInfo_6_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_6.tmp -in2 ${DATA_DIR_TOPP}/FileInfo_6_output.txt )
set_tests_properties("TOPP_FileInfo_6_out1" PROPERTIES DEPENDS "TOPP_FileInfo_6")
add_test("TOPP_FileInfo_7" ${TOPP_BIN_PATH}/FileInfo -in ${DATA_DIR_TOPP}/FileInfo_7_input.consensusXML -s -m -p -no_progress -out FileInfo_7.tmp)
add_test("TOPP_FileInfo_7_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_7.tmp -in2 ${DATA_DIR_TOPP}/FileInfo_7_output.txt )
set_tests_properties("TOPP_FileInfo_7_out1" PROPERTIES DEPENDS "TOPP_FileInfo_7")
add_test("TOPP_FileInfo_9" ${TOPP_BIN_PATH}/FileInfo -in ${DATA_DIR_TOPP}/FileInfo_9_input.mzML -m -p -s -no_progress -out FileInfo_9.tmp)
add_test("TOPP_FileInfo_9_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_9.tmp -in2 ${DATA_DIR_TOPP}/FileInfo_9_output.txt )
set_tests_properties("TOPP_FileInfo_9_out1" PROPERTIES DEPENDS "TOPP_FileInfo_9")
add_test("TOPP_FileInfo_10" ${TOPP_BIN_PATH}/FileInfo -in ${DATA_DIR_TOPP}/FileInfo_10_input.idXML -no_progress -out FileInfo_10.tmp)
add_test("TOPP_FileInfo_10_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_10.tmp -in2 ${DATA_DIR_TOPP}/FileInfo_10_output.txt )
set_tests_properties("TOPP_FileInfo_10_out1" PROPERTIES DEPENDS "TOPP_FileInfo_10")

### FileMerger tests
add_test("TOPP_FileMerger_1" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_1_input1.dta2d ${DATA_DIR_TOPP}/FileMerger_1_input2.dta ${DATA_DIR_TOPP}/FileMerger_1_input3.dta2d -out FileMerger_1.tmp -rt_custom 1 2 3)
add_test("TOPP_FileMerger_1_out1" ${DIFF} -in1 FileMerger_1.tmp -in2 ${DATA_DIR_TOPP}/FileMerger_1_output.mzML )
set_tests_properties("TOPP_FileMerger_1_out1" PROPERTIES DEPENDS "TOPP_FileMerger_1")
add_test("TOPP_FileMerger_2" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_2_input1.txt ${DATA_DIR_TOPP}/FileMerger_2_input2.txt -in_type DTA -out FileMerger_2.tmp -rt_custom 5 10)
add_test("TOPP_FileMerger_2_out1" ${DIFF} -in1 FileMerger_2.tmp -in2 ${DATA_DIR_TOPP}/FileMerger_2_output.mzML )
set_tests_properties("TOPP_FileMerger_2_out1" PROPERTIES DEPENDS "TOPP_FileMerger_2")
add_test("TOPP_FileMerger_3" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_3_input1.dta ${DATA_DIR_TOPP}/FileMerger_3_input2.dta -out FileMerger_3.tmp -rt_auto)
add_test("TOPP_FileMerger_3_out1" ${DIFF} -in1 FileMerger_3.tmp -in2 ${DATA_DIR_TOPP}/FileMerger_3_output.mzML )
set_tests_properties("TOPP_FileMerger_3_out1" PROPERTIES DEPENDS "TOPP_FileMerger_3")
add_test("TOPP_FileMerger_4" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_4_input1.dta2d ${DATA_DIR_TOPP}/FileMerger_4_input2.dta2d -out FileMerger_4.tmp)
add_test("TOPP_FileMerger_4_out1" ${DIFF} -in1 FileMerger_4.tmp -in2 ${DATA_DIR_TOPP}/FileMerger_4_output.mzML )
set_tests_properties("TOPP_FileMerger_4_out1" PROPERTIES DEPENDS "TOPP_FileMerger_4")
add_test("TOPP_FileMerger_5" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_5_input_rt1023.331.dta ${DATA_DIR_TOPP}/FileMerger_5_input_rt2044.334.dta ${DATA_DIR_TOPP}/FileMerger_5_input_rt889.32.dta -rt_filename -user_ms_level -ms_level 2 -out FileMerger_5.tmp)
add_test("TOPP_FileMerger_5_out1" ${DIFF} -in1 FileMerger_5.tmp -in2 ${DATA_DIR_TOPP}/FileMerger_5_output.mzML )
set_tests_properties("TOPP_FileMerger_5_out1" PROPERTIES DEPENDS "TOPP_FileMerger_5")
add_test("TOPP_FileMerger_6" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_6_input1.mzML ${DATA_DIR_TOPP}/FileMerger_6_input2.mzML -out FileMerger_6.tmp)
add_test("TOPP_FileMerger_6_out1" ${DIFF} -in1 FileMerger_6.tmp -in2 ${DATA_DIR_TOPP}/FileMerger_6_output.mzML )
set_tests_properties("TOPP_FileMerger_6_out1" PROPERTIES DEPENDS "TOPP_FileMerger_6")
add_test("TOPP_FileMerger_7" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_7_input1.featureXML ${DATA_DIR_TOPP}/FileMerger_7_input2.featureXML -out FileMerger_7.tmp)
add_test("TOPP_FileMerger_7_out1" ${DIFF} -in1 FileMerger_7.tmp -in2 ${DATA_DIR_TOPP}/FileMerger_7_output.featureXML )
set_tests_properties("TOPP_FileMerger_7_out1" PROPERTIES DEPENDS "TOPP_FileMerger_7")
add_test("TOPP_FileMerger_8" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_8_input1.consensusXML ${DATA_DIR_TOPP}/FileMerger_8_input2.consensusXML -out FileMerger_8.tmp)
add_test("TOPP_FileMerger_8_out1" ${DIFF} -in1 FileMerger_8.tmp -in2 ${DATA_DIR_TOPP}/FileMerger_8_output.consensusXML )
set_tests_properties("TOPP_FileMerger_8_out1" PROPERTIES DEPENDS "TOPP_FileMerger_8")
add_test("TOPP_FileMerger_9" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_9_input1.TraML ${DATA_DIR_TOPP}/FileMerger_9_input2.TraML -out FileMerger_9.tmp)
add_test("TOPP_FileMerger_9_out1" ${DIFF} -in1 FileMerger_9.tmp -in2 ${DATA_DIR_TOPP}/FileMerger_9_output.TraML )
set_tests_properties("TOPP_FileMerger_9_out1" PROPERTIES DEPENDS "TOPP_FileMerger_9")

### IDRTCalibration tests
add_test("TOPP_IDRTCalibration_1" ${TOPP_BIN_PATH}/IDRTCalibration -in ${DATA_DIR_TOPP}/IDRTCalibration_1_input.IdXML -out IDRTCalibration_1_output.tmp -calibrant_1_input 10 -calibrant_2_input 90)
add_test("TOPP_IDRTCalibration_1_out1" ${DIFF} -in1 IDRTCalibration_1_output.tmp -in2 ${DATA_DIR_TOPP}/IDRTCalibration_1_output.IdXML )
set_tests_properties("TOPP_IDRTCalibration_1_out1" PROPERTIES DEPENDS "TOPP_IDRTCalibration_1")
# test empty IDs
add_test("TOPP_IDRTCalibration_2" ${TOPP_BIN_PATH}/IDRTCalibration -in ${DATA_DIR_TOPP}/../data/degenerate_cases/empty.idXML -out IDRTCalibration_2_output.tmp -calibrant_1_input 10 -calibrant_2_input 90)
add_test("TOPP_IDRTCalibration_2_out1" ${DIFF} -in1 IDRTCalibration_2_output.tmp -in2 ${DATA_DIR_TOPP}/IDRTCalibration_2_output.IdXML )
set_tests_properties("TOPP_IDRTCalibration_2_out1" PROPERTIES DEPENDS "TOPP_IDRTCalibration_2")

### ITRAQAnalyzer tests
add_test("TOPP_ITRAQAnalyzer_1" ${TOPP_BIN_PATH}/ITRAQAnalyzer -test -in ${DATA_DIR_TOPP}/ITRAQAnalyzer_input_1.mzML -ini ${DATA_DIR_TOPP}/ITRAQAnalyzer.ini -out ITRAQAnalyzer_output_1.tmp)
add_test("TOPP_ITRAQAnalyzer_1_out1" ${DIFF} -whitelist "<map" "?xml-stylesheet" -in1 ITRAQAnalyzer_output_1.tmp -in2 ${DATA_DIR_TOPP}/ITRAQAnalyzer_output_1.consensusXML )
set_tests_properties("TOPP_ITRAQAnalyzer_1_out1" PROPERTIES DEPENDS "TOPP_ITRAQAnalyzer_1")
# test empty IDs
add_test("TOPP_ITRAQAnalyzer_2" ${TOPP_BIN_PATH}/ITRAQAnalyzer -test -in ${DATA_DIR_TOPP}/ITRAQAnalyzer_input_2.mzML -ini ${DATA_DIR_TOPP}/ITRAQAnalyzer.ini -out ITRAQAnalyzer_output_2.tmp)
add_test("TOPP_ITRAQAnalyzer_2_out1" ${DIFF} -whitelist "<map" "?xml-stylesheet" -in1 ITRAQAnalyzer_output_2.tmp -in2 ${DATA_DIR_TOPP}/ITRAQAnalyzer_output_2.consensusXML )
set_tests_properties("TOPP_ITRAQAnalyzer_2_out1" PROPERTIES DEPENDS "TOPP_ITRAQAnalyzer_2")

### IDConflictResolver tests
add_test("TOPP_IDConflictResolver_1" ${TOPP_BIN_PATH}/IDConflictResolver -test -in ${DATA_DIR_TOPP}/IDConflictResolver_1_input.featureXML -out IDConflictResolver_1_output.tmp)
add_test("TOPP_IDConflictResolver_1_out1" ${DIFF} -in1 IDConflictResolver_1_output.tmp -in2 ${DATA_DIR_TOPP}/IDConflictResolver_1_output.featureXML )
set_tests_properties("TOPP_IDConflictResolver_1_out1" PROPERTIES DEPENDS "TOPP_IDConflictResolver_1")
add_test("TOPP_IDConflictResolver_2" ${TOPP_BIN_PATH}/IDConflictResolver -test -in ${DATA_DIR_TOPP}/IDConflictResolver_2_input.consensusXML -out IDConflictResolver_2_output.tmp)
add_test("TOPP_IDConflictResolver_2_out1" ${DIFF} -in1 IDConflictResolver_2_output.tmp -in2 ${DATA_DIR_TOPP}/IDConflictResolver_2_output.consensusXML )
set_tests_properties("TOPP_IDConflictResolver_2_out1" PROPERTIES DEPENDS "TOPP_IDConflictResolver_2")

### IDFileConverter tests
add_test("TOPP_IDFileConverter_1" ${TOPP_BIN_PATH}/IDFileConverter -in ${DATA_DIR_TOPP}/IDFileConverter_1_input_sequestoutdir -mz_file ${DATA_DIR_TOPP}/IDFileConverter_1_input_mzfile.mzML -out IDFileConverter_1_output.tmp -out_type idXML -ignore_proteins_per_peptide )
add_test("TOPP_IDFileConverter_1_out1" ${DIFF} -in1 IDFileConverter_1_output.tmp -in2 ${DATA_DIR_TOPP}/IDFileConverter_1_output.idXML )
set_tests_properties("TOPP_IDFileConverter_1_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_1")
add_test("TOPP_IDFileConverter_2" ${TOPP_BIN_PATH}/IDFileConverter -in ${PROJECT_SOURCE_DIR}/data/PepXMLFile_test.pepxml -mz_file ${PROJECT_SOURCE_DIR}/data/PepXMLFile_test.mzML -mz_name PepXMLFile_test -use_precursor_data -out IDFileConverter_2_output.tmp -out_type idXML)
add_test("TOPP_IDFileConverter_2_out1" ${DIFF} -in1 IDFileConverter_2_output.tmp -in2 ${DATA_DIR_TOPP}/IDFileConverter_2_output.idXML )
set_tests_properties("TOPP_IDFileConverter_2_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_2")
add_test("TOPP_IDFileConverter_3" ${TOPP_BIN_PATH}/IDFileConverter -in ${DATA_DIR_TOPP}/IDFileConverter_3_input.protXML -out IDFileConverter_3_output.tmp -out_type idXML)
add_test("TOPP_IDFileConverter_3_out1" ${DIFF} -in1 IDFileConverter_3_output.tmp -in2 ${DATA_DIR_TOPP}/IDFileConverter_3_output.idXML )
set_tests_properties("TOPP_IDFileConverter_3_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_3")
add_test("TOPP_IDFileConverter_4" ${TOPP_BIN_PATH}/IDFileConverter -in ${PROJECT_SOURCE_DIR}/data/PepXMLFile_test.pepxml -mz_file ${PROJECT_SOURCE_DIR}/data/PepXMLFile_test.mzML -mz_name PepXMLFile_test -out IDFileConverter_4_output.tmp -out_type idXML)
add_test("TOPP_IDFileConverter_4_out1" ${DIFF} -in1 IDFileConverter_4_output.tmp -in2 ${DATA_DIR_TOPP}/IDFileConverter_4_output.idXML )
set_tests_properties("TOPP_IDFileConverter_4_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_4")

### IDFilter tests
add_test("TOPP_IDFilter_1" ${TOPP_BIN_PATH}/IDFilter -in ${DATA_DIR_TOPP}/IDFilter_1_input.idXML -out IDFilter_1_output.tmp -whitelist:proteins ${DATA_DIR_TOPP}/IDFilter_1_input.fas)
add_test("TOPP_IDFilter_1_out1" ${DIFF} -in1 IDFilter_1_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_1_output.idXML )
set_tests_properties("TOPP_IDFilter_1_out1" PROPERTIES DEPENDS "TOPP_IDFilter_1")
add_test("TOPP_IDFilter_2" ${TOPP_BIN_PATH}/IDFilter -in ${DATA_DIR_TOPP}/IDFilter_2_input.idXML -out IDFilter_2_output.tmp -thresh:pep 1 -thresh:prot 1)
add_test("TOPP_IDFilter_2_out1" ${DIFF} -in1 IDFilter_2_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_2_output.idXML )
set_tests_properties("TOPP_IDFilter_2_out1" PROPERTIES DEPENDS "TOPP_IDFilter_2")
add_test("TOPP_IDFilter_3" ${TOPP_BIN_PATH}/IDFilter -in ${DATA_DIR_TOPP}/IDFilter_3_input.idXML -out IDFilter_3_output.tmp -blacklist:peptides ${DATA_DIR_TOPP}/IDFilter_3_2_input.idXML)
add_test("TOPP_IDFilter_3_out1" ${DIFF} -in1 IDFilter_3_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_3_output.idXML )
set_tests_properties("TOPP_IDFilter_3_out1" PROPERTIES DEPENDS "TOPP_IDFilter_3")
add_test("TOPP_IDFilter_4" ${TOPP_BIN_PATH}/IDFilter -in ${DATA_DIR_TOPP}/IDFilter_4_input.idXML -out IDFilter_4_output.tmp -rt:p_value 0.08)
add_test("TOPP_IDFilter_4_out1" ${DIFF} -in1 IDFilter_4_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_4_output.idXML )
set_tests_properties("TOPP_IDFilter_4_out1" PROPERTIES DEPENDS "TOPP_IDFilter_4")
add_test("TOPP_IDFilter_5" ${TOPP_BIN_PATH}/IDFilter -in ${DATA_DIR_TOPP}/IDFilter_5_input.idXML -out IDFilter_5_output.tmp -score:pep 32 -score:prot 25)
add_test("TOPP_IDFilter_5_out1" ${DIFF} -in1 IDFilter_5_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_5_output.idXML )
set_tests_properties("TOPP_IDFilter_5_out1" PROPERTIES DEPENDS "TOPP_IDFilter_5")
add_test("TOPP_IDFilter_6" ${TOPP_BIN_PATH}/IDFilter -in ${DATA_DIR_TOPP}/IDFilter_6_input.idXML -out IDFilter_6_output.tmp -best:n_peptide_hits 2 -best:n_protein_hits 10)
add_test("TOPP_IDFilter_6_out1" ${DIFF} -in1 IDFilter_6_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_6_output.idXML )
set_tests_properties("TOPP_IDFilter_6_out1" PROPERTIES DEPENDS "TOPP_IDFilter_6")
add_test("TOPP_IDFilter_7" ${TOPP_BIN_PATH}/IDFilter -in ${DATA_DIR_TOPP}/IDFilter_7_input.idXML -out IDFilter_7_output.tmp -unique)
add_test("TOPP_IDFilter_7_out1" ${DIFF} -in1 IDFilter_7_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_7_output.idXML )
set_tests_properties("TOPP_IDFilter_7_out1" PROPERTIES DEPENDS "TOPP_IDFilter_7")

### MapAlignerPoseClustering tests
# featureXML input:
add_test("TOPP_MapAlignerPoseClustering_1" ${TOPP_BIN_PATH}/MapAlignerPoseClustering -test -ini ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_parameters.ini -in ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_input1.featureXML ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_input2.featureXML ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_input3.featureXML -out MapAlignerPoseClustering_1_output1.tmp MapAlignerPoseClustering_1_output2.tmp MapAlignerPoseClustering_1_output3.tmp -trafo_out MapAlignerPoseClustering_1_trafo1.tmp MapAlignerPoseClustering_1_trafo2.tmp MapAlignerPoseClustering_1_trafo3.tmp)
add_test("TOPP_MapAlignerPoseClustering_1_out1" ${DIFF} -in1 MapAlignerPoseClustering_1_output1.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_output1.featureXML )
set_tests_properties("TOPP_MapAlignerPoseClustering_1_out1" PROPERTIES DEPENDS "TOPP_MapAlignerPoseClustering_1")
add_test("TOPP_MapAlignerPoseClustering_1_out2" ${DIFF} -in1 MapAlignerPoseClustering_1_output2.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_output2.featureXML )
set_tests_properties("TOPP_MapAlignerPoseClustering_1_out2" PROPERTIES DEPENDS "TOPP_MapAlignerPoseClustering_1")
add_test("TOPP_MapAlignerPoseClustering_1_out3" ${DIFF} -in1 MapAlignerPoseClustering_1_output3.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_output3.featureXML )
set_tests_properties("TOPP_MapAlignerPoseClustering_1_out3" PROPERTIES DEPENDS "TOPP_MapAlignerPoseClustering_1")
add_test("TOPP_MapAlignerPoseClustering_1_out4" ${DIFF} -in1 MapAlignerPoseClustering_1_trafo1.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_trafo1.trafoXML )
set_tests_properties("TOPP_MapAlignerPoseClustering_1_out4" PROPERTIES DEPENDS "TOPP_MapAlignerPoseClustering_1")
add_test("TOPP_MapAlignerPoseClustering_1_out5" ${DIFF} -in1 MapAlignerPoseClustering_1_trafo2.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_trafo2.trafoXML )
set_tests_properties("TOPP_MapAlignerPoseClustering_1_out5" PROPERTIES DEPENDS "TOPP_MapAlignerPoseClustering_1")
add_test("TOPP_MapAlignerPoseClustering_1_out6" ${DIFF} -in1 MapAlignerPoseClustering_1_trafo3.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_trafo3.trafoXML )
set_tests_properties("TOPP_MapAlignerPoseClustering_1_out6" PROPERTIES DEPENDS "TOPP_MapAlignerPoseClustering_1")
# mzML input:
add_test("TOPP_MapAlignerPoseClustering_2" ${TOPP_BIN_PATH}/MapAlignerPoseClustering -test -ini ${DATA_DIR_TOPP}/MapAlignerPoseClustering_2_parameters.ini -in ${DATA_DIR_TOPP}/MapAlignerPoseClustering_2_input1.mzML ${DATA_DIR_TOPP}/MapAlignerPoseClustering_2_input2.mzML ${DATA_DIR_TOPP}/MapAlignerPoseClustering_2_input3.mzML -out MapAlignerPoseClustering_2_output1.tmp MapAlignerPoseClustering_2_output2.tmp MapAlignerPoseClustering_2_output3.tmp)
add_test("TOPP_MapAlignerPoseClustering_2_out1" ${DIFF} -in1 MapAlignerPoseClustering_2_output1.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerPoseClustering_2_output1.mzML )
set_tests_properties("TOPP_MapAlignerPoseClustering_2_out1" PROPERTIES DEPENDS "TOPP_MapAlignerPoseClustering_2")
add_test("TOPP_MapAlignerPoseClustering_2_out2" ${DIFF} -in1 MapAlignerPoseClustering_2_output2.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerPoseClustering_2_output2.mzML )
set_tests_properties("TOPP_MapAlignerPoseClustering_2_out2" PROPERTIES DEPENDS "TOPP_MapAlignerPoseClustering_2")
add_test("TOPP_MapAlignerPoseClustering_2_out3" ${DIFF} -in1 MapAlignerPoseClustering_2_output3.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerPoseClustering_2_output3.mzML )
set_tests_properties("TOPP_MapAlignerPoseClustering_2_out3" PROPERTIES DEPENDS "TOPP_MapAlignerPoseClustering_2")
# "reference:file" option:
add_test("TOPP_MapAlignerPoseClustering_3" ${TOPP_BIN_PATH}/MapAlignerPoseClustering -test -ini ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_parameters.ini -in ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_input2.featureXML ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_input3.featureXML -out MapAlignerPoseClustering_3_output1.tmp MapAlignerPoseClustering_3_output2.tmp -reference:file ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_input1.featureXML)
add_test("TOPP_MapAlignerPoseClustering_3_out1" ${DIFF} -in1 MapAlignerPoseClustering_3_output1.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_output2.featureXML )
set_tests_properties("TOPP_MapAlignerPoseClustering_3_out1" PROPERTIES DEPENDS "TOPP_MapAlignerPoseClustering_3")
add_test("TOPP_MapAlignerPoseClustering_3_out2" ${DIFF} -in1 MapAlignerPoseClustering_3_output2.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_output3.featureXML )
set_tests_properties("TOPP_MapAlignerPoseClustering_3_out2" PROPERTIES DEPENDS "TOPP_MapAlignerPoseClustering_3")
# "reference:index" option:
add_test("TOPP_MapAlignerPoseClustering_4" ${TOPP_BIN_PATH}/MapAlignerPoseClustering -test -ini ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_parameters.ini -in ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_input1.featureXML ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_input2.featureXML -trafo_out MapAlignerPoseClustering_4_trafo1.tmp MapAlignerPoseClustering_4_trafo2.tmp -reference:index 2)
add_test("TOPP_MapAlignerPoseClustering_4_out1" ${DIFF} -in1 MapAlignerPoseClustering_4_trafo1.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerPoseClustering_4_trafo1.trafoXML )
set_tests_properties("TOPP_MapAlignerPoseClustering_4_out1" PROPERTIES DEPENDS "TOPP_MapAlignerPoseClustering_4")
add_test("TOPP_MapAlignerPoseClustering_4_out2" ${DIFF} -in1 MapAlignerPoseClustering_4_trafo2.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_trafo1.trafoXML )
set_tests_properties("TOPP_MapAlignerPoseClustering_4_out2" PROPERTIES DEPENDS "TOPP_MapAlignerPoseClustering_4")

### MapAlignerIdentification tests:
add_test("TOPP_MapAlignerIdentification_1" ${TOPP_BIN_PATH}/MapAlignerIdentification -test -ini ${DATA_DIR_TOPP}/MapAlignerIdentification_parameters.ini -in ${DATA_DIR_TOPP}/MapAlignerIdentification_1_input1.featureXML ${DATA_DIR_TOPP}/MapAlignerIdentification_1_input2.featureXML -out MapAlignerIdentification_1_output1.tmp MapAlignerIdentification_1_output2.tmp)
add_test("TOPP_MapAlignerIdentification_1_out1" ${DIFF} -in1 MapAlignerIdentification_1_output1.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerIdentification_1_output1.featureXML )
set_tests_properties("TOPP_MapAlignerIdentification_1_out1" PROPERTIES DEPENDS "TOPP_MapAlignerIdentification_1")
add_test("TOPP_MapAlignerIdentification_1_out2" ${DIFF} -in1 MapAlignerIdentification_1_output2.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerIdentification_1_output2.featureXML )
set_tests_properties("TOPP_MapAlignerIdentification_1_out2" PROPERTIES DEPENDS "TOPP_MapAlignerIdentification_1")
# "reference:file" option:
add_test("TOPP_MapAlignerIdentification_2" ${TOPP_BIN_PATH}/MapAlignerIdentification -test -ini ${DATA_DIR_TOPP}/MapAlignerIdentification_parameters.ini -in ${DATA_DIR_TOPP}/MapAlignerIdentification_1_input1.featureXML -out MapAlignerIdentification_2_output1.tmp -reference:file ${DATA_DIR_TOPP}/MapAlignerIdentification_1_input2.featureXML)
add_test("TOPP_MapAlignerIdentification_2_out1" ${DIFF} -in1 MapAlignerIdentification_2_output1.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerIdentification_2_output1.featureXML )
set_tests_properties("TOPP_MapAlignerIdentification_2_out1" PROPERTIES DEPENDS "TOPP_MapAlignerIdentification_2")
# "reference:index" option:
add_test("TOPP_MapAlignerIdentification_3" ${TOPP_BIN_PATH}/MapAlignerIdentification -test -ini ${DATA_DIR_TOPP}/MapAlignerIdentification_parameters.ini -in ${DATA_DIR_TOPP}/MapAlignerIdentification_1_input2.featureXML ${DATA_DIR_TOPP}/MapAlignerIdentification_1_input1.featureXML -out MapAlignerIdentification_3_output1.tmp MapAlignerIdentification_3_output2.tmp -reference:index 1)
add_test("TOPP_MapAlignerIdentification_3_out1" ${DIFF} -in1 MapAlignerIdentification_3_output1.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerIdentification_3_output1.featureXML )
set_tests_properties("TOPP_MapAlignerIdentification_3_out1" PROPERTIES DEPENDS "TOPP_MapAlignerIdentification_3")
add_test("TOPP_MapAlignerIdentification_3_out2" ${DIFF} -in1 MapAlignerIdentification_3_output2.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerIdentification_2_output1.featureXML )
set_tests_properties("TOPP_MapAlignerIdentification_3_out2" PROPERTIES DEPENDS "TOPP_MapAlignerIdentification_3")
add_test("TOPP_MapAlignerIdentification_4" ${TOPP_BIN_PATH}/MapAlignerIdentification -test -ini ${DATA_DIR_TOPP}/MapAlignerIdentification_parameters.ini -in  ${DATA_DIR_TOPP}/MapAlignerIdentification_1_input1.featureXML ${DATA_DIR_TOPP}/MapAlignerIdentification_1_input2.featureXML -out MapAlignerIdentification_4_output1.tmp MapAlignerIdentification_4_output2.tmp -reference:index 2)
add_test("TOPP_MapAlignerIdentification_4_out1" ${DIFF} -in1 MapAlignerIdentification_4_output1.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerIdentification_2_output1.featureXML )
set_tests_properties("TOPP_MapAlignerIdentification_4_out1" PROPERTIES DEPENDS "TOPP_MapAlignerIdentification_4")
add_test("TOPP_MapAlignerIdentification_4_out2" ${DIFF} -in1 MapAlignerIdentification_4_output2.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerIdentification_3_output1.featureXML )
set_tests_properties("TOPP_MapAlignerIdentification_4_out2" PROPERTIES DEPENDS "TOPP_MapAlignerIdentification_4")
# consensusXML input:
add_test("TOPP_MapAlignerIdentification_5" ${TOPP_BIN_PATH}/MapAlignerIdentification -test -ini ${DATA_DIR_TOPP}/MapAlignerIdentification_parameters.ini -in ${DATA_DIR_TOPP}/MapAlignerIdentification_5_input1.consensusXML ${DATA_DIR_TOPP}/MapAlignerIdentification_5_input2.consensusXML -out MapAlignerIdentification_5_output1.tmp MapAlignerIdentification_5_output2.tmp)
add_test("TOPP_MapAlignerIdentification_5_out1" ${DIFF} -in1 MapAlignerIdentification_5_output1.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerIdentification_5_output1.consensusXML )
set_tests_properties("TOPP_MapAlignerIdentification_5_out1" PROPERTIES DEPENDS "TOPP_MapAlignerIdentification_5")
add_test("TOPP_MapAlignerIdentification_5_out2" ${DIFF} -in1 MapAlignerIdentification_5_output2.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerIdentification_5_output2.consensusXML )
set_tests_properties("TOPP_MapAlignerIdentification_5_out2" PROPERTIES DEPENDS "TOPP_MapAlignerIdentification_5")

### MapAlignerSpectrumAlignment tests:
add_test("TOPP_MapAlignerSpectrum_1" ${TOPP_BIN_PATH}/MapAlignerSpectrum -test -ini ${DATA_DIR_TOPP}/MapAlignerSpectrum_parameters.ini -in ${DATA_DIR_TOPP}/MapAlignerSpectrum_1_input1.mzML ${DATA_DIR_TOPP}/MapAlignerSpectrum_1_input2.mzML ${DATA_DIR_TOPP}/MapAlignerSpectrum_1_input3.mzML -out MapAlignerSpectrum_1_output1.tmp MapAlignerSpectrum_1_output2.tmp MapAlignerSpectrum_1_output3.tmp)
add_test("TOPP_MapAlignerSpectrum_1_out1" ${DIFF} -in1 MapAlignerSpectrum_1_output1.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerSpectrum_1_output1.mzML )
set_tests_properties("TOPP_MapAlignerSpectrum_1_out1" PROPERTIES DEPENDS "TOPP_MapAlignerSpectrum_1")
add_test("TOPP_MapAlignerSpectrum_1_out2" ${DIFF} -in1 MapAlignerSpectrum_1_output2.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerSpectrum_1_output2.mzML )
set_tests_properties("TOPP_MapAlignerSpectrum_1_out2" PROPERTIES DEPENDS "TOPP_MapAlignerSpectrum_1")
add_test("TOPP_MapAlignerSpectrum_1_out3" ${DIFF} -in1 MapAlignerSpectrum_1_output3.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerSpectrum_1_output3.mzML )
set_tests_properties("TOPP_MapAlignerSpectrum_1_out3" PROPERTIES DEPENDS "TOPP_MapAlignerSpectrum_1")

### MapRTTransformer tests:
add_test("TOPP_MapRTTransformer_1" ${TOPP_BIN_PATH}/MapRTTransformer -test -in ${DATA_DIR_TOPP}/MapRTTransformer_1_input.featureXML ${DATA_DIR_TOPP}/MapRTTransformer_1_input.featureXML -trafo_in ${DATA_DIR_TOPP}/MapRTTransformer_1_trafo1.trafoXML ${DATA_DIR_TOPP}/MapRTTransformer_1_trafo2.trafoXML -out MapRTTransformer_1_output1.tmp MapRTTransformer_1_output2.tmp)
add_test("TOPP_MapRTTransformer_1_out1" ${DIFF} -in1 MapRTTransformer_1_output1.tmp -in2 ${DATA_DIR_TOPP}/MapRTTransformer_1_output1.featureXML )
set_tests_properties("TOPP_MapRTTransformer_1_out1" PROPERTIES DEPENDS "TOPP_MapRTTransformer_1")
add_test("TOPP_MapRTTransformer_1_out2" ${DIFF} -in1 MapRTTransformer_1_output2.tmp -in2 ${DATA_DIR_TOPP}/MapRTTransformer_1_output2.featureXML )
set_tests_properties("TOPP_MapRTTransformer_1_out2" PROPERTIES DEPENDS "TOPP_MapRTTransformer_1")
add_test("TOPP_MapRTTransformer_2" ${TOPP_BIN_PATH}/MapRTTransformer -test -in ${DATA_DIR_TOPP}/MapRTTransformer_2_input.mzML -trafo_in ${DATA_DIR_TOPP}/MapRTTransformer_2_trafo.trafoXML -out MapRTTransformer_2_output.tmp)
add_test("TOPP_MapRTTransformer_2_out1" ${DIFF} -in1 MapRTTransformer_2_output.tmp -in2 ${DATA_DIR_TOPP}/MapRTTransformer_2_output.mzML )
set_tests_properties("TOPP_MapRTTransformer_2_out1" PROPERTIES DEPENDS "TOPP_MapRTTransformer_2")
# "invert" option:
add_test("TOPP_MapRTTransformer_3" ${TOPP_BIN_PATH}/MapRTTransformer -test -invert -trafo_in ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_trafo2.trafoXML -trafo_out MapRTTransformer_3_output.tmp)
add_test("TOPP_MapRTTransformer_3_out1" ${DIFF} -in1 MapRTTransformer_3_output.tmp -in2 ${DATA_DIR_TOPP}/MapRTTransformer_3_trafo.trafoXML )
set_tests_properties("TOPP_MapRTTransformer_3_out1" PROPERTIES DEPENDS "TOPP_MapRTTransformer_3")
# chromatograms
add_test("TOPP_MapRTTransformer_4" ${TOPP_BIN_PATH}/MapRTTransformer -test -in ${DATA_DIR_TOPP}/MapRTTransformer_4_input.chrom.mzML -trafo_in ${DATA_DIR_TOPP}/MapRTTransformer_4_trafo.trafoXML -out MapRTTransformer_4_output.tmp)
add_test("TOPP_MapRTTransformer_4_out1" ${DIFF} -in1 MapRTTransformer_4_output.tmp -in2 ${DATA_DIR_TOPP}/MapRTTransformer_4_output.chrom.mzML )
set_tests_properties("TOPP_MapRTTransformer_4_out1" PROPERTIES DEPENDS "TOPP_MapRTTransformer_4")

### MascotAdapter tests
add_test("TOPP_MascotAdapter_1" ${TOPP_BIN_PATH}/MascotAdapter -ini ${DATA_DIR_TOPP}/MascotAdapter_1_parameters.ini -mascot_in -in ${DATA_DIR_TOPP}/MascotAdapter_1_input.mzData)
add_test("TOPP_MascotAdapter_1_out1" ${DIFF} -in1 MascotAdapter_1_output.tmp -in2 ${DATA_DIR_TOPP}/MascotAdapter_1_output.mascot_in )
set_tests_properties("TOPP_MascotAdapter_1_out1" PROPERTIES DEPENDS "TOPP_MascotAdapter_1")
add_test("TOPP_MascotAdapter_2" ${TOPP_BIN_PATH}/MascotAdapter -in ${DATA_DIR_TOPP}/MascotAdapter_2_input.mascotXML -out MascotAdapter_2_output.tmp -mascot_out)
add_test("TOPP_MascotAdapter_2_out1" ${DIFF} -in1 MascotAdapter_2_output.tmp -in2 ${DATA_DIR_TOPP}/MascotAdapter_2_output.idXML )
set_tests_properties("TOPP_MascotAdapter_2_out1" PROPERTIES DEPENDS "TOPP_MascotAdapter_2")

### NoiseFilter tests
add_test("TOPP_NoiseFilterGaussian_1" ${TOPP_BIN_PATH}/NoiseFilterGaussian -test -ini ${DATA_DIR_TOPP}/NoiseFilterGaussian_1_parameters.ini -in ${DATA_DIR_TOPP}/NoiseFilterGaussian_1_input.mzML -out NoiseFilterGaussian_1.tmp)
add_test("TOPP_NoiseFilterGaussian_1_out1" ${DIFF} -in1 NoiseFilterGaussian_1.tmp -in2 ${DATA_DIR_TOPP}/NoiseFilterGaussian_1_output.mzML )
set_tests_properties("TOPP_NoiseFilterGaussian_1_out1" PROPERTIES DEPENDS "TOPP_NoiseFilterGaussian_1")
add_test("TOPP_NoiseFilterSGolay_1" ${TOPP_BIN_PATH}/NoiseFilterSGolay -test -ini ${DATA_DIR_TOPP}/NoiseFilterSGolay_1_parameters.ini -in ${DATA_DIR_TOPP}/NoiseFilterSGolay_1_input.mzML -out NoiseFilterSGolay_1.tmp)
add_test("TOPP_NoiseFilterSGolay_1_out1" ${DIFF} -in1 NoiseFilterSGolay_1.tmp -in2 ${DATA_DIR_TOPP}/NoiseFilterSGolay_1_output.mzML )
set_tests_properties("TOPP_NoiseFilterSGolay_1_out1" PROPERTIES DEPENDS "TOPP_NoiseFilterSGolay_1")

### PeakPicker tests
add_test("TOPP_PeakPickerWavelet_1" ${TOPP_BIN_PATH}/PeakPickerWavelet  -test -ini ${DATA_DIR_TOPP}/PeakPickerWavelet_parameters.ini -in ${DATA_DIR_TOPP}/PeakPickerWavelet_input.mzML -out PeakPickerWavelet_1.tmp)
add_test("TOPP_PeakPickerWavelet_1_out1" ${DIFF} -in1 PeakPickerWavelet_1.tmp -in2 ${DATA_DIR_TOPP}/PeakPickerWavelet_output.mzML )
set_tests_properties("TOPP_PeakPickerWavelet_1_out1" PROPERTIES DEPENDS "TOPP_PeakPickerWavelet_1")
add_test("TOPP_PeakPickerWavelet_2" ${TOPP_BIN_PATH}/PeakPickerWavelet  -test -ini ${DATA_DIR_TOPP}/PeakPickerWavelet_deconv_parameters.ini -in ${DATA_DIR_TOPP}/PeakPickerWavelet_deconv_input.mzML -out PeakPickerWavelet_2.tmp)
add_test("TOPP_PeakPickerWavelet_2_out1" ${DIFF} -in1 PeakPickerWavelet_2.tmp -in2 ${DATA_DIR_TOPP}/PeakPickerWavelet_deconv_output.mzML)
set_tests_properties("TOPP_PeakPickerWavelet_2_out1" PROPERTIES DEPENDS "TOPP_PeakPickerWavelet_2")
add_test("TOPP_PeakPickerWavelet_3" ${TOPP_BIN_PATH}/PeakPickerWavelet  -test -ini ${DATA_DIR_TOPP}/PeakPickerWavelet_2Dopt_parameters.ini -in ${DATA_DIR_TOPP}/PeakPickerWavelet_2Dopt_input.mzML -out PeakPickerWavelet_3.tmp)
add_test("TOPP_PeakPickerWavelet_3_out1" ${DIFF} -in1 PeakPickerWavelet_3.tmp -in2 ${DATA_DIR_TOPP}/PeakPickerWavelet_2Dopt_output.mzML)
set_tests_properties("TOPP_PeakPickerWavelet_3_out1" PROPERTIES DEPENDS "TOPP_PeakPickerWavelet_3")
add_test("TOPP_PeakPickerWavelet_4" ${TOPP_BIN_PATH}/PeakPickerWavelet  -test -ini ${DATA_DIR_TOPP}/PeakPickerWavelet_parameters.ini -in ${DATA_DIR_TOPP}/PeakPickerWavelet_input.mzML -out PeakPickerWavelet_4.tmp -threads 2)
add_test("TOPP_PeakPickerWavelet_4_out1" ${DIFF} -in1 PeakPickerWavelet_4.tmp -in2 ${DATA_DIR_TOPP}/PeakPickerWavelet_output.mzML)
set_tests_properties("TOPP_PeakPickerWavelet_4_out1" PROPERTIES DEPENDS "TOPP_PeakPickerWavelet_4")
add_test("TOPP_PeakPickerWavelet_5" ${TOPP_BIN_PATH}/PeakPickerWavelet  -test -ini ${DATA_DIR_TOPP}/PeakPickerWavelet_parameters_noMetaData.ini -in ${DATA_DIR_TOPP}/PeakPickerWavelet_input.mzML -out PeakPickerWavelet_5.tmp -threads 2)
add_test("TOPP_PeakPickerWavelet_5_out1" ${DIFF} -in1 PeakPickerWavelet_5.tmp -in2 ${DATA_DIR_TOPP}/PeakPickerWavelet_output_noMetaData.mzML)
set_tests_properties("TOPP_PeakPickerWavelet_5_out1" PROPERTIES DEPENDS "TOPP_PeakPickerWavelet_5")

# "high_res" algorithm with "ms1_only" option:
add_test("TOPP_PeakPickerHiRes_1" ${TOPP_BIN_PATH}/PeakPickerHiRes -test -ini ${DATA_DIR_TOPP}/PeakPickerHiRes_parameters.ini -in ${DATA_DIR_TOPP}/PeakPickerHiRes_input.mzML -out PeakPickerHiRes_1.tmp)
add_test("TOPP_PeakPickerHiRes_1_out1" ${DIFF} -in1 PeakPickerHiRes_1.tmp -in2 ${DATA_DIR_TOPP}/PeakPickerHiRes_output.mzML)
set_tests_properties("TOPP_PeakPickerHiRes_1_out1" PROPERTIES DEPENDS "TOPP_PeakPickerHiRes_1")

### Resampler tests
add_test("TOPP_Resampler_1" ${TOPP_BIN_PATH}/Resampler -test -in ${DATA_DIR_TOPP}/Resampler_1_input.mzML -out Resampler.mzML -sampling_rate 0.3)
add_test("TOPP_Resampler_1_out1" ${DIFF} -in1 Resampler.mzML -in2 ${DATA_DIR_TOPP}/Resampler_1_output.mzML )
set_tests_properties("TOPP_Resampler_1_out1" PROPERTIES DEPENDS "TOPP_Resampler_1")

### RTModel tests
add_test("TOPP_RTModel_1" ${TOPP_BIN_PATH}/RTModel -in ${DATA_DIR_TOPP}/RTModel_1_input.idXML -out RTModel_1_output.tmp -ini ${DATA_DIR_TOPP}/RTModel_1_parameters.ini)
add_test("TOPP_RTModel_1_out1" ${DIFF} -in1 RTModel_1_output.tmp -in2 ${DATA_DIR_TOPP}/RTModel_1_output.model )
set_tests_properties("TOPP_RTModel_1_out1" PROPERTIES DEPENDS "TOPP_RTModel_1")
add_test("TOPP_RTModel_2" ${TOPP_BIN_PATH}/RTModel -in_positive ${DATA_DIR_TOPP}/RTModel_2_input_positive.idXML -in_negative ${DATA_DIR_TOPP}/RTModel_2_input_negative.idXML -out RTModel_2_output.tmp -ini ${DATA_DIR_TOPP}/RTModel_2_parameters.ini)
add_test("TOPP_RTModel_2_out1" ${DIFF} -in1 RTModel_2_output.tmp -in2 ${DATA_DIR_TOPP}/RTModel_2_output.model )
set_tests_properties("TOPP_RTModel_2_out1" PROPERTIES DEPENDS "TOPP_RTModel_2")
add_test("TOPP_RTModel_2_out2" ${DIFF} -in1 RTModel_2_output.tmp_samples -in2 ${DATA_DIR_TOPP}/RTModel_2_output.samples )
set_tests_properties("TOPP_RTModel_2_out2" PROPERTIES DEPENDS "TOPP_RTModel_2")
add_test("TOPP_RTModel_2_out3" ${DIFF} -in1 RTModel_2_output.tmp_additional_parameters -in2 ${DATA_DIR_TOPP}/RTModel_2_output.additional_parameters )
set_tests_properties("TOPP_RTModel_2_out3" PROPERTIES DEPENDS "TOPP_RTModel_2")
add_test("TOPP_RTModel_3" ${TOPP_BIN_PATH}/RTModel -in ${DATA_DIR_TOPP}/RTModel_3_input.idXML -out RTModel_3_output.tmp -ini ${DATA_DIR_TOPP}/RTModel_3_parameters.ini)
add_test("TOPP_RTModel_3_out1" ${DIFF} -in1 RTModel_3_output.tmp -in2 ${DATA_DIR_TOPP}/RTModel_3_output.model )
set_tests_properties("TOPP_RTModel_3_out1" PROPERTIES DEPENDS "TOPP_RTModel_3")
add_test("TOPP_RTModel_3_out2" ${DIFF} -in1 RTModel_3_output.tmp_samples -in2 ${DATA_DIR_TOPP}/RTModel_3_output.samples )
set_tests_properties("TOPP_RTModel_3_out2" PROPERTIES DEPENDS "TOPP_RTModel_3")
add_test("TOPP_RTModel_4" ${TOPP_BIN_PATH}/RTModel -in ${DATA_DIR_TOPP}/RTModel_4_input.txt -out RTModel_4_output.tmp -ini ${DATA_DIR_TOPP}/RTModel_4_parameters.ini -textfile_input)
add_test("TOPP_RTModel_4_out1" ${DIFF} -in1 RTModel_4_output.tmp -in2 ${DATA_DIR_TOPP}/RTModel_4_output.model )
set_tests_properties("TOPP_RTModel_4_out1" PROPERTIES DEPENDS "TOPP_RTModel_4")
add_test("TOPP_RTModel_4_out2" ${DIFF} -in1 RTModel_4_output.tmp_samples -in2 ${DATA_DIR_TOPP}/RTModel_4_output.samples )
set_tests_properties("TOPP_RTModel_4_out2" PROPERTIES DEPENDS "TOPP_RTModel_4")

### PTModel tests
add_test("TOPP_PTModel_1" ${TOPP_BIN_PATH}/PTModel -in_positive ${DATA_DIR_TOPP}/PTModel_1_input_positive.idXML -in_negative ${DATA_DIR_TOPP}/PTModel_1_input_negative.idXML -out PTModel_1_output.tmp -ini ${DATA_DIR_TOPP}/PTModel_1_parameters.ini)
add_test("TOPP_PTModel_1_out1" ${DIFF} -whitelist "probA" "probB" -in1 PTModel_1_output.tmp -in2 ${DATA_DIR_TOPP}/PTModel_1_output.model )
set_tests_properties("TOPP_PTModel_1_out1" PROPERTIES DEPENDS "TOPP_PTModel_1")
add_test("TOPP_PTModel_1_out2" ${DIFF} -in1 PTModel_1_output.tmp_samples -in2 ${DATA_DIR_TOPP}/PTModel_1_output.samples )
set_tests_properties("TOPP_PTModel_1_out2" PROPERTIES DEPENDS "TOPP_PTModel_1")
add_test("TOPP_PTModel_1_out3" ${DIFF} -in1 PTModel_1_output.tmp_additional_parameters -in2 ${DATA_DIR_TOPP}/PTModel_1_output.additional_parameters )
set_tests_properties("TOPP_PTModel_1_out3" PROPERTIES DEPENDS "TOPP_PTModel_1")

### PTPredict tests
add_test("TOPP_PTPredict_1" ${TOPP_BIN_PATH}/PTPredict -in ${DATA_DIR_TOPP}/PTPredict_1_input.idXML -out PTPredict_1_output.tmp -svm_model ${DATA_DIR_TOPP}/PTPredict_1_input.model)
add_test("TOPP_PTPredict_1_out1" ${DIFF} -in1 PTPredict_1_output.tmp -in2 ${DATA_DIR_TOPP}/PTPredict_1_output.idXML )
set_tests_properties("TOPP_PTPredict_1_out1" PROPERTIES DEPENDS "TOPP_PTPredict_1")

### RTPredict tests
add_test("TOPP_RTPredict_1" ${TOPP_BIN_PATH}/RTPredict -in_id ${DATA_DIR_TOPP}/RTPredict_1_input.idXML -out_id:file RTPredict_1_output.tmp -total_gradient_time 3000 -svm_model ${DATA_DIR_TOPP}/RTPredict_1_input.model)
add_test("TOPP_RTPredict_1_out1" ${DIFF} -in1 RTPredict_1_output.tmp -in2 ${DATA_DIR_TOPP}/RTPredict_1_output.idXML )
set_tests_properties("TOPP_RTPredict_1_out1" PROPERTIES DEPENDS "TOPP_RTPredict_1")
add_test("TOPP_RTPredict_2" ${TOPP_BIN_PATH}/RTPredict -in_id ${DATA_DIR_TOPP}/RTPredict_2_input.idXML -out_id:positive RTPredict_2_output_positive.tmp -out_id:negative RTPredict_2_output_negative.tmp -svm_model ${DATA_DIR_TOPP}/RTPredict_2_input.model)
add_test("TOPP_RTPredict_2_out1" ${DIFF} -in1 RTPredict_2_output_positive.tmp -in2 ${DATA_DIR_TOPP}/RTPredict_2_output_positive.idXML )
set_tests_properties("TOPP_RTPredict_2_out1" PROPERTIES DEPENDS "TOPP_RTPredict_2")
add_test("TOPP_RTPredict_2_out2" ${DIFF} -in1 RTPredict_2_output_negative.tmp -in2 ${DATA_DIR_TOPP}/RTPredict_2_output_negative.idXML )
set_tests_properties("TOPP_RTPredict_2_out2" PROPERTIES DEPENDS "TOPP_RTPredict_2")
add_test("TOPP_RTPredict_3" ${TOPP_BIN_PATH}/RTPredict -in_id ${DATA_DIR_TOPP}/RTPredict_3_input.idXML -out_id:file RTPredict_3_output.tmp -total_gradient_time 1 -svm_model ${DATA_DIR_TOPP}/RTPredict_3_input.model)
add_test("TOPP_RTPredict_3_out1" ${DIFF} -whitelist "?xml-stylesheet" "predicted_RT_p_value" -in1 RTPredict_3_output.tmp -in2 ${DATA_DIR_TOPP}/RTPredict_3_output.idXML )
set_tests_properties("TOPP_RTPredict_3_out1" PROPERTIES DEPENDS "TOPP_RTPredict_3")
add_test("TOPP_RTPredict_4" ${TOPP_BIN_PATH}/RTPredict -in_text ${DATA_DIR_TOPP}/RTPredict_4_input.txt -out_text:file RTPredict_4_output.tmp -total_gradient_time 1 -svm_model ${DATA_DIR_TOPP}/RTPredict_4_input.model)
add_test("TOPP_RTPredict_4_out1" ${DIFF} -in1 RTPredict_4_output.tmp -in2 ${DATA_DIR_TOPP}/RTPredict_4_output.txt )
set_tests_properties("TOPP_RTPredict_4_out1" PROPERTIES DEPENDS "TOPP_RTPredict_4")
add_test("TOPP_RTPredict_5" ${TOPP_BIN_PATH}/RTPredict -in_id ${DATA_DIR_TOPP}/RTPredict_5_input.idXML -out_text:file RTPredict_5_output.tmp -total_gradient_time 1 -svm_model ${DATA_DIR_TOPP}/RTPredict_5_input.model)
add_test("TOPP_RTPredict_5_out1" ${DIFF} -in1 RTPredict_5_output.tmp -in2 ${DATA_DIR_TOPP}/RTPredict_5_output.txt )
set_tests_properties("TOPP_RTPredict_5_out1" PROPERTIES DEPENDS "TOPP_RTPredict_5")

### InspectAdapter tests
add_test("TOPP_InspectAdapter_1" ${TOPP_BIN_PATH}/InspectAdapter -ini ${DATA_DIR_TOPP}/InspectAdapter_1_parameters.ini -trie_dbs ${DATA_DIR_TOPP}/Inspect_FASTAFile_test2.trie -in ${DATA_DIR_TOPP}/Inspect.mzXML -dbs ${DATA_DIR_TOPP}/Inspect_FASTAFile_test.fasta -inspect_in -out InspectAdapter_2_output.tmp)
add_test("TOPP_InspectAdapter_1_out1" ${DIFF} -whitelist "spectra," "db," -in1 InspectAdapter_2_output.tmp -in2 ${DATA_DIR_TOPP}/InspectAdapter_2_output.inspect_in )
set_tests_properties("TOPP_InspectAdapter_1_out1" PROPERTIES DEPENDS "TOPP_InspectAdapter_1")

add_test("TOPP_InspectAdapter_2" ${TOPP_BIN_PATH}/InspectAdapter -ini ${DATA_DIR_TOPP}/InspectAdapter_1_parameters.ini -trie_dbs ${DATA_DIR_TOPP}/Inspect_FASTAFile_test2.trie -in ${DATA_DIR_TOPP}/Inspect.mzData -dbs ${DATA_DIR_TOPP}/Inspect_FASTAFile_test.fasta -inspect_in -out InspectAdapter_3_output.tmp)
add_test("TOPP_InspectAdapter_2_out1" ${DIFF} -whitelist "spectra," "db," -in1 InspectAdapter_3_output.tmp -in2 ${DATA_DIR_TOPP}/InspectAdapter_2_output.inspect_in )
set_tests_properties("TOPP_InspectAdapter_2_out1" PROPERTIES DEPENDS "TOPP_InspectAdapter_2")

# TODO the following tests are waiting for better implementations of InspectAdapter and associated classes
#add_test("TOPP_InspectAdapter_3" ${TOPP_BIN_PATH}/InspectAdapter -ini ${DATA_DIR_TOPP}/InspectAdapter_1_parameters.ini -trie_dbs ${DATA_DIR_TOPP}/Inspect_FASTAFile_test2.trie -in ${DATA_DIR_TOPP}/InspectAdapter.out -dbs ${DATA_DIR_TOPP}/Inspect_FASTAFile_test.fasta -out InspectAdapter_4_output.tmp -inspect_out)
#add_test("TOPP_InspectAdapter_3_out1" ${DIFF} -whitelist "?xml-stylesheet" "IdentificationRun date" -in1 InspectAdapter_4_output.tmp -in2 ${DATA_DIR_TOPP}/InspectAdapter_4_output.idXML )
#set_tests_properties("TOPP_InspectAdapter_3_out1" PROPERTIES DEPENDS "TOPP_InspectAdapter_3")

### SeedListGenerator tests
add_test("TOPP_SeedListGenerator_1" ${TOPP_BIN_PATH}/SeedListGenerator -test -in ${PROJECT_SOURCE_DIR}/data/PepXMLFile_test.mzML -out SeedListGenerator_1_output.tmp)
add_test("TOPP_SeedListGenerator_1_out1" ${DIFF} -in1 SeedListGenerator_1_output.tmp -in2 ${DATA_DIR_TOPP}/SeedListGenerator_1_output.featureXML )
set_tests_properties("TOPP_SeedListGenerator_1_out1" PROPERTIES DEPENDS "TOPP_SeedListGenerator_1")
add_test("TOPP_SeedListGenerator_2" ${TOPP_BIN_PATH}/SeedListGenerator -test -in ${DATA_DIR_TOPP}/IDMapper_1_output.featureXML -out SeedListGenerator_2_output.tmp)
add_test("TOPP_SeedListGenerator_2_out1" ${DIFF} -in1 SeedListGenerator_2_output.tmp -in2 ${DATA_DIR_TOPP}/SeedListGenerator_2_output.featureXML )
set_tests_properties("TOPP_SeedListGenerator_2_out1" PROPERTIES DEPENDS "TOPP_SeedListGenerator_2")
add_test("TOPP_SeedListGenerator_3" ${TOPP_BIN_PATH}/SeedListGenerator -test -in ${PROJECT_SOURCE_DIR}/data/ConsensusXMLFile_1.consensusXML -out SeedListGenerator_3_output1.tmp SeedListGenerator_3_output2.tmp)
add_test("TOPP_SeedListGenerator_3_out1" ${DIFF} -in1 SeedListGenerator_3_output1.tmp -in2 ${DATA_DIR_TOPP}/SeedListGenerator_3_output1.featureXML )
set_tests_properties("TOPP_SeedListGenerator_3_out1" PROPERTIES DEPENDS "TOPP_SeedListGenerator_3")
add_test("TOPP_SeedListGenerator_3_out2" ${DIFF} -in1 SeedListGenerator_3_output2.tmp -in2 ${DATA_DIR_TOPP}/SeedListGenerator_3_output2.featureXML )
set_tests_properties("TOPP_SeedListGenerator_3_out2" PROPERTIES DEPENDS "TOPP_SeedListGenerator_3")

### SequestAdapter tests
# TODO disabled until tool is reactivated
#add_test("TOPP_SequestAdapter_1" ${TOPP_BIN_PATH}/SequestAdapter -ini ${DATA_DIR_TOPP}/SequestAdapter_1_parameters.ini -in ${DATA_DIR_TOPP}/Sequest.mzXML -mz_files ${DATA_DIR_TOPP}/Sequest.mzXML -modifications_xml_file ${DATA_DIR_TOPP}/Sequest_PTMs.xml -out SequestAdapter_2_output.tmp -sequest_in -temp_data_directory ${DATA_DIR_TOPP} -db ${DATA_DIR_TOPP}/Sequest_FASTAFile_test.fasta)
#add_test("TOPP_SequestAdapter_1_out1" ${DIFF} -in1 SequestAdapter_2_output.tmp -in2 ${DATA_DIR_TOPP}/SequestAdapter_2_output.sequest_in)
#add_test("TOPP_SequestAdapter_2" ${TOPP_BIN_PATH}/SequestAdapter -ini ${DATA_DIR_TOPP}/SequestAdapter_1_parameters.ini -in ${DATA_DIR_TOPP}/Sequest.mzData -mz_files ${DATA_DIR_TOPP}/Sequest.mzXML -modifications_xml_file ${DATA_DIR_TOPP}/Sequest_PTMs.xml -out SequestAdapter_3_output.tmp -sequest_in -temp_data_directory ${DATA_DIR_TOPP} -db ${DATA_DIR_TOPP}/Sequest_FASTAFile_test.fasta)
#add_test("TOPP_SequestAdapter_2_out1" ${DIFF} -in1 SequestAdapter_3_output.tmp -in2 ${DATA_DIR_TOPP}/SequestAdapter_2_output.sequest_in)

# TODO the following tests are waiting for better implementations of InspectAdapter and
# associated classes
#add_test("TOPP_SequestAdapter_3" ${TOPP_BIN_PATH}/SequestAdapter -ini ${DATA_DIR_TOPP}/SequestAdapter_2_parameters.ini -mz_files ${DATA_DIR_TOPP}/Sequest.mzXML -modifications_xml_file ${DATA_DIR_TOPP}/Sequest_PTMs.xml -in ${DATA_DIR_TOPP}/Sequest.mzXML -out SequestAdapter_4_output.tmp -temp_data_directory ${DATA_DIR_TOPP} -db ${DATA_DIR_TOPP}/Sequest_FASTAFile_test.fasta)
#add_test("TOPP_SequestAdapter_3_out1" ${DIFF} -in1 SequestAdapter_4_output.tmp -in2 ${DATA_DIR_TOPP}/SequestAdapter_4_output.idXML)

### SpecLibSearcher tests
#add_test("TOPP_SpecLibSearcher_1" ${TOPP_BIN_PATH}/SpecLibSearcher -test -ini ${DATA_DIR_TOPP}/SpecLibSearcher_1_parameters.ini -in ${DATA_DIR_TOPP}/SpecLibSearcher_1.MzData -lib $(DATA_DIR_TOPP)/SpecLibSearcher_1.MSP -out SpecLibSearcher_1.tmp)
#add_test("TOPP_SpecLibSearcher_1_out1" ${DIFF} -in1 SpecLibSearcher_1.tmp  -in2 $(DATA_DIR_TOPP)/SpecLibSearcher_1.idXML)
### PepNovoAdapter tests
#The PepNovoAdapter now only works as a frontend and cannot be run without an installation of PepNovo.Therefore no test possible
#add_test("TOPP_PepNovoAdapter_1" ${TOPP_BIN_PATH}/PepNovoAdapter -ini ${DATA_DIR_TOPP}/PepNovoAdapter_1_parameters.ini -in ${DATA_DIR_TOPP}/PepNovo.mzXML -pepnovo_in -out PepNovoAdapter_3_output.tmp -dta_list ${DATA_DIR_TOPP}/tmp/dta_list.txt -model_directory ${DATA_DIR_TOPP}/tmp/ -temp_data_directory ${DATA_DIR_TOPP}/tmp/)
#add_test("TOPP_PepNovoAdapter_1_out1" ${DIFF} -in1 ${DATA_DIR_TOPP}/tmp/PepNovo_PTMs_.txt -in2 ${DATA_DIR_TOPP}/tmp/PepNovo_PTMs.txt)
#TODO ANDREAS - We have to clean up the /tmp/ directory to run this test multiple times
#add_test("TOPP_PepNovoAdapter_2" ${TOPP_BIN_PATH}/PepNovoAdapter -ini ${DATA_DIR_TOPP}/PepNovoAdapter_1_parameters.ini -in ${DATA_DIR_TOPP}/PepNovo.mzData -pepnovo_in -out PepNovoAdapter_4_output.tmp -temp_data_directory ${DATA_DIR_TOPP})
#add_test("TOPP_PepNovoAdapter_2_out1" ${DIFF} -in1 ${DATA_DIR_TOPP}/PepNovo_PTMs_.txt -in2 ${DATA_DIR_TOPP}/PepNovo_PTMs.txt)
#add_test("TOPP_PepNovoAdapter_3" ${TOPP_BIN_PATH}/PepNovoAdapter -ini ${DATA_DIR_TOPP}/PepNovoAdapter_5_parameters.ini -in ${DATA_DIR_TOPP}/PepNovoAdapter_5_output.pepnovo_out -out PepNovoAdapter_5_output.tmp -pepnovo_out -dta_list ${DATA_DIR_TOPP}/tmp/dta_list.txt -model_directory ${DATA_DIR_TOPP}/tmp/ -temp_data_directory ${DATA_DIR_TOPP}/tmp/ -modifications_xml_file ${DATA_DIR_TOPP}/PepNovo_PTMs.xml -mz_files ${DATA_DIR_TOPP}/PepNovo.mzXML)
#add_test("TOPP_PepNovoAdapter_3_out1" ${DIFF} -whitelist "?xml-stylesheet" "date_group_1" -in1 PepNovoAdapter_5_output.tmp -in2 ${DATA_DIR_TOPP}/PepNovoAdapter_5_output.idXML)

### SpectraFilter tests
add_test("TOPP_SpectraFilterSqrtMower_1" ${TOPP_BIN_PATH}/SpectraFilterSqrtMower -test -in ${DATA_DIR_TOPP}/SpectraFilterSqrtMower_1_input.mzML -out SpectraFilterSqrtMower.tmp)
add_test("TOPP_SpectraFilterSqrtMower_1_out1" ${DIFF} -in1 SpectraFilterSqrtMower.tmp -in2 ${DATA_DIR_TOPP}/SpectraFilterSqrtMower_1_output.mzML )
set_tests_properties("TOPP_SpectraFilterSqrtMower_1_out1" PROPERTIES DEPENDS "TOPP_SpectraFilterSqrtMower_1")
add_test("TOPP_SpectraFilterWindowMower_1" ${TOPP_BIN_PATH}/SpectraFilterWindowMower -test -in ${DATA_DIR_TOPP}/SpectraFilterWindowMower_1_input.mzML -out SpectraFilterWindowMower_1.tmp)
add_test("TOPP_SpectraFilterWindowMower_1_out1" ${DIFF} -in1 SpectraFilterWindowMower_1.tmp -in2 ${DATA_DIR_TOPP}/SpectraFilterWindowMower_1_output.mzML )
set_tests_properties("TOPP_SpectraFilterWindowMower_1_out1" PROPERTIES DEPENDS "TOPP_SpectraFilterWindowMower_1")
add_test("TOPP_SpectraFilterWindowMower_2" ${TOPP_BIN_PATH}/SpectraFilterWindowMower -test -in ${DATA_DIR_TOPP}/SpectraFilterWindowMower_2_input.mzML -out SpectraFilterWindowMower_2.tmp -ini ${DATA_DIR_TOPP}/SpectraFilterWindowMower_2_parameters.ini)
add_test("TOPP_SpectraFilterWindowMower_2_out1" ${DIFF} -in1 SpectraFilterWindowMower_2.tmp -in2 ${DATA_DIR_TOPP}/SpectraFilterWindowMower_2_output.mzML )
set_tests_properties("TOPP_SpectraFilterWindowMower_2_out1" PROPERTIES DEPENDS "TOPP_SpectraFilterWindowMower_2")

### InternalCalibration tests
add_test("TOPP_InternalCalibration_1" ${TOPP_BIN_PATH}/InternalCalibration -test -in ${DATA_DIR_TOPP}/InternalCalibration_2_input.mzML -out InternalCalibration.tmp -ref_peaks ${DATA_DIR_TOPP}/InternalCalibration_ref_masses -type spectrumwise)
#add_test("TOPP_InternalCalibration_1_out1" ${DIFF} -in1 InternalCalibration.tmp -in2 ${DATA_DIR_TOPP}/InternalCalibration_2_output.mzML )
#set_tests_properties("TOPP_InternalCalibration_1_out1" PROPERTIES DEPENDS "TOPP_InternalCalibration_1")

### TMTAnalyzer
add_test("TOPP_TMTAnalyzer_1" ${TOPP_BIN_PATH}/TMTAnalyzer -test -ini ${DATA_DIR_TOPP}/TMTAnalyzer.ini -in ${DATA_DIR_TOPP}/TMTAnalyzer_input_1.mzML -out TMTAnalyzer_1.tmp)
add_test("TOPP_TMTAnalyzer_1_out1" ${DIFF} -whitelist "<map" "?xml-stylesheet" -in1 TMTAnalyzer_1.tmp -in2 ${DATA_DIR_TOPP}/TMTAnalyzer_output_1.consensusXML )
set_tests_properties("TOPP_TMTAnalyzer_1_out1" PROPERTIES DEPENDS "TOPP_TMTAnalyzer_1")

### TOFCalibration tests
add_test("TOPP_TOFCalibration_1" ${TOPP_BIN_PATH}/TOFCalibration -test -in ${DATA_DIR_TOPP}/TOFCalibration_1_input.mzML -out TOFCalibration.tmp -ref_masses ${DATA_DIR_TOPP}/TOFCalibration_ref_masses -ini ${DATA_DIR_TOPP}/TOFCalibration_parameters.ini -tof_const ${DATA_DIR_TOPP}/TOFCalibration_const -ext_calibrants ${DATA_DIR_TOPP}/TOFCalibration_1_calibrants.mzML)
#add_test("TOPP_TOFCalibration_1_out1" ${DIFF} -in1 TOFCalibration.tmp -in2 ${DATA_DIR_TOPP}/TOFCalibration_1_output.mzML )
#set_tests_properties("TOPP_TOFCalibration_1_out1" PROPERTIES DEPENDS "TOPP_TOFCalibration_1")
add_test("TOPP_TOFCalibration_2" ${TOPP_BIN_PATH}/TOFCalibration -test -in ${DATA_DIR_TOPP}/TOFCalibration_2_input.mzML -out TOFCalibration.tmp -ref_masses ${DATA_DIR_TOPP}/TOFCalibration_ref_masses -ini ${DATA_DIR_TOPP}/TOFCalibration_parameters.ini -tof_const ${DATA_DIR_TOPP}/TOFCalibration_const -ext_calibrants ${DATA_DIR_TOPP}/TOFCalibration_2_calibrants.mzML -peak_data)
#add_test("TOPP_TOFCalibration_2_out1" ${DIFF} -in1 TOFCalibration.tmp -in2 ${DATA_DIR_TOPP}/TOFCalibration_2_output.mzML )
#set_tests_properties("TOPP_TOFCalibration_2_out1" PROPERTIES DEPENDS "TOPP_TOFCalibration_2")

### TextExporter tests
add_test("TOPP_TextExporter_1" ${TOPP_BIN_PATH}/TextExporter -in ${DATA_DIR_TOPP}/TextExporter_1_input.featureXML -no_progress -out TextExporter_1_output.tmp)
add_test("TOPP_TextExporter_1_out1" ${DIFF} -in1 TextExporter_1_output.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_1_output.txt )
set_tests_properties("TOPP_TextExporter_1_out1" PROPERTIES DEPENDS "TOPP_TextExporter_1")

add_test("TOPP_TextExporter_2" ${TOPP_BIN_PATH}/TextExporter -in ${DATA_DIR_TOPP}/TextExporter_2_input.consensusXML -ini ${DATA_DIR_TOPP}/TextExporter_2_parameters.ini)
add_test("TOPP_TextExporter_2_out1" ${DIFF} -whitelist "extracted from" -in1 TextExporter_2_consensus_centroids.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_2_consensus_centroids.wsv )
set_tests_properties("TOPP_TextExporter_2_out1" PROPERTIES DEPENDS "TOPP_TextExporter_2")
add_test("TOPP_TextExporter_2_out2" ${DIFF} -whitelist "extracted from" -in1 TextExporter_2_consensus_elements.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_2_consensus_elements.wsv )
set_tests_properties("TOPP_TextExporter_2_out2" PROPERTIES DEPENDS "TOPP_TextExporter_2")
add_test("TOPP_TextExporter_2_out3" ${DIFF} -whitelist "extracted from" -in1 TextExporter_2_consensus_features.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_2_consensus_features.wsv )
set_tests_properties("TOPP_TextExporter_2_out3" PROPERTIES DEPENDS "TOPP_TextExporter_2")
add_test("TOPP_TextExporter_2_out4" ${DIFF} -whitelist "extracted from" -in1 TextExporter_2_consensus_tsv.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_2_consensus_tsv.tsv )
set_tests_properties("TOPP_TextExporter_2_out4" PROPERTIES DEPENDS "TOPP_TextExporter_2")

add_test("TOPP_TextExporter_3" ${TOPP_BIN_PATH}/TextExporter -in ${DATA_DIR_TOPP}/TextExporter_3_input.idXML -no_progress -out TextExporter_3_output.tmp)
add_test("TOPP_TextExporter_3_out1" ${DIFF} -in1 TextExporter_3_output.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_3_output.txt )
set_tests_properties("TOPP_TextExporter_3_out1" PROPERTIES DEPENDS "TOPP_TextExporter_3")

add_test("TOPP_TextExporter_4" ${TOPP_BIN_PATH}/TextExporter -in ${DATA_DIR_TOPP}/TextExporter_3_input.idXML -no_progress -out TextExporter_4_output.tmp -proteins_only)
add_test("TOPP_TextExporter_4_out1" ${DIFF} -in1 TextExporter_4_output.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_4_output_proteins.txt )
set_tests_properties("TOPP_TextExporter_4_out1" PROPERTIES DEPENDS "TOPP_TextExporter_4")

add_test("TOPP_TextExporter_5" ${TOPP_BIN_PATH}/TextExporter -in ${DATA_DIR_TOPP}/TextExporter_5_input.idXML -no_progress -out TextExporter_5_output.tmp -peptides_only -first_dim_rt)
add_test("TOPP_TextExporter_5_out1" ${DIFF} -in1 TextExporter_5_output.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_5_output_peptides.txt )
set_tests_properties("TOPP_TextExporter_5_out1" PROPERTIES DEPENDS "TOPP_TextExporter_5")

add_test("TOPP_TextExporter_6" ${TOPP_BIN_PATH}/TextExporter -in ${DATA_DIR_TOPP}/TextExporter_6_input.featureXML -no_progress -out TextExporter_6_output.tmp -no_ids)
add_test("TOPP_TextExporter_6_out1" ${DIFF} -in1 TextExporter_6_output.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_6_output.txt )
set_tests_properties("TOPP_TextExporter_6_out1" PROPERTIES DEPENDS "TOPP_TextExporter_6")

add_test("TOPP_TextExporter_7" ${TOPP_BIN_PATH}/TextExporter -in ${DATA_DIR_TOPP}/TextExporter_7_input.consensusXML -ini ${DATA_DIR_TOPP}/TextExporter_7_parameters.ini)
add_test("TOPP_TextExporter_7_out1" ${DIFF} -whitelist "extracted from" -in1 TextExporter_7_consensus_tsv.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_7_consensus_tsv.tsv )
set_tests_properties("TOPP_TextExporter_7_out1" PROPERTIES DEPENDS "TOPP_TextExporter_7")

add_test("TOPP_TextExporter_8" ${TOPP_BIN_PATH}/TextExporter -in ${DATA_DIR_TOPP}/TextExporter_1_input.featureXML -no_ids -no_progress -out TextExporter_8_output.tmp)
add_test("TOPP_TextExporter_8_out1" ${DIFF} -in1 TextExporter_8_output.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_8_output.txt )
set_tests_properties("TOPP_TextExporter_8_out1" PROPERTIES DEPENDS "TOPP_TextExporter_8")

### SILACAnalyzer tests
# mzML
add_test("TOPP_SILACAnalyzer_1" ${TOPP_BIN_PATH}/SILACAnalyzer -test -ini ${DATA_DIR_TOPP}/SILACAnalyzer.ini -in ${DATA_DIR_TOPP}/SILACAnalyzer.mzML -out SILACAnalyzer.consensusXML.tmp -out_features SILACAnalyzer.featureXML.tmp)
add_test("TOPP_SILACAnalyzer_1_out1" ${DIFF} -in1 SILACAnalyzer.consensusXML.tmp -in2 ${DATA_DIR_TOPP}/SILACAnalyzer.consensusXML )
set_tests_properties("TOPP_SILACAnalyzer_1_out1" PROPERTIES DEPENDS "TOPP_SILACAnalyzer_1")
add_test("TOPP_SILACAnalyzer_1_out2" ${DIFF} -in1 SILACAnalyzer.featureXML.tmp -in2 ${DATA_DIR_TOPP}/SILACAnalyzer.featureXML )
set_tests_properties("TOPP_SILACAnalyzer_1_out2" PROPERTIES DEPENDS "TOPP_SILACAnalyzer_1")
# empty mzML
add_test("TOPP_SILACAnalyzer_2" ${TOPP_BIN_PATH}/SILACAnalyzer -test -ini ${DATA_DIR_TOPP}/SILACAnalyzer.ini -in ${DATA_DIR_TOPP}/SILACAnalyzer_empty.mzML)
set_tests_properties("TOPP_SILACAnalyzer_2" PROPERTIES WILL_FAIL 1)
# mzML with negativ RT
add_test("TOPP_SILACAnalyzer_3" ${TOPP_BIN_PATH}/SILACAnalyzer -test -ini ${DATA_DIR_TOPP}/SILACAnalyzer.ini -in ${DATA_DIR_TOPP}/SILACAnalyzer_negativ.mzML -out SILACAnalyzer_negativ.consensusXML.tmp)
add_test("TOPP_SILACAnalyzer_3_out1" ${DIFF} -in1 SILACAnalyzer_negativ.consensusXML.tmp -in2 ${DATA_DIR_TOPP}/SILACAnalyzer_negativ.consensusXML )
set_tests_properties("TOPP_SILACAnalyzer_3_out1" PROPERTIES DEPENDS "TOPP_SILACAnalyzer_3")
add_test("TOPP_FeatureFinderRaw_1" ${TOPP_BIN_PATH}/FeatureFinderRaw -test -ini ${DATA_DIR_TOPP}/FeatureFinderRaw.ini -in ${DATA_DIR_TOPP}/FeatureFinderRaw_1.mzML -out FeatureFinderRaw_1.featureXML.tmp)
add_test("TOPP_FeatureFinderRaw_1_out" ${DIFF} -in1 FeatureFinderRaw_1.featureXML.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinderRaw_1.featureXML )
set_tests_properties("TOPP_FeatureFinderRaw_1_out" PROPERTIES DEPENDS "TOPP_FeatureFinderRaw_1")


### FeatureLinker tests
# "labeled" algorithm:
add_test("TOPP_FeatureLinkerLabeled_1" ${TOPP_BIN_PATH}/FeatureLinkerLabeled -test -ini ${DATA_DIR_TOPP}/FeatureLinkerLabeled_1_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerLabeled_1_input.featureXML -out FeatureLinkerLabeled_1_output.tmp)
add_test("TOPP_FeatureLinkerLabeled_1_out1" ${DIFF} -in1 FeatureLinkerLabeled_1_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinkerLabeled_1_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerLabeled_1_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerLabeled_1")
add_test("TOPP_FeatureLinkerLabeled_2" ${TOPP_BIN_PATH}/FeatureLinkerLabeled -test -ini ${DATA_DIR_TOPP}/FeatureLinkerLabeled_2_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerLabeled_2_input.featureXML -out FeatureLinkerLabeled_2_output.tmp)
add_test("TOPP_FeatureLinkerLabeled_2_out1" ${DIFF} -in1 FeatureLinkerLabeled_2_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinkerLabeled_2_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerLabeled_2_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerLabeled_2")
# "unlabeled" algorithm:
add_test("TOPP_FeatureLinkerUnlabeled_1" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeled -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input2.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input3.featureXML -out FeatureLinkerUnlabeled_1_output.tmp)
add_test("TOPP_FeatureLinkerUnlabeled_1_out1" ${DIFF} -in1 FeatureLinkerUnlabeled_1_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerUnlabeled_1_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeled_1")
add_test("TOPP_FeatureLinkerUnlabeled_2" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeled -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_2_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_2_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_2_input2.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_2_input3.featureXML -out FeatureLinkerUnlabeled_2_output.tmp)
add_test("TOPP_FeatureLinkerUnlabeled_2_out1" ${DIFF} -in1 FeatureLinkerUnlabeled_2_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_2_output.consensusXML)
set_tests_properties("TOPP_FeatureLinkerUnlabeled_2_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeled_2")
add_test("TOPP_FeatureLinkerUnlabeled_3" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeled -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_3_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_3_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_3_input2.featureXML -out FeatureLinkerUnlabeled_3_output.tmp)
add_test("TOPP_FeatureLinkerUnlabeled_3_out1" ${DIFF} -in1 FeatureLinkerUnlabeled_3_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_3_output.consensusXML)
set_tests_properties("TOPP_FeatureLinkerUnlabeled_3_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeled_3")
# "unlabeled" algorithm with "use_identifications" option:
add_test("TOPP_FeatureLinkerUnlabeled_4" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeled -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_4_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input2.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input3.featureXML -out FeatureLinkerUnlabeled_4_output.tmp)
add_test("TOPP_FeatureLinkerUnlabeled_4_out1" ${DIFF} -in1 FeatureLinkerUnlabeled_4_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_4_output.consensusXML)
set_tests_properties("TOPP_FeatureLinkerUnlabeled_4_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeled_4")
# "unlabeled_qt" algorithm:
add_test("TOPP_FeatureLinkerUnlabeledQT_1" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeledQT -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_1_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input2.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input3.featureXML -out FeatureLinkerUnlabeledQT_1_output.tmp)
add_test("TOPP_FeatureLinkerUnlabeledQT_1_out1" ${DIFF} -in1 FeatureLinkerUnlabeledQT_1_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_1_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerUnlabeledQT_1_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeledQT_1")
# "unlabeled_qt" algorithm with consensusXML input and "keep_subelements" flag:
add_test("TOPP_FeatureLinkerUnlabeledQT_2" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeledQT -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_2_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_1_output.consensusXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_1_output.consensusXML -out FeatureLinkerUnlabeledQT_2_output.tmp)
add_test("TOPP_FeatureLinkerUnlabeledQT_2_out1" ${DIFF} -in1 FeatureLinkerUnlabeledQT_2_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_2_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerUnlabeledQT_2_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeledQT_2")

### IDMapper tests
add_test("TOPP_IDMapper_1" ${TOPP_BIN_PATH}/IDMapper -test -in ${DATA_DIR_TOPP}/IDMapper_1_input.featureXML -out IDMapper_1_output.tmp -id ${DATA_DIR_TOPP}/IDMapper_1_input.idXML -mz_measure Da -mz_tolerance 1 -ignore_charge)
add_test("TOPP_IDMapper_1_out1" ${DIFF} -in1 IDMapper_1_output.tmp -in2 ${DATA_DIR_TOPP}/IDMapper_1_output.featureXML )
set_tests_properties("TOPP_IDMapper_1_out1" PROPERTIES DEPENDS "TOPP_IDMapper_1")
add_test("TOPP_IDMapper_2" ${TOPP_BIN_PATH}/IDMapper -test -in ${DATA_DIR_TOPP}/IDMapper_2_input.consensusXML -out IDMapper_2_output.tmp -id ${DATA_DIR_TOPP}/IDMapper_2_input.idXML  -mz_measure Da -mz_tolerance 1)
add_test("TOPP_IDMapper_2_out1" ${DIFF} -in1 IDMapper_2_output.tmp -in2 ${DATA_DIR_TOPP}/IDMapper_2_output.consensusXML )
set_tests_properties("TOPP_IDMapper_2_out1" PROPERTIES DEPENDS "TOPP_IDMapper_2")
add_test("TOPP_IDMapper_3" ${TOPP_BIN_PATH}/IDMapper -test -in ${DATA_DIR_TOPP}/IDMapper_3_input.featureXML -out IDMapper_3_output.tmp -id ${DATA_DIR_TOPP}/IDMapper_3_input.idXML  -mz_measure ppm -mz_tolerance 3 -rt_tolerance 4 -ignore_charge)
add_test("TOPP_IDMapper_3_out1" ${DIFF} -in1 IDMapper_3_output.tmp -in2 ${DATA_DIR_TOPP}/IDMapper_3_output.featureXML )
set_tests_properties("TOPP_IDMapper_3_out1" PROPERTIES DEPENDS "TOPP_IDMapper_3")

### IDRipper tests
add_test("TOPP_IDRipper_1" ${TOPP_BIN_PATH}/IDRipper -test -in ${DATA_DIR_TOPP}/IDRipper_1_input.idXML -out dummy.tmp )
add_test("TOPP_IDRipper_1_out1" ${DIFF} -in1 IDRipper_1_output_1.tmp -in2 ${DATA_DIR_TOPP}/IDRipper_1_output1.idXML )
set_tests_properties("TOPP_IDRipper_1_out1" PROPERTIES DEPENDS "TOPP_IDRipper_1")
add_test("TOPP_IDRipper_1_out2" ${DIFF} -in1 IDRipper_1_output_2.tmp -in2 ${DATA_DIR_TOPP}/IDRipper_1_output2.idXML )
set_tests_properties("TOPP_IDRipper_1_out2" PROPERTIES DEPENDS "TOPP_IDRipper_1")
add_test("TOPP_IDRipper_2" ${TOPP_BIN_PATH}/IDRipper -test -in ${DATA_DIR_TOPP}/IDRipper_2_input.idXML -out_path dummy.tmp)
add_test("TOPP_IDRipper_2_out1" ${DIFF} -in1 IDRipper_2_output_1.tmp -in2 ${DATA_DIR_TOPP}/IDRipper_2_output1.idXML )
set_tests_properties("TOPP_IDRipper_2_out1" PROPERTIES DEPENDS "TOPP_IDRipper_2")
add_test("TOPP_IDRipper_2_out2" ${DIFF} -in1 IDRipper_2_output_2.tmp -in2 ${DATA_DIR_TOPP}/IDRipper_2_output2.idXML )
set_tests_properties("TOPP_IDRipper_2_out2" PROPERTIES DEPENDS "TOPP_IDRipper_2")
# use IDMerger results to test. Input -> IDMerger -> IDRipper -> Output == Input
add_test("TOPP_IDRipper_3_prep" ${TOPP_BIN_PATH}/IDMerger -test -in ${DATA_DIR_TOPP}/IDRipper_3_input1.idXML ${DATA_DIR_TOPP}/IDRipper_3_input2.idXML -out IDRipper_3_output.tmp -annotate_file_origin)
add_test("TOPP_IDRipper_3_prep_out1" ${DIFF} -whitelist "?xml-stylesheet" "file_origin" -in1 IDRipper_3_output.tmp -in2 ${DATA_DIR_TOPP}/IDRipper_3_output.idXML )
set_tests_properties("TOPP_IDRipper_3_prep_out1" PROPERTIES DEPENDS "TOPP_IDRipper_3_prep")
add_test("TOPP_IDRipper_3" ${TOPP_BIN_PATH}/IDRipper -test -in IDRipper_3_output.tmp -out dummy.tmp)
add_test("TOPP_IDRipper_3_out1" ${DIFF} -in1 IDRipper_3_input1.idXML -in2 ${DATA_DIR_TOPP}/IDRipper_3_input1.idXML )
set_tests_properties("TOPP_IDRipper_3_out1" PROPERTIES DEPENDS "TOPP_IDRipper_3")
add_test("TOPP_IDRipper_3_out2" ${DIFF} -in1 IDRipper_3_input2.idXML -in2 ${DATA_DIR_TOPP}/IDRipper_3_input2.idXML )
set_tests_properties("TOPP_IDRipper_3_out2" PROPERTIES DEPENDS "TOPP_IDRipper_3")


### ConsensusID tests
add_test("TOPP_ConsensusID_1" ${TOPP_BIN_PATH}/ConsensusID -ini ${DATA_DIR_TOPP}/ConsensusID_1_parameters.ini -in ${DATA_DIR_TOPP}/ConsensusID_1_input.idXML -out ConsensusID_1_output.tmp)
add_test("TOPP_ConsensusID_1_out1" ${DIFF} -whitelist "?xml-stylesheet" "IdentificationRun date" -in1 ConsensusID_1_output.tmp -in2 ${DATA_DIR_TOPP}/ConsensusID_1_output.idXML )
set_tests_properties("TOPP_ConsensusID_1_out1" PROPERTIES DEPENDS "TOPP_ConsensusID_1")
add_test("TOPP_ConsensusID_2" ${TOPP_BIN_PATH}/ConsensusID -ini ${DATA_DIR_TOPP}/ConsensusID_2_parameters.ini -out ConsensusID_2_output.tmp -in ${DATA_DIR_TOPP}/ConsensusID_2_input.featureXML)
add_test("TOPP_ConsensusID_2_out1" ${DIFF} -whitelist "IdentificationRun id" -in1 ConsensusID_2_output.tmp -in2 ${DATA_DIR_TOPP}/ConsensusID_2_output.featureXML )
set_tests_properties("TOPP_ConsensusID_2_out1" PROPERTIES DEPENDS "TOPP_ConsensusID_2")
add_test("TOPP_ConsensusID_3" ${TOPP_BIN_PATH}/ConsensusID -ini ${DATA_DIR_TOPP}/ConsensusID_3_parameters.ini -out ConsensusID_3_output.tmp -in ${DATA_DIR_TOPP}/ConsensusID_3_input.consensusXML)
add_test("TOPP_ConsensusID_3_out1" ${DIFF} -whitelist "?xml-stylesheet" "IdentificationRun id" -in1 ConsensusID_3_output.tmp -in2 ${DATA_DIR_TOPP}/ConsensusID_3_output.consensusXML )
set_tests_properties("TOPP_ConsensusID_3_out1" PROPERTIES DEPENDS "TOPP_ConsensusID_3")


### PrecursorIonSelector tests
add_test("TOPP_PrecursorIonSelector_1" ${TOPP_BIN_PATH}/PrecursorIonSelector -in ${DATA_DIR_TOPP}/PrecursorIonSelector_features.featureXML -sim_results PrecursorIonSelector_1_output.tmp -ini ${DATA_DIR_TOPP}/PrecursorIonSelector_1_parameters.ini -db_path ${DATA_DIR_TOPP}/PrecursorIonSelector_db.fasta -ids ${DATA_DIR_TOPP}/PrecursorIonSelector_ids.IdXML)
add_test("TOPP_PrecursorIonSelector_1_out1" ${DIFF} -in1 PrecursorIonSelector_1_output.tmp -in2 ${DATA_DIR_TOPP}/PrecursorIonSelector_1_output.txt )
set_tests_properties("TOPP_PrecursorIonSelector_1_out1" PROPERTIES DEPENDS "TOPP_PrecursorIonSelector_1")
add_test("TOPP_PrecursorIonSelector_2" ${TOPP_BIN_PATH}/PrecursorIonSelector -in ${DATA_DIR_TOPP}/PrecursorIonSelector_features.featureXML -out PrecursorIonSelector_2_output.tmp -ini ${DATA_DIR_TOPP}/PrecursorIonSelector_2_parameters.ini -db_path ${DATA_DIR_TOPP}/PrecursorIonSelector_db.fasta -ids ${DATA_DIR_TOPP}/PrecursorIonSelector_ids.IdXML)
add_test("TOPP_PrecursorIonSelector_2_out1" ${DIFF} -in1 PrecursorIonSelector_2_output.tmp -in2 ${DATA_DIR_TOPP}/PrecursorIonSelector_2_output.featureXML )
set_tests_properties("TOPP_PrecursorIonSelector_2_out1" PROPERTIES DEPENDS "TOPP_PrecursorIonSelector_2")

### CompNovo tests
add_test("TOPP_CompNovo_1" ${TOPP_BIN_PATH}/CompNovo -in ${DATA_DIR_TOPP}/CompNovo_1_input.mzML -ini ${DATA_DIR_TOPP}/CompNovo_1.ini -out CompNovo_1_output.tmp)
add_test("TOPP_CompNovo_1_out1" ${DIFF} -whitelist "?xml-stylesheet" "IdentificationRun date" -in1 CompNovo_1_output.tmp -in2 ${DATA_DIR_TOPP}/CompNovo_1_output.idXML )
set_tests_properties("TOPP_CompNovo_1_out1" PROPERTIES DEPENDS "TOPP_CompNovo_1")
add_test("TOPP_CompNovoCID_1" ${TOPP_BIN_PATH}/CompNovoCID -in ${DATA_DIR_TOPP}/CompNovoCID_1_input.mzML -ini ${DATA_DIR_TOPP}/CompNovoCID_1.ini -out CompNovoCID_1_output.tmp)
add_test("TOPP_CompNovoCID_1_out1" ${DIFF} -whitelist "?xml-stylesheet" "IdentificationRun date" -in1 CompNovoCID_1_output.tmp -in2 ${DATA_DIR_TOPP}/CompNovoCID_1_output.idXML )
set_tests_properties("TOPP_CompNovoCID_1_out1" PROPERTIES DEPENDS "TOPP_CompNovoCID_1")

### PrecursorMassCorrector tests
add_test("TOPP_PrecursorMassCorrector_1" ${TOPP_BIN_PATH}/PrecursorMassCorrector -in ${DATA_DIR_TOPP}/PrecursorMassCorrector_1_input.mzML -ini ${DATA_DIR_TOPP}/PrecursorMassCorrector_1.ini -out PrecursorMassCorrector_1_output.tmp)
add_test("TOPP_PrecursorMassCorrector_1_out1" ${DIFF} -in1 PrecursorMassCorrector_1_output.tmp -in2 ${DATA_DIR_TOPP}/PrecursorMassCorrector_1_output.mzML )
set_tests_properties("TOPP_PrecursorMassCorrector_1_out1" PROPERTIES DEPENDS "TOPP_PrecursorMassCorrector_1")

### ProteinInference tests
add_test("TOPP_ProteinInference_1" ${TOPP_BIN_PATH}/ProteinInference -test -in ${DATA_DIR_TOPP}/ProteinInference.idXML -out ProteinInference_output.tmp)
add_test("TOPP_ProteinInference_1_out1" ${DIFF} -in1 ProteinInference_output.tmp -in2 ${DATA_DIR_TOPP}/ProteinInference_proteins.idXML )
set_tests_properties("TOPP_ProteinInference_1_out1" PROPERTIES DEPENDS "TOPP_ProteinInference_1")

### FalseDiscoveryRate_tests (1: default param, 2: with runs separated, 3: with charge variants individual, 4: with charge variants and no charge 1 decoys -> q_value 0)
add_test("TOPP_FalseDiscoveryRate_1" ${TOPP_BIN_PATH}/FalseDiscoveryRate -test  -ini ${DATA_DIR_TOPP}/FalseDiscoveryRate_1.ini -in ${DATA_DIR_TOPP}/FalseDiscoveryRate_OMSSA.idXML -out FalseDiscoveryRate_output_1.tmp)
add_test("TOPP_FalseDiscoveryRate_1_out1" ${DIFF} -whitelist "?xml-stylesheet" -in1 FalseDiscoveryRate_output_1.tmp -in2 ${DATA_DIR_TOPP}/FalseDiscoveryRate_output_1.idXML )
set_tests_properties("TOPP_FalseDiscoveryRate_1_out1" PROPERTIES DEPENDS "TOPP_FalseDiscoveryRate_1")
add_test("TOPP_FalseDiscoveryRate_2" ${TOPP_BIN_PATH}/FalseDiscoveryRate -test  -ini ${DATA_DIR_TOPP}/FalseDiscoveryRate_2.ini -in ${DATA_DIR_TOPP}/FalseDiscoveryRate_OMSSA.idXML -out FalseDiscoveryRate_output_2.tmp)
add_test("TOPP_FalseDiscoveryRate_2_out1" ${DIFF} -whitelist "?xml-stylesheet" -in1 FalseDiscoveryRate_output_2.tmp -in2 ${DATA_DIR_TOPP}/FalseDiscoveryRate_output_2.idXML )
set_tests_properties("TOPP_FalseDiscoveryRate_2_out1" PROPERTIES DEPENDS "TOPP_FalseDiscoveryRate_2")
add_test("TOPP_FalseDiscoveryRate_3" ${TOPP_BIN_PATH}/FalseDiscoveryRate -test  -ini ${DATA_DIR_TOPP}/FalseDiscoveryRate_3.ini -in ${DATA_DIR_TOPP}/FalseDiscoveryRate_OMSSA.idXML -out FalseDiscoveryRate_output_3.tmp)
add_test("TOPP_FalseDiscoveryRate_3_out1" ${DIFF} -whitelist "?xml-stylesheet" -in1 FalseDiscoveryRate_output_3.tmp -in2 ${DATA_DIR_TOPP}/FalseDiscoveryRate_output_3.idXML )
set_tests_properties("TOPP_FalseDiscoveryRate_3_out1" PROPERTIES DEPENDS "TOPP_FalseDiscoveryRate_3")
add_test("TOPP_FalseDiscoveryRate_4" ${TOPP_BIN_PATH}/FalseDiscoveryRate -test  -ini ${DATA_DIR_TOPP}/FalseDiscoveryRate_3.ini -in ${DATA_DIR_TOPP}/FalseDiscoveryRate_OMSSA_4.idXML -out FalseDiscoveryRate_output_4.tmp)
add_test("TOPP_FalseDiscoveryRate_4_out1" ${DIFF} -whitelist "?xml-stylesheet" -in1 FalseDiscoveryRate_output_4.tmp -in2 ${DATA_DIR_TOPP}/FalseDiscoveryRate_output_4.idXML )
set_tests_properties("TOPP_FalseDiscoveryRate_4_out1" PROPERTIES DEPENDS "TOPP_FalseDiscoveryRate_4")

###IDPosteriorErrorProbability
add_test("TOPP_IDPosteriorErrorProbability_1" ${TOPP_BIN_PATH}/IDPosteriorErrorProbability -test -in ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_Mascot_input.idXML -out 	IDPosteriorErrorProbability_output_1.tmp)
add_test("TOPP_IDPosteriorErrorProbability_1_out1" ${DIFF} -in1 IDPosteriorErrorProbability_output_1.tmp -in2 ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_Mascot_output.idXML )
set_tests_properties("TOPP_IDPosteriorErrorProbability_1_out1" PROPERTIES DEPENDS "TOPP_IDPosteriorErrorProbability_1")
add_test("TOPP_IDPosteriorErrorProbability_2" ${TOPP_BIN_PATH}/IDPosteriorErrorProbability -test -in ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_XTandem_input.idXML -out IDPosteriorErrorProbability_output_2.tmp)
add_test("TOPP_IDPosteriorErrorProbability_2_out1" ${DIFF} -in1 IDPosteriorErrorProbability_output_2.tmp -in2 ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_XTandem_output.idXML )
set_tests_properties("TOPP_IDPosteriorErrorProbability_2_out1" PROPERTIES DEPENDS "TOPP_IDPosteriorErrorProbability_2")
add_test("TOPP_IDPosteriorErrorProbability_3" ${TOPP_BIN_PATH}/IDPosteriorErrorProbability -test -in ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_OMSSA_input.idXML -out IDPosteriorErrorProbability_output_3.tmp)
add_test("TOPP_IDPosteriorErrorProbability_3_out1" ${DIFF} -in1 IDPosteriorErrorProbability_output_3.tmp -in2 ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_OMSSA_output.idXML )
set_tests_properties("TOPP_IDPosteriorErrorProbability_3_out1" PROPERTIES DEPENDS "TOPP_IDPosteriorErrorProbability_3")
add_test("TOPP_IDPosteriorErrorProbability_4" ${TOPP_BIN_PATH}/IDPosteriorErrorProbability -test -in ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_OMSSA_input2.idXML -out IDPosteriorErrorProbability_output_4.tmp -split_charge)
add_test("TOPP_IDPosteriorErrorProbability_4_out1" ${DIFF} -in1 IDPosteriorErrorProbability_output_4.tmp -in2 ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_OMSSA_output2.idXML )
set_tests_properties("TOPP_IDPosteriorErrorProbability_4_out1" PROPERTIES DEPENDS "TOPP_IDPosteriorErrorProbability_4")
add_test("TOPP_IDPosteriorErrorProbability_5" ${TOPP_BIN_PATH}/IDPosteriorErrorProbability -test -in ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_XTandem_input2.idXML -out IDPosteriorErrorProbability_output_5.tmp -split_charge)
add_test("TOPP_IDPosteriorErrorProbability_5_out1" ${DIFF} -in1 IDPosteriorErrorProbability_output_5.tmp -in2 ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_XTandem_output2.idXML )
set_tests_properties("TOPP_IDPosteriorErrorProbability_5_out1" PROPERTIES DEPENDS "TOPP_IDPosteriorErrorProbability_5")
add_test("TOPP_IDPosteriorErrorProbability_6" ${TOPP_BIN_PATH}/IDPosteriorErrorProbability -test -in ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_Mascot_input2.idXML -out IDPosteriorErrorProbability_output_6.tmp -split_charge)
add_test("TOPP_IDPosteriorErrorProbability_6_out1" ${DIFF} -in1 IDPosteriorErrorProbability_output_6.tmp -in2 ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_Mascot_output2.idXML )
set_tests_properties("TOPP_IDPosteriorErrorProbability_6_out1" PROPERTIES DEPENDS "TOPP_IDPosteriorErrorProbability_6")
add_test("TOPP_IDPosteriorErrorProbability_7" ${TOPP_BIN_PATH}/IDPosteriorErrorProbability -test -in ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_bad_data.idXML -out IDPosteriorErrorPorbability_bad_data.idXML -ignore_bad_data)

### ProteinResolver tests
### TODO test for idXML
add_test("TOPP_ProteinResolver_1" ${TOPP_BIN_PATH}/ProteinResolver -test -in ${DATA_DIR_TOPP}/ProteinResolver_1_input.consensusXML -fasta ${DATA_DIR_TOPP}/ProteinResolver_1_input.fasta -protein_groups ProteinResolver_1_output1.tmp -peptide_table ProteinResolver_1_output2.tmp -protein_table ProteinResolver_1_output3.tmp)
add_test("TOPP_ProteinResolver_1_out1" ${DIFF} -in1 ProteinResolver_1_output1.tmp -in2 ${DATA_DIR_TOPP}/ProteinResolver_1_output1.txt)
set_tests_properties("TOPP_ProteinResolver_1_out1" PROPERTIES DEPENDS "TOPP_ProteinResolver_1")
add_test("TOPP_ProteinResolver_1_out2" ${DIFF} -in1 ProteinResolver_1_output2.tmp -in2 ${DATA_DIR_TOPP}/ProteinResolver_1_output2.txt)
set_tests_properties("TOPP_ProteinResolver_1_out2" PROPERTIES DEPENDS "TOPP_ProteinResolver_1")
add_test("TOPP_ProteinResolver_1_out3" ${DIFF} -in1 ProteinResolver_1_output3.tmp -in2 ${DATA_DIR_TOPP}/ProteinResolver_1_output3.txt)
set_tests_properties("TOPP_ProteinResolver_1_out3" PROPERTIES DEPENDS "TOPP_ProteinResolver_1")
#add_test("TOPP_ProteinResolver_2" ${TOPP_BIN_PATH}/ProteinResolver -test -in ${DATA_DIR_TOPP}/ProteinResolver_2_input.consensusXML -fasta ${DATA_DIR_TOPP}/ProteinResolver_2_input.fasta -design ${DATA_DIR_TOPP}/ProteinResolver_2_input.txt -protein_groups ProteinResolver_2_output1.tmp -peptide_table ProteinResolver_2_output2.tmp -protein_table ProteinResolver_2_output3.tmp)
#add_test("TOPP_ProteinResolver_2_out1" ${DIFF} -in1 ProteinResolver_2_output1.tmp -in2 ${DATA_DIR_TOPP}/ProteinResolver_2_output1.txt)
#set_tests_properties("TOPP_ProteinResolver_1_out1" PROPERTIES DEPENDS "TOPP_ProteinResolver_1")
#add_test("TOPP_ProteinResolver_2_out2" ${DIFF} -in1 ProteinResolver_2_output2.tmp -in2 ${DATA_DIR_TOPP}/ProteinResolver_2_output2.txt)
#set_tests_properties("TOPP_ProteinResolver_1_out1" PROPERTIES DEPENDS "TOPP_ProteinResolver_1")
#add_test("TOPP_ProteinResolver_2_out3" ${DIFF} -in1 ProteinResolver_2_output3.tmp -in2 ${DATA_DIR_TOPP}/ProteinResolver_2_output3.txt)
#set_tests_properties("TOPP_ProteinResolver_1_out1" PROPERTIES DEPENDS "TOPP_ProteinResolver_1")

### ProteinQuantifier tests
add_test("TOPP_ProteinQuantifier_1" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.featureXML -out ProteinQuantifier_1_output1.tmp -peptide_out ProteinQuantifier_1_output2.tmp)
add_test("TOPP_ProteinQuantifier_1_out1" ${DIFF} -in1 ProteinQuantifier_1_output1.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_1_output1.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_1_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_1")
add_test("TOPP_ProteinQuantifier_1_out2" ${DIFF} -in1 ProteinQuantifier_1_output2.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_1_output2.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_1_out2" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_1")
add_test("TOPP_ProteinQuantifier_2" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.featureXML -out ProteinQuantifier_2_output1.tmp -peptide_out ProteinQuantifier_2_output2.tmp -top 2 -include_all -filter_charge -average sum)
add_test("TOPP_ProteinQuantifier_2_out1" ${DIFF} -in1 ProteinQuantifier_2_output1.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_2_output1.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_2_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_2")
add_test("TOPP_ProteinQuantifier_2_out2" ${DIFF} -in1 ProteinQuantifier_2_output2.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_2_output2.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_2_out2" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_2")
# handling of modified peptides:
add_test("TOPP_ProteinQuantifier_3" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_3_input.featureXML -out ProteinQuantifier_3_output1.tmp -peptide_out ProteinQuantifier_3_output2.tmp -top 2 -include_all -average mean)
add_test("TOPP_ProteinQuantifier_3_out1" ${DIFF} -in1 ProteinQuantifier_3_output1.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_3_output1.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_3_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_3")
add_test("TOPP_ProteinQuantifier_3_out2" ${DIFF} -in1 ProteinQuantifier_3_output2.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_3_output2.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_3_out2" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_3")
# consensusXML input, different options for protein quantification:
add_test("TOPP_ProteinQuantifier_4" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.consensusXML -out ProteinQuantifier_4_output.tmp -top 0 -average sum)
add_test("TOPP_ProteinQuantifier_4_out1" ${DIFF} -in1 ProteinQuantifier_4_output.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_4_output.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_4_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_4")
add_test("TOPP_ProteinQuantifier_5" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.consensusXML -out ProteinQuantifier_5_output.tmp -top 3 -average sum)
add_test("TOPP_ProteinQuantifier_5_out1" ${DIFF} -in1 ProteinQuantifier_5_output.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_5_output.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_5_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_5")
add_test("TOPP_ProteinQuantifier_6" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.consensusXML -out ProteinQuantifier_6_output.tmp -top 3 -include_all -average sum)
add_test("TOPP_ProteinQuantifier_6_out1" ${DIFF} -in1 ProteinQuantifier_6_output.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_6_output.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_6_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_6")
add_test("TOPP_ProteinQuantifier_7" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.consensusXML -out ProteinQuantifier_7_output.tmp -top 0 -consensus:fix_peptides -average sum)
add_test("TOPP_ProteinQuantifier_7_out1" ${DIFF} -in1 ProteinQuantifier_7_output.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_7_output.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_7_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_7")
add_test("TOPP_ProteinQuantifier_8" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.consensusXML -out ProteinQuantifier_8_output.tmp -top 3 -consensus:fix_peptides -average sum)
add_test("TOPP_ProteinQuantifier_8_out1" ${DIFF} -in1 ProteinQuantifier_8_output.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_8_output.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_8_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_8")
add_test("TOPP_ProteinQuantifier_9" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.consensusXML -out ProteinQuantifier_9_output.tmp -top 3 -include_all -consensus:fix_peptides -average sum)
add_test("TOPP_ProteinQuantifier_9_out1" ${DIFF} -in1 ProteinQuantifier_9_output.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_9_output.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_9_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_9")
# mzTab output:
add_test("TOPP_ProteinQuantifier_10" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.consensusXML -mzTab_out ProteinQuantifier_10_output.tmp -top 0 -average sum)
add_test("TOPP_ProteinQuantifier_10_out1" ${DIFF} -in1 ProteinQuantifier_10_output.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_10_output.mzTab)
set_tests_properties("TOPP_ProteinQuantifier_10_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_10")
add_test("TOPP_ProteinQuantifier_11" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.featureXML -mzTab_out ProteinQuantifier_11_output.tmp -top 0 -average sum -filter_charge)
add_test("TOPP_ProteinQuantifier_11_out1" ${DIFF} -in1 ProteinQuantifier_11_output.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_11_output.mzTab)
set_tests_properties("TOPP_ProteinQuantifier_11_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_11")
#ratios
add_test("TOPP_ProteinQuantifier_12" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.consensusXML -out ProteinQuantifier_12_output.tmp -top 3 -include_all -consensus:fix_peptides -average sum -ratios)
add_test("TOPP_ProteinQuantifier_12_out1" ${DIFF} -in1 ProteinQuantifier_12_output.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_12_output.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_12_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_12")
add_test("TOPP_ProteinQuantifier_13" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.consensusXML -out ProteinQuantifier_13_output.tmp -top 3 -include_all -consensus:fix_peptides -average sum -ratiosSILAC)
add_test("TOPP_ProteinQuantifier_13_out1" ${DIFF} -in1 ProteinQuantifier_13_output.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_13_output.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_13_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_13")

### InclusionExclusionListCreator tests
add_test("TOPP_InclusionExclusionListCreator_1" ${TOPP_BIN_PATH}/InclusionExclusionListCreator -test -include ${DATA_DIR_TOPP}/InclusionExclusionListCreator.featureXML -out InclusionExclusionListCreator_1_output.tmp)
add_test("TOPP_InclusionExclusionListCreator_1_out1" ${DIFF} -in1 InclusionExclusionListCreator_1_output.tmp -in2 ${DATA_DIR_TOPP}/InclusionExclusionListCreator_1_output.txt )
set_tests_properties("TOPP_InclusionExclusionListCreator_1_out1" PROPERTIES DEPENDS "TOPP_InclusionExclusionListCreator_1")
add_test("TOPP_InclusionExclusionListCreator_2" ${TOPP_BIN_PATH}/InclusionExclusionListCreator -test -exclude ${DATA_DIR_TOPP}/InclusionExclusionListCreator.featureXML -out InclusionExclusionListCreator_2_output.tmp -ini ${DATA_DIR_TOPP}/InclusionExclusionListCreator_2.ini)
add_test("TOPP_InclusionExclusionListCreator_2_out1" ${DIFF} -in1 InclusionExclusionListCreator_2_output.tmp -in2 ${DATA_DIR_TOPP}/InclusionExclusionListCreator_2_output.txt )
set_tests_properties("TOPP_InclusionExclusionListCreator_2_out1" PROPERTIES DEPENDS "TOPP_InclusionExclusionListCreator_2")
add_test("TOPP_InclusionExclusionListCreator_3" ${TOPP_BIN_PATH}/InclusionExclusionListCreator -test -include ${DATA_DIR_TOPP}/InclusionExclusionListCreator.fasta -out InclusionExclusionListCreator_3_output.tmp -inclusion_charges 1 2 -rt_model ${DATA_DIR_TOPP}/InclusionExclusionListCreator_rt.model)
add_test("TOPP_InclusionExclusionListCreator_3_out1" ${DIFF} -in1 InclusionExclusionListCreator_3_output.tmp -in2 ${DATA_DIR_TOPP}/InclusionExclusionListCreator_3_output.txt )
set_tests_properties("TOPP_InclusionExclusionListCreator_3_out1" PROPERTIES DEPENDS "TOPP_InclusionExclusionListCreator_3")
add_test("TOPP_InclusionExclusionListCreator_4" ${TOPP_BIN_PATH}/InclusionExclusionListCreator -test -exclude ${DATA_DIR_TOPP}/InclusionExclusionListCreator.fasta -out InclusionExclusionListCreator_4_output.tmp -exclusion_charges 1 2 -rt_model ${DATA_DIR_TOPP}/InclusionExclusionListCreator_rt.model  -ini ${DATA_DIR_TOPP}/InclusionExclusionListCreator_4.ini)
add_test("TOPP_InclusionExclusionListCreator_4_out1" ${DIFF} -in1 InclusionExclusionListCreator_4_output.tmp -in2 ${DATA_DIR_TOPP}/InclusionExclusionListCreator_4_output.txt )
set_tests_properties("TOPP_InclusionExclusionListCreator_4_out1" PROPERTIES DEPENDS "TOPP_InclusionExclusionListCreator_4")
add_test("TOPP_InclusionExclusionListCreator_5" ${TOPP_BIN_PATH}/InclusionExclusionListCreator -test -exclude ${DATA_DIR_TOPP}/InclusionExclusionListCreator.IdXML -out InclusionExclusionListCreator_5_output.tmp -exclusion_charges 1 2 -rt_model ${DATA_DIR_TOPP}/InclusionExclusionListCreator_rt.model)
add_test("TOPP_InclusionExclusionListCreator_5_out1" ${DIFF} -in1 InclusionExclusionListCreator_5_output.tmp -in2 ${DATA_DIR_TOPP}/InclusionExclusionListCreator_5_output.txt )
set_tests_properties("TOPP_InclusionExclusionListCreator_5_out1" PROPERTIES DEPENDS "TOPP_InclusionExclusionListCreator_5")
add_test("TOPP_InclusionExclusionListCreator_6" ${TOPP_BIN_PATH}/InclusionExclusionListCreator -test -exclude ${DATA_DIR_TOPP}/InclusionExclusionListCreator.IdXML -out InclusionExclusionListCreator_6_output.tmp -exclusion_charges 1 2 -rt_model ${DATA_DIR_TOPP}/InclusionExclusionListCreator_rt.model  -ini ${DATA_DIR_TOPP}/InclusionExclusionListCreator_6.ini)
add_test("TOPP_InclusionExclusionListCreator_6_out1" ${DIFF} -in1 InclusionExclusionListCreator_6_output.tmp -in2 ${DATA_DIR_TOPP}/InclusionExclusionListCreator_6_output.txt )
set_tests_properties("TOPP_InclusionExclusionListCreator_6_out1" PROPERTIES DEPENDS "TOPP_InclusionExclusionListCreator_6")

### PeptideIndexer tests
add_test("TOPP_PeptideIndexer_1" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_1.idXML -out PeptideIndexer_1_out.tmp -allow_unmatched)
add_test("TOPP_PeptideIndexer_1_out" ${DIFF} -in1 PeptideIndexer_1_out.tmp -in2 ${DATA_DIR_TOPP}/PeptideIndexer_1_out.idXML )
set_tests_properties("TOPP_PeptideIndexer_1_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_1")
add_test("TOPP_PeptideIndexer_2" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_1.idXML -out PeptideIndexer_2_out.tmp -allow_unmatched -write_protein_sequence)
add_test("TOPP_PeptideIndexer_2_out" ${DIFF} -in1 PeptideIndexer_2_out.tmp -in2 ${DATA_DIR_TOPP}/PeptideIndexer_2_out.idXML )
set_tests_properties("TOPP_PeptideIndexer_2_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_2")
add_test("TOPP_PeptideIndexer_3" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_1.idXML -out PeptideIndexer_3_out.tmp -allow_unmatched -keep_unreferenced_proteins)
add_test("TOPP_PeptideIndexer_3_out" ${DIFF} -in1 PeptideIndexer_3_out.tmp -in2 ${DATA_DIR_TOPP}/PeptideIndexer_3_out.idXML )
set_tests_properties("TOPP_PeptideIndexer_3_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_3")
add_test("TOPP_PeptideIndexer_4" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_1.idXML -out PeptideIndexer_4_out.tmp -allow_unmatched -aaa_max 0 -write_protein_sequence)
add_test("TOPP_PeptideIndexer_4_out" ${DIFF} -in1 PeptideIndexer_4_out.tmp -in2 ${DATA_DIR_TOPP}/PeptideIndexer_4_out.idXML )
set_tests_properties("TOPP_PeptideIndexer_4_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_4")

### ExecutePipeline tests (as substitute for TOPPAS) - the ResourceFiles are in binary tree, as they have been configured from a .in file (see above)!
add_test("TOPP_ExecutePipeline_1" ${TOPP_BIN_PATH}/ExecutePipeline -test -in ${DATA_DIR_TOPP}/ExecutePipeline_1.toppas -resource_file ${DATA_DIR_TOPP_BIN}/ExecutePipeline_1.trf -out_dir .)
set_tests_properties("TOPP_ExecutePipeline_1" PROPERTIES WILL_FAIL 1)

### UTILS tests (if this becomes too much, we should create a separate directory)
add_test("UTILS_ImageCreator_1" ${TOPP_BIN_PATH}/ImageCreator -test -in ${DATA_DIR_TOPP}/ImageCreator_1_input.mzML -out ImageCreator1.bmp -width 20 -height 15)
add_test("UTILS_ImageCreator_1_out1" ${DIFF} -in1 ImageCreator1.bmp -in2 ${DATA_DIR_TOPP}/ImageCreator_1_output.bmp )
set_tests_properties("UTILS_ImageCreator_1_out1" PROPERTIES DEPENDS "UTILS_ImageCreator_1")
# testing forced out_type
add_test("UTILS_ImageCreator_2" ${TOPP_BIN_PATH}/ImageCreator -test -in ${DATA_DIR_TOPP}/ImageCreator_2_input.mzML -out ImageCreator2.png -out_type BMP -width 20 -height 15 -precursors -precursor_size 3 -precursor_color green -log_intensity)
add_test("UTILS_ImageCreator_2_out1" ${DIFF} -in1 ImageCreator2.png -in2 ${DATA_DIR_TOPP}/ImageCreator_2_output.bmp )
set_tests_properties("UTILS_ImageCreator_2_out1" PROPERTIES DEPENDS "UTILS_ImageCreator_2")

add_test("UTILS_IDSplitter_1" ${TOPP_BIN_PATH}/IDSplitter -test -in ${DATA_DIR_TOPP}/IDMapper_1_output.featureXML -out IDSplitter_1_output1.tmp -id_out IDSplitter_1_output2.tmp)
add_test("UTILS_IDSplitter_1_out1" ${DIFF} -in1 IDSplitter_1_output1.tmp -in2 ${DATA_DIR_TOPP}/IDSplitter_1_output.featureXML )
set_tests_properties("UTILS_IDSplitter_1_out1" PROPERTIES DEPENDS "UTILS_IDSplitter_1")
add_test("UTILS_IDSplitter_1_out2" ${DIFF} -in1 IDSplitter_1_output2.tmp -in2 ${DATA_DIR_TOPP}/IDSplitter_1_output.idXML )
set_tests_properties("UTILS_IDSplitter_1_out2" PROPERTIES DEPENDS "UTILS_IDSplitter_1")

add_test("UTILS_MassCalculator_1" ${TOPP_BIN_PATH}/MassCalculator -test -in ${DATA_DIR_TOPP}/MassCalculator_1_input.txt -out MassCalculator_1_output.tmp -charge 0 1 -format table -separator ,)
add_test("UTILS_MassCalculator_1_out1" ${DIFF} -in1 MassCalculator_1_output.tmp -in2 ${DATA_DIR_TOPP}/MassCalculator_1_output.csv )
set_tests_properties("UTILS_MassCalculator_1_out1" PROPERTIES DEPENDS "UTILS_MassCalculator_1")
add_test("UTILS_MassCalculator_2" ${TOPP_BIN_PATH}/MassCalculator -test -in_seq "LDQWLC(Carbamidomethyl)EKL" "(Glu->pyro-Glu)EAM(Oxidation)APKHK" "RANVM(Oxidation)DYR" "FGVEQDVDMVFASFIR" -out MassCalculator_2_output.tmp -charge 1 2 3)
add_test("UTILS_MassCalculator_2_out1" ${DIFF} -in1 MassCalculator_2_output.tmp -in2 ${DATA_DIR_TOPP}/MassCalculator_2_output.txt )
set_tests_properties("UTILS_MassCalculator_2_out1" PROPERTIES DEPENDS "UTILS_MassCalculator_2")


add_test("UTILS_INIUpdater_1" ${TOPP_BIN_PATH}/INIUpdater -test -in ${DATA_DIR_TOPP}/INIUpdater_1_noupdate.toppas -out INIUpdater_1_noupdate.toppas.tmp)
add_test("UTILS_INIUpdater_1_out" ${DIFF} -in1 INIUpdater_1_noupdate.toppas.tmp -in2 ${DATA_DIR_TOPP}/INIUpdater_1_noupdate.toppas )
set_tests_properties("UTILS_INIUpdater_1_out" PROPERTIES DEPENDS "UTILS_INIUpdater_1")

add_test("UTILS_INIUpdater_2" ${TOPP_BIN_PATH}/INIUpdater -test -in ${DATA_DIR_TOPP}/INIUpdater_2_broken.toppas -out INIUpdater_2_broken.toppas.tmp)
set_tests_properties("UTILS_INIUpdater_2" PROPERTIES WILL_FAIL 1)

add_test("UTILS_INIUpdater_3" ${TOPP_BIN_PATH}/INIUpdater -test -in ${DATA_DIR_TOPP}/INIUpdater_3_old.toppas -out INIUpdater_3_old.toppas.tmp)
add_test("UTILS_INIUpdater_3_out" ${DIFF} -in1 INIUpdater_3_old.toppas.tmp -in2 ${DATA_DIR_TOPP}/INIUpdater_3_old_out.toppas )
set_tests_properties("UTILS_INIUpdater_3_out" PROPERTIES DEPENDS "UTILS_INIUpdater_3")


add_test("UTILS_DecoyDatabase_1" ${TOPP_BIN_PATH}/DecoyDatabase -test -in ${DATA_DIR_TOPP}/DecoyDatabase_1.fasta -out DecoyDatabase_1.fasta.tmp)
add_test("UTILS_DecoyDatabase_1_out" ${DIFF} -in1 DecoyDatabase_1.fasta.tmp -in2 ${DATA_DIR_TOPP}/DecoyDatabase_1_out.fasta )
set_tests_properties("UTILS_DecoyDatabase_1_out" PROPERTIES DEPENDS "UTILS_DecoyDatabase_1")

add_test("UTILS_DecoyDatabase_2" ${TOPP_BIN_PATH}/DecoyDatabase -test -in ${DATA_DIR_TOPP}/DecoyDatabase_1.fasta -out DecoyDatabase_2.fasta.tmp -decoy_string "blabla" -decoy_string_position "prefix" -append)
add_test("UTILS_DecoyDatabase_2_out" ${DIFF} -in1 DecoyDatabase_2.fasta.tmp -in2 ${DATA_DIR_TOPP}/DecoyDatabase_2_out.fasta )
set_tests_properties("UTILS_DecoyDatabase_2_out" PROPERTIES DEPENDS "UTILS_DecoyDatabase_21")
