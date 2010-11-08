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
#
# Warning: When passing parameters avoid usage of "<" and ">" signs, as this will confuse CDASH when 
#          interpreting the test results file (which happens to be XML)
#
# Hint: Most TOPP tools add DataProcessing information to their output.
#       The -test flag ensures that this output is independent of local path settings etc, so use it when *testing* TOPP tools!

### General settings
set(TOPP_BIN_PATH ${OPENMS_BINARY_DIR})
set(DATA_DIR_TOPP ${PROJECT_SOURCE_DIR}/TOPP)
set(DATA_DIR_TOPP_BIN ${PROJECT_BINARY_DIR}/TOPP)
set(DIFF ${TOPP_BIN_PATH}/FuzzyDiff -ini ${DATA_DIR_TOPP}/FuzzyDiff.ini)

### test if -write_wdsl works
add_test("TOPP_WSDL_1" ${TOPP_BIN_PATH}/AdditiveSeries -write_wsdl WSDL_out.tmp)
add_test("TOPP_WSDL_1_out1" ${DIFF} -in1 WSDL_out.tmp -in2 ${DATA_DIR_TOPP}/WSDL_out.wsdl)

### test if -write_ini works (TODO: configure names of TOPP tools and types automatically)
set(ALL_WRITEINI_TOOLS
AdditiveSeries
BaselineFilter
"CompNovo -type CompNovo"
ConsensusID
DBExporter
DBImporter
DTAExtractor
Decharger
ExecutePipeline
FalseDiscoveryRate
"FeatureFinder -type centroided"
"FeatureLinker -type unlabeled"
"FeatureLinker -type unlabeled_qt"
FileConverter
FileFilter
FileInfo
FileMerger
#GenericWrapper ##disabled for now, as we need an automatic way to determine its types, changing frequently
IDDecoyProbability
IDPosteriorErrorProbability
IDFileConverter
IDFilter
IDMapper
IDMerger
IDRTCalibration
#INIFileEditor
"ITRAQAnalyzer -type 4plex"
InspectAdapter
InternalCalibration
"MapAligner -type pose_clustering_affine"
MapNormalizer
MascotAdapter
MascotAdapterOnline
"NoiseFilter -type sgolay"
OMSSAAdapter
PILISIdentification
"PILISModel -type training"
PTModel
PTPredict
"PeakPicker -type wavelet"
PepNovoAdapter
PeptideIndexer
PrecursorIonSelector
PrecursorMassCorrector
ProteinInference
ProteinQuantifier
RTModel
RTPredict
Resampler
"SILACAnalyzer -type double"
SeedListGenerator
#SequestAdapter #(reenable when fixed)
SpecLibSearcher
"SpectraFilter -type BernNorm"
TOFCalibration
#TOPPAS
#TOPPView
TextExporter
XTandemAdapter
InclusionExclusionListCreator
SpectraMerger

### UTILS ####
CVInspector
CaapConvert
DeMeanderize
DecoyDatabase
Digestor
DigestorMotif
ERPairFinder
FFEval
FuzzyDiff
#HistView (not derived from TOPPBase)
IDExtractor
IDMassAccuracy
IDSplitter
IdXMLEvaluation
ImageCreator
LabeledEval
MassCalculator
MRMPairFinder
"MSSimulator -type o18"
MapAlignmentEvaluation
OpenMSInfo
SemanticValidator
SequenceCoverageCalculator
SpecLibCreator
SpectrumGeneratorNetworkTrainer
#UniqueIdAssigner (not derived from TOPPBase)
XMLValidator

)
#####
## one test for each TOPP tool, checking -write_ini functionality:
#####

foreach (i ${ALL_WRITEINI_TOOLS})
  string(REGEX MATCH "^(.*) (.*) (.*)" has_type ${i}) ## search for +3 components in TOPP name (when it has a -type)
  if (has_type STREQUAL "") ##regex did not match
    add_test("TOPPWRITEINI_${i}" ${TOPP_BIN_PATH}/${i} -write_ini "${i}.ini.tmp")  
  else() ## regex matched
    set(toolname ${CMAKE_MATCH_1})
    set(args ${CMAKE_MATCH_2} ${CMAKE_MATCH_3})
    add_test("TOPPWRITEINI_${toolname}" ${TOPP_BIN_PATH}/${toolname} ${args} -write_ini "${toolname}.ini.tmp")
  endif()
endforeach(i)


##############################################################################################################################################################################
###### TOPP tests  ############
##############################################################################################################################################################################

### AdditiveSeries tests
add_test("TOPP_AdditiveSeries_1" ${TOPP_BIN_PATH}/AdditiveSeries -in ${DATA_DIR_TOPP}/AdditiveSeries_1_feat11.featureXML ${DATA_DIR_TOPP}/AdditiveSeries_1_feat20.featureXML ${DATA_DIR_TOPP}/AdditiveSeries_1_feat22.featureXML ${DATA_DIR_TOPP}/AdditiveSeries_1_feat34.featureXML -ini ${DATA_DIR_TOPP}/AdditiveSeries_1_parameters.ini)
add_test("TOPP_AdditiveSeries_1_out1" ${DIFF} -in1 AdditiveSeries_1_output.tmp -in2 ${DATA_DIR_TOPP}/AdditiveSeries_1_output.txt)
add_test("TOPP_AdditiveSeries_1_out2" ${DIFF} -in1 AdditiveSeries_1_gnuplot_tmp.cmd -in2 ${DATA_DIR_TOPP}/AdditiveSeries_1_gnuplot.cmd)
add_test("TOPP_AdditiveSeries_1_out3" ${DIFF} -in1 AdditiveSeries_1_gnuplot_tmp.dat -in2 ${DATA_DIR_TOPP}/AdditiveSeries_1_gnuplot.dat)
add_test("TOPP_AdditiveSeries_1_out4" ${DIFF} -in1 AdditiveSeries_1_gnuplot_tmp.err -in2 ${DATA_DIR_TOPP}/AdditiveSeries_1_gnuplot.err)

### IDMerger tests
add_test("TOPP_IDMerger_1" ${TOPP_BIN_PATH}/IDMerger -test -in ${DATA_DIR_TOPP}/IDMerger_1_input1.idXML ${DATA_DIR_TOPP}/IDMerger_1_input2.idXML -out IDMerger_1_output.tmp)
add_test("TOPP_IDMerger_1_out1" ${DIFF} -in1 IDMerger_1_output.tmp -in2 ${DATA_DIR_TOPP}/IDMerger_1_output.idXML)
add_test("TOPP_IDMerger_2" ${TOPP_BIN_PATH}/IDMerger -test -pepxml_protxml -in ${DATA_DIR_TOPP}/IDMerger_2_input1.idXML ${DATA_DIR_TOPP}/IDMerger_2_input2.idXML -out IDMerger_2_output.tmp)
add_test("TOPP_IDMerger_2_out1" ${DIFF} -in1 IDMerger_2_output.tmp -in2 ${DATA_DIR_TOPP}/IDMerger_2_output.idXML)
# conflicting IDs:
add_test("TOPP_IDMerger_3" ${TOPP_BIN_PATH}/IDMerger -test -in ${DATA_DIR_TOPP}/IDMerger_1_input1.idXML ${DATA_DIR_TOPP}/IDMerger_1_input1.idXML -out IDMerger_3_output.tmp)
add_test("TOPP_IDMerger_3_out1" ${DIFF} -in1 IDMerger_3_output.tmp -in2 ${DATA_DIR_TOPP}/IDMerger_3_output.idXML)

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
add_test("TOPP_Decharger_1" ${TOPP_BIN_PATH}/Decharger -test -in ${DATA_DIR_TOPP}/Decharger_input.featureXML -ini ${DATA_DIR_TOPP}/Decharger_input.ini -out_fm Decharger_output_fm.tmp -out_cm Decharger_output.tmp -outpairs Decharger_p_output.tmp)
add_test("TOPP_Decharger_1_out1" ${DIFF} -in1 Decharger_output.tmp -in2 ${DATA_DIR_TOPP}/Decharger_output.consensusXML)

### DTAExtractor tests
add_test("TOPP_DTAExtractor_1" ${TOPP_BIN_PATH}/DTAExtractor -in ${DATA_DIR_TOPP}/DTAExtractor_1_input.mzML -out DTAExtractor -rt :61)
add_test("TOPP_DTAExtractor_1_out1" ${DIFF} -in1 DTAExtractor_RT60.dta -in2 ${DATA_DIR_TOPP}/DTAExtractor_1_output.dta)
add_test("TOPP_DTAExtractor_2" ${TOPP_BIN_PATH}/DTAExtractor -in ${DATA_DIR_TOPP}/DTAExtractor_1_input.mzML -out DTAExtractor -level 1)
add_test("TOPP_DTAExtractor_2_out1" ${DIFF} -in1 DTAExtractor_RT60.dta -in2 ${DATA_DIR_TOPP}/DTAExtractor_2_output.dta)
add_test("TOPP_DTAExtractor_3" ${TOPP_BIN_PATH}/DTAExtractor -in ${DATA_DIR_TOPP}/DTAExtractor_1_input.mzML -out DTAExtractor -level 2 -mz :1000)
add_test("TOPP_DTAExtractor_3_out1" ${DIFF} -in1 DTAExtractor_RT140_MZ5.dta -in2 ${DATA_DIR_TOPP}/DTAExtractor_3_output.dta)

### FeatureFinder tests
# "simple" algorithm: deprecated
# add_test("TOPP_FeatureFinder_1" ${TOPP_BIN_PATH}/FeatureFinder -test -ini ${DATA_DIR_TOPP}/FeatureFinder_1_parameters.ini -in ${DATA_DIR_TOPP}/FeatureFinder_1_input.mzML -out FeatureFinder.tmp)
# add_test("TOPP_FeatureFinder_1_out1" ${DIFF} -in1 FeatureFinder.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinder_1_output.featureXML)
# "simplest" algorithm: deprecated
# add_test("TOPP_FeatureFinder_2" ${TOPP_BIN_PATH}/FeatureFinder -test -ini ${DATA_DIR_TOPP}/FeatureFinder_2_parameters.ini -in ${DATA_DIR_TOPP}/FeatureFinder_1_input.mzML -out FeatureFinder.tmp)
# add_test("TOPP_FeatureFinder_2_out1" ${DIFF} -in1 FeatureFinder.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinder_2_output.featureXML)
# "centroided" algorithm:
add_test("TOPP_FeatureFinder_4" ${TOPP_BIN_PATH}/FeatureFinder -test -ini ${DATA_DIR_TOPP}/FeatureFinder_4_parameters.ini -in ${DATA_DIR_TOPP}/FeatureFinder_4_input.mzML -out FeatureFinder.tmp)
add_test("TOPP_FeatureFinder_4_out1" ${DIFF} -in1 FeatureFinder.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinder_4_output.featureXML)
# "mrm" algorithm:
add_test("TOPP_FeatureFinder_5" ${TOPP_BIN_PATH}/FeatureFinder -test -in ${DATA_DIR_TOPP}/FeatureFinder_5_input.mzML -ini ${DATA_DIR_TOPP}/FeatureFinder_5_parameters.ini -out FeatureFinder.tmp)
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
add_test("TOPP_FileConverter_10" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_10_input.edta -no_progress -out FileConverter.tmp -out_type featureXML)
add_test("TOPP_FileConverter_10_out1" ${DIFF} -in1 FileConverter.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_10_output.featureXML)
#SpecArray
add_test("TOPP_FileConverter_11" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_11_input.pepList -no_progress -out FileConverter.tmp -out_type featureXML)
add_test("TOPP_FileConverter_11_out1" ${DIFF} -in1 FileConverter.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_11_output.featureXML)
#msInspect
add_test("TOPP_FileConverter_12" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_12_input.peptides.tsv -no_progress -out FileConverter.tmp -out_type featureXML)
add_test("TOPP_FileConverter_12_out1" ${DIFF} -in1 FileConverter.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_12_output.featureXML)
#kroenik
add_test("TOPP_FileConverter_13" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_13_input.peptides.kroenik -no_progress -out FileConverter.tmp -out_type featureXML)
add_test("TOPP_FileConverter_13_out1" ${DIFF} -in1 FileConverter.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_13_output.featureXML)
# featureXML to consensusXML
add_test("TOPP_FileConverter_14" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_9_output.featureXML -no_progress -out FileConverter.tmp -out_type consensusXML)
add_test("TOPP_FileConverter_14_out1" ${DIFF} -in1 FileConverter.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_14_output.consensusXML)


### FileFilter tests
add_test("TOPP_FileFilter_1" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out FileFilter.tmp -rt :30 -mz :1000 -int :20000 -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_1_out1" ${DIFF} -in1 FileFilter.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_1_output.mzML)
add_test("TOPP_FileFilter_2" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out FileFilter.tmp -rt 30: -mz 1000: -int 100: -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_2_out1" ${DIFF} -in1 FileFilter.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_2_output.mzML)
add_test("TOPP_FileFilter_3" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out FileFilter.tmp -level 2 -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_3_out1" ${DIFF} -in1 FileFilter.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_3_output.mzML)
add_test("TOPP_FileFilter_4" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_4_input.mzML -out FileFilter.tmp -remove_zoom -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_4_out1" ${DIFF} -in1 FileFilter.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_4_output.mzML)
add_test("TOPP_FileFilter_5" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_5_input.featureXML -out FileFilter.tmp -rt :1000 -mz :480 -int :79000 -charge :3 -q :0.6 -in_type featureXML -out_type featureXML)
add_test("TOPP_FileFilter_5_out1" ${DIFF} -in1 FileFilter.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_5_out.featureXML)
add_test("TOPP_FileFilter_6" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_5_input.featureXML -out FileFilter.tmp -rt 1000: -mz 440: -int 70000: -charge 3: -q 0.51: -in_type featureXML -out_type featureXML)
add_test("TOPP_FileFilter_6_out1" ${DIFF} -in1 FileFilter.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_6_out.featureXML)
add_test("TOPP_FileFilter_7" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_7_input.mzML -out FileFilter.tmp -int 7000: -level 1 2 3 -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_7_out1" ${DIFF} -in1 FileFilter.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_7_output.mzML)
add_test("TOPP_FileFilter_8" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_8_input.consensusXML -out FileFilter.tmp -rt 600:1400 -mz 700:2300 -int 1100:6000 -in_type consensusXML -out_type consensusXML)
add_test("TOPP_FileFilter_8_out1" ${DIFF} -in1 FileFilter.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_8_output.consensusXML)
add_test("TOPP_FileFilter_9" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_9_input.mzML -out FileFilter.tmp -remove_mode SelectedIonMonitoring -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_9_out1" ${DIFF} -in1 FileFilter.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_9_output.mzML)
add_test("TOPP_FileFilter_10" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_10_input.mzML -out FileFilter.tmp -remove_activation "Collision-induced dissociation" -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_10_out1" ${DIFF} -in1 FileFilter.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_10_output.mzML)
add_test("TOPP_FileFilter_11" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_11_input.mzML -out FileFilter.tmp -remove_activation "Plasma desorption" -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_11_out1" ${DIFF} -in1 FileFilter.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_11_output.mzML)
add_test("TOPP_FileFilter_12" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_12_input.mzML -out FileFilter.tmp -remove_chromatograms -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_12_out1" ${DIFF} -in1 FileFilter.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_12_output.mzML)
add_test("TOPP_FileFilter_13" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_13_input.consensusXML -out FileFilter.tmp -in_type consensusXML -out_type featureXML -map 2)
add_test("TOPP_FileFilter_13_out1" ${DIFF} -in1 FileFilter.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_13_output.featureXML)
add_test("TOPP_FileFilter_14" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_14_input.consensusXML -out FileFilter.tmp -in_type consensusXML -out_type consensusXML -map 0 2)
add_test("TOPP_FileFilter_14_out1" ${DIFF} -in1 FileFilter.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_14_output.consensusXML)

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
add_test("TOPP_FileMerger_7" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_7_input1.featureXML ${DATA_DIR_TOPP}/FileMerger_7_input2.featureXML -out FileMerger.tmp)
add_test("TOPP_FileMerger_7_out1" ${DIFF} -in1 FileMerger.tmp -in2 ${DATA_DIR_TOPP}/FileMerger_7_output.featureXML)


### IDRTCalibration tests
add_test("TOPP_IDRTCalibration_1" ${TOPP_BIN_PATH}/IDRTCalibration -in ${DATA_DIR_TOPP}/IDRTCalibration_1_input.IdXML -out IDRTCalibration_1_output.tmp -calibrant_1_input 10 -calibrant_2_input 90)
add_test("TOPP_IDRTCalibration_1_out1" ${DIFF} -in1 IDRTCalibration_1_output.tmp -in2 ${DATA_DIR_TOPP}/IDRTCalibration_1_output.IdXML)

### IDFileConverter tests
add_test("TOPP_IDFileConverter_1" ${TOPP_BIN_PATH}/IDFileConverter -in ${DATA_DIR_TOPP}/IDFileConverter_1_input_sequestoutdir -mz_file ${DATA_DIR_TOPP}/IDFileConverter_1_input_mzfile.mzML -out IDFileConverter_1_output.tmp -out_type idXML -ignore_proteins_per_peptide )
add_test("TOPP_IDFileConverter_1_out1" ${DIFF} -in1 IDFileConverter_1_output.tmp -in2 ${DATA_DIR_TOPP}/IDFileConverter_1_output.idXML)
add_test("TOPP_IDFileConverter_2" ${TOPP_BIN_PATH}/IDFileConverter -in ${PROJECT_SOURCE_DIR}/data/PepXMLFile_test.pepxml -mz_file ${PROJECT_SOURCE_DIR}/data/PepXMLFile_test.mzML -mz_name PepXMLFile_test -out IDFileConverter_2_output.tmp -out_type idXML)
add_test("TOPP_IDFileConverter_2_out1" ${DIFF} -in1 IDFileConverter_2_output.tmp -in2 ${DATA_DIR_TOPP}/IDFileConverter_2_output.idXML)
add_test("TOPP_IDFileConverter_3" ${TOPP_BIN_PATH}/IDFileConverter -in ${DATA_DIR_TOPP}/IDFileConverter_3_input.protXML -out IDFileConverter_3_output.tmp -out_type idXML)
add_test("TOPP_IDFileConverter_3_out1" ${DIFF} -in1 IDFileConverter_3_output.tmp -in2 ${DATA_DIR_TOPP}/IDFileConverter_3_output.idXML)

### IDFilter tests
add_test("TOPP_IDFilter_1" ${TOPP_BIN_PATH}/IDFilter -in ${DATA_DIR_TOPP}/IDFilter_1_input.idXML -out IDFilter_1_output.tmp -whitelist:proteins ${DATA_DIR_TOPP}/IDFilter_1_input.fas)
add_test("TOPP_IDFilter_1_out1" ${DIFF} -in1 IDFilter_1_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_1_output.idXML)
add_test("TOPP_IDFilter_2" ${TOPP_BIN_PATH}/IDFilter -in ${DATA_DIR_TOPP}/IDFilter_2_input.idXML -out IDFilter_2_output.tmp -thresh:pep 1 -thresh:prot 1)
add_test("TOPP_IDFilter_2_out1" ${DIFF} -in1 IDFilter_2_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_2_output.idXML)
add_test("TOPP_IDFilter_3" ${TOPP_BIN_PATH}/IDFilter -in ${DATA_DIR_TOPP}/IDFilter_3_input.idXML -out IDFilter_3_output.tmp -blacklist:peptides ${DATA_DIR_TOPP}/IDFilter_3_2_input.idXML)
add_test("TOPP_IDFilter_3_out1" ${DIFF} -in1 IDFilter_3_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_3_output.idXML)
add_test("TOPP_IDFilter_4" ${TOPP_BIN_PATH}/IDFilter -in ${DATA_DIR_TOPP}/IDFilter_4_input.idXML -out IDFilter_4_output.tmp -rt:p_value 0.08)
add_test("TOPP_IDFilter_4_out1" ${DIFF} -in1 IDFilter_4_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_4_output.idXML)
add_test("TOPP_IDFilter_5" ${TOPP_BIN_PATH}/IDFilter -in ${DATA_DIR_TOPP}/IDFilter_5_input.idXML -out IDFilter_5_output.tmp -score:pep 32 -score:prot 25)
add_test("TOPP_IDFilter_5_out1" ${DIFF} -in1 IDFilter_5_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_5_output.idXML)
add_test("TOPP_IDFilter_6" ${TOPP_BIN_PATH}/IDFilter -in ${DATA_DIR_TOPP}/IDFilter_6_input.idXML -out IDFilter_6_output.tmp -best:n_peptide_hits 2 -best:n_protein_hits 10)
add_test("TOPP_IDFilter_6_out1" ${DIFF} -in1 IDFilter_6_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_6_output.idXML)
add_test("TOPP_IDFilter_7" ${TOPP_BIN_PATH}/IDFilter -in ${DATA_DIR_TOPP}/IDFilter_7_input.idXML -out IDFilter_7_output.tmp -unique)
add_test("TOPP_IDFilter_7_out1" ${DIFF} -in1 IDFilter_7_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_7_output.idXML)

### MapAligner tests
# "pose_clustering_affine" algorithm:
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
# "spectrum_alignment" algorithm:
add_test("TOPP_MapAligner_3" ${TOPP_BIN_PATH}/MapAligner -test -ini ${DATA_DIR_TOPP}/MapAligner_3_parameters.ini -in ${DATA_DIR_TOPP}/MapAligner_3_input1.mzML ${DATA_DIR_TOPP}/MapAligner_3_input2.mzML ${DATA_DIR_TOPP}/MapAligner_3_input3.mzML)
add_test("TOPP_MapAligner_3_out1" ${DIFF} -in1 MapAligner_3_output1.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_3_output1.mzML)
add_test("TOPP_MapAligner_3_out2" ${DIFF} -in1 MapAligner_3_output2.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_3_output2.mzML)
add_test("TOPP_MapAligner_3_out3" ${DIFF} -in1 MapAligner_3_output3.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_3_output3.mzML)
# "apply_given_trafo" algorithm:
add_test("TOPP_MapAligner_4" ${TOPP_BIN_PATH}/MapAligner -test -type apply_given_trafo -in ${DATA_DIR_TOPP}/MapAligner_4_input.featureXML ${DATA_DIR_TOPP}/MapAligner_4_input.featureXML -given_transformations ${DATA_DIR_TOPP}/MapAligner_4_trafo1.trafoXML ${DATA_DIR_TOPP}/MapAligner_4_trafo2.trafoXML -out MapAligner_4_output1.tmp MapAligner_4_output2.tmp)
add_test("TOPP_MapAligner_4_out1" ${DIFF} -in1 MapAligner_4_output1.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_4_output1.featureXML)
add_test("TOPP_MapAligner_4_out2" ${DIFF} -in1 MapAligner_4_output2.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_4_output2.featureXML)
add_test("TOPP_MapAligner_5" ${TOPP_BIN_PATH}/MapAligner -test -type apply_given_trafo -in ${DATA_DIR_TOPP}/MapAligner_5_input.mzML -given_transformations ${DATA_DIR_TOPP}/MapAligner_5_trafo.trafoXML -out MapAligner_5_output.tmp)
add_test("TOPP_MapAligner_5_out1" ${DIFF} -in1 MapAligner_5_output.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_5_output.mzML)
# "pose_clustering_affine" algorithm with "symmetric_regression" option:
add_test("TOPP_MapAligner_6" ${TOPP_BIN_PATH}/MapAligner -test -ini ${DATA_DIR_TOPP}/MapAligner_6_parameters.ini -in ${DATA_DIR_TOPP}/MapAligner_6_input1.featureXML ${DATA_DIR_TOPP}/MapAligner_6_input2.featureXML -out MapAligner_6_output1.tmp MapAligner_6_output2.tmp -transformations MapAligner_6_trafo1.tmp MapAligner_6_trafo2.tmp)
add_test("TOPP_MapAligner_6_out1" ${DIFF} -in1 MapAligner_6_trafo1.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_6_trafo1.trafoXML)
add_test("TOPP_MapAligner_6_out2" ${DIFF} -in1 MapAligner_6_trafo2.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_6_trafo2.trafoXML)
# "identification" algorithm:
add_test("TOPP_MapAligner_7" ${TOPP_BIN_PATH}/MapAligner -test -ini ${DATA_DIR_TOPP}/MapAligner_7_parameters.ini -in ${DATA_DIR_TOPP}/MapAligner_7_input1.featureXML ${DATA_DIR_TOPP}/ProteinQuantifier_input.featureXML -out MapAligner_7_output1.tmp MapAligner_7_output2.tmp)
add_test("TOPP_MapAligner_7_out1" ${DIFF} -in1 MapAligner_7_output1.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_7_output1.featureXML)
add_test("TOPP_MapAligner_7_out2" ${DIFF} -in1 MapAligner_7_output2.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_7_output2.featureXML)
# "identification" algorithm with "reference:file" option:
add_test("TOPP_MapAligner_8" ${TOPP_BIN_PATH}/MapAligner -test -ini ${DATA_DIR_TOPP}/MapAligner_7_parameters.ini -in ${DATA_DIR_TOPP}/MapAligner_7_input1.featureXML -out MapAligner_8_output1.tmp -reference:file ${DATA_DIR_TOPP}/ProteinQuantifier_input.featureXML)
add_test("TOPP_MapAligner_8_out1" ${DIFF} -in1 MapAligner_8_output1.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_8_output1.featureXML)
# "identification" algorithm with "reference:index" option:
add_test("TOPP_MapAligner_9" ${TOPP_BIN_PATH}/MapAligner -test -ini ${DATA_DIR_TOPP}/MapAligner_7_parameters.ini -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.featureXML ${DATA_DIR_TOPP}/MapAligner_7_input1.featureXML -out MapAligner_9_output1.tmp MapAligner_9_output2.tmp -reference:index 1)
add_test("TOPP_MapAligner_9_out1" ${DIFF} -in1 MapAligner_9_output1.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_9_output1.featureXML)
add_test("TOPP_MapAligner_9_out2" ${DIFF} -in1 MapAligner_9_output2.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_8_output1.featureXML)
add_test("TOPP_MapAligner_10" ${TOPP_BIN_PATH}/MapAligner -test -ini ${DATA_DIR_TOPP}/MapAligner_7_parameters.ini -in ${DATA_DIR_TOPP}/MapAligner_7_input1.featureXML ${DATA_DIR_TOPP}/ProteinQuantifier_input.featureXML -out MapAligner_10_output1.tmp MapAligner_10_output2.tmp -reference:index 2)
add_test("TOPP_MapAligner_10_out1" ${DIFF} -in1 MapAligner_10_output1.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_8_output1.featureXML)
add_test("TOPP_MapAligner_10_out2" ${DIFF} -in1 MapAligner_10_output2.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_9_output1.featureXML)
# "pose_clustering_affine" algorithm with "reference:file" option:
add_test("TOPP_MapAligner_11" ${TOPP_BIN_PATH}/MapAligner -test -ini ${DATA_DIR_TOPP}/MapAligner_1_parameters.ini -in ${DATA_DIR_TOPP}/MapAligner_1_input2.featureXML ${DATA_DIR_TOPP}/MapAligner_1_input3.featureXML -out MapAligner_11_output1.tmp MapAligner_11_output2.tmp -transformations -reference:file ${DATA_DIR_TOPP}/MapAligner_1_input1.featureXML)
add_test("TOPP_MapAligner_11_out1" ${DIFF} -in1 MapAligner_11_output1.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_1_output2.featureXML)
add_test("TOPP_MapAligner_11_out2" ${DIFF} -in1 MapAligner_11_output2.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_1_output3.featureXML)
# "identification" algorithm with consensusXML input:
add_test("TOPP_MapAligner_12" ${TOPP_BIN_PATH}/MapAligner -test -ini ${DATA_DIR_TOPP}/MapAligner_7_parameters.ini -in ${DATA_DIR_TOPP}/MapAligner_12_input1.consensusXML ${DATA_DIR_TOPP}/MapAligner_12_input2.consensusXML -out MapAligner_12_output1.tmp MapAligner_12_output2.tmp)
add_test("TOPP_MapAligner_12_out1" ${DIFF} -in1 MapAligner_12_output1.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_12_output1.consensusXML)
add_test("TOPP_MapAligner_12_out2" ${DIFF} -in1 MapAligner_12_output2.tmp -in2 ${DATA_DIR_TOPP}/MapAligner_12_output2.consensusXML)

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
# "high_res" algorithm with "ms1_only" option:
add_test("TOPP_PeakPicker_5" ${TOPP_BIN_PATH}/PeakPicker -test -ini ${DATA_DIR_TOPP}/PeakPicker_highres_parameters.ini -in ${DATA_DIR_TOPP}/PeakPicker_highres_input.mzML -out PeakPicker.tmp)
add_test("TOPP_PeakPicker_5_out1" ${DIFF} -in1 PeakPicker.tmp -in2 ${DATA_DIR_TOPP}/PeakPicker_highres_output.mzML)

### Resampler tests
add_test("TOPP_Resampler_1" ${TOPP_BIN_PATH}/Resampler -test -in ${DATA_DIR_TOPP}/Resampler_1_input.mzML -out Resampler.mzML -sampling_rate 0.3)
add_test("TOPP_Resampler_1_out1" ${DIFF} -in1 Resampler.mzML -in2 ${DATA_DIR_TOPP}/Resampler_1_output.mzML)

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
add_test("TOPP_RTPredict_1" ${TOPP_BIN_PATH}/RTPredict -in_id ${DATA_DIR_TOPP}/RTPredict_1_input.idXML -out_id:file RTPredict_1_output.tmp -total_gradient_time 3000 -svm_model ${DATA_DIR_TOPP}/RTPredict_1_input.model)
add_test("TOPP_RTPredict_1_out1" ${DIFF} -in1 RTPredict_1_output.tmp -in2 ${DATA_DIR_TOPP}/RTPredict_1_output.idXML)
add_test("TOPP_RTPredict_2" ${TOPP_BIN_PATH}/RTPredict -in_id ${DATA_DIR_TOPP}/RTPredict_2_input.idXML -out_id:positive RTPredict_2_output_positive.tmp -out_id:negative RTPredict_2_output_negative.tmp -svm_model ${DATA_DIR_TOPP}/RTPredict_2_input.model)
add_test("TOPP_RTPredict_2_out1" ${DIFF} -in1 RTPredict_2_output_positive.tmp -in2 ${DATA_DIR_TOPP}/RTPredict_2_output_positive.idXML)
add_test("TOPP_RTPredict_2_out2" ${DIFF} -in1 RTPredict_2_output_negative.tmp -in2 ${DATA_DIR_TOPP}/RTPredict_2_output_negative.idXML)
add_test("TOPP_RTPredict_3" ${TOPP_BIN_PATH}/RTPredict -in_id ${DATA_DIR_TOPP}/RTPredict_3_input.idXML -out_id:file RTPredict_3_output.tmp -total_gradient_time 1 -svm_model ${DATA_DIR_TOPP}/RTPredict_3_input.model)
add_test("TOPP_RTPredict_3_out1" ${DIFF} -whitelist "?xml-stylesheet" "predicted_RT_p_value" -in1 RTPredict_3_output.tmp -in2 ${DATA_DIR_TOPP}/RTPredict_3_output.idXML)
add_test("TOPP_RTPredict_4" ${TOPP_BIN_PATH}/RTPredict -in_text ${DATA_DIR_TOPP}/RTPredict_4_input.txt -out_text:file RTPredict_4_output.tmp -total_gradient_time 1 -svm_model ${DATA_DIR_TOPP}/RTPredict_4_input.model)
add_test("TOPP_RTPredict_4_out1" ${DIFF} -in1 RTPredict_4_output.tmp -in2 ${DATA_DIR_TOPP}/RTPredict_4_output.txt)
add_test("TOPP_RTPredict_5" ${TOPP_BIN_PATH}/RTPredict -in_id ${DATA_DIR_TOPP}/RTPredict_5_input.idXML -out_text:file RTPredict_5_output.tmp -total_gradient_time 1 -svm_model ${DATA_DIR_TOPP}/RTPredict_5_input.model)
add_test("TOPP_RTPredict_5_out1" ${DIFF} -in1 RTPredict_5_output.tmp -in2 ${DATA_DIR_TOPP}/RTPredict_5_output.txt)

### InspectAdapter tests
add_test("TOPP_InspectAdapter_1" ${TOPP_BIN_PATH}/InspectAdapter -ini ${DATA_DIR_TOPP}/InspectAdapter_1_parameters.ini -trie_dbs ${DATA_DIR_TOPP}/Inspect_FASTAFile_test2.trie -in ${DATA_DIR_TOPP}/Inspect.mzXML -dbs ${DATA_DIR_TOPP}/Inspect_FASTAFile_test.fasta -inspect_in -out InspectAdapter_2_output.tmp)

add_test("TOPP_InspectAdapter_1_out1" ${DIFF} -whitelist "spectra," "db," -in1 InspectAdapter_2_output.tmp -in2 ${DATA_DIR_TOPP}/InspectAdapter_2_output.inspect_in)

add_test("TOPP_InspectAdapter_2" ${TOPP_BIN_PATH}/InspectAdapter -ini ${DATA_DIR_TOPP}/InspectAdapter_1_parameters.ini -trie_dbs ${DATA_DIR_TOPP}/Inspect_FASTAFile_test2.trie -in ${DATA_DIR_TOPP}/Inspect.mzData -dbs ${DATA_DIR_TOPP}/Inspect_FASTAFile_test.fasta -inspect_in -out InspectAdapter_3_output.tmp)

add_test("TOPP_InspectAdapter_2_out1" ${DIFF} -whitelist "spectra," "db," -in1 InspectAdapter_3_output.tmp -in2 ${DATA_DIR_TOPP}/InspectAdapter_2_output.inspect_in)

# TODO the following tests are waiting for better implementations of InspectAdapter and associated classes
#add_test("TOPP_InspectAdapter_3" ${TOPP_BIN_PATH}/InspectAdapter -ini ${DATA_DIR_TOPP}/InspectAdapter_1_parameters.ini -trie_dbs ${DATA_DIR_TOPP}/Inspect_FASTAFile_test2.trie -in ${DATA_DIR_TOPP}/InspectAdapter.out -dbs ${DATA_DIR_TOPP}/Inspect_FASTAFile_test.fasta -out InspectAdapter_4_output.tmp -inspect_out)
#add_test("TOPP_InspectAdapter_3_out1" ${DIFF} -whitelist "?xml-stylesheet" "IdentificationRun date" -in1 InspectAdapter_4_output.tmp -in2 ${DATA_DIR_TOPP}/InspectAdapter_4_output.idXML)

### SeedListGenerator tests
add_test("TOPP_SeedListGenerator_1" ${TOPP_BIN_PATH}/SeedListGenerator -test -in ${PROJECT_SOURCE_DIR}/data/PepXMLFile_test.mzML -out SeedListGenerator_1_output.tmp)
add_test("TOPP_SeedListGenerator_1_out1" ${DIFF} -in1 SeedListGenerator_1_output.tmp -in2 ${DATA_DIR_TOPP}/SeedListGenerator_1_output.featureXML)
add_test("TOPP_SeedListGenerator_2" ${TOPP_BIN_PATH}/SeedListGenerator -test -in ${DATA_DIR_TOPP}/IDMapper_1_output.featureXML -out SeedListGenerator_2_output.tmp)
add_test("TOPP_SeedListGenerator_2_out1" ${DIFF} -in1 SeedListGenerator_2_output.tmp -in2 ${DATA_DIR_TOPP}/SeedListGenerator_2_output.featureXML)
add_test("TOPP_SeedListGenerator_3" ${TOPP_BIN_PATH}/SeedListGenerator -test -in ${PROJECT_SOURCE_DIR}/data/ConsensusXMLFile_1.consensusXML -out SeedListGenerator_3_output1.tmp SeedListGenerator_3_output2.tmp)
add_test("TOPP_SeedListGenerator_3_out1" ${DIFF} -in1 SeedListGenerator_3_output1.tmp -in2 ${DATA_DIR_TOPP}/SeedListGenerator_3_output1.featureXML)
add_test("TOPP_SeedListGenerator_3_out2" ${DIFF} -in1 SeedListGenerator_3_output2.tmp -in2 ${DATA_DIR_TOPP}/SeedListGenerator_3_output2.featureXML)

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
add_test("TOPP_SpectraFilter_1" ${TOPP_BIN_PATH}/SpectraFilter -test -in ${DATA_DIR_TOPP}/SpectraFilter_1_input.mzML -out SpectraFilter.tmp -type SqrtMower)
add_test("TOPP_SpectraFilter_1_out1" ${DIFF} -in1 SpectraFilter.tmp -in2 ${DATA_DIR_TOPP}/SpectraFilter_1_output.mzML)
add_test("TOPP_SpectraFilter_2" ${TOPP_BIN_PATH}/SpectraFilter -test -in ${DATA_DIR_TOPP}/SpectraFilter_2_input.mzML -out SpectraFilter.tmp -type WindowMower)
add_test("TOPP_SpectraFilter_2_out1" ${DIFF} -in1 SpectraFilter.tmp -in2 ${DATA_DIR_TOPP}/SpectraFilter_2_output.mzML)
add_test("TOPP_SpectraFilter_3" ${TOPP_BIN_PATH}/SpectraFilter -test -in ${DATA_DIR_TOPP}/SpectraFilter_3_input.mzML -out SpectraFilter.tmp -ini ${DATA_DIR_TOPP}/SpectraFilter_3_parameters.ini -type WindowMower)
add_test("TOPP_SpectraFilter_3_out1" ${DIFF} -in1 SpectraFilter.tmp -in2 ${DATA_DIR_TOPP}/SpectraFilter_3_output.mzML)

### InternalCalibration tests
add_test("TOPP_InternalCalibration_1" ${TOPP_BIN_PATH}/InternalCalibration -test -in ${DATA_DIR_TOPP}/InternalCalibration_2_input.mzML -out InternalCalibration.tmp -ref_peaks ${DATA_DIR_TOPP}/InternalCalibration_ref_masses -type spectrumwise)
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


### SILACAnalyzer tests
add_test("TOPP_SILACAnalyzer_1" ${TOPP_BIN_PATH}/SILACAnalyzer -test -ini ${DATA_DIR_TOPP}/SILACAnalyzer_double.ini -in ${DATA_DIR_TOPP}/SILACAnalyzer_double.mzML -out SILACAnalyzer.tmp -out_visual SILACAnalyzer_vis.tmp)
add_test("TOPP_SILACAnalyzer_1_out1" ${DIFF} -in1 SILACAnalyzer.tmp -in2 ${DATA_DIR_TOPP}/SILACAnalyzer_double.consensusXML)
add_test("TOPP_SILACAnalyzer_1_out2" ${DIFF} -in1 SILACAnalyzer_vis.tmp -in2 ${DATA_DIR_TOPP}/SILACAnalyzer_double.featureXML)

add_test("TOPP_SILACAnalyzer_2" ${TOPP_BIN_PATH}/SILACAnalyzer -test -ini ${DATA_DIR_TOPP}/SILACAnalyzer_triple.ini -in ${DATA_DIR_TOPP}/SILACAnalyzer_triple.mzML -out SILACAnalyzer.tmp -out_visual SILACAnalyzer_vis.tmp)
add_test("TOPP_SILACAnalyzer_2_out1" ${DIFF} -in1 SILACAnalyzer.tmp -in2 ${DATA_DIR_TOPP}/SILACAnalyzer_triple.consensusXML)
add_test("TOPP_SILACAnalyzer_2_out2" ${DIFF} -in1 SILACAnalyzer_vis.tmp -in2 ${DATA_DIR_TOPP}/SILACAnalyzer_triple.featureXML)

### FeatureLinker tests
# "labeled" algorithm:
add_test("TOPP_FeatureLinker_1" ${TOPP_BIN_PATH}/FeatureLinker -test -ini ${DATA_DIR_TOPP}/FeatureLinker_1_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinker_1_input.featureXML)
add_test("TOPP_FeatureLinker_1_out1" ${DIFF} -in1 FeatureLinker_1_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinker_1_output.consensusXML)
# "unlabeled" algorithm:
add_test("TOPP_FeatureLinker_2" ${TOPP_BIN_PATH}/FeatureLinker -test -ini ${DATA_DIR_TOPP}/FeatureLinker_2_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinker_2_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinker_2_input2.featureXML ${DATA_DIR_TOPP}/FeatureLinker_2_input3.featureXML)
add_test("TOPP_FeatureLinker_2_out1" ${DIFF} -in1 FeatureLinker_2_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinker_2_output.consensusXML)
add_test("TOPP_FeatureLinker_3" ${TOPP_BIN_PATH}/FeatureLinker -test -ini ${DATA_DIR_TOPP}/FeatureLinker_3_parameters.ini  -in ${DATA_DIR_TOPP}/FeatureLinker_3_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinker_3_input2.featureXML ${DATA_DIR_TOPP}/FeatureLinker_3_input3.featureXML)
add_test("TOPP_FeatureLinker_3_out1" ${DIFF} -in1 FeatureLinker_3_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinker_3_output.consensusXML)
add_test("TOPP_FeatureLinker_4" ${TOPP_BIN_PATH}/FeatureLinker -test -ini ${DATA_DIR_TOPP}/FeatureLinker_4_parameters.ini  -in ${DATA_DIR_TOPP}/FeatureLinker_4_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinker_4_input2.featureXML)
add_test("TOPP_FeatureLinker_4_out1" ${DIFF} -in1 FeatureLinker_4_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinker_4_output.consensusXML)
# "labeled" algorithm:
add_test("TOPP_FeatureLinker_5" ${TOPP_BIN_PATH}/FeatureLinker -test -ini ${DATA_DIR_TOPP}/FeatureLinker_5_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinker_5_input.featureXML)
add_test("TOPP_FeatureLinker_5_out1" ${DIFF} -in1 FeatureLinker_5_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinker_5_output.consensusXML)
# "unlabeled" algorithm with "use_identifications" option:
add_test("TOPP_FeatureLinker_6" ${TOPP_BIN_PATH}/FeatureLinker -test -ini ${DATA_DIR_TOPP}/FeatureLinker_6_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinker_2_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinker_2_input2.featureXML ${DATA_DIR_TOPP}/FeatureLinker_2_input3.featureXML -out FeatureLinker_6_output.tmp)
add_test("TOPP_FeatureLinker_6_out1" ${DIFF} -in1 FeatureLinker_6_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinker_6_output.consensusXML)
# "unlabeled_qt" algorithm:
add_test("TOPP_FeatureLinker_7" ${TOPP_BIN_PATH}/FeatureLinker -test -ini ${DATA_DIR_TOPP}/FeatureLinker_7_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinker_2_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinker_2_input2.featureXML ${DATA_DIR_TOPP}/FeatureLinker_2_input3.featureXML -out FeatureLinker_7_output.tmp)
add_test("TOPP_FeatureLinker_7_out1" ${DIFF} -in1 FeatureLinker_7_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinker_7_output.consensusXML)
add_test("TOPP_FeatureLinker_8" ${TOPP_BIN_PATH}/FeatureLinker -test -ini ${DATA_DIR_TOPP}/FeatureLinker_8_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinker_7_output.consensusXML ${DATA_DIR_TOPP}/FeatureLinker_7_output.consensusXML -out FeatureLinker_8_output.tmp)
add_test("TOPP_FeatureLinker_8_out1" ${DIFF} -in1 FeatureLinker_8_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinker_8_output.consensusXML)

### IDMapper tests
add_test("TOPP_IDMapper_1" ${TOPP_BIN_PATH}/IDMapper -test -in ${DATA_DIR_TOPP}/IDMapper_1_input.featureXML -out IDMapper_1_output.tmp -id ${DATA_DIR_TOPP}/IDMapper_1_input.idXML -mz_measure Da -mz_tolerance 1)
add_test("TOPP_IDMapper_1_out1" ${DIFF} -in1 IDMapper_1_output.tmp -in2 ${DATA_DIR_TOPP}/IDMapper_1_output.featureXML)
add_test("TOPP_IDMapper_2" ${TOPP_BIN_PATH}/IDMapper -test -in ${DATA_DIR_TOPP}/IDMapper_2_input.consensusXML -out IDMapper_2_output.tmp -id ${DATA_DIR_TOPP}/IDMapper_2_input.idXML  -mz_measure Da -mz_tolerance 1)
add_test("TOPP_IDMapper_2_out1" ${DIFF} -in1 IDMapper_2_output.tmp -in2 ${DATA_DIR_TOPP}/IDMapper_2_output.consensusXML)
add_test("TOPP_IDMapper_3" ${TOPP_BIN_PATH}/IDMapper -test -in ${DATA_DIR_TOPP}/IDMapper_3_input.featureXML -out IDMapper_3_output.tmp -id ${DATA_DIR_TOPP}/IDMapper_3_input.idXML  -mz_measure ppm -mz_tolerance 3 -rt_tolerance 4)
add_test("TOPP_IDMapper_3_out1" ${DIFF} -in1 IDMapper_3_output.tmp -in2 ${DATA_DIR_TOPP}/IDMapper_3_output.featureXML)

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

### PrecursorMassCorrector tests
add_test("TOPP_PrecursorMassCorrector_1" ${TOPP_BIN_PATH}/PrecursorMassCorrector -in ${DATA_DIR_TOPP}/PrecursorMassCorrector_1_input.mzML -ini ${DATA_DIR_TOPP}/PrecursorMassCorrector_1.ini -out PrecursorMassCorrector_1_output.tmp)
add_test("TOPP_PrecursorMassCorrector_1_out1" ${DIFF} -in1 PrecursorMassCorrector_1_output.tmp -in2 ${DATA_DIR_TOPP}/PrecursorMassCorrector_1_output.mzML)

### ProteinInference tests
add_test("TOPP_ProteinInference_1" ${TOPP_BIN_PATH}/ProteinInference -test -in ${DATA_DIR_TOPP}/ProteinInference.idXML -out ProteinInference_output.tmp)
add_test("TOPP_ProteinInference_1_out1" ${DIFF} -in1 ProteinInference_output.tmp -in2 ${DATA_DIR_TOPP}/ProteinInference_proteins.idXML)

### FalseDiscoveryRate_tests (1: default param, 2: with runs separated, 3: with charge variants individual, 4: with charge variants and no charge 1 decoys -> q_value 0)
add_test("TOPP_FalseDiscoveryRate_1" ${TOPP_BIN_PATH}/FalseDiscoveryRate -test  -ini ${DATA_DIR_TOPP}/FalseDiscoveryRate_1.ini -in ${DATA_DIR_TOPP}/FalseDiscoveryRate_OMSSA.idXML -out FalseDiscoveryRate_output.tmp)
add_test("TOPP_FalseDiscoveryRate_1_out1" ${DIFF} -whitelist "?xml-stylesheet" -in1 FalseDiscoveryRate_output.tmp -in2 ${DATA_DIR_TOPP}/FalseDiscoveryRate_output_1.idXML)
add_test("TOPP_FalseDiscoveryRate_2" ${TOPP_BIN_PATH}/FalseDiscoveryRate -test  -ini ${DATA_DIR_TOPP}/FalseDiscoveryRate_2.ini -in ${DATA_DIR_TOPP}/FalseDiscoveryRate_OMSSA.idXML -out FalseDiscoveryRate_output.tmp)
add_test("TOPP_FalseDiscoveryRate_2_out1" ${DIFF} -whitelist "?xml-stylesheet" -in1 FalseDiscoveryRate_output.tmp -in2 ${DATA_DIR_TOPP}/FalseDiscoveryRate_output_2.idXML)
add_test("TOPP_FalseDiscoveryRate_3" ${TOPP_BIN_PATH}/FalseDiscoveryRate -test  -ini ${DATA_DIR_TOPP}/FalseDiscoveryRate_3.ini -in ${DATA_DIR_TOPP}/FalseDiscoveryRate_OMSSA.idXML -out FalseDiscoveryRate_output.tmp)
add_test("TOPP_FalseDiscoveryRate_3_out1" ${DIFF} -whitelist "?xml-stylesheet" -in1 FalseDiscoveryRate_output.tmp -in2 ${DATA_DIR_TOPP}/FalseDiscoveryRate_output_3.idXML)
add_test("TOPP_FalseDiscoveryRate_4" ${TOPP_BIN_PATH}/FalseDiscoveryRate -test  -ini ${DATA_DIR_TOPP}/FalseDiscoveryRate_3.ini -in ${DATA_DIR_TOPP}/FalseDiscoveryRate_OMSSA_4.idXML -out FalseDiscoveryRate_output.tmp)
add_test("TOPP_FalseDiscoveryRate_4_out1" ${DIFF} -whitelist "?xml-stylesheet" -in1 FalseDiscoveryRate_output.tmp -in2 ${DATA_DIR_TOPP}/FalseDiscoveryRate_output_4.idXML)

###IDPosteriorErrorProbability
add_test("TOPP_IDPosteriorErrorProbability_1" ${TOPP_BIN_PATH}/IDPosteriorErrorProbability -test -in ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_Mascot_input.idXML -out 	IDPosteriorErrorProbability_output.tmp)
add_test("TOPP_IDPosteriorErrorProbability_1_out1" ${DIFF} -in1 IDPosteriorErrorProbability_output.tmp -in2 ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_Mascot_output.idXML)
add_test("TOPP_IDPosteriorErrorProbability_2" ${TOPP_BIN_PATH}/IDPosteriorErrorProbability -test -in ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_XTandem_input.idXML -out IDPosteriorErrorProbability_output.tmp)
add_test("TOPP_IDPosteriorErrorProbability_2_out1" ${DIFF} -in1 IDPosteriorErrorProbability_output.tmp -in2 ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_XTandem_output.idXML)
add_test("TOPP_IDPosteriorErrorProbability_3" ${TOPP_BIN_PATH}/IDPosteriorErrorProbability -test -in ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_OMSSA_input.idXML -out IDPosteriorErrorProbability_output.tmp)
add_test("TOPP_IDPosteriorErrorProbability_3_out1" ${DIFF} -in1 IDPosteriorErrorProbability_output.tmp -in2 ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_OMSSA_output.idXML)
add_test("TOPP_IDPosteriorErrorProbability_4" ${TOPP_BIN_PATH}/IDPosteriorErrorProbability -test -in ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_OMSSA_input2.idXML -out IDPosteriorErrorProbability_output.tmp -split_charge)
add_test("TOPP_IDPosteriorErrorProbability_4_out1" ${DIFF} -in1 IDPosteriorErrorProbability_output.tmp -in2 ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_OMSSA_output2.idXML)
add_test("TOPP_IDPosteriorErrorProbability_5" ${TOPP_BIN_PATH}/IDPosteriorErrorProbability -test -in ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_XTandem_input2.idXML -out IDPosteriorErrorProbability_output.tmp -split_charge)
add_test("TOPP_IDPosteriorErrorProbability_5_out1" ${DIFF} -in1 IDPosteriorErrorProbability_output.tmp -in2 ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_XTandem_output2.idXML)
add_test("TOPP_IDPosteriorErrorProbability_6" ${TOPP_BIN_PATH}/IDPosteriorErrorProbability -test -in ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_Mascot_input2.idXML -out IDPosteriorErrorProbability_output.tmp -split_charge)
add_test("TOPP_IDPosteriorErrorProbability_6_out1" ${DIFF} -in1 IDPosteriorErrorProbability_output.tmp -in2 ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_Mascot_output2.idXML)

### ProteinQuantifier tests
add_test("TOPP_ProteinQuantifier_1" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.featureXML -out ProteinQuantifier_1_output1.tmp -peptide_out ProteinQuantifier_1_output2.tmp)
add_test("TOPP_ProteinQuantifier_1_out1" ${DIFF} -in1 ProteinQuantifier_1_output1.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_1_output1.txt -whitelist "computed from file")
add_test("TOPP_ProteinQuantifier_1_out2" ${DIFF} -in1 ProteinQuantifier_1_output2.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_1_output2.txt -whitelist "computed from file")
add_test("TOPP_ProteinQuantifier_2" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.featureXML -out ProteinQuantifier_2_output1.tmp -peptide_out ProteinQuantifier_2_output2.tmp -top 2 -include_all -filter_charge -average sum)
add_test("TOPP_ProteinQuantifier_2_out1" ${DIFF} -in1 ProteinQuantifier_2_output1.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_2_output1.txt -whitelist "computed from file")
add_test("TOPP_ProteinQuantifier_2_out2" ${DIFF} -in1 ProteinQuantifier_2_output2.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_2_output2.txt -whitelist "computed from file")
# handling of modified peptides:
add_test("TOPP_ProteinQuantifier_3" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_3_input.featureXML -out ProteinQuantifier_3_output1.tmp -peptide_out ProteinQuantifier_3_output2.tmp -top 2 -include_all -average mean)
add_test("TOPP_ProteinQuantifier_3_out1" ${DIFF} -in1 ProteinQuantifier_3_output1.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_3_output1.txt -whitelist "computed from file")
add_test("TOPP_ProteinQuantifier_3_out2" ${DIFF} -in1 ProteinQuantifier_3_output2.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_3_output2.txt -whitelist "computed from file")
# consensusXML input, different options for protein quantification:
add_test("TOPP_ProteinQuantifier_4" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.consensusXML -out ProteinQuantifier_4_output.tmp -top 0 -average sum)
add_test("TOPP_ProteinQuantifier_4_out1" ${DIFF} -in1 ProteinQuantifier_4_output.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_4_output.txt -whitelist "computed from file")
add_test("TOPP_ProteinQuantifier_5" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.consensusXML -out ProteinQuantifier_5_output.tmp -top 3 -average sum)
add_test("TOPP_ProteinQuantifier_5_out1" ${DIFF} -in1 ProteinQuantifier_5_output.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_5_output.txt -whitelist "computed from file")
add_test("TOPP_ProteinQuantifier_6" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.consensusXML -out ProteinQuantifier_6_output.tmp -top 3 -include_all -average sum)
add_test("TOPP_ProteinQuantifier_6_out1" ${DIFF} -in1 ProteinQuantifier_6_output.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_6_output.txt -whitelist "computed from file")
add_test("TOPP_ProteinQuantifier_7" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.consensusXML -out ProteinQuantifier_7_output.tmp -top 0 -consensus:fix_peptides -average sum)
add_test("TOPP_ProteinQuantifier_7_out1" ${DIFF} -in1 ProteinQuantifier_7_output.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_7_output.txt -whitelist "computed from file")
add_test("TOPP_ProteinQuantifier_8" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.consensusXML -out ProteinQuantifier_8_output.tmp -top 3 -consensus:fix_peptides -average sum)
add_test("TOPP_ProteinQuantifier_8_out1" ${DIFF} -in1 ProteinQuantifier_8_output.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_8_output.txt -whitelist "computed from file")
add_test("TOPP_ProteinQuantifier_9" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.consensusXML -out ProteinQuantifier_9_output.tmp -top 3 -include_all -consensus:fix_peptides -average sum)
add_test("TOPP_ProteinQuantifier_9_out1" ${DIFF} -in1 ProteinQuantifier_9_output.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_9_output.txt -whitelist "computed from file")

### InclusionExclusionListCreator tests
add_test("TOPP_InclusionExclusionListCreator_1" ${TOPP_BIN_PATH}/InclusionExclusionListCreator -test -include ${DATA_DIR_TOPP}/InclusionExclusionListCreator.featureXML -out InclusionExclusionListCreator_1_output.tmp)
add_test("TOPP_InclusionExclusionListCreator_1_out1" ${DIFF} -in1 InclusionExclusionListCreator_1_output.tmp -in2 ${DATA_DIR_TOPP}/InclusionExclusionListCreator_1_output.txt )
add_test("TOPP_InclusionExclusionListCreator_2" ${TOPP_BIN_PATH}/InclusionExclusionListCreator -test -exclude ${DATA_DIR_TOPP}/InclusionExclusionListCreator.featureXML -out InclusionExclusionListCreator_2_output.tmp -rel_rt_window_size 0.5)
add_test("TOPP_InclusionExclusionListCreator_2_out1" ${DIFF} -in1 InclusionExclusionListCreator_2_output.tmp -in2 ${DATA_DIR_TOPP}/InclusionExclusionListCreator_2_output.txt)
add_test("TOPP_InclusionExclusionListCreator_3" ${TOPP_BIN_PATH}/InclusionExclusionListCreator -test -include ${DATA_DIR_TOPP}/InclusionExclusionListCreator.fasta -out InclusionExclusionListCreator_3_output.tmp -inclusion_charges 1 2 -rt_model ${DATA_DIR_TOPP}/InclusionExclusionListCreator_rt.model)
add_test("TOPP_InclusionExclusionListCreator_3_out1" ${DIFF} -in1 InclusionExclusionListCreator_3_output.tmp -in2 ${DATA_DIR_TOPP}/InclusionExclusionListCreator_3_output.txt)
add_test("TOPP_InclusionExclusionListCreator_4" ${TOPP_BIN_PATH}/InclusionExclusionListCreator -test -exclude ${DATA_DIR_TOPP}/InclusionExclusionListCreator.fasta -out InclusionExclusionListCreator_4_output.tmp -exclusion_charges 1 2 -rt_model ${DATA_DIR_TOPP}/InclusionExclusionListCreator_rt.model -missed_cleavages 1 -rt_in_seconds)
add_test("TOPP_InclusionExclusionListCreator_4_out1" ${DIFF} -in1 InclusionExclusionListCreator_4_output.tmp -in2 ${DATA_DIR_TOPP}/InclusionExclusionListCreator_4_output.txt)
add_test("TOPP_InclusionExclusionListCreator_5" ${TOPP_BIN_PATH}/InclusionExclusionListCreator -test -exclude ${DATA_DIR_TOPP}/InclusionExclusionListCreator.IdXML -out InclusionExclusionListCreator_5_output.tmp -exclusion_charges 1 2 -rt_model ${DATA_DIR_TOPP}/InclusionExclusionListCreator_rt.model)
add_test("TOPP_InclusionExclusionListCreator_5_out1" ${DIFF} -in1 InclusionExclusionListCreator_5_output.tmp -in2 ${DATA_DIR_TOPP}/InclusionExclusionListCreator_5_output.txt)



### UTILS tests (if this becomes too much, we should create a separate directory)
add_test("UTILS_ImageCreator_1" ${TOPP_BIN_PATH}/ImageCreator -test -in ${DATA_DIR_TOPP}/ImageCreator_1_input.mzML -out ImageCreator1.png -width 20 -height 15)
add_test("UTILS_ImageCreator_1_out1" ${DIFF} -in1 ImageCreator1.png -in2 ${DATA_DIR_TOPP}/ImageCreator_1_output.png)
add_test("UTILS_ImageCreator_2" ${TOPP_BIN_PATH}/ImageCreator -test -in ${DATA_DIR_TOPP}/ImageCreator_2_input.mzML -out ImageCreator2.png -width 20 -height 15 -precursors -precursor_size 3 -precursor_color green -log_intensity)
add_test("UTILS_ImageCreator_2_out1" ${DIFF} -in1 ImageCreator2.png -in2 ${DATA_DIR_TOPP}/ImageCreator_2_output.png)

add_test("UTILS_IDSplitter_1" ${TOPP_BIN_PATH}/IDSplitter -test -in ${DATA_DIR_TOPP}/IDMapper_1_output.featureXML -out IDSplitter_1_output1.tmp -id_out IDSplitter_1_output2.tmp)
add_test("UTILS_IDSplitter_1_out1" ${DIFF} -in1 IDSplitter_1_output1.tmp -in2 ${DATA_DIR_TOPP}/IDSplitter_1_output.featureXML)
add_test("UTILS_IDSplitter_1_out2" ${DIFF} -in1 IDSplitter_1_output2.tmp -in2 ${DATA_DIR_TOPP}/IDSplitter_1_output.idXML)

add_test("UTILS_MassCalculator_1" ${TOPP_BIN_PATH}/MassCalculator -test -in ${DATA_DIR_TOPP}/MassCalculator_1_input.txt -out MassCalculator_1_output.tmp -charge 0 1 -format table -separator ,)
add_test("UTILS_MassCalculator_1_out1" ${DIFF} -in1 MassCalculator_1_output.tmp -in2 ${DATA_DIR_TOPP}/MassCalculator_1_output.csv)
add_test("UTILS_MassCalculator_2" ${TOPP_BIN_PATH}/MassCalculator -test -in_seq "LDQWLC(Carbamidomethyl)EKL" "(Glu->pyro-Glu)EAM(Oxidation)APKHK" "RANVM(Oxidation)DYR" "FGVEQDVDMVFASFIR" -out MassCalculator_2_output.tmp -charge 1 2 3)
add_test("UTILS_MassCalculator_2_out1" ${DIFF} -in1 MassCalculator_2_output.tmp -in2 ${DATA_DIR_TOPP}/MassCalculator_2_output.txt)
