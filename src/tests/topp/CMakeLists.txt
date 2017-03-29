# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2016.
#
# This software is released under a three-clause BSD license:
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#  * Neither the name of any author or any participating institution
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.
# For a full list of authors, refer to the file AUTHORS.
# --------------------------------------------------------------------------
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
# INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche, Chris Bielow $
# $Authors: Andreas Bertsch, Chris Bielow, Stephan Aiche $
# --------------------------------------------------------------------------

project("OpenMS_TOPP_tests")
cmake_minimum_required(VERSION 2.8.3 FATAL_ERROR)

#
# README:
#
# This file contains the declarations of all TOPP/UTILS tests that are executed by CTest
#
# In order to add a test do the following:
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
# <num>   Index that is increased for each call of a TOPP tool during a test.
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
#       The -test flag ensures that this output is independent of local path settings etc., so use it when *testing* TOPP tools!
# Hint: Provide individual names for each test output file e.g., the first test for the FileConverter should have FileConverter_1.tmp
#       as output file.

#------------------------------------------------------------------------------
# General settings
#------------------------------------------------------------------------------
set(TOPP_BIN_PATH ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
set(DATA_DIR_TOPP ${PROJECT_SOURCE_DIR})
set(DATA_DIR_TOPP_BIN ${PROJECT_BINARY_DIR})
set(INDEX_WHITELIST "offset" "indexListOffset")
set(DIFF ${TOPP_BIN_PATH}/FuzzyDiff -test -ini ${DATA_DIR_TOPP}/FuzzyDiff.ini)

## configure resource files (for TOPPAS/ExecutePipeline) - there is no other way to specify the input
configure_file(${DATA_DIR_TOPP}/ExecutePipeline_1.trf.in ${DATA_DIR_TOPP_BIN}/ExecutePipeline_1.trf)

#------------------------------------------------------------------------------
# test if -write_wdsl works
#------------------------------------------------------------------------------
add_test("TOPP_WSDL_1" ${TOPP_BIN_PATH}/AdditiveSeries -test -write_wsdl WSDL_out.tmp)
add_test("TOPP_WSDL_1_out1" ${DIFF} -in1 WSDL_out.tmp -in2 ${DATA_DIR_TOPP}/WSDL_out.wsdl)
set_tests_properties("TOPP_WSDL_1_out1" PROPERTIES DEPENDS "TOPP_WSDL_1")

#------------------------------------------------------------------------------
# WRITE_INI test for each TOPP/UTILS tool, checking -write_ini functionality
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# TOPP
set(TOPP_executables_INI ${TOPP_TOOLS})
list(REMOVE_ITEM TOPP_executables_INI GenericWrapper) # GenericWrapper (has a -type flag)
add_test("TOPPWRITEINI_GenericWrapper" ${TOPP_BIN_PATH}/GenericWrapper -test -type RAWFileConvert -write_ini "GenericWrapper.ini.tmp")
foreach (i ${TOPP_executables_INI})
  add_test("TOPPWRITEINI_${i}" "${TOPP_BIN_PATH}/${i}" -test -write_ini "${i}.ini.tmp")
  add_test("TOPPWRITEINI_${i}_SectionName" ${CMAKE_COMMAND} -test -DNAME=${i} -DFILE=${i}.ini.tmp -P "${DATA_DIR_TOPP}/check_ini.cmake")
  set_tests_properties("TOPPWRITEINI_${i}_SectionName" PROPERTIES DEPENDS "TOPPWRITEINI_${i}")
endforeach(i)

#------------------------------------------------------------------------------
# UTILS
set(UTILS_executables_INI ${UTILS_TOOLS})
list(REMOVE_ITEM UTILS_executables_INI OpenMSInfo)

if(NOT HAS_XSERVER) # remove IDEvaluator if we don't have a running Xserver
  list(REMOVE_ITEM UTILS_executables_INI IDEvaluator)
endif()

foreach (i ${UTILS_executables_INI})
  add_test("UTILSWRITEINI_${i}" "${TOPP_BIN_PATH}/${i}" -test -write_ini "${i}.ini.tmp")
  add_test("UTILSWRITEINI_${i}_SectionName" ${CMAKE_COMMAND} -test -DNAME=${i} -DFILE=${i}.ini.tmp -P "${DATA_DIR_TOPP}/check_ini.cmake")
  set_tests_properties("UTILSWRITEINI_${i}_SectionName" PROPERTIES DEPENDS "UTILSWRITEINI_${i}")
endforeach(i)

#------------------------------------------------------------------------------
# WRITECTD test to check whether -write_ctd works
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# UTILS
set(UTILS_executables_CTD ${UTILS_TOOLS})
list(REMOVE_ITEM UTILS_executables_CTD OpenMSInfo)
if(NOT HAS_XSERVER) # remove IDEvaluator if we don't have a running Xserver
	list(REMOVE_ITEM UTILS_executables_CTD IDEvaluator)
endif()
foreach (i ${UTILS_executables_CTD})
  add_test("UTILSWRITECTD_${i}" "${TOPP_BIN_PATH}/${i}" -test -write_ctd)
endforeach()

#------------------------------------------------------------------------------
# TOPP
set(TOPP_executables_CTD ${TOPP_TOOLS})
list(REMOVE_ITEM TOPP_executables_CTD GenericWrapper)
foreach (i ${TOPP_executables_CTD})
  add_test("TOPPWRITECTD_${i}" "${TOPP_BIN_PATH}/${i}" -test -write_ctd)
endforeach()

#------------------------------------------------------------------------------
# test if write_ini does the correct updates when given an old INI file (using
# one arbitrary TOPP tool)
#------------------------------------------------------------------------------
add_test("TOPPWRITEINI_OVERWRITE" ${TOPP_BIN_PATH}/PeakPickerWavelet -test -write_ini "WRITE_INI.ini.tmp" -ini "${DATA_DIR_TOPP}/WRITE_INI_IN.ini")
add_test("TOPPWRITEINI_OVERWRITE_out" ${DIFF} -in1 WRITE_INI.ini.tmp -in2 "${DATA_DIR_TOPP}/WRITE_INI_OUT.ini")
set_tests_properties("TOPPWRITEINI_OVERWRITE_out" PROPERTIES DEPENDS "TOPPWRITEINI_OVERWRITE")
##
## Test if invalid arguments or values (from INI or CLI) lead to immediate exit
##
## PARAM-VALUES:
add_test("TOPP_INI_INVALIDVALUE" ${TOPP_BIN_PATH}/PeakPickerHiRes -ini ${DATA_DIR_TOPP}/degenerate_cases/PPHiRes_invalidValue.ini )
# > Invalid string parameter value 'THISVALUEISINVALID' for parameter 'processOption' given! Valid values are: 'inmemory,lowmemory'. Updating failed!
# > Parameters passed to 'PeakPickerHiRes' are invalid. To prevent usage of wrong defaults, please update/fix the parameters!
set_tests_properties("TOPP_INI_INVALIDVALUE" PROPERTIES WILL_FAIL 1)
## ... the same as explicit command line option (instead of INI)
add_test("TOPP_CLI_INVALIDVALUE" ${TOPP_BIN_PATH}/PeakPickerHiRes -processOption THISVALUEISINVALID )
set_tests_properties("TOPP_CLI_INVALIDVALUE" PROPERTIES WILL_FAIL 1)
 
## test algorithm section parameters (not directly registered in TOPPBase); they should fail just the same
add_test("TOPP_INI_INVALIDVALUE_SECTION" ${TOPP_BIN_PATH}/PeakPickerHiRes -ini ${DATA_DIR_TOPP}/degenerate_cases/PPHiRes_invalidValueSection.ini )
set_tests_properties("TOPP_INI_INVALIDVALUE_SECTION" PROPERTIES WILL_FAIL 1)
## ... the same as explicit command line option (instead of INI)
add_test("TOPP_CLI_INVALIDVALUE_SECTION" ${TOPP_BIN_PATH}/PeakPickerHiRes -algorithm:report_FWHM_unit THISVALUEISINVALID )
set_tests_properties("TOPP_CLI_INVALIDVALUE_SECTION" PROPERTIES WILL_FAIL 1)

## PARAM-NAMES:
add_test("TOPP_INI_INVALIDNAME" ${TOPP_BIN_PATH}/PeakPickerHiRes -ini ${DATA_DIR_TOPP}/degenerate_cases/PPHiRes_invalidParamName.ini )
# > Unknown (or deprecated) Parameter 'INVALIDPARAMNAME_BASE' given in old parameter file! Update failed!
# > Unknown (or deprecated) Parameter 'algorithm:invalidParamName' given in old parameter file! Update failed!
# > Parameters passed to 'PeakPickerHiRes' are invalid. To prevent usage of wrong defaults, please update/fix the parameters!
set_tests_properties("TOPP_INI_INVALIDNAME" PROPERTIES WILL_FAIL 1)
## ... the same as explicit command line option (instead of INI)
add_test("TOPP_CLI_INVALIDNAME" ${TOPP_BIN_PATH}/PeakPickerHiRes -algorithm:invalidParamName somevalue )
set_tests_properties("TOPP_CLI_INVALIDNAME" PROPERTIES WILL_FAIL 1) 
 
#------------------------------------------------------------------------------
# TOPP tests
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# AdditiveSeries tests
add_test("TOPP_AdditiveSeries_1" ${TOPP_BIN_PATH}/AdditiveSeries -test -in ${DATA_DIR_TOPP}/AdditiveSeries_1_feat11.featureXML ${DATA_DIR_TOPP}/AdditiveSeries_1_feat20.featureXML ${DATA_DIR_TOPP}/AdditiveSeries_1_feat22.featureXML ${DATA_DIR_TOPP}/AdditiveSeries_1_feat34.featureXML -ini ${DATA_DIR_TOPP}/AdditiveSeries_1_parameters.ini)
add_test("TOPP_AdditiveSeries_1_out1" ${DIFF} -in1 AdditiveSeries_1_output.tmp -in2 ${DATA_DIR_TOPP}/AdditiveSeries_1_output.txt )
set_tests_properties("TOPP_AdditiveSeries_1_out1" PROPERTIES DEPENDS "TOPP_AdditiveSeries_1")
add_test("TOPP_AdditiveSeries_1_out2" ${DIFF} -in1 AdditiveSeries_1_gnuplot_tmp.cmd -in2 ${DATA_DIR_TOPP}/AdditiveSeries_1_gnuplot.cmd )
set_tests_properties("TOPP_AdditiveSeries_1_out2" PROPERTIES DEPENDS "TOPP_AdditiveSeries_1")
add_test("TOPP_AdditiveSeries_1_out3" ${DIFF} -in1 AdditiveSeries_1_gnuplot_tmp.dat -in2 ${DATA_DIR_TOPP}/AdditiveSeries_1_gnuplot.dat )
set_tests_properties("TOPP_AdditiveSeries_1_out3" PROPERTIES DEPENDS "TOPP_AdditiveSeries_1")
add_test("TOPP_AdditiveSeries_1_out4" ${DIFF} -in1 AdditiveSeries_1_gnuplot_tmp.err -in2 ${DATA_DIR_TOPP}/AdditiveSeries_1_gnuplot.err )
set_tests_properties("TOPP_AdditiveSeries_1_out4" PROPERTIES DEPENDS "TOPP_AdditiveSeries_1")

#------------------------------------------------------------------------------
# IDMerger tests
add_test("TOPP_IDMerger_1" ${TOPP_BIN_PATH}/IDMerger -test -in ${DATA_DIR_TOPP}/IDMerger_1_input1.idXML ${DATA_DIR_TOPP}/IDMerger_1_input2.idXML -out IDMerger_1_output.tmp -annotate_file_origin)
add_test("TOPP_IDMerger_1_out1" ${DIFF} -in1 IDMerger_1_output.tmp -in2 ${DATA_DIR_TOPP}/IDMerger_1_output.idXML)
set_tests_properties("TOPP_IDMerger_1_out1" PROPERTIES DEPENDS "TOPP_IDMerger_1")
# merge pepXML and protXML:
add_test("TOPP_IDMerger_2" ${TOPP_BIN_PATH}/IDMerger -test -pepxml_protxml -in ${DATA_DIR_TOPP}/IDMerger_2_input1.idXML ${DATA_DIR_TOPP}/IDMerger_2_input2.idXML -out IDMerger_2_output.tmp)
add_test("TOPP_IDMerger_2_out1" ${DIFF} -in1 IDMerger_2_output.tmp -in2 ${DATA_DIR_TOPP}/IDMerger_2_output.idXML )
set_tests_properties("TOPP_IDMerger_2_out1" PROPERTIES DEPENDS "TOPP_IDMerger_2")
# conflicting IDs:
add_test("TOPP_IDMerger_3" ${TOPP_BIN_PATH}/IDMerger -test -in ${DATA_DIR_TOPP}/IDMerger_1_input1.idXML ${DATA_DIR_TOPP}/IDMerger_1_input1.idXML -out IDMerger_3_output.tmp)
add_test("TOPP_IDMerger_3_out1" ${DIFF} -in1 IDMerger_3_output.tmp -in2 ${DATA_DIR_TOPP}/IDMerger_3_output.idXML )
set_tests_properties("TOPP_IDMerger_3_out1" PROPERTIES DEPENDS "TOPP_IDMerger_3")
# empty IDs:
add_test("TOPP_IDMerger_4" ${TOPP_BIN_PATH}/IDMerger -test -in ${DATA_DIR_TOPP}/degenerated_empty.idXML ${DATA_DIR_TOPP}/degenerated_empty.idXML -out IDMerger_4_output.tmp)
add_test("TOPP_IDMerger_4_out1" ${DIFF} -in1 IDMerger_4_output.tmp -in2 ${DATA_DIR_TOPP}/IDMerger_4_output.idXML )
set_tests_properties("TOPP_IDMerger_4_out1" PROPERTIES DEPENDS "TOPP_IDMerger_4")
# "add_to" option:
add_test("TOPP_IDMerger_5" ${TOPP_BIN_PATH}/IDMerger -test -in ${DATA_DIR_TOPP}/IDMerger_2_input1.idXML -add_to ${DATA_DIR_TOPP}/IDMerger_5_input1.idXML -out IDMerger_5_output.tmp)
add_test("TOPP_IDMerger_5_out1" ${DIFF} -in1 IDMerger_5_output.tmp -in2 ${DATA_DIR_TOPP}/IDMerger_5_output.idXML )
set_tests_properties("TOPP_IDMerger_5_out1" PROPERTIES DEPENDS "TOPP_IDMerger_5")

#------------------------------------------------------------------------------
# BaselineFilter tests
add_test("TOPP_BaselineFilter_1" ${TOPP_BIN_PATH}/BaselineFilter -test -in ${DATA_DIR_TOPP}/BaselineFilter_input.mzML -out BaselineFilter.tmp -struc_elem_length 1.5)
add_test("TOPP_BaselineFilter_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 BaselineFilter.tmp -in2 ${DATA_DIR_TOPP}/BaselineFilter_output.mzML )
set_tests_properties("TOPP_BaselineFilter_1_out1" PROPERTIES DEPENDS "TOPP_BaselineFilter_1")

#------------------------------------------------------------------------------
# ConsensusMapNormalizer tests
add_test("TOPP_ConsensusMapNormalizer_1" ${TOPP_BIN_PATH}/ConsensusMapNormalizer -test -in ${DATA_DIR_TOPP}/ConsensusMapNormalizer_input.consensusXML -out ConsensusMapNormalizer.tmp)
add_test("TOPP_ConsensusMapNormalizer_1_out1" ${DIFF} -in1 ConsensusMapNormalizer.tmp -in2 ${DATA_DIR_TOPP}/ConsensusMapNormalizer_output.consensusXML )
set_tests_properties("TOPP_ConsensusMapNormalizer_1_out1" PROPERTIES DEPENDS "TOPP_ConsensusMapNormalizer_1")

#------------------------------------------------------------------------------
# MapNormalizer tests
add_test("TOPP_MapNormalizer_1" ${TOPP_BIN_PATH}/MapNormalizer -test -in ${DATA_DIR_TOPP}/MapNormalizer_input.mzML -out MapNormalizer.tmp)
add_test("TOPP_MapNormalizer_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 MapNormalizer.tmp -in2 ${DATA_DIR_TOPP}/MapNormalizer_output.mzML )
set_tests_properties("TOPP_MapNormalizer_1_out1" PROPERTIES DEPENDS "TOPP_MapNormalizer_1")

#------------------------------------------------------------------------------
# Decharger tests
add_test("TOPP_Decharger_1" ${TOPP_BIN_PATH}/Decharger -test -in ${DATA_DIR_TOPP}/Decharger_input.featureXML -ini ${DATA_DIR_TOPP}/Decharger_input.ini -out_fm Decharger_output_fm.tmp -out_cm Decharger_output.tmp -outpairs Decharger_p_output.tmp)
add_test("TOPP_Decharger_1_out1" ${DIFF} -whitelist "id=" "href=" -in1 Decharger_output.tmp -in2 ${DATA_DIR_TOPP}/Decharger_output.consensusXML )
set_tests_properties("TOPP_Decharger_1_out1" PROPERTIES DEPENDS "TOPP_Decharger_1" )

#------------------------------------------------------------------------------
# DTAExtractor tests
add_test("TOPP_DTAExtractor_1" ${TOPP_BIN_PATH}/DTAExtractor -test -in ${DATA_DIR_TOPP}/DTAExtractor_1_input.mzML -out DTAExtractor -rt :61)
add_test("TOPP_DTAExtractor_1_out1" ${DIFF} -in1 DTAExtractor_RT60.dta -in2 ${DATA_DIR_TOPP}/DTAExtractor_1_output.dta )
set_tests_properties("TOPP_DTAExtractor_1_out1" PROPERTIES DEPENDS "TOPP_DTAExtractor_1")
add_test("TOPP_DTAExtractor_2" ${TOPP_BIN_PATH}/DTAExtractor -test -in ${DATA_DIR_TOPP}/DTAExtractor_1_input.mzML -out DTAExtractor -level 1)
add_test("TOPP_DTAExtractor_2_out1" ${DIFF} -in1 DTAExtractor_RT60.dta -in2 ${DATA_DIR_TOPP}/DTAExtractor_2_output.dta )
set_tests_properties("TOPP_DTAExtractor_2_out1" PROPERTIES DEPENDS "TOPP_DTAExtractor_2")
add_test("TOPP_DTAExtractor_3" ${TOPP_BIN_PATH}/DTAExtractor -test -in ${DATA_DIR_TOPP}/DTAExtractor_1_input.mzML -out DTAExtractor -level 2 -mz :1000)
add_test("TOPP_DTAExtractor_3_out1" ${DIFF} -in1 DTAExtractor_RT140_MZ5.dta -in2 ${DATA_DIR_TOPP}/DTAExtractor_3_output.dta )
set_tests_properties("TOPP_DTAExtractor_3_out1" PROPERTIES DEPENDS "TOPP_DTAExtractor_3")

#------------------------------------------------------------------------------
# MassTraceExtractor test
add_test("TOPP_MassTraceExtractor_1" ${TOPP_BIN_PATH}/MassTraceExtractor -test -ini ${DATA_DIR_TOPP}/MassTraceExtractor.ini -in ${DATA_DIR_TOPP}/MassTraceExtractor_1_input.mzML -out MassTraceExtractor_1.tmp)
add_test("TOPP_MassTraceExtractor_1_out" ${DIFF} -whitelist "id=" -in1 MassTraceExtractor_1.tmp -in2 ${DATA_DIR_TOPP}/MassTraceExtractor_1_output.featureXML )
set_tests_properties("TOPP_MassTraceExtractor_1_out" PROPERTIES DEPENDS "TOPP_MassTraceExtractor_1")
add_test("TOPP_MassTraceExtractor_2" ${TOPP_BIN_PATH}/MassTraceExtractor -test -ini ${DATA_DIR_TOPP}/MassTraceExtractor_2.ini -in ${DATA_DIR_TOPP}/MassTraceExtractor_1_input.mzML -out MassTraceExtractor_2.tmp)
add_test("TOPP_MassTraceExtractor_2_out" ${DIFF} -whitelist "id=" -in1 MassTraceExtractor_2.tmp -in2 ${DATA_DIR_TOPP}/MassTraceExtractor_2_output.featureXML )
set_tests_properties("TOPP_MassTraceExtractor_2_out" PROPERTIES DEPENDS "TOPP_MassTraceExtractor_2")

#------------------------------------------------------------------------------
# FeatureFinderMetabo test
add_test("TOPP_FeatureFinderMetabo_1" ${TOPP_BIN_PATH}/FeatureFinderMetabo -test -ini ${DATA_DIR_TOPP}/FeatureFinderMetabo.ini -in ${DATA_DIR_TOPP}/FeatureFinderMetabo_1_input.mzML -out FeatureFinderMetabo_1.tmp)
add_test("TOPP_FeatureFinderMetabo_1_out" ${DIFF} -whitelist "id=" -in1 FeatureFinderMetabo_1.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinderMetabo_1_output.featureXML)
set_tests_properties("TOPP_FeatureFinderMetabo_1_out" PROPERTIES DEPENDS "TOPP_FeatureFinderMetabo_1")
add_test("TOPP_FeatureFinderMetabo_2" ${TOPP_BIN_PATH}/FeatureFinderMetabo -test -ini ${DATA_DIR_TOPP}/FeatureFinderMetabo_2_noEPD.ini -in ${DATA_DIR_TOPP}/FeatureFinderMetabo_2_input.mzML -out FeatureFinderMetabo_2.tmp)
add_test("TOPP_FeatureFinderMetabo_2_out" ${DIFF} -whitelist "id=" "completion_time" -in1 FeatureFinderMetabo_2.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinderMetabo_2_noEPD_output.featureXML)
set_tests_properties("TOPP_FeatureFinderMetabo_2_out" PROPERTIES DEPENDS "TOPP_FeatureFinderMetabo_2")

#------------------------------------------------------------------------------
# FeatureFinderCentroided test
add_test("TOPP_FeatureFinderCentroided_1" ${TOPP_BIN_PATH}/FeatureFinderCentroided -test -ini ${DATA_DIR_TOPP}/FeatureFinderCentroided_1_parameters.ini -in ${DATA_DIR_TOPP}/FeatureFinderCentroided_1_input.mzML -out FeatureFinderCentroided_1.tmp)
add_test("TOPP_FeatureFinderCentroided_1_out1" ${DIFF} -whitelist "id=" -in1 FeatureFinderCentroided_1.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinderCentroided_1_output.featureXML )
set_tests_properties("TOPP_FeatureFinderCentroided_1_out1" PROPERTIES DEPENDS "TOPP_FeatureFinderCentroided_1")

#------------------------------------------------------------------------------
# FeatureFinderIdentification test
## internal IDs only:
add_test("TOPP_FeatureFinderIdentification_1" ${TOPP_BIN_PATH}/FeatureFinderIdentification -test -in ${DATA_DIR_TOPP}/FeatureFinderIdentification_1_input.mzML -id ${DATA_DIR_TOPP}/FeatureFinderIdentification_1_input.idXML -out FeatureFinderIdentification_1.tmp -extract:mz_window 0.1 -detect:peak_width 60 -model:type none)
add_test("TOPP_FeatureFinderIdentification_1_out1" ${DIFF} -whitelist "id=" -in1 FeatureFinderIdentification_1.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinderIdentification_1_output.featureXML)
set_tests_properties("TOPP_FeatureFinderIdentification_1_out1" PROPERTIES DEPENDS "TOPP_FeatureFinderIdentification_1")
## with (faked) external IDs; fix SVM parameters to avoid randomness:
add_test("TOPP_FeatureFinderIdentification_2" ${TOPP_BIN_PATH}/FeatureFinderIdentification -test -in ${DATA_DIR_TOPP}/FeatureFinderIdentification_1_input.mzML -id ${DATA_DIR_TOPP}/FeatureFinderIdentification_2_input.idXML -id_ext ${DATA_DIR_TOPP}/FeatureFinderIdentification_1_input.idXML -out FeatureFinderIdentification_2.tmp -extract:mz_window 0.1 -detect:peak_width 60 -svm:no_selection -svm:log2_C 1 -svm:log2_gamma 1 -model:type none)
add_test("TOPP_FeatureFinderIdentification_2_out1" ${DIFF} -whitelist "id=" -in1 FeatureFinderIdentification_2.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinderIdentification_2_output.featureXML)
set_tests_properties("TOPP_FeatureFinderIdentification_2_out1" PROPERTIES DEPENDS "TOPP_FeatureFinderIdentification_2")
## with elution model fitting:
add_test("TOPP_FeatureFinderIdentification_3" ${TOPP_BIN_PATH}/FeatureFinderIdentification -test -in ${DATA_DIR_TOPP}/FeatureFinderIdentification_1_input.mzML -id ${DATA_DIR_TOPP}/FeatureFinderIdentification_1_input.idXML -out FeatureFinderIdentification_3.tmp -extract:mz_window 0.1 -detect:peak_width 60 -model:type symmetric)
add_test("TOPP_FeatureFinderIdentification_3_out1" ${DIFF} -whitelist "id=" -in1 FeatureFinderIdentification_3.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinderIdentification_3_output.featureXML)
set_tests_properties("TOPP_FeatureFinderIdentification_3_out1" PROPERTIES DEPENDS "TOPP_FeatureFinderIdentification_3")

#------------------------------------------------------------------------------
# FeatureFinderMRM test
add_test("TOPP_FeatureFinderMRM_1" ${TOPP_BIN_PATH}/FeatureFinderMRM -test -in ${DATA_DIR_TOPP}/FeatureFinderMRM_1_input.mzML -ini ${DATA_DIR_TOPP}/FeatureFinderMRM_1_parameters.ini -out FeatureFinderMRM_1.tmp)
add_test("TOPP_FeatureFinderMRM_1_out1" ${DIFF} -whitelist "id=" -in1 FeatureFinderMRM_1.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinderMRM_1_output.featureXML )
set_tests_properties("TOPP_FeatureFinderMRM_1_out1" PROPERTIES DEPENDS "TOPP_FeatureFinderMRM_1")

#------------------------------------------------------------------------------
# FeatureFinderMultiplex test
add_test("TOPP_FeatureFinderMultiplex_1" ${TOPP_BIN_PATH}/FeatureFinderMultiplex -test -in ${DATA_DIR_TOPP}/FeatureFinderMultiplex_1_input.mzML -ini ${DATA_DIR_TOPP}/FeatureFinderMultiplex_1_parameters.ini -out FeatureFinderMultiplex_1.tmp -out_features FeatureFinderMultiplex_2.tmp -out_mzq FeatureFinderMultiplex_3.tmp)
add_test("TOPP_FeatureFinderMultiplex_1_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureFinderMultiplex_1.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_1_output.consensusXML )
set_tests_properties("TOPP_FeatureFinderMultiplex_1_out1" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_1")
add_test("TOPP_FeatureFinderMultiplex_1_out2" ${DIFF} -whitelist "id=" "href=" -in1 FeatureFinderMultiplex_2.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_1_output.featureXML )
set_tests_properties("TOPP_FeatureFinderMultiplex_1_out2" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_1")
add_test("TOPP_FeatureFinderMultiplex_1_out3" ${DIFF} -whitelist "id=" "href=" -in1 FeatureFinderMultiplex_3.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_1_output.mzq )
set_tests_properties("TOPP_FeatureFinderMultiplex_1_out3" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_1")
add_test("TOPP_FeatureFinderMultiplex_2" ${TOPP_BIN_PATH}/FeatureFinderMultiplex -test -in ${DATA_DIR_TOPP}/FeatureFinderMultiplex_2_input.mzML -ini ${DATA_DIR_TOPP}/FeatureFinderMultiplex_2_parameters.ini -out FeatureFinderMultiplex_4.tmp -out_features FeatureFinderMultiplex_5.tmp -out_mzq FeatureFinderMultiplex_6.tmp)
add_test("TOPP_FeatureFinderMultiplex_2_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureFinderMultiplex_4.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_2_output.consensusXML )
set_tests_properties("TOPP_FeatureFinderMultiplex_2_out1" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_2")
add_test("TOPP_FeatureFinderMultiplex_2_out2" ${DIFF} -whitelist "id=" "href=" -in1 FeatureFinderMultiplex_5.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_2_output.featureXML )
set_tests_properties("TOPP_FeatureFinderMultiplex_2_out2" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_2")
add_test("TOPP_FeatureFinderMultiplex_2_out3" ${DIFF} -whitelist "id=" "href=" -in1 FeatureFinderMultiplex_6.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_2_output.mzq )
set_tests_properties("TOPP_FeatureFinderMultiplex_2_out3" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_2")
add_test("TOPP_FeatureFinderMultiplex_3" ${TOPP_BIN_PATH}/FeatureFinderMultiplex -test -in ${DATA_DIR_TOPP}/FeatureFinderMultiplex_3_input.mzML -ini ${DATA_DIR_TOPP}/FeatureFinderMultiplex_3_parameters.ini -out FeatureFinderMultiplex_7.tmp -out_features FeatureFinderMultiplex_8.tmp -out_mzq FeatureFinderMultiplex_9.tmp)
add_test("TOPP_FeatureFinderMultiplex_3_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureFinderMultiplex_7.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_3_output.consensusXML )
set_tests_properties("TOPP_FeatureFinderMultiplex_3_out1" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_3")
add_test("TOPP_FeatureFinderMultiplex_3_out2" ${DIFF} -whitelist "id=" "href=" -in1 FeatureFinderMultiplex_8.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_3_output.featureXML )
set_tests_properties("TOPP_FeatureFinderMultiplex_3_out2" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_3")
add_test("TOPP_FeatureFinderMultiplex_3_out3" ${DIFF} -whitelist "id=" "href=" -in1 FeatureFinderMultiplex_9.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_3_output.mzq )
set_tests_properties("TOPP_FeatureFinderMultiplex_3_out3" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_3")
add_test("TOPP_FeatureFinderMultiplex_4" ${TOPP_BIN_PATH}/FeatureFinderMultiplex -test -in ${DATA_DIR_TOPP}/FeatureFinderMultiplex_3_input.mzML -ini ${DATA_DIR_TOPP}/FeatureFinderMultiplex_3_parameters.ini -out FeatureFinderMultiplex_10.tmp -out_features FeatureFinderMultiplex_11.tmp -out_mzq FeatureFinderMultiplex_12.tmp)
add_test("TOPP_FeatureFinderMultiplex_4_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureFinderMultiplex_10.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_3_output.consensusXML )
set_tests_properties("TOPP_FeatureFinderMultiplex_4_out1" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_4")
add_test("TOPP_FeatureFinderMultiplex_4_out2" ${DIFF} -whitelist "id=" "href=" -in1 FeatureFinderMultiplex_11.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_3_output.featureXML )
set_tests_properties("TOPP_FeatureFinderMultiplex_4_out2" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_4")
add_test("TOPP_FeatureFinderMultiplex_4_out3" ${DIFF} -whitelist "id=" "href=" -in1 FeatureFinderMultiplex_12.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_3_output.mzq )
set_tests_properties("TOPP_FeatureFinderMultiplex_4_out3" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_4")

#------------------------------------------------------------------------------
# FileConverter tests
add_test("TOPP_FileConverter_1" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_1_input.mzData -out FileConverter_1.tmp -out_type mzML)
add_test("TOPP_FileConverter_1_out1" ${DIFF} -whitelist "id=" "location=" "href=" -in1 FileConverter_1.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_1_output.mzML )
set_tests_properties("TOPP_FileConverter_1_out1" PROPERTIES DEPENDS "TOPP_FileConverter_1")
add_test("TOPP_FileConverter_2" ${TOPP_BIN_PATH}/FileConverter -test  -in ${DATA_DIR_TOPP}/FileConverter_2_input.someInputDTA2D -in_type dta2d -out FileConverter_2.tmp -out_type mzML)
add_test("TOPP_FileConverter_2_out1" ${DIFF} -whitelist "id=" "location=" "href=" -in1 FileConverter_2.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_2_output.mzML )
set_tests_properties("TOPP_FileConverter_2_out1" PROPERTIES DEPENDS "TOPP_FileConverter_2")
add_test("TOPP_FileConverter_3" ${TOPP_BIN_PATH}/FileConverter -test  -in ${DATA_DIR_TOPP}/FileConverter_3_input.featureXML -out FileConverter_3.tmp -out_type mzML)
add_test("TOPP_FileConverter_3_out1" ${DIFF} -in1 FileConverter_3.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_3_output.mzML )
set_tests_properties("TOPP_FileConverter_3_out1" PROPERTIES DEPENDS "TOPP_FileConverter_3")
add_test("TOPP_FileConverter_4" ${TOPP_BIN_PATH}/FileConverter -test  -in ${DATA_DIR_TOPP}/FileConverter_4_input.mzXML -out FileConverter_4.tmp -out_type mzML)
add_test("TOPP_FileConverter_4_out1" ${DIFF} -whitelist "location=" -in1 FileConverter_4.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_4_output.mzML )
set_tests_properties("TOPP_FileConverter_4_out1" PROPERTIES DEPENDS "TOPP_FileConverter_4")
add_test("TOPP_FileConverter_5" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_5_input.mzML -out FileConverter_5.tmp -out_type mzXML)
add_test("TOPP_FileConverter_5_out1" ${DIFF} -in1 FileConverter_5.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_5_output.mzXML )
set_tests_properties("TOPP_FileConverter_5_out1" PROPERTIES DEPENDS "TOPP_FileConverter_5")
add_test("TOPP_FileConverter_6" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_6_input.mzML -out FileConverter_6.tmp -out_type mzXML)
add_test("TOPP_FileConverter_6_out1" ${DIFF} -in1 FileConverter_6.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_6_output.mzXML )
set_tests_properties("TOPP_FileConverter_6_out1" PROPERTIES DEPENDS "TOPP_FileConverter_6")
add_test("TOPP_FileConverter_7" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_7_input.consensusXML -out FileConverter_7.tmp -out_type featureXML)
add_test("TOPP_FileConverter_7_out1" ${DIFF} -in1 FileConverter_7.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_7_output.featureXML )
set_tests_properties("TOPP_FileConverter_7_out1" PROPERTIES DEPENDS "TOPP_FileConverter_7")
add_test("TOPP_FileConverter_8" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_8_input.mzML -out FileConverter_8.tmp -out_type mzData)
add_test("TOPP_FileConverter_8_out1" ${DIFF} -whitelist "<pathToFile>" -in1 FileConverter_8.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_8_output.mzData )
set_tests_properties("TOPP_FileConverter_8_out1" PROPERTIES DEPENDS "TOPP_FileConverter_8")
add_test("TOPP_FileConverter_9" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_9_input.consensusXML -out FileConverter_9.tmp -out_type featureXML)
add_test("TOPP_FileConverter_9_out1" ${DIFF} -in1 FileConverter_9.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_9_output.featureXML )
set_tests_properties("TOPP_FileConverter_9_out1" PROPERTIES DEPENDS "TOPP_FileConverter_9")
add_test("TOPP_FileConverter_10" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_10_input.edta -no_progress -out FileConverter_10.tmp -out_type featureXML)
add_test("TOPP_FileConverter_10_out1" ${DIFF} -whitelist "id=" -in1 FileConverter_10.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_10_output.featureXML )
set_tests_properties("TOPP_FileConverter_10_out1" PROPERTIES DEPENDS "TOPP_FileConverter_10")
#SpecArray
add_test("TOPP_FileConverter_11" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_11_input.peplist -no_progress -out FileConverter_11.tmp -out_type featureXML)
add_test("TOPP_FileConverter_11_out1" ${DIFF} -whitelist "id=" -in1 FileConverter_11.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_11_output.featureXML )
set_tests_properties("TOPP_FileConverter_11_out1" PROPERTIES DEPENDS "TOPP_FileConverter_11")
#msInspect
add_test("TOPP_FileConverter_12" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_12_input.peptides.tsv -no_progress -out FileConverter_12.tmp -out_type featureXML)
add_test("TOPP_FileConverter_12_out1" ${DIFF} -whitelist "id=" -in1 FileConverter_12.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_12_output.featureXML )
set_tests_properties("TOPP_FileConverter_12_out1" PROPERTIES DEPENDS "TOPP_FileConverter_12")
#kroenik
add_test("TOPP_FileConverter_13" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_13_input.peptides.kroenik -no_progress -out FileConverter_13.tmp -out_type featureXML)
add_test("TOPP_FileConverter_13_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileConverter_13.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_13_output.featureXML )
set_tests_properties("TOPP_FileConverter_13_out1" PROPERTIES DEPENDS "TOPP_FileConverter_13")
# featureXML to consensusXML
add_test("TOPP_FileConverter_14" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_9_output.featureXML -no_progress -out FileConverter_14.tmp -out_type consensusXML)
add_test("TOPP_FileConverter_14_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileConverter_14.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_14_output.consensusXML )
set_tests_properties("TOPP_FileConverter_14_out1" PROPERTIES DEPENDS "TOPP_FileConverter_14")
# EDTA to consensusXML
add_test("TOPP_FileConverter_15" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_10_input.edta -no_progress -out FileConverter_15.tmp -out_type consensusXML)
add_test("TOPP_FileConverter_15_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileConverter_15.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_15_output.consensusXML )
set_tests_properties("TOPP_FileConverter_15_out1" PROPERTIES DEPENDS "TOPP_FileConverter_15")
# consensus EDTA to consensusXML
add_test("TOPP_FileConverter_16" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_16_input.edta -no_progress -out FileConverter_16.tmp -out_type consensusXML)
add_test("TOPP_FileConverter_16_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileConverter_16.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_16_output.consensusXML )
set_tests_properties("TOPP_FileConverter_16_out1" PROPERTIES DEPENDS "TOPP_FileConverter_16")
# consensusXML to ibspectra without out_type
add_test("TOPP_FileConverter_17" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_17_input.consensusXML -no_progress -out FileConverter_17.csv)
add_test("TOPP_FileConverter_17_out1" ${DIFF} -in1 FileConverter_17.csv -in2 ${DATA_DIR_TOPP}/FileConverter_17_output.ibspectra.csv )
set_tests_properties("TOPP_FileConverter_17_out1" PROPERTIES DEPENDS "TOPP_FileConverter_17")
# consensusXML to ibspectra with out_type -> note that input and expected output are the same as in FileConverter_17
add_test("TOPP_FileConverter_18" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_17_input.consensusXML -no_progress -out FileConverter_18.tmp -out_type csv)
add_test("TOPP_FileConverter_18_out1" ${DIFF} -in1 FileConverter_18.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_17_output.ibspectra.csv )
set_tests_properties("TOPP_FileConverter_18_out1" PROPERTIES DEPENDS "TOPP_FileConverter_18")
# Purpose: test the -write_mzML_index and the -process_lowmemory options (for mzML and mzXML)
add_test("TOPP_FileConverter_19" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out FileConverter_19.tmp  -write_mzML_index -process_lowmemory -in_type mzML -out_type mzML)
add_test("TOPP_FileConverter_19_out" ${DIFF} -in1 FileConverter_19.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_19_output.mzML )
set_tests_properties("TOPP_FileConverter_19_out" PROPERTIES DEPENDS "TOPP_FileConverter_19")
# Purpose: test the writing and reading back the same featureXML 
add_test("TOPP_FileConverter_20" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_20_input.featureXML -out FileConverter_20.tmp -in_type featureXML -out_type featureXML)
add_test("TOPP_FileConverter_20_out" ${DIFF} -in1 FileConverter_20.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_20_output.featureXML )
set_tests_properties("TOPP_FileConverter_20_out" PROPERTIES DEPENDS "TOPP_FileConverter_20")
add_test("TOPP_FileConverter_21" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_4_input.mzXML -out FileConverter_21.tmp -out_type mzML -process_lowmemory)
# TODO sourceFile does not seem to be handled correctly, we should fix this! 
add_test("TOPP_FileConverter_21_out" ${DIFF} -whitelist "location=" "MS:100056" "SHA-1" -in1 FileConverter_21.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_4_output.mzML)
set_tests_properties("TOPP_FileConverter_21_out" PROPERTIES DEPENDS "TOPP_FileConverter_21")
add_test("TOPP_FileConverter_22" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_8_input.mzML -out FileConverter_22.tmp -out_type mzData -process_lowmemory)
set_tests_properties("TOPP_FileConverter_22" PROPERTIES WILL_FAIL 1) ## low memory not implemented yet
# Corrupt file with empty spectra and chromatograms, it should be handled gracefully
add_test("TOPP_FileConverter_23" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_23_input.mzML -out FileConverter_23.tmp -out_type mzML)
add_test("TOPP_FileConverter_23_out" ${DIFF} -in1 FileConverter_23.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_23_output.mzML -whitelist "location=")
set_tests_properties("TOPP_FileConverter_23_out" PROPERTIES DEPENDS "TOPP_FileConverter_23")
# Purpose: test the numpress conversion options (using full accuracy and 1ppm accuracy @ 100 m/z)
add_test("TOPP_FileConverter_24" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out FileConverter_24.tmp  -write_mzML_index -process_lowmemory -in_type mzML -out_type mzML -lossy_compression -lossy_mass_accuracy 0.0001)
add_test("TOPP_FileConverter_24_out" ${DIFF} -in1 FileConverter_24.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_24_output.mzML )
set_tests_properties("TOPP_FileConverter_24_out" PROPERTIES DEPENDS "TOPP_FileConverter_24")
add_test("TOPP_FileConverter_25" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out FileConverter_25.tmp  -write_mzML_index -process_lowmemory -in_type mzML -out_type mzML -lossy_compression)
add_test("TOPP_FileConverter_25_out" ${DIFF} -in1 FileConverter_25.tmp -in2 ${DATA_DIR_TOPP}/FileConverter_25_output.mzML )
set_tests_properties("TOPP_FileConverter_25_out" PROPERTIES DEPENDS "TOPP_FileConverter_25")

#------------------------------------------------------------------------------
# FileFilter tests
add_test("TOPP_FileFilter_1" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out FileFilter_1.tmp -rt :30 -mz :1000 -int :20000 -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_1_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_1.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_1_output.mzML )
set_tests_properties("TOPP_FileFilter_1_out1" PROPERTIES DEPENDS "TOPP_FileFilter_1")
add_test("TOPP_FileFilter_2" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out FileFilter_2.tmp -rt 30: -mz 1000: -int 100: -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_2_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_2.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_2_output.mzML )
set_tests_properties("TOPP_FileFilter_2_out1" PROPERTIES DEPENDS "TOPP_FileFilter_2")
add_test("TOPP_FileFilter_3" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out FileFilter_3.tmp -peak_options:level 2 -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_3_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_3.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_3_output.mzML )
set_tests_properties("TOPP_FileFilter_3_out1" PROPERTIES DEPENDS "TOPP_FileFilter_3")
add_test("TOPP_FileFilter_4" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_4_input.mzML -out FileFilter_4.tmp -spectra:remove_zoom -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_4_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_4.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_4_output.mzML )
set_tests_properties("TOPP_FileFilter_4_out1" PROPERTIES DEPENDS "TOPP_FileFilter_4")
add_test("TOPP_FileFilter_5" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_5_input.featureXML -out FileFilter_5.tmp -rt :1000 -mz :480 -int :79000 -f_and_c:charge :3 -feature:q :0.6 -in_type featureXML -out_type featureXML)
add_test("TOPP_FileFilter_5_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_5.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_5_out.featureXML )
set_tests_properties("TOPP_FileFilter_5_out1" PROPERTIES DEPENDS "TOPP_FileFilter_5")
add_test("TOPP_FileFilter_6" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_5_input.featureXML -out FileFilter_6.tmp -rt 1000: -mz 440: -int 70000: -f_and_c:charge 3: -feature:q 0.51: -in_type featureXML -out_type featureXML)
add_test("TOPP_FileFilter_6_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_6.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_6_out.featureXML )
set_tests_properties("TOPP_FileFilter_6_out1" PROPERTIES DEPENDS "TOPP_FileFilter_6")
add_test("TOPP_FileFilter_7" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_7_input.mzML -out FileFilter_7.tmp -int 7000: -peak_options:level 1 2 3 -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_7_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_7.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_7_output.mzML )
set_tests_properties("TOPP_FileFilter_7_out1" PROPERTIES DEPENDS "TOPP_FileFilter_7")
add_test("TOPP_FileFilter_8" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_8_input.consensusXML -out FileFilter_8.tmp -rt 600:1400 -mz 700:2300 -int 1100:6000 -in_type consensusXML -out_type consensusXML)
add_test("TOPP_FileFilter_8_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_8.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_8_output.consensusXML )
set_tests_properties("TOPP_FileFilter_8_out1" PROPERTIES DEPENDS "TOPP_FileFilter_8")
add_test("TOPP_FileFilter_9" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_9_input.mzML -out FileFilter_9.tmp -spectra:remove_mode SelectedIonMonitoring -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_9_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_9.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_9_output.mzML)
set_tests_properties("TOPP_FileFilter_9_out1" PROPERTIES DEPENDS "TOPP_FileFilter_9")
add_test("TOPP_FileFilter_10" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_10_input.mzML -out FileFilter_10.tmp -spectra:remove_activation "Collision-induced dissociation" -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_10_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_10.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_10_output.mzML )
set_tests_properties("TOPP_FileFilter_10_out1" PROPERTIES DEPENDS "TOPP_FileFilter_10")
add_test("TOPP_FileFilter_11" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_11_input.mzML -out FileFilter_11.tmp -spectra:remove_activation "Plasma desorption" -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_11_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_11.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_11_output.mzML )
set_tests_properties("TOPP_FileFilter_11_out1" PROPERTIES DEPENDS "TOPP_FileFilter_11")
add_test("TOPP_FileFilter_12" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_12_input.mzML -out FileFilter_12.tmp -peak_options:remove_chromatograms -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_12_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_12.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_12_output.mzML )
set_tests_properties("TOPP_FileFilter_12_out1" PROPERTIES DEPENDS "TOPP_FileFilter_12")
add_test("TOPP_FileFilter_13" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_13_input.consensusXML -out FileFilter_13.tmp -in_type consensusXML -out_type featureXML -consensus:map 2)
add_test("TOPP_FileFilter_13_out1" ${DIFF} -whitelist "id=" "href="  -in1 FileFilter_13.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_13_output.featureXML )
set_tests_properties("TOPP_FileFilter_13_out1" PROPERTIES DEPENDS "TOPP_FileFilter_13")
add_test("TOPP_FileFilter_14" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_14_input.consensusXML -out FileFilter_14.tmp -in_type consensusXML -out_type consensusXML -consensus:map 0 2)
add_test("TOPP_FileFilter_14_out1" ${DIFF} -whitelist "id=" "href="  -in1 FileFilter_14.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_14_output.consensusXML )
set_tests_properties("TOPP_FileFilter_14_out1" PROPERTIES DEPENDS "TOPP_FileFilter_14")
add_test("TOPP_FileFilter_15" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_15_input.featureXML -out FileFilter_15.tmp -id:sequences_whitelist Oxidation -id:remove_unassigned_ids)
add_test("TOPP_FileFilter_15_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_15.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_15_output.featureXML )
set_tests_properties("TOPP_FileFilter_15_out1" PROPERTIES DEPENDS "TOPP_FileFilter_15")
add_test("TOPP_FileFilter_16" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_15_input.featureXML -out FileFilter_16.tmp -id:sequences_whitelist Oxidation -id:remove_unassigned_ids -mz 400:600 -rt 3000:4000)
add_test("TOPP_FileFilter_16_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_16.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_16_output.featureXML )
set_tests_properties("TOPP_FileFilter_16_out1" PROPERTIES DEPENDS "TOPP_FileFilter_16")
add_test("TOPP_FileFilter_17" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_15_input.featureXML -out FileFilter_17.tmp -id:remove_annotated_features -mz 400:600 -rt 3000:4000)
add_test("TOPP_FileFilter_17_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_17.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_17_output.featureXML )
set_tests_properties("TOPP_FileFilter_17_out1" PROPERTIES DEPENDS "TOPP_FileFilter_17")
add_test("TOPP_FileFilter_18" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_18_input.consensusXML -out FileFilter_18.tmp -id:sequences_whitelist Oxidation -id:remove_unassigned_ids)
add_test("TOPP_FileFilter_18_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_18.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_18_output.consensusXML )
set_tests_properties("TOPP_FileFilter_18_out1" PROPERTIES DEPENDS "TOPP_FileFilter_18")
add_test("TOPP_FileFilter_19" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_18_input.consensusXML -out FileFilter_19.tmp -id:remove_unannotated_features -mz 400:600 -rt 3000:4000)
add_test("TOPP_FileFilter_19_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_19.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_19_output.consensusXML )
set_tests_properties("TOPP_FileFilter_19_out1" PROPERTIES DEPENDS "TOPP_FileFilter_19")
add_test("TOPP_FileFilter_20" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_15_input.featureXML -out FileFilter_20.tmp -id:accessions_whitelist YDL217C -id:remove_unassigned_ids)
add_test("TOPP_FileFilter_20_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_20.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_20_output.featureXML )
set_tests_properties("TOPP_FileFilter_20_out1" PROPERTIES DEPENDS "TOPP_FileFilter_20")
add_test("TOPP_FileFilter_21" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_15_input.featureXML -out FileFilter_21.tmp -id:remove_unassigned_ids -id:remove_unannotated_features -id:keep_best_score_id)
add_test("TOPP_FileFilter_21_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_21.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_21_output.featureXML )
set_tests_properties("TOPP_FileFilter_21_out1" PROPERTIES DEPENDS "TOPP_FileFilter_21")
add_test("TOPP_FileFilter_22" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_22_input.consensusXML -out FileFilter_22.tmp -f_and_c:remove_meta distinct_charges gt "1,2")
add_test("TOPP_FileFilter_22_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_22.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_22_output.consensusXML )
set_tests_properties("TOPP_FileFilter_22_out1" PROPERTIES DEPENDS "TOPP_FileFilter_22")
add_test("TOPP_FileFilter_23" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_22_input.consensusXML -out FileFilter_23.tmp -f_and_c:remove_meta distinct_charges_size gt 2)
add_test("TOPP_FileFilter_23_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_23.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_22_output.consensusXML )
set_tests_properties("TOPP_FileFilter_23_out1" PROPERTIES DEPENDS "TOPP_FileFilter_23")
add_test("TOPP_FileFilter_24" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_22_input.consensusXML -out FileFilter_24.tmp -f_and_c:remove_meta DOESNOTEXIST lt "whatever")
add_test("TOPP_FileFilter_24_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_24.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_24_output.consensusXML )
set_tests_properties("TOPP_FileFilter_24_out1" PROPERTIES DEPENDS "TOPP_FileFilter_24")

add_test("TOPP_FileFilter_25" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_25_input.mzML.gz -id:blacklist ${DATA_DIR_TOPP}/FileFilter_25_input.idXML -out FileFilter_25.tmp -id:mz 0.05 -id:rt 1 -id:blacklist_imperfect)
add_test("TOPP_FileFilter_25_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_25.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_25_output.mzML )
set_tests_properties("TOPP_FileFilter_25_out1" PROPERTIES DEPENDS "TOPP_FileFilter_25")

add_test("TOPP_FileFilter_26" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_25_input.mzML.gz -id:blacklist ${DATA_DIR_TOPP}/FileFilter_25_input.idXML -out FileFilter_26.tmp -id:blacklist_imperfect)
add_test("TOPP_FileFilter_26_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_26.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_25_output.mzML ) ## with missing id:mz and id:rt -- to test the default values
set_tests_properties("TOPP_FileFilter_26_out1" PROPERTIES DEPENDS "TOPP_FileFilter_26")

add_test("TOPP_FileFilter_27" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_25_input.mzML.gz -id:blacklist ${DATA_DIR_TOPP}/FileFilter_25_input.idXML -out FileFilter_27.tmp  -id:mz 0.05 -id:rt 1 )
set_tests_properties("TOPP_FileFilter_27" PROPERTIES WILL_FAIL 1) ## has 2 imperfect matches

add_test("TOPP_FileFilter_28" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_28_input.mzML.gz -peak_options:pc_mz_range 832:836 -out FileFilter_28.tmp -peak_options:level 2) ## only one MS2 remaining
add_test("TOPP_FileFilter_28_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_28.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_28_output.mzML )
set_tests_properties("TOPP_FileFilter_28_out1" PROPERTIES DEPENDS "TOPP_FileFilter_28")

add_test("TOPP_FileFilter_29" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_28_input.mzML.gz -peak_options:pc_mz_range 832:836 -out FileFilter_29.tmp -peak_options:level 1 2) ## only one MS2 remaining, and lots of MS1
add_test("TOPP_FileFilter_29_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_29.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_29_output.mzML )
set_tests_properties("TOPP_FileFilter_29_out1" PROPERTIES DEPENDS "TOPP_FileFilter_29")

add_test("TOPP_FileFilter_30" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_28_input.mzML.gz -peak_options:pc_mz_range 832:836 -out FileFilter_30.tmp -peak_options:level 1 2) ## no MS2 remaining since -mz filter empties it and it will be deleted, and lots of MS1
add_test("TOPP_FileFilter_30_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_30.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_30_output.mzML )
set_tests_properties("TOPP_FileFilter_30_out1" PROPERTIES DEPENDS "TOPP_FileFilter_30")

add_test("TOPP_FileFilter_31" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_31_34_input.mzML -spectra:remove_isolation_window_width :2 -out FileFilter_31.tmp )
add_test("TOPP_FileFilter_31_out" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_31.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_31_remove_collision_energy.mzML )
set_tests_properties("TOPP_FileFilter_31_out" PROPERTIES DEPENDS "TOPP_FileFilter_31")

add_test("TOPP_FileFilter_32" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_31_34_input.mzML -spectra:select_isolation_window_width :2 -out FileFilter_32.tmp )
add_test("TOPP_FileFilter_32_out" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_32.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_32_select_collision_energy.mzML )
set_tests_properties("TOPP_FileFilter_32_out" PROPERTIES DEPENDS "TOPP_FileFilter_32")

add_test("TOPP_FileFilter_33" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_31_34_input.mzML -spectra:remove_collision_energy :35 -out FileFilter_33.tmp )
add_test("TOPP_FileFilter_33_out" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_33.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_33_remove_isolation_window.mzML )
set_tests_properties("TOPP_FileFilter_33_out" PROPERTIES DEPENDS "TOPP_FileFilter_33")

add_test("TOPP_FileFilter_34" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_31_34_input.mzML -spectra:select_collision_energy :35 -out FileFilter_34.tmp )
add_test("TOPP_FileFilter_34_out" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_34.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_34_select_isolation_window.mzML )
set_tests_properties("TOPP_FileFilter_34_out" PROPERTIES DEPENDS "TOPP_FileFilter_34")

# Purpose: test the -peak_options:indexed_file to write mzML files with/without indices
add_test("TOPP_FileFilter_35" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out FileFilter_35.tmp -peak_options:indexed_file true -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_35_out" ${DIFF} -in1 FileFilter_35.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_35_output.mzML )
set_tests_properties("TOPP_FileFilter_35_out" PROPERTIES DEPENDS "TOPP_FileFilter_35")
add_test("TOPP_FileFilter_36" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out FileFilter_36.tmp -peak_options:indexed_file false -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_36_out" ${DIFF} -in1 FileFilter_36.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_36_output.mzML )
set_tests_properties("TOPP_FileFilter_36_out" PROPERTIES DEPENDS "TOPP_FileFilter_36")

# Purpose: test the -numpress options for mzML files
add_test("TOPP_FileFilter_37" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out FileFilter_37.tmp -test -peak_options:numpress:intensity slof -peak_options:numpress:masstime linear -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_37_out" ${DIFF} -in1 FileFilter_37.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_37_output.mzML )
set_tests_properties("TOPP_FileFilter_37_out" PROPERTIES DEPENDS "TOPP_FileFilter_37")
add_test("TOPP_FileFilter_38" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out FileFilter_38.tmp -test -peak_options:numpress:intensity pic -peak_options:numpress:masstime linear -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_38_out" ${DIFF} -in1 FileFilter_38.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_38_output.mzML )
set_tests_properties("TOPP_FileFilter_38_out" PROPERTIES DEPENDS "TOPP_FileFilter_38")
add_test("TOPP_FileFilter_39" ${TOPP_BIN_PATH}/FileFilter -in ${DATA_DIR_TOPP}/FileFilter_39_input_corrupt.mzML -out FileFilter_39.tmp -test)
set_tests_properties("TOPP_FileFilter_39" PROPERTIES WILL_FAIL 1) ## this mzML has a corrupt numpress encoded string and should fail

add_test("TOPP_FileFilter_40" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_40_input.mzML -out FileFilter_40.tmp -spectra:select_polarity positive -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_40_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_40.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_40_output.mzML)
set_tests_properties("TOPP_FileFilter_40_out1" PROPERTIES DEPENDS "TOPP_FileFilter_40")
add_test("TOPP_FileFilter_41" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_40_input.mzML -out FileFilter_41.tmp -spectra:select_polarity negative -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_41_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_41.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_41_output.mzML)
set_tests_properties("TOPP_FileFilter_41_out1" PROPERTIES DEPENDS "TOPP_FileFilter_41")
add_test("TOPP_FileFilter_42" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_40_input.mzML -out FileFilter_42.tmp -spectra:select_polarity "" -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_42_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_42.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_42_output.mzML)
set_tests_properties("TOPP_FileFilter_42_out1" PROPERTIES DEPENDS "TOPP_FileFilter_42")
# Test XML escaping of our XML handlers
add_test("TOPP_FileFilter_43" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_43_input.mzML -out FileFilter_43.tmp.mzML -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_43_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_43.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_43_output.mzML)
# We check whether the XML we just produced can be read again
add_test("TOPP_FileFilter_43_read_again" ${TOPP_BIN_PATH}/FileFilter -test -in FileFilter_43.tmp.mzML -out FileFilter_43.dummy.tmp  -in_type mzML -out_type mzML)
set_tests_properties("TOPP_FileFilter_43_out1" PROPERTIES DEPENDS "TOPP_FileFilter_43")
set_tests_properties("TOPP_FileFilter_43_read_again" PROPERTIES DEPENDS "TOPP_FileFilter_43")
# test corrupt numpress file, see #1641
add_test("TOPP_FileFilter_44" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_44_input.mzML -out FileFilter_44.tmp -test -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_44_out" ${DIFF} -in1 FileFilter_44.tmp -in2 ${DATA_DIR_TOPP}/FileFilter_44_output.mzML )
set_tests_properties("TOPP_FileFilter_44_out" PROPERTIES DEPENDS "TOPP_FileFilter_44")

#------------------------------------------------------------------------------
# FileInfo tests
add_test("TOPP_FileInfo_1" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_1_input.dta -in_type dta -no_progress -out FileInfo_1.tmp)
add_test("TOPP_FileInfo_1_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_1.tmp -in2 ${DATA_DIR_TOPP}/FileInfo_1_output.txt )
set_tests_properties("TOPP_FileInfo_1_out1" PROPERTIES DEPENDS "TOPP_FileInfo_1")
add_test("TOPP_FileInfo_2" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_2_input.dta2d -no_progress -out FileInfo_2.tmp)
add_test("TOPP_FileInfo_2_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_2.tmp -in2 ${DATA_DIR_TOPP}/FileInfo_2_output.txt )
set_tests_properties("TOPP_FileInfo_2_out1" PROPERTIES DEPENDS "TOPP_FileInfo_2")
add_test("TOPP_FileInfo_3" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_3_input.featureXML -m -s -p -no_progress -out FileInfo_3.tmp)
add_test("TOPP_FileInfo_3_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_3.tmp -in2 ${DATA_DIR_TOPP}/FileInfo_3_output.txt )
set_tests_properties("TOPP_FileInfo_3_out1" PROPERTIES DEPENDS "TOPP_FileInfo_3")
add_test("TOPP_FileInfo_4" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_4_input.mzXML -m -no_progress -out FileInfo_4.tmp)
add_test("TOPP_FileInfo_4_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_4.tmp -in2 ${DATA_DIR_TOPP}/FileInfo_4_output.txt )
set_tests_properties("TOPP_FileInfo_4_out1" PROPERTIES DEPENDS "TOPP_FileInfo_4")
add_test("TOPP_FileInfo_5" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_5_input.mzDat -in_type mzData -m -s -no_progress -out FileInfo_5.tmp)
add_test("TOPP_FileInfo_5_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_5.tmp -in2 ${DATA_DIR_TOPP}/FileInfo_5_output.txt )
set_tests_properties("TOPP_FileInfo_5_out1" PROPERTIES DEPENDS "TOPP_FileInfo_5")
add_test("TOPP_FileInfo_6" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_6_input.mzData -d -s -no_progress -out FileInfo_6.tmp)
add_test("TOPP_FileInfo_6_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_6.tmp -in2 ${DATA_DIR_TOPP}/FileInfo_6_output.txt )
set_tests_properties("TOPP_FileInfo_6_out1" PROPERTIES DEPENDS "TOPP_FileInfo_6")
add_test("TOPP_FileInfo_7" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_7_input.consensusXML -s -m -p -no_progress -out FileInfo_7.tmp)
add_test("TOPP_FileInfo_7_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_7.tmp -in2 ${DATA_DIR_TOPP}/FileInfo_7_output.txt )
set_tests_properties("TOPP_FileInfo_7_out1" PROPERTIES DEPENDS "TOPP_FileInfo_7")
add_test("TOPP_FileInfo_9" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_9_input.mzML -m -p -s -no_progress -out FileInfo_9.tmp)
add_test("TOPP_FileInfo_9_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_9.tmp -in2 ${DATA_DIR_TOPP}/FileInfo_9_output.txt )
set_tests_properties("TOPP_FileInfo_9_out1" PROPERTIES DEPENDS "TOPP_FileInfo_9")
add_test("TOPP_FileInfo_10" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_10_input.idXML -no_progress -out FileInfo_10.tmp)
add_test("TOPP_FileInfo_10_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_10.tmp -in2 ${DATA_DIR_TOPP}/FileInfo_10_output.txt )
set_tests_properties("TOPP_FileInfo_10_out1" PROPERTIES DEPENDS "TOPP_FileInfo_10")
add_test("TOPP_FileInfo_11" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_11_input.mzML -i -no_progress)
set_tests_properties("TOPP_FileInfo_11" PROPERTIES WILL_FAIL 1) ## the mzML has no index
add_test("TOPP_FileInfo_12" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_12_input.mzML -i  -no_progress)
set_tests_properties("TOPP_FileInfo_12" PROPERTIES WILL_FAIL 0) ## this mzML has an index, should succeed
add_test("TOPP_FileInfo_13" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_13_input.consensusXML -no_progress) # empty file should not crash
add_test("TOPP_FileInfo_14" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_14_input.mzid -v -no_progress -out FileInfo_14.tmp)
add_test("TOPP_FileInfo_14_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_14.tmp -in2 ${DATA_DIR_TOPP}/FileInfo_14_output.txt )
set_tests_properties("TOPP_FileInfo_14_out1" PROPERTIES DEPENDS "TOPP_FileInfo_14")
add_test("TOPP_FileInfo_15" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_15_input.mzid -v -no_progress -out FileInfo_15.tmp)
add_test("TOPP_FileInfo_15_out1" ${DIFF} -whitelist "File name" "line 327 column 183" -in1 FileInfo_15.tmp -in2 ${DATA_DIR_TOPP}/FileInfo_15_output.txt )
set_tests_properties("TOPP_FileInfo_15_out1" PROPERTIES DEPENDS "TOPP_FileInfo_15")
add_test("TOPP_FileInfo_16" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_16_input.trafoXML -out FileInfo_16_output.tmp)
add_test("TOPP_FileInfo_16_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_16_output.tmp -in2 ${DATA_DIR_TOPP}/FileInfo_16_output.txt )
set_tests_properties("TOPP_FileInfo_16_out1" PROPERTIES DEPENDS "TOPP_FileInfo_16")

#------------------------------------------------------------------------------
# FileMerger tests
add_test("TOPP_FileMerger_1" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_1_input1.dta2d ${DATA_DIR_TOPP}/FileMerger_1_input2.dta ${DATA_DIR_TOPP}/FileMerger_1_input3.dta2d -out FileMerger_1.tmp -raw:rt_custom 1 2 3)
add_test("TOPP_FileMerger_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} "id=" "location=" "href=" -in1 FileMerger_1.tmp -in2 ${DATA_DIR_TOPP}/FileMerger_1_output.mzML )
set_tests_properties("TOPP_FileMerger_1_out1" PROPERTIES DEPENDS "TOPP_FileMerger_1")
add_test("TOPP_FileMerger_2" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_2_input1.dta ${DATA_DIR_TOPP}/FileMerger_2_input2.dta -in_type dta -out FileMerger_2.tmp -raw:rt_custom 5 10)
add_test("TOPP_FileMerger_2_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} "id=" "location=" "href=" -in1 FileMerger_2.tmp -in2 ${DATA_DIR_TOPP}/FileMerger_2_output.mzML )
set_tests_properties("TOPP_FileMerger_2_out1" PROPERTIES DEPENDS "TOPP_FileMerger_2")
add_test("TOPP_FileMerger_3" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_3_input1.dta ${DATA_DIR_TOPP}/FileMerger_3_input2.dta -out FileMerger_3.tmp -raw:rt_auto)
add_test("TOPP_FileMerger_3_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} "id=" "location=" "href=" -in1 FileMerger_3.tmp -in2 ${DATA_DIR_TOPP}/FileMerger_3_output.mzML )
set_tests_properties("TOPP_FileMerger_3_out1" PROPERTIES DEPENDS "TOPP_FileMerger_3")
add_test("TOPP_FileMerger_4" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_4_input1.dta2d ${DATA_DIR_TOPP}/FileMerger_4_input2.dta2d -out FileMerger_4.tmp)
add_test("TOPP_FileMerger_4_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} "id=" "location=" "href=" -in1 FileMerger_4.tmp -in2 ${DATA_DIR_TOPP}/FileMerger_4_output.mzML )
set_tests_properties("TOPP_FileMerger_4_out1" PROPERTIES DEPENDS "TOPP_FileMerger_4")
add_test("TOPP_FileMerger_5" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_5_input_rt1023.331.dta ${DATA_DIR_TOPP}/FileMerger_5_input_rt2044.334.dta ${DATA_DIR_TOPP}/FileMerger_5_input_rt889.32.dta -raw:rt_filename -raw:ms_level 2 -out FileMerger_5.tmp)
add_test("TOPP_FileMerger_5_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} "id=" "location=" "href=" -in1 FileMerger_5.tmp -in2 ${DATA_DIR_TOPP}/FileMerger_5_output.mzML )
set_tests_properties("TOPP_FileMerger_5_out1" PROPERTIES DEPENDS "TOPP_FileMerger_5")
add_test("TOPP_FileMerger_6" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_6_input1.mzML ${DATA_DIR_TOPP}/FileMerger_6_input2.mzML -out FileMerger_6.tmp)
add_test("TOPP_FileMerger_6_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} "id=" "location=" "href=" -in1 FileMerger_6.tmp -in2 ${DATA_DIR_TOPP}/FileMerger_6_output.mzML )
set_tests_properties("TOPP_FileMerger_6_out1" PROPERTIES DEPENDS "TOPP_FileMerger_6")
add_test("TOPP_FileMerger_7" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_7_input1.featureXML ${DATA_DIR_TOPP}/FileMerger_7_input2.featureXML -out FileMerger_7.tmp)
add_test("TOPP_FileMerger_7_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} "id=" "location=" "href=" -in1 FileMerger_7.tmp -in2 ${DATA_DIR_TOPP}/FileMerger_7_output.featureXML )
set_tests_properties("TOPP_FileMerger_7_out1" PROPERTIES DEPENDS "TOPP_FileMerger_7")
add_test("TOPP_FileMerger_8" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_8_input1.consensusXML ${DATA_DIR_TOPP}/FileMerger_8_input2.consensusXML -out FileMerger_8.tmp)
add_test("TOPP_FileMerger_8_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} "id=" "location=" "href=" -in1 FileMerger_8.tmp -in2 ${DATA_DIR_TOPP}/FileMerger_8_output.consensusXML )
set_tests_properties("TOPP_FileMerger_8_out1" PROPERTIES DEPENDS "TOPP_FileMerger_8")
add_test("TOPP_FileMerger_9" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_9_input1.traML ${DATA_DIR_TOPP}/FileMerger_9_input2.traML -out FileMerger_9.tmp)
add_test("TOPP_FileMerger_9_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} "id=" "location=" "href=" -in1 FileMerger_9.tmp -in2 ${DATA_DIR_TOPP}/FileMerger_9_output.traML )
set_tests_properties("TOPP_FileMerger_9_out1" PROPERTIES DEPENDS "TOPP_FileMerger_9")
# "rt_concat" options:
add_test("TOPP_FileMerger_10" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_6_input2.mzML ${DATA_DIR_TOPP}/FileMerger_6_input2.mzML -out FileMerger_10_output.tmp -rt_concat:gap 10.0 -rt_concat:trafo_out FileMerger_10_trafo1.tmp FileMerger_10_trafo2.tmp)
add_test("TOPP_FileMerger_10_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} "id=" "location=" "href=" "offset" "indexListOffset" -in1 FileMerger_10_output.tmp -in2 ${DATA_DIR_TOPP}/FileMerger_10_output.mzML )
set_tests_properties("TOPP_FileMerger_10_out1" PROPERTIES DEPENDS "TOPP_FileMerger_10")
add_test("TOPP_FileMerger_10_out2" ${DIFF} -in1 FileMerger_10_trafo1.tmp -in2 ${DATA_DIR_TOPP}/FileMerger_10_trafo1.trafoXML)
set_tests_properties("TOPP_FileMerger_10_out2" PROPERTIES DEPENDS "TOPP_FileMerger_10")
add_test("TOPP_FileMerger_10_out3" ${DIFF} -in1 FileMerger_10_trafo2.tmp -in2 ${DATA_DIR_TOPP}/FileMerger_10_trafo2.trafoXML)
set_tests_properties("TOPP_FileMerger_10_out3" PROPERTIES DEPENDS "TOPP_FileMerger_10")

#------------------------------------------------------------------------------
# HighResPrecursorMassCorrector tests
add_test("TOPP_HighResPrecursorMassCorrector_1"
${TOPP_BIN_PATH}/HighResPrecursorMassCorrector -test -in ${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_2860_1103_3.mzML -feature:in ${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_2860_1103_3.featureXML -out HighResPrecursorMassCorrector_2860_1103_3_out.tmp)
add_test("TOPP_HighResPrecursorMassCorrector_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 HighResPrecursorMassCorrector_2860_1103_3_out.tmp -in2 ${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_2860_1103_3_out.mzML )
set_tests_properties("TOPP_HighResPrecursorMassCorrector_1_out1" PROPERTIES DEPENDS "TOPP_HighResPrecursorMassCorrector_1")
# precursor is on the 4th isotopic trace (with 0 beeing the mono-isotopic one)
add_test("TOPP_HighResPrecursorMassCorrector_2"
	${TOPP_BIN_PATH}/HighResPrecursorMassCorrector -test -in
	${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_1035_1178_4.mzML -feature:in
	${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_1035_1178_4.featureXML
	-feature:max_trace 4
	-feature:mz_tolerance 10
	-out HighResPrecursorMassCorrector_1035_1178_4_out.tmp)
add_test("TOPP_HighResPrecursorMassCorrector_2_out1" ${DIFF} -whitelist ${INDEX_WHITELIST}  -in1 HighResPrecursorMassCorrector_1035_1178_4_out.tmp -in2 ${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_1035_1178_4_out.mzML )
set_tests_properties("TOPP_HighResPrecursorMassCorrector_2_out1" PROPERTIES DEPENDS "TOPP_HighResPrecursorMassCorrector_2")
add_test("TOPP_HighResPrecursorMassCorrector_3" ${TOPP_BIN_PATH}/HighResPrecursorMassCorrector -test -in ${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_2538_1091_2.mzML -feature:in ${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_2538_1091_2.featureXML -out HighResPrecursorMassCorrector_2538_1091_2_out.tmp)
add_test("TOPP_HighResPrecursorMassCorrector_3_out1" ${DIFF} -whitelist ${INDEX_WHITELIST}  -in1 HighResPrecursorMassCorrector_2538_1091_2_out.tmp -in2 ${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_2538_1091_2_out.mzML )
set_tests_properties("TOPP_HighResPrecursorMassCorrector_3_out1" PROPERTIES DEPENDS "TOPP_HighResPrecursorMassCorrector_3")
add_test("TOPP_HighResPrecursorMassCorrector_4" ${TOPP_BIN_PATH}/HighResPrecursorMassCorrector -test -in ${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_2810_1091_3.mzML -feature:in ${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_2810_1091_3.featureXML -out HighResPrecursorMassCorrector_2810_1091_3_out.tmp)
add_test("TOPP_HighResPrecursorMassCorrector_4_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 HighResPrecursorMassCorrector_2810_1091_3_out.tmp -in2 ${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_2810_1091_3_out.mzML )
set_tests_properties("TOPP_HighResPrecursorMassCorrector_4_out1" PROPERTIES DEPENDS "TOPP_HighResPrecursorMassCorrector_4")
add_test("TOPP_HighResPrecursorMassCorrector_5" ${TOPP_BIN_PATH}/HighResPrecursorMassCorrector -test -in ${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_3070_1191_3.mzML -feature:in ${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_3070_1191_3.featureXML -out HighResPrecursorMassCorrector_3070_1191_3_out.tmp)
add_test("TOPP_HighResPrecursorMassCorrector_5_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 HighResPrecursorMassCorrector_3070_1191_3_out.tmp -in2 ${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_3070_1191_3_out.mzML )
set_tests_properties("TOPP_HighResPrecursorMassCorrector_5_out1" PROPERTIES DEPENDS "TOPP_HighResPrecursorMassCorrector_5")
add_test("TOPP_HighResPrecursorMassCorrector_6" ${TOPP_BIN_PATH}/HighResPrecursorMassCorrector -test -in ${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_6.mzML -highest_intensity_peak:mz_tolerance 0.2 -out HighResPrecursorMassCorrector_6_out.tmp)
add_test("TOPP_HighResPrecursorMassCorrector_6_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 HighResPrecursorMassCorrector_6_out.tmp -in2 ${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_6_out.mzML )
set_tests_properties("TOPP_HighResPrecursorMassCorrector_6_out1" PROPERTIES DEPENDS "TOPP_HighResPrecursorMassCorrector_6")

#------------------------------------------------------------------------------
# IDRTCalibration tests
add_test("TOPP_IDRTCalibration_1" ${TOPP_BIN_PATH}/IDRTCalibration -test -in ${DATA_DIR_TOPP}/IDRTCalibration_1_input.idXML -out IDRTCalibration_1_output.tmp -calibrant_1_input 10 -calibrant_2_input 90)
add_test("TOPP_IDRTCalibration_1_out1" ${DIFF} -in1 IDRTCalibration_1_output.tmp -in2 ${DATA_DIR_TOPP}/IDRTCalibration_1_output.idXML )
set_tests_properties("TOPP_IDRTCalibration_1_out1" PROPERTIES DEPENDS "TOPP_IDRTCalibration_1")
# test empty IDs
add_test("TOPP_IDRTCalibration_2" ${TOPP_BIN_PATH}/IDRTCalibration -test -in ${DATA_DIR_TOPP}/degenerated_empty.idXML -out IDRTCalibration_2_output.tmp -calibrant_1_input 10 -calibrant_2_input 90)
add_test("TOPP_IDRTCalibration_2_out1" ${DIFF} -in1 IDRTCalibration_2_output.tmp -in2 ${DATA_DIR_TOPP}/IDRTCalibration_2_output.idXML )
set_tests_properties("TOPP_IDRTCalibration_2_out1" PROPERTIES DEPENDS "TOPP_IDRTCalibration_2")

#------------------------------------------------------------------------------
# ITRAQAnalyzer tests
add_test("TOPP_ITRAQAnalyzer_1" ${TOPP_BIN_PATH}/ITRAQAnalyzer -test -in ${DATA_DIR_TOPP}/ITRAQAnalyzer_input_1.mzML -ini ${DATA_DIR_TOPP}/ITRAQAnalyzer.ini -out ITRAQAnalyzer_output_1.tmp)
add_test("TOPP_ITRAQAnalyzer_1_out1" ${DIFF} -whitelist "id=" "<map" "?xml-stylesheet" -in1 ITRAQAnalyzer_output_1.tmp -in2 ${DATA_DIR_TOPP}/ITRAQAnalyzer_output_1.consensusXML )
set_tests_properties("TOPP_ITRAQAnalyzer_1_out1" PROPERTIES DEPENDS "TOPP_ITRAQAnalyzer_1")
# test empty IDs
add_test("TOPP_ITRAQAnalyzer_2" ${TOPP_BIN_PATH}/ITRAQAnalyzer -test -in ${DATA_DIR_TOPP}/ITRAQAnalyzer_input_2.mzML -ini ${DATA_DIR_TOPP}/ITRAQAnalyzer.ini -out ITRAQAnalyzer_output_2.tmp)
add_test("TOPP_ITRAQAnalyzer_2_out1" ${DIFF} -whitelist "id=" "<map" "?xml-stylesheet" -in1 ITRAQAnalyzer_output_2.tmp -in2 ${DATA_DIR_TOPP}/ITRAQAnalyzer_output_2.consensusXML )
set_tests_properties("TOPP_ITRAQAnalyzer_2_out1" PROPERTIES DEPENDS "TOPP_ITRAQAnalyzer_2")

#------------------------------------------------------------------------------
# IsobaricAnalyzer tests
add_test("TOPP_IsobaricAnalyzer_1" ${TOPP_BIN_PATH}/IsobaricAnalyzer -test -in ${DATA_DIR_TOPP}/IsobaricAnalyzer_input_1.mzML -ini ${DATA_DIR_TOPP}/IsobaricAnalyzer.ini -out IsobaricAnalyzer_output_1.tmp)
add_test("TOPP_IsobaricAnalyzer_1_out1" ${DIFF} -whitelist "id=" "<map" "?xml-stylesheet" -in1 IsobaricAnalyzer_output_1.tmp -in2 ${DATA_DIR_TOPP}/IsobaricAnalyzer_output_1.consensusXML )
set_tests_properties("TOPP_IsobaricAnalyzer_1_out1" PROPERTIES DEPENDS "TOPP_IsobaricAnalyzer_1")

add_test("TOPP_IsobaricAnalyzer_TMTTenPlexMethod_1" ${TOPP_BIN_PATH}/IsobaricAnalyzer -test -in ${DATA_DIR_TOPP}/TMTTenPlexMethod_test.mzML -ini ${DATA_DIR_TOPP}/TMTTenPlexMethod_test.ini -out TMTTenPlexMethod_ouput.tmp)
add_test("TOPP_IsobaricAnalyzer_TMTTenPlexMethod_1_out1" ${DIFF} -whitelist "id=" "<map" "?xml-stylesheet" -in1 TMTTenPlexMethod_ouput.tmp -in2 ${DATA_DIR_TOPP}/TMTTenPlexMethod_test.consensusXML )
set_tests_properties("TOPP_IsobaricAnalyzer_TMTTenPlexMethod_1_out1" PROPERTIES DEPENDS "TOPP_IsobaricAnalyzer_TMTTenPlexMethod_1")

#------------------------------------------------------------------------------
# IDConflictResolver tests
add_test("TOPP_IDConflictResolver_1" ${TOPP_BIN_PATH}/IDConflictResolver -test -in ${DATA_DIR_TOPP}/IDConflictResolver_1_input.featureXML -out IDConflictResolver_1_output.tmp)
add_test("TOPP_IDConflictResolver_1_out1" ${DIFF} -in1 IDConflictResolver_1_output.tmp -in2 ${DATA_DIR_TOPP}/IDConflictResolver_1_output.featureXML )
set_tests_properties("TOPP_IDConflictResolver_1_out1" PROPERTIES DEPENDS "TOPP_IDConflictResolver_1")
add_test("TOPP_IDConflictResolver_2" ${TOPP_BIN_PATH}/IDConflictResolver -test -in ${DATA_DIR_TOPP}/IDConflictResolver_2_input.consensusXML -out IDConflictResolver_2_output.tmp)
add_test("TOPP_IDConflictResolver_2_out1" ${DIFF} -in1 IDConflictResolver_2_output.tmp -in2 ${DATA_DIR_TOPP}/IDConflictResolver_2_output.consensusXML )
set_tests_properties("TOPP_IDConflictResolver_2_out1" PROPERTIES DEPENDS "TOPP_IDConflictResolver_2")

#------------------------------------------------------------------------------
# IDFileConverter tests
# Mascot XML to idXML:
add_test("TOPP_IDFileConverter_1" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_1_input1.mascotXML -mz_file ${DATA_DIR_TOPP}/IDFileConverter_1_input2.mzML -out IDFileConverter_1_output.tmp -out_type idXML)
add_test("TOPP_IDFileConverter_1_out1" ${DIFF} -in1 IDFileConverter_1_output.tmp -in2 ${DATA_DIR_TOPP}/IDFileConverter_1_output.idXML )
set_tests_properties("TOPP_IDFileConverter_1_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_1")
# pepXML to idXML (using precursor RTs):
add_test("TOPP_IDFileConverter_2" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/PepXMLFile_test.pepxml -mz_file ${DATA_DIR_TOPP}/PepXMLFile_test.mzML -mz_name PepXMLFile_test -out IDFileConverter_2_output.tmp -out_type idXML)
add_test("TOPP_IDFileConverter_2_out1" ${DIFF} -in1 IDFileConverter_2_output.tmp -in2 ${DATA_DIR_TOPP}/IDFileConverter_2_output.idXML )
set_tests_properties("TOPP_IDFileConverter_2_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_2")
# protXML to idXML:
add_test("TOPP_IDFileConverter_3" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_3_input.protXML -out IDFileConverter_3_output.tmp -out_type idXML)
add_test("TOPP_IDFileConverter_3_out1" ${DIFF} -in1 IDFileConverter_3_output.tmp -in2 ${DATA_DIR_TOPP}/IDFileConverter_3_output.idXML )
set_tests_properties("TOPP_IDFileConverter_3_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_3")
# pepXML to idXML (using MS2 RTs):
add_test("TOPP_IDFileConverter_4" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/PepXMLFile_test.pepxml -mz_file ${DATA_DIR_TOPP}/PepXMLFile_test.mzML -mz_name PepXMLFile_test -out IDFileConverter_4_output.tmp -out_type idXML)
add_test("TOPP_IDFileConverter_4_out1" ${DIFF} -in1 IDFileConverter_4_output.tmp -in2 ${DATA_DIR_TOPP}/IDFileConverter_4_output.idXML )
set_tests_properties("TOPP_IDFileConverter_4_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_4")
# Mascot XML to idXML (with user-defined regular expression):
add_test("TOPP_IDFileConverter_5" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_1_input1.mascotXML -mz_file ${DATA_DIR_TOPP}/IDFileConverter_1_input2.mzML -out IDFileConverter_5_output.tmp -ini ${DATA_DIR_TOPP}/IDFileConverter_5_parameters.ini)
add_test("TOPP_IDFileConverter_5_out1" ${DIFF} -in1 IDFileConverter_5_output.tmp -in2 ${DATA_DIR_TOPP}/IDFileConverter_5_output.idXML )
set_tests_properties("TOPP_IDFileConverter_5_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_5")
# pepXML (exported from Mascot) to idXML:
add_test("TOPP_IDFileConverter_6" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_6_input1.pepXML -mz_file ${DATA_DIR_TOPP}/IDFileConverter_1_input2.mzML -mz_name F025589.dat.mzML -out IDFileConverter_6_output.tmp -out_type idXML)
add_test("TOPP_IDFileConverter_6_out1" ${DIFF} -in1 IDFileConverter_6_output.tmp -in2 ${DATA_DIR_TOPP}/IDFileConverter_6_output.idXML )
set_tests_properties("TOPP_IDFileConverter_6_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_6")
# X! Tandem XML to idXML:
add_test("TOPP_IDFileConverter_7" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_7_input1.xml -out IDFileConverter_7_output.tmp -out_type idXML)
add_test("TOPP_IDFileConverter_7_out1" ${DIFF} -in1 IDFileConverter_7_output.tmp -in2 ${DATA_DIR_TOPP}/IDFileConverter_7_output.idXML -whitelist "IdentificationRun date" )
set_tests_properties("TOPP_IDFileConverter_7_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_7")
# mzIdentML to idXML (input file is "MzIdentMLFile_whole.mzid" from the class test):
add_test("TOPP_IDFileConverter_8" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_8_input.mzid -out IDFileConverter_8_output.tmp -out_type idXML)
add_test("TOPP_IDFileConverter_8_out1" ${DIFF} -in1 IDFileConverter_8_output.tmp -in2 ${DATA_DIR_TOPP}/IDFileConverter_8_output.idXML -whitelist "IdentificationRun date" )
set_tests_properties("TOPP_IDFileConverter_8_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_8")
#  idXML to mzIdentML (input file is "IdXMLFile_whole.idXML" from the class test):
add_test("TOPP_IDFileConverter_9" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_9_input.idXML -out IDFileConverter_9_output.tmp -out_type mzid)
add_test("TOPP_IDFileConverter_9_out1" ${DIFF} -in1 IDFileConverter_9_output.tmp -in2 ${DATA_DIR_TOPP}/IDFileConverter_9_output.mzid -whitelist "creationDate" "id" "spectraData_ref" "searchDatabase_ref")
set_tests_properties("TOPP_IDFileConverter_9_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_9")
# pepXML with very small delta mass as modification (most probably an artefact and removed during parsing)
add_test("TOPP_IDFileConverter_10" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_10_input.pepXML -out IDFileConverter_10_output.tmp -out_type idXML)
add_test("TOPP_IDFileConverter_10_out1" ${DIFF} -in1 IDFileConverter_10_output.tmp -in2 ${DATA_DIR_TOPP}/IDFileConverter_10_output.idXML )
set_tests_properties("TOPP_IDFileConverter_10_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_10")
# tsv to idXML
add_test("TOPP_IDFileConverter_11" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_11_input.tsv -out IDFileConverter_11_output.tmp -out_type idXML)
add_test("TOPP_IDFileConverter_11_out1" ${DIFF} -in1 IDFileConverter_11_output.tmp -in2 ${DATA_DIR_TOPP}/IDFileConverter_11_output.idXML)
set_tests_properties("TOPP_IDFileConverter_11_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_11")
# Percolator .psms to idXML (score: q-value)
add_test("TOPP_IDFileConverter_12" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_12_input.psms -out IDFileConverter_12_output.tmp -out_type idXML -score_type qvalue)
add_test("TOPP_IDFileConverter_12_out1" ${DIFF} -in1 IDFileConverter_12_output.tmp -in2 ${DATA_DIR_TOPP}/IDFileConverter_12_output.idXML -whitelist "IdentificationRun date")
set_tests_properties("TOPP_IDFileConverter_12_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_12")
# Percolator .psms to idXML (score: PEP)
add_test("TOPP_IDFileConverter_13" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_12_input.psms -out IDFileConverter_13_output.tmp -out_type idXML -score_type PEP)
add_test("TOPP_IDFileConverter_13_out1" ${DIFF} -in1 IDFileConverter_13_output.tmp -in2 ${DATA_DIR_TOPP}/IDFileConverter_13_output.idXML -whitelist "IdentificationRun date")
set_tests_properties("TOPP_IDFileConverter_13_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_13")
# Percolator .psms to idXML (score: raw score)
add_test("TOPP_IDFileConverter_14" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_12_input.psms -out IDFileConverter_14_output.tmp -out_type idXML -score_type score)
add_test("TOPP_IDFileConverter_14_out1" ${DIFF} -in1 IDFileConverter_14_output.tmp -in2 ${DATA_DIR_TOPP}/IDFileConverter_14_output.idXML -whitelist "IdentificationRun date")
set_tests_properties("TOPP_IDFileConverter_14_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_14")
# mzIdentML to idXML with RT look-up in mzML:
add_test("TOPP_IDFileConverter_15" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/THIRDPARTY/MSGFPlusAdapter_1_out.mzid -out IDFileConverter_15_output.tmp -out_type idXML -mz_file ${DATA_DIR_TOPP}/THIRDPARTY/spectra.mzML)
add_test("TOPP_IDFileConverter_15_out1" ${DIFF} -in1 IDFileConverter_15_output.tmp -in2 ${DATA_DIR_TOPP}/IDFileConverter_15_output.idXML -whitelist "IdentificationRun date" "db=" "spectra_data")
set_tests_properties("TOPP_IDFileConverter_15_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_15")
# pepXML (from iProphet) to idXML:
add_test("TOPP_IDFileConverter_16" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_16_input.pepXML -out IDFileConverter_16_output.tmp -out_type idXML)
add_test("TOPP_IDFileConverter_16_out1" ${DIFF} -in1 IDFileConverter_16_output.tmp -in2 ${DATA_DIR_TOPP}/IDFileConverter_16_output.idXML)
set_tests_properties("TOPP_IDFileConverter_16_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_16")
# idXML to pepXML (like in LuciphorAdapter)
# Tests  E-Value, hyperscore
add_test("TOPP_IDFileConverter_17" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_17_input.idXML -out IDFileConverter_17_output.tmp -out_type pepXML)
add_test("TOPP_IDFileConverter_17_out1" ${DIFF} -in1 IDFileConverter_17_output.tmp -in2 ${DATA_DIR_TOPP}/IDFileConverter_17_output.pepXML)
set_tests_properties("TOPP_IDFileConverter_17_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_17")
# Tests nextscore
add_test("TOPP_IDFileConverter_18" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_18_input.idXML -out IDFileConverter_18_output.tmp -out_type pepXML)
add_test("TOPP_IDFileConverter_18_out1" ${DIFF} -in1 IDFileConverter_17_output.tmp -in2 ${DATA_DIR_TOPP}/IDFileConverter_17_output.pepXML)
set_tests_properties("TOPP_IDFileConverter_18_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_18")
# Tests Percolator_score
add_test("TOPP_IDFileConverter_19" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_19_input.idXML -out IDFileConverter_19_output.tmp -out_type pepXML)
add_test("TOPP_IDFileConverter_19_out1" ${DIFF} -in1 IDFileConverter_17_output.tmp -in2 ${DATA_DIR_TOPP}/IDFileConverter_17_output.pepXML)
set_tests_properties("TOPP_IDFileConverter_19_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_19")
# Tests hyperscore
add_test("TOPP_IDFileConverter_20" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_20_input.idXML -out IDFileConverter_20_output.tmp -out_type pepXML)
add_test("TOPP_IDFileConverter_20_out1" ${DIFF} -in1 IDFileConverter_17_output.tmp -in2 ${DATA_DIR_TOPP}/IDFileConverter_17_output.pepXML)
set_tests_properties("TOPP_IDFileConverter_20_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_20")

#------------------------------------------------------------------------------
# IDFilter tests
add_test("TOPP_IDFilter_1" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_1_input.idXML -out IDFilter_1_output.tmp -whitelist:proteins ${DATA_DIR_TOPP}/IDFilter_1_input.fas)
add_test("TOPP_IDFilter_1_out1" ${DIFF} -in1 IDFilter_1_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_1_output.idXML )
set_tests_properties("TOPP_IDFilter_1_out1" PROPERTIES DEPENDS "TOPP_IDFilter_1")
add_test("TOPP_IDFilter_2" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_2_input.idXML -out IDFilter_2_output.tmp -thresh:pep 1 -thresh:prot 1)
add_test("TOPP_IDFilter_2_out1" ${DIFF} -in1 IDFilter_2_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_2_output.idXML )
set_tests_properties("TOPP_IDFilter_2_out1" PROPERTIES DEPENDS "TOPP_IDFilter_2")
add_test("TOPP_IDFilter_3" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_3_input.idXML -out IDFilter_3_output.tmp -blacklist:peptides ${DATA_DIR_TOPP}/IDFilter_3_2_input.idXML)
add_test("TOPP_IDFilter_3_out1" ${DIFF} -in1 IDFilter_3_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_3_output.idXML )
set_tests_properties("TOPP_IDFilter_3_out1" PROPERTIES DEPENDS "TOPP_IDFilter_3")
add_test("TOPP_IDFilter_4" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_4_input.idXML -out IDFilter_4_output.tmp -rt:p_value 0.08)
add_test("TOPP_IDFilter_4_out1" ${DIFF} -in1 IDFilter_4_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_4_output.idXML )
set_tests_properties("TOPP_IDFilter_4_out1" PROPERTIES DEPENDS "TOPP_IDFilter_4")

add_test("TOPP_IDFilter_5" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_5_input.idXML -out IDFilter_5_output.tmp -score:pep 32 -score:prot 25)
add_test("TOPP_IDFilter_5_out1" ${DIFF} -in1 IDFilter_5_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_5_output.idXML )
set_tests_properties("TOPP_IDFilter_5_out1" PROPERTIES DEPENDS "TOPP_IDFilter_5")
add_test("TOPP_IDFilter_5a" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_5_input.idXML -out IDFilter_5a_output.tmp -score:pep 32 )
add_test("TOPP_IDFilter_5a_out1" ${DIFF} -in1 IDFilter_5a_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_5_output.idXML )
set_tests_properties("TOPP_IDFilter_5a_out1" PROPERTIES DEPENDS "TOPP_IDFilter_5a")

add_test("TOPP_IDFilter_5b" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_5_input.idXML -out IDFilter_5b_output.tmp -score:prot 25)
add_test("TOPP_IDFilter_5b_out1" ${DIFF} -in1 IDFilter_5b_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_5b_output.idXML )
set_tests_properties("TOPP_IDFilter_5b_out1" PROPERTIES DEPENDS "TOPP_IDFilter_5b")

add_test("TOPP_IDFilter_5c" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_5_input.idXML -out IDFilter_5c_output.tmp -score:prot 25 -delete_unreferenced_peptide_hits)
add_test("TOPP_IDFilter_5c_out1" ${DIFF} -in1 IDFilter_5c_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_5c_output.idXML )
set_tests_properties("TOPP_IDFilter_5c_out1" PROPERTIES DEPENDS "TOPP_IDFilter_5c")

add_test("TOPP_IDFilter_6" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_6_input.idXML -out IDFilter_6_output.tmp -best:n_peptide_hits 2 -best:n_protein_hits 10)
add_test("TOPP_IDFilter_6_out1" ${DIFF} -in1 IDFilter_6_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_6_output.idXML )
set_tests_properties("TOPP_IDFilter_6_out1" PROPERTIES DEPENDS "TOPP_IDFilter_6")
add_test("TOPP_IDFilter_7" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_7_input.idXML -out IDFilter_7_output.tmp -unique)
add_test("TOPP_IDFilter_7_out1" ${DIFF} -in1 IDFilter_7_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_7_output.idXML )
set_tests_properties("TOPP_IDFilter_7_out1" PROPERTIES DEPENDS "TOPP_IDFilter_7")
# filter by precursor RT and m/z:
add_test("TOPP_IDFilter_8" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_8_input.idXML -out IDFilter_8_output.tmp -precursor:rt 200:350 -precursor:mz 999:1000)
add_test("TOPP_IDFilter_8_out1" ${DIFF} -in1 IDFilter_8_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_8_output.idXML )
set_tests_properties("TOPP_IDFilter_8_out1" PROPERTIES DEPENDS "TOPP_IDFilter_8")
# check for issue 887 regression (multiple runs)
add_test("TOPP_IDFilter_9" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_9_input.idXML -out IDFilter_9_output.tmp -score:pep 0.05)
add_test("TOPP_IDFilter_9_out1" ${DIFF} -in1 IDFilter_9_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_9_output.idXML )
set_tests_properties("TOPP_IDFilter_9_out1" PROPERTIES DEPENDS "TOPP_IDFilter_9")
# update protein groups ("IDFilter_10_input.idXML" is same as "IDFileConverter_3_output.idXML", i.e. converted protXML):
add_test("TOPP_IDFilter_10" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_10_input.idXML -out IDFilter_10_output.tmp -score:prot 0.3 -delete_unreferenced_peptide_hits)
add_test("TOPP_IDFilter_10_out1" ${DIFF} -in1 IDFilter_10_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_10_output.idXML )
set_tests_properties("TOPP_IDFilter_10_out1" PROPERTIES DEPENDS "TOPP_IDFilter_10")
# remove decoys and update protein groups:
add_test("TOPP_IDFilter_11" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_11_input.idXML -out IDFilter_11_output.tmp -remove_decoys -delete_unreferenced_peptide_hits)
add_test("TOPP_IDFilter_11_out1" ${DIFF} -in1 IDFilter_11_output.tmp -in2 ${DATA_DIR_TOPP}/IDFilter_11_output.idXML )
set_tests_properties("TOPP_IDFilter_11_out1" PROPERTIES DEPENDS "TOPP_IDFilter_11")

#------------------------------------------------------------------------------
# MapAlignerPoseClustering tests
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
add_test("TOPP_MapAlignerPoseClustering_2_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 MapAlignerPoseClustering_2_output1.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerPoseClustering_2_output1.mzML )
set_tests_properties("TOPP_MapAlignerPoseClustering_2_out1" PROPERTIES DEPENDS "TOPP_MapAlignerPoseClustering_2")
add_test("TOPP_MapAlignerPoseClustering_2_out2" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 MapAlignerPoseClustering_2_output2.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerPoseClustering_2_output2.mzML )
set_tests_properties("TOPP_MapAlignerPoseClustering_2_out2" PROPERTIES DEPENDS "TOPP_MapAlignerPoseClustering_2")
add_test("TOPP_MapAlignerPoseClustering_2_out3" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 MapAlignerPoseClustering_2_output3.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerPoseClustering_2_output3.mzML )
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

#------------------------------------------------------------------------------
# MapAlignerIdentification tests:
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

#------------------------------------------------------------------------------
# MapAlignerSpectrum tests:
add_test("TOPP_MapAlignerSpectrum_1" ${TOPP_BIN_PATH}/MapAlignerSpectrum -test -ini ${DATA_DIR_TOPP}/MapAlignerSpectrum_parameters.ini -in ${DATA_DIR_TOPP}/MapAlignerSpectrum_1_input1.mzML ${DATA_DIR_TOPP}/MapAlignerSpectrum_1_input2.mzML ${DATA_DIR_TOPP}/MapAlignerSpectrum_1_input3.mzML -out MapAlignerSpectrum_1_output1.tmp MapAlignerSpectrum_1_output2.tmp MapAlignerSpectrum_1_output3.tmp)
add_test("TOPP_MapAlignerSpectrum_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 MapAlignerSpectrum_1_output1.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerSpectrum_1_output1.mzML )
set_tests_properties("TOPP_MapAlignerSpectrum_1_out1" PROPERTIES DEPENDS "TOPP_MapAlignerSpectrum_1")
add_test("TOPP_MapAlignerSpectrum_1_out2" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 MapAlignerSpectrum_1_output2.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerSpectrum_1_output2.mzML )
set_tests_properties("TOPP_MapAlignerSpectrum_1_out2" PROPERTIES DEPENDS "TOPP_MapAlignerSpectrum_1")
add_test("TOPP_MapAlignerSpectrum_1_out3" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 MapAlignerSpectrum_1_output3.tmp -in2 ${DATA_DIR_TOPP}/MapAlignerSpectrum_1_output3.mzML )
set_tests_properties("TOPP_MapAlignerSpectrum_1_out3" PROPERTIES DEPENDS "TOPP_MapAlignerSpectrum_1")

#------------------------------------------------------------------------------
# MapRTTransformer tests
# featureXML input:
add_test("TOPP_MapRTTransformer_1" ${TOPP_BIN_PATH}/MapRTTransformer -test -in ${DATA_DIR_TOPP}/MapRTTransformer_1_input.featureXML -trafo_in ${DATA_DIR_TOPP}/MapRTTransformer_trafo_linear.trafoXML -out MapRTTransformer_1_output.tmp)
add_test("TOPP_MapRTTransformer_1_out1" ${DIFF} -in1 MapRTTransformer_1_output.tmp -in2 ${DATA_DIR_TOPP}/MapRTTransformer_1_output.featureXML )
set_tests_properties("TOPP_MapRTTransformer_1_out1" PROPERTIES DEPENDS "TOPP_MapRTTransformer_1")
# mzML input:
add_test("TOPP_MapRTTransformer_2" ${TOPP_BIN_PATH}/MapRTTransformer -test -in ${DATA_DIR_TOPP}/MapRTTransformer_2_input.mzML -trafo_in ${DATA_DIR_TOPP}/MapRTTransformer_trafo_linear.trafoXML -out MapRTTransformer_2_output.tmp)
add_test("TOPP_MapRTTransformer_2_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 MapRTTransformer_2_output.tmp -in2 ${DATA_DIR_TOPP}/MapRTTransformer_2_output.mzML )
set_tests_properties("TOPP_MapRTTransformer_2_out1" PROPERTIES DEPENDS "TOPP_MapRTTransformer_2")
# "invert" option:
add_test("TOPP_MapRTTransformer_3" ${TOPP_BIN_PATH}/MapRTTransformer -test -invert -trafo_in ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_trafo2.trafoXML -trafo_out MapRTTransformer_3_output.tmp)
add_test("TOPP_MapRTTransformer_3_out1" ${DIFF} -in1 MapRTTransformer_3_output.tmp -in2 ${DATA_DIR_TOPP}/MapRTTransformer_3_trafo.trafoXML )
set_tests_properties("TOPP_MapRTTransformer_3_out1" PROPERTIES DEPENDS "TOPP_MapRTTransformer_3")
# chromatograms
add_test("TOPP_MapRTTransformer_4" ${TOPP_BIN_PATH}/MapRTTransformer -test -in ${DATA_DIR_TOPP}/MapRTTransformer_4_input.chrom.mzML -trafo_in ${DATA_DIR_TOPP}/MapRTTransformer_trafo_linear.trafoXML -out MapRTTransformer_4_output.tmp)
add_test("TOPP_MapRTTransformer_4_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 MapRTTransformer_4_output.tmp -in2 ${DATA_DIR_TOPP}/MapRTTransformer_4_output.chrom.mzML )
set_tests_properties("TOPP_MapRTTransformer_4_out1" PROPERTIES DEPENDS "TOPP_MapRTTransformer_4")
# no transformation:
add_test("TOPP_MapRTTransformer_5" ${TOPP_BIN_PATH}/MapRTTransformer -test -in ${DATA_DIR_TOPP}/MapRTTransformer_1_input.featureXML -trafo_in ${DATA_DIR_TOPP}/MapRTTransformer_trafo_none.trafoXML -out MapRTTransformer_5_output.tmp)
add_test("TOPP_MapRTTransformer_5_out1" ${DIFF} -in1 MapRTTransformer_5_output.tmp -in2 ${DATA_DIR_TOPP}/MapRTTransformer_5_output.featureXML )
set_tests_properties("TOPP_MapRTTransformer_5_out1" PROPERTIES DEPENDS "TOPP_MapRTTransformer_5")
# "store_original_rt" option:
add_test("TOPP_MapRTTransformer_6" ${TOPP_BIN_PATH}/MapRTTransformer -test -in ${DATA_DIR_TOPP}/MapRTTransformer_1_input.featureXML -trafo_in ${DATA_DIR_TOPP}/MapRTTransformer_trafo_linear.trafoXML -out MapRTTransformer_6_output.tmp -store_original_rt)
add_test("TOPP_MapRTTransformer_6_out1" ${DIFF} -in1 MapRTTransformer_6_output.tmp -in2 ${DATA_DIR_TOPP}/MapRTTransformer_6_output.featureXML )
set_tests_properties("TOPP_MapRTTransformer_6_out1" PROPERTIES DEPENDS "TOPP_MapRTTransformer_6")

#------------------------------------------------------------------------------
# MascotAdapter tests
add_test("TOPP_MascotAdapter_1" ${TOPP_BIN_PATH}/MascotAdapter -test -ini ${DATA_DIR_TOPP}/MascotAdapter_1_parameters.ini -mascot_in -in ${DATA_DIR_TOPP}/MascotAdapter_1_input.mzData)
add_test("TOPP_MascotAdapter_1_out1" ${DIFF} -in1 MascotAdapter_1_output.tmp -in2 ${DATA_DIR_TOPP}/MascotAdapter_1_output.mascot_in )
set_tests_properties("TOPP_MascotAdapter_1_out1" PROPERTIES DEPENDS "TOPP_MascotAdapter_1")
add_test("TOPP_MascotAdapter_2" ${TOPP_BIN_PATH}/MascotAdapter -test -in ${DATA_DIR_TOPP}/MascotAdapter_2_input.mascotXML -out MascotAdapter_2_output.tmp -mascot_out)
add_test("TOPP_MascotAdapter_2_out1" ${DIFF} -in1 MascotAdapter_2_output.tmp -in2 ${DATA_DIR_TOPP}/MascotAdapter_2_output.idXML )
set_tests_properties("TOPP_MascotAdapter_2_out1" PROPERTIES DEPENDS "TOPP_MascotAdapter_2")

if(NOT DISABLE_OPENSWATH)
  #------------------------------------------------------------------------------
  # MRM / SWATH tests (from OpenSWATH)
  add_test("TOPP_OpenSwathAssayGenerator_test_1"
${TOPP_BIN_PATH}/OpenSwathAssayGenerator -test -in ${DATA_DIR_TOPP}/OpenSwathAssayGenerator_input.TraML -out OpenSwathAssayGenerator_output.TraML.tmp -min_transitions 6 -max_transitions 6 -allowed_fragment_types b,y -allowed_fragment_charges 2,3 -enable_detection_specific_losses -enable_reannotation)
  add_test("TOPP_OpenSwathAssayGenerator_test_1_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathAssayGenerator_output.TraML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathAssayGenerator_output.TraML)
  set_tests_properties("TOPP_OpenSwathAssayGenerator_test_1_out1" PROPERTIES DEPENDS "TOPP_OpenSwathAssayGenerator_test_1")
  
  add_test("TOPP_OpenSwathDecoyGenerator_test_1"
${TOPP_BIN_PATH}/OpenSwathDecoyGenerator -test -in ${DATA_DIR_TOPP}/OpenSwathDecoyGenerator_input.TraML -out OpenSwathDecoyGenerator.TraML.tmp -method pseudo-reverse -remove_CNterm_mods)
  add_test("TOPP_OpenSwathDecoyGenerator_test_1_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathDecoyGenerator.TraML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathDecoyGenerator_output.TraML)
  set_tests_properties("TOPP_OpenSwathDecoyGenerator_test_1_out1" PROPERTIES DEPENDS "TOPP_OpenSwathDecoyGenerator_test_1")
  add_test("TOPP_OpenSwathDecoyGenerator_test_2"
${TOPP_BIN_PATH}/OpenSwathDecoyGenerator -test -in ${DATA_DIR_TOPP}/OpenSwathDecoyGenerator_input_2.TraML -out OpenSwathDecoyGenerator_2.TraML.tmp -append -method pseudo-reverse -remove_CNterm_mods -remove_unannotated -mz_threshold 0.8)
  add_test("TOPP_OpenSwathDecoyGenerator_test_2_out1" ${DIFF} -in1 OpenSwathDecoyGenerator_2.TraML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathDecoyGenerator_output_2.TraML)
  set_tests_properties("TOPP_OpenSwathDecoyGenerator_test_2_out1" PROPERTIES DEPENDS "TOPP_OpenSwathDecoyGenerator_test_2")
  add_test("TOPP_OpenSwathDecoyGenerator_test_3"
${TOPP_BIN_PATH}/OpenSwathDecoyGenerator -test -in ${DATA_DIR_TOPP}/OpenSwathDecoyGenerator_input_3.TraML -out OpenSwathDecoyGenerator_3.TraML.tmp -method pseudo-reverse -remove_CNterm_mods)
  add_test("TOPP_OpenSwathDecoyGenerator_test_3_out1" ${DIFF} -in1 OpenSwathDecoyGenerator_3.TraML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathDecoyGenerator_output_3.TraML)
  set_tests_properties("TOPP_OpenSwathDecoyGenerator_test_3_out1" PROPERTIES DEPENDS "TOPP_OpenSwathDecoyGenerator_test_3")

  # TargetedFileConverter : convert from and to different targeted files
  # Contains heavy and light transitions in 2 different notations
  add_test("TOPP_ConvertTSVToTraML_test_1" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ${DATA_DIR_TOPP}/ConvertTSVToTraML_1_input.tsv -out ConvertTSVToTraML_output.TraML.tmp -out_type TraML)
  add_test("TOPP_ConvertTSVToTraML_test_1_out1" ${DIFF} -in1 ConvertTSVToTraML_output.TraML.tmp -in2 ${DATA_DIR_TOPP}/ConvertTSVToTraML_output.TraML)
  set_tests_properties("TOPP_ConvertTSVToTraML_test_1_out1" PROPERTIES DEPENDS "TOPP_ConvertTSVToTraML_test_1")
  add_test("TOPP_ConvertTSVToTraML_test_2" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ${DATA_DIR_TOPP}/ConvertTSVToTraML_2_input.tsv -out ConvertTSVToTraML_2_output.TraML.tmp -algorithm:retentionTimeInterpretation minutes -out_type TraML)
  add_test("TOPP_ConvertTSVToTraML_test_2_out1" ${DIFF} -in1 ConvertTSVToTraML_2_output.TraML.tmp -in2 ${DATA_DIR_TOPP}/ConvertTSVToTraML_2_output.TraML)
  set_tests_properties("TOPP_ConvertTSVToTraML_test_2_out1" PROPERTIES DEPENDS "TOPP_ConvertTSVToTraML_test_2")
  add_test("TOPP_ConvertTSVToTraML_test_3" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ${DATA_DIR_TOPP}/ConvertTSVToTraML_2_input.tsv -out ConvertTSVToTraML_3_output.TraML.tmp -algorithm:retentionTimeInterpretation seconds -out_type TraML)
  add_test("TOPP_ConvertTSVToTraML_test_3_out1" ${DIFF} -in1 ConvertTSVToTraML_3_output.TraML.tmp -in2 ${DATA_DIR_TOPP}/ConvertTSVToTraML_3_output.TraML)
  set_tests_properties("TOPP_ConvertTSVToTraML_test_3_out1" PROPERTIES DEPENDS "TOPP_ConvertTSVToTraML_test_3")
  add_test("TOPP_ConvertTSVToTraML_test_4" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ${DATA_DIR_TOPP}/ConvertTSVToTraML_4_input.mrm -out ConvertTSVToTraML_4_output.TraML.tmp -out_type TraML)
  add_test("TOPP_ConvertTSVToTraML_test_4_out1" ${DIFF} -in1 ConvertTSVToTraML_4_output.TraML.tmp -in2 ${DATA_DIR_TOPP}/ConvertTSVToTraML_4_output.TraML)
  set_tests_properties("TOPP_ConvertTSVToTraML_test_4_out1" PROPERTIES DEPENDS "TOPP_ConvertTSVToTraML_test_4")

  # Proteomics PQP - TraML Conversion
  add_test("TOPP_TargetedFileConverter_test_1_prepare" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ${DATA_DIR_TOPP}/TargetedFileConverter_1_input.TraML -out ${DATA_DIR_TOPP}/TargetedFileConverter_1_output.pqp.tmp -out_type pqp)
  add_test("TOPP_TargetedFileConverter_test_1" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ${DATA_DIR_TOPP}/TargetedFileConverter_1_output.pqp.tmp -in_type pqp -out ${DATA_DIR_TOPP}/TargetedFileConverter_1_output.TraML.tmp -out_type TraML)
  add_test("TOPP_TargetedFileConverter_test_1_out1" ${DIFF} -in1 ${DATA_DIR_TOPP}/TargetedFileConverter_1_output.TraML.tmp -in2 ${DATA_DIR_TOPP}/TargetedFileConverter_1_output.TraML)
  set_tests_properties("TOPP_TargetedFileConverter_test_1" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_1_prepare")
  set_tests_properties("TOPP_TargetedFileConverter_test_1_out1" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_1")
  add_test("TOPP_TargetedFileConverter_test_2" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ${DATA_DIR_TOPP}/TargetedFileConverter_1_output.pqp.tmp -in_type pqp -out ${DATA_DIR_TOPP}/TargetedFileConverter_2_output.TraML.tmp -out_type TraML -legacy_traml_id)
  add_test("TOPP_TargetedFileConverter_test_2_out1" ${DIFF} -in1 ${DATA_DIR_TOPP}/TargetedFileConverter_2_output.TraML.tmp -in2 ${DATA_DIR_TOPP}/TargetedFileConverter_2_output.TraML)
  set_tests_properties("TOPP_TargetedFileConverter_test_2" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_1_prepare")
  set_tests_properties("TOPP_TargetedFileConverter_test_2_out1" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_2")

  # Metabolomics
  add_test("TOPP_ConvertTSVToTraML_test_5" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ${DATA_DIR_TOPP}/ConvertTSVToTraML_5_input.tsv -out ConvertTSVToTraML_5_output.TraML.tmp -out_type TraML)
  add_test("TOPP_ConvertTSVToTraML_test_5_out1" ${DIFF} -in1 ConvertTSVToTraML_5_output.TraML.tmp -in2 ${DATA_DIR_TOPP}/ConvertTSVToTraML_5_output.TraML)
  set_tests_properties("TOPP_ConvertTSVToTraML_test_5_out1" PROPERTIES DEPENDS "TOPP_ConvertTSVToTraML_test_5")

  add_test("TOPP_ConvertTraMLToTSV_test_1" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ${DATA_DIR_TOPP}/ConvertTSVToTraML_output.TraML -out ConvertTraMLToTSV_output.tsv.tmp -out_type tsv)
  add_test("TOPP_ConvertTraMLToTSV_test_1_out1" ${DIFF} -in1 ConvertTraMLToTSV_output.tsv.tmp -in2 ${DATA_DIR_TOPP}/ConvertTraMLToTSV_output.tsv)
  set_tests_properties("TOPP_ConvertTraMLToTSV_test_1_out1" PROPERTIES DEPENDS "TOPP_ConvertTraMLToTSV_test_1")
  # Metabolomics
  add_test("TOPP_ConvertTraMLToTSV_test_2" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ${DATA_DIR_TOPP}/ConvertTSVToTraML_5_output.TraML -out ConvertTraMLToTSV_output_2.tsv.tmp -out_type tsv)
  add_test("TOPP_ConvertTraMLToTSV_test_2_out1" ${DIFF} -in1 ConvertTraMLToTSV_output_2.tsv.tmp -in2 ${DATA_DIR_TOPP}/ConvertTraMLToTSV_output_2.tsv)
  set_tests_properties("TOPP_ConvertTraMLToTSV_test_2_out1" PROPERTIES DEPENDS "TOPP_ConvertTraMLToTSV_test_2")

  # Metabolomics PQP - TraML Conversion
  add_test("TOPP_TargetedFileConverter_test_3_prepare" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ${DATA_DIR_TOPP}/TargetedFileConverter_3_input.TraML -out TargetedFileConverter_3_output.pqp.tmp -out_type pqp)
  add_test("TOPP_TargetedFileConverter_test_3" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in TargetedFileConverter_3_output.pqp.tmp -in_type pqp -out TargetedFileConverter_3_output.TraML.tmp -out_type TraML)
  add_test("TOPP_TargetedFileConverter_test_3_out1" ${DIFF} -in1 TargetedFileConverter_3_output.TraML.tmp -in2 ${DATA_DIR_TOPP}/TargetedFileConverter_3_output.TraML)
  set_tests_properties("TOPP_TargetedFileConverter_test_3" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_3_prepare")
  set_tests_properties("TOPP_TargetedFileConverter_test_3_out1" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_3")
  add_test("TOPP_TargetedFileConverter_test_4" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in TargetedFileConverter_3_output.pqp.tmp -in_type pqp -out TargetedFileConverter_4_output.TraML.tmp -out_type TraML -legacy_traml_id)
  add_test("TOPP_TargetedFileConverter_test_4_out1" ${DIFF} -in1 TargetedFileConverter_4_output.TraML.tmp -in2 ${DATA_DIR_TOPP}/TargetedFileConverter_4_output.TraML)
  set_tests_properties("TOPP_TargetedFileConverter_test_4_out1" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_4")
  set_tests_properties("TOPP_TargetedFileConverter_test_4" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_3")

  add_test("TOPP_MRMMapping_test_1" ${TOPP_BIN_PATH}/MRMMapper -in ${DATA_DIR_TOPP}/MRMMapping_input.chrom.mzML -tr ${DATA_DIR_TOPP}/MRMMapping_input.TraML -out MRMMapping_output.chrom.mzML.tmp -test)
  add_test("TOPP_MRMMapping_test_1_out_1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 MRMMapping_output.chrom.mzML.tmp -in2 ${DATA_DIR_TOPP}/MRMMapping_output.chrom.mzML)
  set_tests_properties("TOPP_MRMMapping_test_1_out_1" PROPERTIES DEPENDS "TOPP_MRMMapping_test_1")

  add_test("TOPP_OpenSwathFeatureXMLToTSV_test_1" ${TOPP_BIN_PATH}/OpenSwathFeatureXMLToTSV -in ${DATA_DIR_TOPP}/OpenSwathFeatureXMLToTSV_input.featureXML -tr ${DATA_DIR_TOPP}/OpenSwathFeatureXMLToTSV_input.TraML -out OpenSwathFeatureXMLToTSV_output.short.csv.tmp -short_format -test)
  add_test("TOPP_OpenSwathFeatureXMLToTSV_test_1_out1" ${DIFF} -in1 ${DATA_DIR_TOPP}/OpenSwathFeatureXMLToTSV_output.short.csv -in2 OpenSwathFeatureXMLToTSV_output.short.csv.tmp)
  set_tests_properties("TOPP_OpenSwathFeatureXMLToTSV_test_1_out1" PROPERTIES DEPENDS "TOPP_OpenSwathFeatureXMLToTSV_test_1")
  add_test("TOPP_OpenSwathFeatureXMLToTSV_test_2" ${TOPP_BIN_PATH}/OpenSwathFeatureXMLToTSV -in ${DATA_DIR_TOPP}/OpenSwathFeatureXMLToTSV_input.featureXML -tr ${DATA_DIR_TOPP}/OpenSwathFeatureXMLToTSV_input.TraML -out OpenSwathFeatureXMLToTSV_output.long.csv.tmp -test)
  add_test("TOPP_OpenSwathFeatureXMLToTSV_test_2_out1" ${DIFF} -in1 ${DATA_DIR_TOPP}/OpenSwathFeatureXMLToTSV_output.long.csv -in2 OpenSwathFeatureXMLToTSV_output.long.csv.tmp)
  set_tests_properties("TOPP_OpenSwathFeatureXMLToTSV_test_2_out1" PROPERTIES DEPENDS "TOPP_OpenSwathFeatureXMLToTSV_test_2")
  add_test("TOPP_OpenSwathFeatureXMLToTSV_test_3" ${TOPP_BIN_PATH}/OpenSwathFeatureXMLToTSV -in ${DATA_DIR_TOPP}/OpenSwathFeatureXMLToTSV_input.featureXML -tr ${DATA_DIR_TOPP}/OpenSwathFeatureXMLToTSV_input.TraML -out OpenSwathFeatureXMLToTSV_3_output.short.csv.tmp -short_format -best_scoring_peptide main_var_xx_lda_prelim_score -test)
  add_test("TOPP_OpenSwathFeatureXMLToTSV_test_3_out1" ${DIFF} -in1 ${DATA_DIR_TOPP}/OpenSwathFeatureXMLToTSV_3_output.short.csv -in2 OpenSwathFeatureXMLToTSV_3_output.short.csv.tmp)
  set_tests_properties("TOPP_OpenSwathFeatureXMLToTSV_test_3_out1" PROPERTIES DEPENDS "TOPP_OpenSwathFeatureXMLToTSV_test_3")

  add_test("TOPP_OpenSwathChromatogramExtractor_test_1" ${TOPP_BIN_PATH}/OpenSwathChromatogramExtractor -in ${DATA_DIR_TOPP}/OpenSwathChromatogramExtractor_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathChromatogramExtractor_input.TraML -out OpenSwathChromatogramExtractor_output.mzML.tmp -test)
  add_test("TOPP_OpenSwathChromatogramExtractor_test_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 OpenSwathChromatogramExtractor_output.mzML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathChromatogramExtractor_output.mzML)
  set_tests_properties("TOPP_OpenSwathChromatogramExtractor_test_1_out1" PROPERTIES DEPENDS "TOPP_OpenSwathChromatogramExtractor_test_1")
  add_test("TOPP_OpenSwathChromatogramExtractor_test_2" ${TOPP_BIN_PATH}/OpenSwathChromatogramExtractor -in ${DATA_DIR_TOPP}/OpenSwathChromatogramExtractor_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathChromatogramExtractor_input.TraML -rt_norm ${DATA_DIR_TOPP}/OpenSwathChromatogramExtractor_input.trafoXML -out OpenSwathChromatogramExtractor_output_2.mzML.tmp -test -rt_window 50)
  add_test("TOPP_OpenSwathChromatogramExtractor_test_2_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 OpenSwathChromatogramExtractor_output_2.mzML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathChromatogramExtractor_output_2.mzML)
  set_tests_properties("TOPP_OpenSwathChromatogramExtractor_test_2_out1" PROPERTIES DEPENDS "TOPP_OpenSwathChromatogramExtractor_test_2")
  add_test("TOPP_OpenSwathChromatogramExtractor_test_3" ${TOPP_BIN_PATH}/OpenSwathChromatogramExtractor -in ${DATA_DIR_TOPP}/OpenSwathChromatogramExtractor_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathChromatogramExtractor_input.TraML -out OpenSwathChromatogramExtractor_output_3.mzML.tmp -test -extract_MS1)
  add_test("TOPP_OpenSwathChromatogramExtractor_test_3_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 OpenSwathChromatogramExtractor_output_3.mzML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathChromatogramExtractor_output_3.mzML)
  set_tests_properties("TOPP_OpenSwathChromatogramExtractor_test_3_out1" PROPERTIES DEPENDS "TOPP_OpenSwathChromatogramExtractor_test_3")

  add_test("TOPP_OpenSwathAnalyzer_test_1" ${TOPP_BIN_PATH}/OpenSwathAnalyzer -in ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input_chrom.mzML -tr ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input.TraML -out OpenSwathAnalyzer_1_output.featureXML.tmp -rt_norm ${DATA_DIR_TOPP}/OpenSwathAnalyzer_input.trafoXML -test)
  add_test("TOPP_OpenSwathAnalyzer_test_1_out1" ${DIFF} -whitelist "id="  -in1 OpenSwathAnalyzer_1_output.featureXML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_output.featureXML)
  set_tests_properties("TOPP_OpenSwathAnalyzer_test_1_out1" PROPERTIES DEPENDS "TOPP_OpenSwathAnalyzer_test_1")
  add_test("TOPP_OpenSwathAnalyzer_test_2" ${TOPP_BIN_PATH}/OpenSwathAnalyzer -in ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input_chrom.mzML -tr ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input.TraML -out OpenSwathAnalyzer_2_output.featureXML.tmp -swath_files ${DATA_DIR_TOPP}/OpenSwathAnalyzer_2_swathfile.mzML -test)
  add_test("TOPP_OpenSwathAnalyzer_test_2_out1" ${DIFF} -whitelist "id="  -in1 OpenSwathAnalyzer_2_output.featureXML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathAnalyzer_2_output.featureXML)
  set_tests_properties("TOPP_OpenSwathAnalyzer_test_2_out1" PROPERTIES DEPENDS "TOPP_OpenSwathAnalyzer_test_2")
  add_test("TOPP_OpenSwathAnalyzer_test_5_mod" ${TOPP_BIN_PATH}/OpenSwathAnalyzer -in ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input_chrom.mzML -tr ${DATA_DIR_TOPP}/OpenSwathAnalyzer_mod_input.TraML -out OpenSwathAnalyzer_5_output.featureXML.tmp -swath_files ${DATA_DIR_TOPP}/OpenSwathAnalyzer_2_swathfile.mzML -ini ${DATA_DIR_TOPP}/OpenSwathAnalyzer_5.ini -test)
  add_test("TOPP_OpenSwathAnalyzer_test_5_mod_out1" ${DIFF} -whitelist "id="  -in1 OpenSwathAnalyzer_5_output.featureXML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathAnalyzer_5_output.featureXML)
  set_tests_properties("TOPP_OpenSwathAnalyzer_test_5_mod_out1" PROPERTIES DEPENDS "TOPP_OpenSwathAnalyzer_test_5_mod")
  add_test("TOPP_OpenSwathAnalyzer_test_6_nomod" ${TOPP_BIN_PATH}/OpenSwathAnalyzer -in ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input_chrom.mzML -tr ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input.TraML -out OpenSwathAnalyzer_6_output.featureXML.tmp -swath_files ${DATA_DIR_TOPP}/OpenSwathAnalyzer_2_swathfile.mzML -ini ${DATA_DIR_TOPP}/OpenSwathAnalyzer_5.ini -test)
  add_test("TOPP_OpenSwathAnalyzer_test_6_nomod_out1" ${DIFF} -whitelist "id="  -in1 OpenSwathAnalyzer_6_output.featureXML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathAnalyzer_6_output.featureXML)
  set_tests_properties("TOPP_OpenSwathAnalyzer_test_6_nomod_out1" PROPERTIES DEPENDS "TOPP_OpenSwathAnalyzer_test_6_nomod")
  #// TODO write test with TrafoXML


  add_test("TOPP_OpenSwathAnalyzer_test_7_backgroundSubtraction" ${TOPP_BIN_PATH}/OpenSwathAnalyzer -in ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input_chrom.mzML -tr ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input.TraML -out OpenSwathAnalyzer_7_output.featureXML.tmp -ini ${DATA_DIR_TOPP}/OpenSwathAnalyzer_7_backgroundSubtraction.ini -test)
  add_test("TOPP_OpenSwathAnalyzer_test_7_backgroundSubtraction_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathAnalyzer_7_output.featureXML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathAnalyzer_7_output.featureXML)
  set_tests_properties("TOPP_OpenSwathAnalyzer_test_7_backgroundSubtraction_out1" PROPERTIES DEPENDS "TOPP_OpenSwathAnalyzer_test_7_backgroundSubtraction")

  add_test("TOPP_OpenSwathRTNormalizer_test_1" ${TOPP_BIN_PATH}/OpenSwathRTNormalizer -in ${DATA_DIR_TOPP}/OpenSwathRTNormalizer_1_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathRTNormalizer_1_input.TraML -out OpenSwathRTNormalizer_1_output.trafoXML.tmp -test)
  add_test("TOPP_OpenSwathRTNormalizer_test_1_out1" ${DIFF} -in1 OpenSwathRTNormalizer_1_output.trafoXML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathRTNormalizer_1_output.trafoXML)
  set_tests_properties("TOPP_OpenSwathRTNormalizer_test_1_out1" PROPERTIES DEPENDS "TOPP_OpenSwathRTNormalizer_test_1")
  add_test("TOPP_OpenSwathRTNormalizer_test_2" ${TOPP_BIN_PATH}/OpenSwathRTNormalizer -in ${DATA_DIR_TOPP}/OpenSwathRTNormalizer_1_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathRTNormalizer_1_input.TraML -out OpenSwathRTNormalizer_1_output.trafoXML.tmp -test -estimateBestPeptides -peptideEstimation:NrRTBins 30 -peptideEstimation:MinBinsFilled 30)
  set_tests_properties("TOPP_OpenSwathRTNormalizer_test_2" PROPERTIES WILL_FAIL 1) ## will not reach bin coverage
  add_test("TOPP_OpenSwathRTNormalizer_test_3" ${TOPP_BIN_PATH}/OpenSwathRTNormalizer -in ${DATA_DIR_TOPP}/OpenSwathRTNormalizer_1_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathRTNormalizer_1_input.TraML -out OpenSwathRTNormalizer_3_output.trafoXML.tmp -test -estimateBestPeptides -peptideEstimation:NrRTBins 3 -peptideEstimation:MinBinsFilled 3)
  add_test("TOPP_OpenSwathRTNormalizer_test_3_out1" ${DIFF} -in1 OpenSwathRTNormalizer_3_output.trafoXML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathRTNormalizer_3_output.trafoXML)
  set_tests_properties("TOPP_OpenSwathRTNormalizer_test_3_out1" PROPERTIES DEPENDS "TOPP_OpenSwathRTNormalizer_test_3")
  add_test("TOPP_OpenSwathRTNormalizer_test_4" ${TOPP_BIN_PATH}/OpenSwathRTNormalizer -in ${DATA_DIR_TOPP}/OpenSwathRTNormalizer_1_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathRTNormalizer_1_input.TraML -out OpenSwathRTNormalizer_4_output.trafoXML.tmp -test -outlierDetection:outlierMethod iter_residual)
  add_test("TOPP_OpenSwathRTNormalizer_test_4_out1" ${DIFF} -in1 OpenSwathRTNormalizer_4_output.trafoXML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathRTNormalizer_4_output.trafoXML)
  set_tests_properties("TOPP_OpenSwathRTNormalizer_test_4_out1" PROPERTIES DEPENDS "TOPP_OpenSwathRTNormalizer_test_4")

  add_test("TOPP_OpenSwathConfidenceScoring_1" ${TOPP_BIN_PATH}/OpenSwathConfidenceScoring -test -in ${DATA_DIR_TOPP}/OpenSwathFeatureXMLToTSV_input.featureXML -lib ${DATA_DIR_TOPP}/OpenSwathFeatureXMLToTSV_input.TraML -trafo ${DATA_DIR_TOPP}/OpenSwathConfidenceScoring_1_input.trafoXML -transitions 2 -decoys 1 -out OpenSwathConfidenceScoring_1_output.tmp)
  add_test("TOPP_OpenSwathConfidenceScoring_1_out1" ${DIFF} -in1 OpenSwathConfidenceScoring_1_output.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathConfidenceScoring_1_output.featureXML)
  set_tests_properties("TOPP_OpenSwathConfidenceScoring_1_out1" PROPERTIES DEPENDS "TOPP_OpenSwathConfidenceScoring_1")

  ### OpenSwath with FileCacher tests (first spectra, then chromatograms)
  add_test("TOPP_OpenSwathMzMLFileCacher_test_1_step1" ${TOPP_BIN_PATH}/OpenSwathMzMLFileCacher -in ${DATA_DIR_TOPP}/OpenSwathMzMLFileCacher_1_input.mzML -out OpenSwathMzMLFileCacher_1_input.cached.tmp.mzML -test)
  add_test("TOPP_OpenSwathMzMLFileCacher_test_1_step2" ${TOPP_BIN_PATH}/OpenSwathMzMLFileCacher -in OpenSwathMzMLFileCacher_1_input.cached.tmp.mzML -out OpenSwathMzMLFileCacher_1_output.tmp.mzML -convert_back -test)
  add_test("TOPP_OpenSwathMzMLFileCacher_test_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 OpenSwathMzMLFileCacher_1_input.cached.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathMzMLFileCacher_1_output_1.mzML)
  add_test("TOPP_OpenSwathMzMLFileCacher_test_1_out2" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 OpenSwathMzMLFileCacher_1_output.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathMzMLFileCacher_1_output_2.mzML)
  set_tests_properties("TOPP_OpenSwathMzMLFileCacher_test_1_step2" PROPERTIES DEPENDS "TOPP_OpenSwathMzMLFileCacher_test_1_step1")
  set_tests_properties("TOPP_OpenSwathMzMLFileCacher_test_1_out1" PROPERTIES DEPENDS "TOPP_OpenSwathMzMLFileCacher_test_1_step1")
  set_tests_properties("TOPP_OpenSwathMzMLFileCacher_test_1_out2" PROPERTIES DEPENDS "TOPP_OpenSwathMzMLFileCacher_test_1_step2")

  add_test("TOPP_OpenSwathMzMLFileCacher_test_2_step1" ${TOPP_BIN_PATH}/OpenSwathMzMLFileCacher -in ${DATA_DIR_TOPP}/OpenSwathMzMLFileCacher_2_input.chrom.mzML -out OpenSwathMzMLFileCacher_2_input.chrom.cached.tmp.mzML -test)
  add_test("TOPP_OpenSwathMzMLFileCacher_test_2_step2" ${TOPP_BIN_PATH}/OpenSwathMzMLFileCacher -in OpenSwathMzMLFileCacher_2_input.chrom.cached.tmp.mzML -out OpenSwathMzMLFileCacher_2_output.chrom.tmp.mzML -convert_back -test)
  add_test("TOPP_OpenSwathMzMLFileCacher_test_2_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 OpenSwathMzMLFileCacher_2_input.chrom.cached.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathMzMLFileCacher_2_output_1.chrom.mzML)
  add_test("TOPP_OpenSwathMzMLFileCacher_test_2_out2" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 OpenSwathMzMLFileCacher_2_output.chrom.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathMzMLFileCacher_2_output_2.chrom.mzML)
  set_tests_properties("TOPP_OpenSwathMzMLFileCacher_test_2_step2" PROPERTIES DEPENDS "TOPP_OpenSwathMzMLFileCacher_test_2_step1")
  set_tests_properties("TOPP_OpenSwathMzMLFileCacher_test_2_out1" PROPERTIES DEPENDS "TOPP_OpenSwathMzMLFileCacher_test_2_step1")
  set_tests_properties("TOPP_OpenSwathMzMLFileCacher_test_2_out2" PROPERTIES DEPENDS "TOPP_OpenSwathMzMLFileCacher_test_2_step2")

  add_test("TOPP_OpenSwathMzMLFileCacher_test_3_step1" ${TOPP_BIN_PATH}/OpenSwathMzMLFileCacher -in ${DATA_DIR_TOPP}/OpenSwathMzMLFileCacher_1_input.mzML -out OpenSwathMzMLFileCacher_3_input.tmp.sqMass -test)
  add_test("TOPP_OpenSwathMzMLFileCacher_test_3_step2" ${TOPP_BIN_PATH}/OpenSwathMzMLFileCacher -in OpenSwathMzMLFileCacher_3_input.tmp.sqMass -out OpenSwathMzMLFileCacher_3_output.tmp.mzML -test)
  add_test("TOPP_OpenSwathMzMLFileCacher_test_3_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 OpenSwathMzMLFileCacher_3_output.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathMzMLFileCacher_3_output.mzML)
  set_tests_properties("TOPP_OpenSwathMzMLFileCacher_test_3_step2" PROPERTIES DEPENDS "TOPP_OpenSwathMzMLFileCacher_test_3_step1")
  set_tests_properties("TOPP_OpenSwathMzMLFileCacher_test_3_out1" PROPERTIES DEPENDS "TOPP_OpenSwathMzMLFileCacher_test_3_step2")

  add_test("TOPP_OpenSwathMzMLFileCacher_test_4_step1" ${TOPP_BIN_PATH}/OpenSwathMzMLFileCacher -in ${DATA_DIR_TOPP}/OpenSwathMzMLFileCacher_2_input.chrom.mzML -out OpenSwathMzMLFileCacher_4_input.tmp.sqMass -test)
  add_test("TOPP_OpenSwathMzMLFileCacher_test_4_step2" ${TOPP_BIN_PATH}/OpenSwathMzMLFileCacher -in OpenSwathMzMLFileCacher_4_input.tmp.sqMass -out OpenSwathMzMLFileCacher_4_output.tmp.mzML -test)
  add_test("TOPP_OpenSwathMzMLFileCacher_test_4_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 OpenSwathMzMLFileCacher_4_output.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathMzMLFileCacher_4_output.mzML)
  set_tests_properties("TOPP_OpenSwathMzMLFileCacher_test_4_step2" PROPERTIES DEPENDS "TOPP_OpenSwathMzMLFileCacher_test_4_step1")
  set_tests_properties("TOPP_OpenSwathMzMLFileCacher_test_4_out1" PROPERTIES DEPENDS "TOPP_OpenSwathMzMLFileCacher_test_4_step2")

  # Testing the OpenSwathAnalyzer together with the FileCacher (for spectra and chromatograms)
  add_test("TOPP_OpenSwathAnalyzer_test_3_prepare" ${TOPP_BIN_PATH}/OpenSwathMzMLFileCacher -in ${DATA_DIR_TOPP}/OpenSwathAnalyzer_2_swathfile.mzML -out OpenSwathAnalyzer_3_swathfile.mzML.cached.tmp -out_type mzML -test)
  add_test("TOPP_OpenSwathAnalyzer_test_3" ${TOPP_BIN_PATH}/OpenSwathAnalyzer -in ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input_chrom.mzML -tr ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input.TraML -out MRMFeatureFinderScore_output_3.featureXML.tmp -swath_files OpenSwathAnalyzer_3_swathfile.mzML.cached.tmp -test)
  add_test("TOPP_OpenSwathAnalyzer_test_3_out1" ${DIFF} -whitelist "id=" -in1 MRMFeatureFinderScore_output_3.featureXML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathAnalyzer_2_output.featureXML)
  set_tests_properties("TOPP_OpenSwathAnalyzer_test_3" PROPERTIES DEPENDS "TOPP_OpenSwathAnalyzer_test_3_prepare")
  set_tests_properties("TOPP_OpenSwathAnalyzer_test_3_out1" PROPERTIES DEPENDS "TOPP_OpenSwathAnalyzer_test_3")
  add_test("TOPP_OpenSwathAnalyzer_test_4_prepare1" ${TOPP_BIN_PATH}/OpenSwathMzMLFileCacher -in ${DATA_DIR_TOPP}/OpenSwathAnalyzer_2_swathfile.mzML -out OpenSwathAnalyzer_4_swathfile.mzML.cached.tmp -out_type mzML -test)
  add_test("TOPP_OpenSwathAnalyzer_test_4_prepare2" ${TOPP_BIN_PATH}/OpenSwathMzMLFileCacher -in ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input_chrom.mzML -out OpenSwathAnalyzer_4_input_chrom.mzML.cached.tmp -out_type mzML -test)
  add_test("TOPP_OpenSwathAnalyzer_test_4" ${TOPP_BIN_PATH}/OpenSwathAnalyzer -in OpenSwathAnalyzer_4_input_chrom.mzML.cached.tmp -tr ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input.TraML -out MRMFeatureFinderScore_output_4.featureXML.tmp -swath_files OpenSwathAnalyzer_4_swathfile.mzML.cached.tmp -test)
  add_test("TOPP_OpenSwathAnalyzer_test_4_out1" ${DIFF} -whitelist "id=" -in1 MRMFeatureFinderScore_output_4.featureXML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathAnalyzer_2_output.featureXML)
  set_tests_properties("TOPP_OpenSwathAnalyzer_test_4" PROPERTIES DEPENDS "TOPP_OpenSwathAnalyzer_test_4_prepare1;TOPP_OpenSwathAnalyzer_test_4_prepare2")
  set_tests_properties("TOPP_OpenSwathAnalyzer_test_4_out1" PROPERTIES DEPENDS "TOPP_OpenSwathAnalyzer_test_4")

  add_test("UTILS_MRMTransitionGroupPicker_test_1" ${TOPP_BIN_PATH}/MRMTransitionGroupPicker -in ${DATA_DIR_TOPP}/MRMTransitionGroupPicker_1_input.mzML -tr ${DATA_DIR_TOPP}/MRMTransitionGroupPicker_1_input.TraML -out MRMTransitionGroupPicker_test_1.featureXML.tmp -test -algorithm:PeakPickerMRM:remove_overlapping_peaks true -algorithm:PeakPickerMRM:method legacy)
  add_test("UTILS_MRMTransitionGroupPicker_test_1_out" ${DIFF} -whitelist "id=" -in1 MRMTransitionGroupPicker_test_1.featureXML.tmp -in2 ${DATA_DIR_TOPP}/MRMTransitionGroupPicker_1_output.featureXML)
  set_tests_properties("UTILS_MRMTransitionGroupPicker_test_1_out" PROPERTIES DEPENDS "UTILS_MRMTransitionGroupPicker_test_1")

  #------------------------------------------------------------------------------
  # OpenSwath Workflow tests
  add_test("TOPP_OpenSwathWorkflow_1" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.TraML -rt_norm ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.trafoXML -out_chrom OpenSwathWorkflow_1.chrom.mzML.tmp -out_features OpenSwathWorkflow_1.featureXML.tmp -test)
  add_test("TOPP_OpenSwathWorkflow_1_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_1.featureXML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_output.featureXML)
  add_test("TOPP_OpenSwathWorkflow_1_out2" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_1.chrom.mzML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_output.chrom.mzML)
  set_tests_properties("TOPP_OpenSwathWorkflow_1_out1" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_1")
  set_tests_properties("TOPP_OpenSwathWorkflow_1_out2" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_1")

  add_test("TOPP_OpenSwathWorkflow_2" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_2_input.mzXML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_2_input.TraML -rt_norm ${DATA_DIR_TOPP}/OpenSwathWorkflow_2_input.trafoXML -out_chrom OpenSwathWorkflow_2.chrom.mzML.tmp -out_features OpenSwathWorkflow_2.featureXML.tmp  -test) 
  add_test("TOPP_OpenSwathWorkflow_2_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_2.featureXML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_2_output.featureXML)
  add_test("TOPP_OpenSwathWorkflow_2_out2" ${DIFF} -in1 OpenSwathWorkflow_2.chrom.mzML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_2_output.chrom.mzML)
  set_tests_properties("TOPP_OpenSwathWorkflow_2_out1" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_2")
  set_tests_properties("TOPP_OpenSwathWorkflow_2_out2" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_2")

  # Also test with ms1 traces (precursor ion trace)
  add_test("TOPP_OpenSwathWorkflow_3" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.TraML -rt_norm ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.trafoXML -out_chrom OpenSwathWorkflow_3.chrom.mzML.tmp -out_features OpenSwathWorkflow_3.featureXML.tmp -test -use_ms1_traces)
  add_test("TOPP_OpenSwathWorkflow_3_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_3.featureXML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_3_output.featureXML)
  add_test("TOPP_OpenSwathWorkflow_3_out2" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_3.chrom.mzML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_3_output.chrom.mzML)
  set_tests_properties("TOPP_OpenSwathWorkflow_3_out1" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_3")
  set_tests_properties("TOPP_OpenSwathWorkflow_3_out2" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_3")

  # Also test whether it is able to write csv output
  add_test("TOPP_OpenSwathWorkflow_4" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.TraML -rt_norm ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.trafoXML -out_chrom OpenSwathWorkflow_4.chrom.mzML.tmp -out_tsv OpenSwathWorkflow_4.tsv.tmp -test -use_ms1_traces)
  add_test("TOPP_OpenSwathWorkflow_4_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_4.chrom.mzML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_3_output.chrom.mzML)
  set_tests_properties("TOPP_OpenSwathWorkflow_4_out1" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_4")
  # We cannot currently test the correctness of the csv output

  # Also test with readoptions cache
  add_test("TOPP_OpenSwathWorkflow_5" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.TraML -rt_norm ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.trafoXML -out_chrom OpenSwathWorkflow_5.chrom.mzML.tmp -out_features OpenSwathWorkflow_5.featureXML.tmp -test -use_ms1_traces -readOptions cache -tempDirectory "osw_file5")
  add_test("TOPP_OpenSwathWorkflow_5_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_5.featureXML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_3_output.featureXML)
  add_test("TOPP_OpenSwathWorkflow_5_out2" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_5.chrom.mzML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_3_output.chrom.mzML)
  set_tests_properties("TOPP_OpenSwathWorkflow_5_out1" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_5")
  set_tests_properties("TOPP_OpenSwathWorkflow_5_out2" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_5")

  # Also test with readoptions cacheWorkingInMemory
  add_test("TOPP_OpenSwathWorkflow_6" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.TraML -rt_norm ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.trafoXML -out_chrom OpenSwathWorkflow_6.chrom.mzML.tmp -out_features OpenSwathWorkflow_6.featureXML.tmp -test -use_ms1_traces -readOptions cacheWorkingInMemory -tempDirectory "osw_file6")
  add_test("TOPP_OpenSwathWorkflow_6_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_6.featureXML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_3_output.featureXML)
  add_test("TOPP_OpenSwathWorkflow_6_out2" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_6.chrom.mzML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_3_output.chrom.mzML)
  set_tests_properties("TOPP_OpenSwathWorkflow_6_out1" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_6")
  set_tests_properties("TOPP_OpenSwathWorkflow_6_out2" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_6")

  # Test with input swath windows file
  add_test("TOPP_OpenSwathWorkflow_7" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.TraML -rt_norm ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.trafoXML -out_chrom OpenSwathWorkflow_7.chrom.mzML.tmp -out_features OpenSwathWorkflow_7.featureXML.tmp -test -use_ms1_traces -swath_windows_file ${DATA_DIR_TOPP}/swath_windows.txt)
  add_test("TOPP_OpenSwathWorkflow_7_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_7.featureXML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_3_output.featureXML)
  add_test("TOPP_OpenSwathWorkflow_7_out2" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_7.chrom.mzML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_3_output.chrom.mzML)
  set_tests_properties("TOPP_OpenSwathWorkflow_7_out1" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_7")
  set_tests_properties("TOPP_OpenSwathWorkflow_7_out2" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_7")

  # Test with faulty swath windows file
  add_test("TOPP_OpenSwathWorkflow_8" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.TraML -rt_norm ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.trafoXML -out_chrom OpenSwathWorkflow_8.chrom.mzML.tmp -out_features OpenSwathWorkflow_8.featureXML.tmp -test -use_ms1_traces -swath_windows_file ${DATA_DIR_TOPP}/swath_windows_overlap.txt)
  set_tests_properties("TOPP_OpenSwathWorkflow_8" PROPERTIES WILL_FAIL 1) # file contains overlaps
  add_test("TOPP_OpenSwathWorkflow_9" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.TraML -rt_norm ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.trafoXML -out_chrom OpenSwathWorkflow_9.chrom.mzML.tmp -out_features OpenSwathWorkflow_9.featureXML.tmp -test -use_ms1_traces -swath_windows_file ${DATA_DIR_TOPP}/swath_windows_tooshort.txt)
  set_tests_properties("TOPP_OpenSwathWorkflow_9" PROPERTIES WILL_FAIL 1) # file does not have enough entries
  # works with force (does not check for overlaps and missing swath windows any more)
  add_test("TOPP_OpenSwathWorkflow_10" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.TraML -rt_norm ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.trafoXML -out_chrom OpenSwathWorkflow_10.chrom.mzML.tmp -out_features OpenSwathWorkflow_10.featureXML.tmp -test -use_ms1_traces -swath_windows_file ${DATA_DIR_TOPP}/swath_windows_overlap.txt -force)

  # test SONAR (also test iRT correction ... )
  # - use iRT correction (but no outlier removal since only 2 datapoints)
  # - use quadratic fitting for m/z correction
  # - use 550 ppm for extraction (0.05 Da @ 100 m/z)
  add_test("TOPP_OpenSwathWorkflow_11" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_11_input.mzML -tr_irt ${DATA_DIR_TOPP}/OpenSwathWorkflow_11_input.TraML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_11_input.TraML -mz_extraction_window 0.2 -rt_extraction_window -1 -Scoring:Scores:use_sonar_scores -sonar -out_chrom OpenSwathWorkflow_11.chrom.mzML.tmp -out_features OpenSwathWorkflow_11.featureXML.tmp -outlierDetection:outlierMethod none -mz_correction_function quadratic_regression_delta_ppm -irt_mz_extraction_window 550 -ppm_irtwindow -test)
  add_test("TOPP_OpenSwathWorkflow_11_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_11.chrom.mzML.tmp  -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_11_output.chrom.mzML)
  add_test("TOPP_OpenSwathWorkflow_11_out2" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_11.featureXML.tmp  -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_11_output.featureXML)
  set_tests_properties("TOPP_OpenSwathWorkflow_11_out1" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_11")
  set_tests_properties("TOPP_OpenSwathWorkflow_11_out2" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_11")

  # 20 ppm is too small for extraction and no iRT peptides are found (thus the tool fails)
  add_test("TOPP_OpenSwathWorkflow_12" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_11_input.mzML -tr_irt ${DATA_DIR_TOPP}/OpenSwathWorkflow_11_input.TraML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_11_input.TraML -mz_extraction_window 0.2 -rt_extraction_window -1 -Scoring:Scores:use_sonar_scores -sonar -out_chrom OpenSwathWorkflow_11.chrom.mzML.tmp -out_features OpenSwathWorkflow_11.featureXML.tmp -outlierDetection:outlierMethod none -mz_correction_function quadratic_regression_delta_ppm -irt_mz_extraction_window 20 -ppm_irtwindow -test)
  set_tests_properties("TOPP_OpenSwathWorkflow_12" PROPERTIES WILL_FAIL 1) # there are no points for iRT regression

  # Test whether PQP files can be read and OSW output can be written
  add_test("TOPP_OpenSwathWorkflow_13_prepare" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.TraML -out OpenSwathWorkflow_1_input.pqp.tmp -out_type pqp)
  add_test("TOPP_OpenSwathWorkflow_13" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.mzML -tr OpenSwathWorkflow_1_input.pqp.tmp -tr_type pqp -rt_norm ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.trafoXML -out_chrom OpenSwathWorkflow_13.chrom.mzML.tmp -out_osw OpenSwathWorkflow_13.osw -test -use_ms1_traces)
  add_test("TOPP_OpenSwathWorkflow_13_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_4.chrom.mzML.tmp -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_3_output.chrom.mzML)
  set_tests_properties("TOPP_OpenSwathWorkflow_13" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_13_prepare")
  set_tests_properties("TOPP_OpenSwathWorkflow_13_out1" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_13")
  # We cannot currently test the correctness of the PQP and OSW files

  add_test("TOPP_OpenSwathFileSplitter_1" ${TOPP_BIN_PATH}/OpenSwathFileSplitter -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.mzML -outputDirectory "" -test)
  add_test("TOPP_OpenSwathFileSplitter_1_out1" ${DIFF} -whitelist "id=" -in1 openswath_tmpfile_ms1.mzML -in2 ${DATA_DIR_TOPP}/openswath_tmpfile_ms1.mzML)
  add_test("TOPP_OpenSwathFileSplitter_1_out2" ${DIFF} -whitelist "id=" -in1 openswath_tmpfile_4.mzML -in2 ${DATA_DIR_TOPP}/openswath_tmpfile_4.mzML)
  set_tests_properties("TOPP_OpenSwathFileSplitter_1_out1" PROPERTIES DEPENDS "TOPP_OpenSwathFileSplitter_1")
  set_tests_properties("TOPP_OpenSwathFileSplitter_1_out2" PROPERTIES DEPENDS "TOPP_OpenSwathFileSplitter_1")

endif(NOT DISABLE_OPENSWATH)

#------------------------------------------------------------------------------
# NoiseFilter tests
add_test("TOPP_NoiseFilterGaussian_1" ${TOPP_BIN_PATH}/NoiseFilterGaussian -test -ini ${DATA_DIR_TOPP}/NoiseFilterGaussian_1_parameters.ini -in ${DATA_DIR_TOPP}/NoiseFilterGaussian_1_input.mzML -out NoiseFilterGaussian_1.tmp)
add_test("TOPP_NoiseFilterGaussian_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 NoiseFilterGaussian_1.tmp -in2 ${DATA_DIR_TOPP}/NoiseFilterGaussian_1_output.mzML )
set_tests_properties("TOPP_NoiseFilterGaussian_1_out1" PROPERTIES DEPENDS "TOPP_NoiseFilterGaussian_1")
add_test("TOPP_NoiseFilterGaussian_2" ${TOPP_BIN_PATH}/NoiseFilterGaussian -test -ini ${DATA_DIR_TOPP}/NoiseFilterGaussian_2_parameters.ini -in ${DATA_DIR_TOPP}/NoiseFilterGaussian_2_input.chrom.mzML -out NoiseFilterGaussian_2.tmp)
add_test("TOPP_NoiseFilterGaussian_2_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 NoiseFilterGaussian_2.tmp -in2 ${DATA_DIR_TOPP}/NoiseFilterGaussian_2_output.chrom.mzML )
set_tests_properties("TOPP_NoiseFilterGaussian_2_out1" PROPERTIES DEPENDS "TOPP_NoiseFilterGaussian_2")

# lowmem variation of the algorithm
add_test("TOPP_NoiseFilterGaussian_3" ${TOPP_BIN_PATH}/NoiseFilterGaussian -test -ini ${DATA_DIR_TOPP}/NoiseFilterGaussian_1_parameters.ini -in ${DATA_DIR_TOPP}/NoiseFilterGaussian_1_input.mzML -out NoiseFilterGaussian_3.tmp -processOption lowmemory)
add_test("TOPP_NoiseFilterGaussian_3_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 NoiseFilterGaussian_3.tmp -in2 ${DATA_DIR_TOPP}/NoiseFilterGaussian_1_output.mzML )
set_tests_properties("TOPP_NoiseFilterGaussian_3_out1" PROPERTIES DEPENDS "TOPP_NoiseFilterGaussian_3")
add_test("TOPP_NoiseFilterGaussian_4" ${TOPP_BIN_PATH}/NoiseFilterGaussian -test -ini ${DATA_DIR_TOPP}/NoiseFilterGaussian_2_parameters.ini -in ${DATA_DIR_TOPP}/NoiseFilterGaussian_2_input.chrom.mzML -out NoiseFilterGaussian_4.tmp -processOption lowmemory)
add_test("TOPP_NoiseFilterGaussian_4_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 NoiseFilterGaussian_4.tmp -in2 ${DATA_DIR_TOPP}/NoiseFilterGaussian_2_output.chrom.mzML )
set_tests_properties("TOPP_NoiseFilterGaussian_4_out1" PROPERTIES DEPENDS "TOPP_NoiseFilterGaussian_4")

add_test("TOPP_NoiseFilterSGolay_1" ${TOPP_BIN_PATH}/NoiseFilterSGolay -test -ini ${DATA_DIR_TOPP}/NoiseFilterSGolay_1_parameters.ini -in ${DATA_DIR_TOPP}/NoiseFilterSGolay_1_input.mzML -out NoiseFilterSGolay_1.tmp)
add_test("TOPP_NoiseFilterSGolay_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 NoiseFilterSGolay_1.tmp -in2 ${DATA_DIR_TOPP}/NoiseFilterSGolay_1_output.mzML )
set_tests_properties("TOPP_NoiseFilterSGolay_1_out1" PROPERTIES DEPENDS "TOPP_NoiseFilterSGolay_1")
add_test("TOPP_NoiseFilterSGolay_2" ${TOPP_BIN_PATH}/NoiseFilterSGolay -test -ini ${DATA_DIR_TOPP}/NoiseFilterSGolay_2_parameters.ini -in ${DATA_DIR_TOPP}/NoiseFilterSGolay_2_input.chrom.mzML -out NoiseFilterSGolay_2.tmp)
add_test("TOPP_NoiseFilterSGolay_2_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 NoiseFilterSGolay_2.tmp -in2 ${DATA_DIR_TOPP}/NoiseFilterSGolay_2_output.chrom.mzML )
set_tests_properties("TOPP_NoiseFilterSGolay_2_out1" PROPERTIES DEPENDS "TOPP_NoiseFilterSGolay_2")

# lowmem variation of the algorithm
add_test("TOPP_NoiseFilterSGolay_3" ${TOPP_BIN_PATH}/NoiseFilterSGolay -test -ini ${DATA_DIR_TOPP}/NoiseFilterSGolay_1_parameters.ini -in ${DATA_DIR_TOPP}/NoiseFilterSGolay_1_input.mzML -out NoiseFilterSGolay_3.tmp -processOption lowmemory)
add_test("TOPP_NoiseFilterSGolay_3_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 NoiseFilterSGolay_3.tmp -in2 ${DATA_DIR_TOPP}/NoiseFilterSGolay_1_output.mzML )
set_tests_properties("TOPP_NoiseFilterSGolay_3_out1" PROPERTIES DEPENDS "TOPP_NoiseFilterSGolay_3")
add_test("TOPP_NoiseFilterSGolay_4" ${TOPP_BIN_PATH}/NoiseFilterSGolay -test -ini ${DATA_DIR_TOPP}/NoiseFilterSGolay_2_parameters.ini -in ${DATA_DIR_TOPP}/NoiseFilterSGolay_2_input.chrom.mzML -out NoiseFilterSGolay_4.tmp -processOption lowmemory)
add_test("TOPP_NoiseFilterSGolay_4_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 NoiseFilterSGolay_4.tmp -in2 ${DATA_DIR_TOPP}/NoiseFilterSGolay_2_output.chrom.mzML )
set_tests_properties("TOPP_NoiseFilterSGolay_4_out1" PROPERTIES DEPENDS "TOPP_NoiseFilterSGolay_4")

#------------------------------------------------------------------------------
# PeakPicker tests
add_test("TOPP_PeakPickerWavelet_1" ${TOPP_BIN_PATH}/PeakPickerWavelet  -test -ini ${DATA_DIR_TOPP}/PeakPickerWavelet_parameters.ini -in ${DATA_DIR_TOPP}/PeakPickerWavelet_input.mzML -out PeakPickerWavelet_1.tmp)
add_test("TOPP_PeakPickerWavelet_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 PeakPickerWavelet_1.tmp -in2 ${DATA_DIR_TOPP}/PeakPickerWavelet_output.mzML )
set_tests_properties("TOPP_PeakPickerWavelet_1_out1" PROPERTIES DEPENDS "TOPP_PeakPickerWavelet_1")
add_test("TOPP_PeakPickerWavelet_2" ${TOPP_BIN_PATH}/PeakPickerWavelet  -test -ini ${DATA_DIR_TOPP}/PeakPickerWavelet_deconv_parameters.ini -in ${DATA_DIR_TOPP}/PeakPickerWavelet_deconv_input.mzML -out PeakPickerWavelet_2.tmp)
add_test("TOPP_PeakPickerWavelet_2_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 PeakPickerWavelet_2.tmp -in2 ${DATA_DIR_TOPP}/PeakPickerWavelet_deconv_output.mzML)
set_tests_properties("TOPP_PeakPickerWavelet_2_out1" PROPERTIES DEPENDS "TOPP_PeakPickerWavelet_2")
add_test("TOPP_PeakPickerWavelet_3" ${TOPP_BIN_PATH}/PeakPickerWavelet  -test -ini ${DATA_DIR_TOPP}/PeakPickerWavelet_2Dopt_parameters.ini -in ${DATA_DIR_TOPP}/PeakPickerWavelet_2Dopt_input.mzML -out PeakPickerWavelet_3.tmp)
if(NOT DISABLE_WAVELET2DTEST)
  #  
  # See https://github.com/OpenMS/OpenMS/issues/1202
  # On Ubuntu 14.04 in Release mode only (not Debug mode), the output of this
  # file is slightly different due to a small numerical error. Since FuzzyDiff is
  # not able to compare the numbers in base64, we have to "whitelist" the
  # corresponding base64 string.
  # See https://github.com/OpenMS/OpenMS/issues/1555
  # For travis, PPWavelet_3 does not work, so we disable it
  add_test("TOPP_PeakPickerWavelet_3_out1" ${DIFF} -in1 PeakPickerWavelet_3.tmp -in2 ${DATA_DIR_TOPP}/PeakPickerWavelet_2Dopt_output.mzML -whitelist ${INDEX_WHITELIST} "3bAxk01KjUAAAACgmVGNQIvf6JpEWo1Ah")
  set_tests_properties("TOPP_PeakPickerWavelet_3_out1" PROPERTIES DEPENDS "TOPP_PeakPickerWavelet_3")
endif(NOT DISABLE_WAVELET2DTEST)
add_test("TOPP_PeakPickerWavelet_4" ${TOPP_BIN_PATH}/PeakPickerWavelet  -test -ini ${DATA_DIR_TOPP}/PeakPickerWavelet_parameters.ini -in ${DATA_DIR_TOPP}/PeakPickerWavelet_input.mzML -out PeakPickerWavelet_4.tmp -threads 2)
add_test("TOPP_PeakPickerWavelet_4_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 PeakPickerWavelet_4.tmp -in2 ${DATA_DIR_TOPP}/PeakPickerWavelet_output.mzML)
set_tests_properties("TOPP_PeakPickerWavelet_4_out1" PROPERTIES DEPENDS "TOPP_PeakPickerWavelet_4")
add_test("TOPP_PeakPickerWavelet_5" ${TOPP_BIN_PATH}/PeakPickerWavelet  -test -ini ${DATA_DIR_TOPP}/PeakPickerWavelet_parameters_noMetaData.ini -in ${DATA_DIR_TOPP}/PeakPickerWavelet_input.mzML -out PeakPickerWavelet_5.tmp -threads 2)
add_test("TOPP_PeakPickerWavelet_5_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 PeakPickerWavelet_5.tmp -in2 ${DATA_DIR_TOPP}/PeakPickerWavelet_output_noMetaData.mzML)
set_tests_properties("TOPP_PeakPickerWavelet_5_out1" PROPERTIES DEPENDS "TOPP_PeakPickerWavelet_5")

# "high_res" algorithm with "ms_levels" option set to only 1:
# .. spectra
add_test("TOPP_PeakPickerHiRes_1" ${TOPP_BIN_PATH}/PeakPickerHiRes -test -ini ${DATA_DIR_TOPP}/PeakPickerHiRes_parameters.ini -in ${DATA_DIR_TOPP}/PeakPickerHiRes_input.mzML -out PeakPickerHiRes_1.tmp)
add_test("TOPP_PeakPickerHiRes_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 PeakPickerHiRes_1.tmp -in2 ${DATA_DIR_TOPP}/PeakPickerHiRes_output.mzML)
set_tests_properties("TOPP_PeakPickerHiRes_1_out1" PROPERTIES DEPENDS "TOPP_PeakPickerHiRes_1")
# .. chromatogram
add_test("TOPP_PeakPickerHiRes_2" ${TOPP_BIN_PATH}/PeakPickerHiRes -test -ini ${DATA_DIR_TOPP}/PeakPickerHiRes_parameters.ini -in ${DATA_DIR_TOPP}/PeakPickerHiRes_2_input.mzML -out PeakPickerHiRes_2.tmp)
add_test("TOPP_PeakPickerHiRes_2_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 PeakPickerHiRes_2.tmp -in2 ${DATA_DIR_TOPP}/PeakPickerHiRes_2_output.mzML)
set_tests_properties("TOPP_PeakPickerHiRes_2_out1" PROPERTIES DEPENDS "TOPP_PeakPickerHiRes_2")

# "high_res" peakpicker with the lowmemory option
# .. spectra
add_test("TOPP_PeakPickerHiRes_3" ${TOPP_BIN_PATH}/PeakPickerHiRes -test -ini ${DATA_DIR_TOPP}/PeakPickerHiRes_parameters.ini -in ${DATA_DIR_TOPP}/PeakPickerHiRes_input.mzML -out PeakPickerHiRes_3.tmp -processOption lowmemory)
## Note: the output SHOULD be identical to 'PeakPickerHiRes_output.mzML', but due to a missing 'Dataprocessing' entry (which is not known when writing the mzML header), we need an extra output file
add_test("TOPP_PeakPickerHiRes_3_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 PeakPickerHiRes_3.tmp -in2 ${DATA_DIR_TOPP}/PeakPickerHiRes_output_lowMem.mzML)
set_tests_properties("TOPP_PeakPickerHiRes_3_out1" PROPERTIES DEPENDS "TOPP_PeakPickerHiRes_3")
# .. chromatogram
add_test("TOPP_PeakPickerHiRes_4" ${TOPP_BIN_PATH}/PeakPickerHiRes -test -ini ${DATA_DIR_TOPP}/PeakPickerHiRes_parameters.ini -in ${DATA_DIR_TOPP}/PeakPickerHiRes_2_input.mzML -out PeakPickerHiRes_4.tmp -processOption lowmemory)
add_test("TOPP_PeakPickerHiRes_4_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 PeakPickerHiRes_4.tmp -in2 ${DATA_DIR_TOPP}/PeakPickerHiRes_2_output.mzML)
set_tests_properties("TOPP_PeakPickerHiRes_4_out1" PROPERTIES DEPENDS "TOPP_PeakPickerHiRes_4")

# peak picking of a smoothed chromatogram (like in PeakPickerMRM; requires disabling of spacing constraints):
add_test("TOPP_PeakPickerHiRes_5" ${TOPP_BIN_PATH}/PeakPickerHiRes -test -in ${DATA_DIR_TOPP}/PeakPickerHiRes_5_input.mzML -out PeakPickerHiRes_5.tmp)
add_test("TOPP_PeakPickerHiRes_5_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 PeakPickerHiRes_5.tmp -in2 ${DATA_DIR_TOPP}/PeakPickerHiRes_5_output.mzML)
set_tests_properties("TOPP_PeakPickerHiRes_5_out1" PROPERTIES DEPENDS "TOPP_PeakPickerHiRes_5")

#------------------------------------------------------------------------------
# UTILS_PeakPickerIterative
add_test("UTILS_PeakPickerIterative_1" ${TOPP_BIN_PATH}/PeakPickerIterative -in ${DATA_DIR_TOPP}/PeakPickerIterative_1_input.mzML -ini ${DATA_DIR_TOPP}/PeakPickerIterative_1.ini -out PeakPickerIterative.mzML.tmp -test)
add_test("UTILS_PeakPickerIterative_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 PeakPickerIterative.mzML.tmp -in2 ${DATA_DIR_TOPP}/PeakPickerIterative_1_output.mzML)
set_tests_properties("UTILS_PeakPickerIterative_1_out1" PROPERTIES DEPENDS "UTILS_PeakPickerIterative_1")
add_test("UTILS_PeakPickerIterative_2" ${TOPP_BIN_PATH}/PeakPickerIterative -in ${DATA_DIR_TOPP}/PeakPickerIterative_2_input.mzML -ini ${DATA_DIR_TOPP}/PeakPickerIterative_2.ini -out PeakPickerIterative_2.mzML.tmp -test)
add_test("UTILS_PeakPickerIterative_2_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 PeakPickerIterative_2.mzML.tmp -in2 ${DATA_DIR_TOPP}/PeakPickerIterative_2_output.mzML)
set_tests_properties("UTILS_PeakPickerIterative_2_out1" PROPERTIES DEPENDS "UTILS_PeakPickerIterative_2")

# Low Memory PeakPicker
add_test("TOPP_LowMemPeakPickerHiRes_1" ${TOPP_BIN_PATH}/LowMemPeakPickerHiRes -test -ini ${DATA_DIR_TOPP}/LowMemPeakPickerHiRes_parameters.ini -in ${DATA_DIR_TOPP}/PeakPickerHiRes_input.mzML -out LowMemPeakPickerHiRes_1.tmp)
add_test("TOPP_LowMemPeakPickerHiRes_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 LowMemPeakPickerHiRes_1.tmp -in2 ${DATA_DIR_TOPP}/LowMemPeakPickerHiRes_output.mzML)
set_tests_properties("TOPP_LowMemPeakPickerHiRes_1_out1" PROPERTIES DEPENDS "TOPP_LowMemPeakPickerHiRes_1")
add_test("TOPP_LowMemPeakPickerHiRes_RandomAccess_1" ${TOPP_BIN_PATH}/LowMemPeakPickerHiRes_RandomAccess -test -ini ${DATA_DIR_TOPP}/LowMemPeakPickerHiRes_RandomAccess_parameters.ini -in ${DATA_DIR_TOPP}/PeakPickerHiRes_RandomAccess_input.mzML -out LowMemPeakPickerHiRes_RandomAccess_1.tmp)
add_test("TOPP_LowMemPeakPickerHiRes_RandomAccess_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 LowMemPeakPickerHiRes_RandomAccess_1.tmp -in2 ${DATA_DIR_TOPP}/LowMemPeakPickerHiRes_RandomAccess_output.mzML)
set_tests_properties("TOPP_LowMemPeakPickerHiRes_RandomAccess_1_out1" PROPERTIES DEPENDS "TOPP_LowMemPeakPickerHiRes_RandomAccess_1")

if(WITH_GUI)
  #------------------------------------------------------------------------------
  # Resampler tests
  add_test("TOPP_Resampler_1" ${TOPP_BIN_PATH}/Resampler -test -in ${DATA_DIR_TOPP}/Resampler_1_input.mzML -out Resampler.mzML -sampling_rate 0.3)
  add_test("TOPP_Resampler_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 Resampler.mzML -in2 ${DATA_DIR_TOPP}/Resampler_1_output.mzML )
  set_tests_properties("TOPP_Resampler_1_out1" PROPERTIES DEPENDS "TOPP_Resampler_1")
endif()

#------------------------------------------------------------------------------
# RNPxlSearch
add_test("TOPP_RNPxlSearch_1" ${TOPP_BIN_PATH}/RNPxlSearch -test -in
${DATA_DIR_TOPP}/RNPxlSearch_1_input.mzML -database ${DATA_DIR_TOPP}/RNPxlSearch_1_input.fasta -out RNPxlSearch_1_output.tmp -ini ${DATA_DIR_TOPP}/RNPxlSearch_1_parameters.ini)
add_test("TOPP_RNPxlSearch_1_out1" ${DIFF} -whitelist "IdentificationRun date"
"db=" -in1 RNPxlSearch_1_output.tmp -in2 ${DATA_DIR_TOPP}/RNPxlSearch_1_output.idXML )
set_tests_properties("TOPP_RNPxlSearch_1_out1" PROPERTIES DEPENDS "TOPP_RNPxlSearch_1")

#------------------------------------------------------------------------------
# RTModel tests
add_test("TOPP_RTModel_1" ${TOPP_BIN_PATH}/RTModel -test -in ${DATA_DIR_TOPP}/RTModel_1_input.idXML -out RTModel_1_output.tmp -ini ${DATA_DIR_TOPP}/RTModel_1_parameters.ini)
add_test("TOPP_RTModel_1_out1" ${DIFF} -in1 RTModel_1_output.tmp -in2 ${DATA_DIR_TOPP}/RTModel_1_output.model )
set_tests_properties("TOPP_RTModel_1_out1" PROPERTIES DEPENDS "TOPP_RTModel_1")
add_test("TOPP_RTModel_2" ${TOPP_BIN_PATH}/RTModel -test -in_positive ${DATA_DIR_TOPP}/RTModel_2_input_positive.idXML -in_negative ${DATA_DIR_TOPP}/RTModel_2_input_negative.idXML -out RTModel_2_output.tmp -ini ${DATA_DIR_TOPP}/RTModel_2_parameters.ini)
add_test("TOPP_RTModel_2_out1" ${DIFF} -in1 RTModel_2_output.tmp -in2 ${DATA_DIR_TOPP}/RTModel_2_output.model )
set_tests_properties("TOPP_RTModel_2_out1" PROPERTIES DEPENDS "TOPP_RTModel_2")
add_test("TOPP_RTModel_2_out2" ${DIFF} -in1 RTModel_2_output.tmp_samples -in2 ${DATA_DIR_TOPP}/RTModel_2_output.samples )
set_tests_properties("TOPP_RTModel_2_out2" PROPERTIES DEPENDS "TOPP_RTModel_2")
add_test("TOPP_RTModel_2_out3" ${DIFF} -in1 RTModel_2_output.tmp_additional_parameters -in2 ${DATA_DIR_TOPP}/RTModel_2_output.additional_parameters )
set_tests_properties("TOPP_RTModel_2_out3" PROPERTIES DEPENDS "TOPP_RTModel_2")
add_test("TOPP_RTModel_3" ${TOPP_BIN_PATH}/RTModel -test -in ${DATA_DIR_TOPP}/RTModel_3_input.idXML -out RTModel_3_output.tmp -ini ${DATA_DIR_TOPP}/RTModel_3_parameters.ini)
add_test("TOPP_RTModel_3_out1" ${DIFF} -in1 RTModel_3_output.tmp -in2 ${DATA_DIR_TOPP}/RTModel_3_output.model )
set_tests_properties("TOPP_RTModel_3_out1" PROPERTIES DEPENDS "TOPP_RTModel_3")
add_test("TOPP_RTModel_3_out2" ${DIFF} -in1 RTModel_3_output.tmp_samples -in2 ${DATA_DIR_TOPP}/RTModel_3_output.samples )
set_tests_properties("TOPP_RTModel_3_out2" PROPERTIES DEPENDS "TOPP_RTModel_3")
add_test("TOPP_RTModel_4" ${TOPP_BIN_PATH}/RTModel -test -in ${DATA_DIR_TOPP}/RTModel_4_input.txt -out RTModel_4_output.tmp -ini ${DATA_DIR_TOPP}/RTModel_4_parameters.ini)
add_test("TOPP_RTModel_4_out1" ${DIFF} -in1 RTModel_4_output.tmp -in2 ${DATA_DIR_TOPP}/RTModel_4_output.model )
set_tests_properties("TOPP_RTModel_4_out1" PROPERTIES DEPENDS "TOPP_RTModel_4")
add_test("TOPP_RTModel_4_out2" ${DIFF} -in1 RTModel_4_output.tmp_samples -in2 ${DATA_DIR_TOPP}/RTModel_4_output.samples )
set_tests_properties("TOPP_RTModel_4_out2" PROPERTIES DEPENDS "TOPP_RTModel_4")

#------------------------------------------------------------------------------
# PTModel tests
add_test("TOPP_PTModel_1" ${TOPP_BIN_PATH}/PTModel -test -in_positive ${DATA_DIR_TOPP}/PTModel_1_input_positive.idXML -in_negative ${DATA_DIR_TOPP}/PTModel_1_input_negative.idXML -out PTModel_1_output.tmp -ini ${DATA_DIR_TOPP}/PTModel_1_parameters.ini)
add_test("TOPP_PTModel_1_out1" ${DIFF} -whitelist "probA" "probB" -in1 PTModel_1_output.tmp -in2 ${DATA_DIR_TOPP}/PTModel_1_output.model )
set_tests_properties("TOPP_PTModel_1_out1" PROPERTIES DEPENDS "TOPP_PTModel_1")
add_test("TOPP_PTModel_1_out2" ${DIFF} -in1 PTModel_1_output.tmp_samples -in2 ${DATA_DIR_TOPP}/PTModel_1_output.samples )
set_tests_properties("TOPP_PTModel_1_out2" PROPERTIES DEPENDS "TOPP_PTModel_1")
add_test("TOPP_PTModel_1_out3" ${DIFF} -in1 PTModel_1_output.tmp_additional_parameters -in2 ${DATA_DIR_TOPP}/PTModel_1_output.additional_parameters )
set_tests_properties("TOPP_PTModel_1_out3" PROPERTIES DEPENDS "TOPP_PTModel_1")

#------------------------------------------------------------------------------
# PTPredict tests
add_test("TOPP_PTPredict_1" ${TOPP_BIN_PATH}/PTPredict -test -in ${DATA_DIR_TOPP}/PTPredict_1_input.idXML -out PTPredict_1_output.tmp -svm_model ${DATA_DIR_TOPP}/PTPredict_1_input.model)
add_test("TOPP_PTPredict_1_out1" ${DIFF} -in1 PTPredict_1_output.tmp -in2 ${DATA_DIR_TOPP}/PTPredict_1_output.idXML )
set_tests_properties("TOPP_PTPredict_1_out1" PROPERTIES DEPENDS "TOPP_PTPredict_1")

#------------------------------------------------------------------------------
# RTPredict tests
add_test("TOPP_RTPredict_1" ${TOPP_BIN_PATH}/RTPredict -test -in_id ${DATA_DIR_TOPP}/RTPredict_1_input.idXML -out_id:file RTPredict_1_output.tmp -total_gradient_time 3000 -svm_model ${DATA_DIR_TOPP}/RTPredict_1_input.model)
add_test("TOPP_RTPredict_1_out1" ${DIFF} -in1 RTPredict_1_output.tmp -in2 ${DATA_DIR_TOPP}/RTPredict_1_output.idXML )
set_tests_properties("TOPP_RTPredict_1_out1" PROPERTIES DEPENDS "TOPP_RTPredict_1")
add_test("TOPP_RTPredict_2" ${TOPP_BIN_PATH}/RTPredict -test -in_id ${DATA_DIR_TOPP}/RTPredict_2_input.idXML -out_id:positive RTPredict_2_output_positive.tmp -out_id:negative RTPredict_2_output_negative.tmp -svm_model ${DATA_DIR_TOPP}/RTPredict_2_input.model)
add_test("TOPP_RTPredict_2_out1" ${DIFF} -in1 RTPredict_2_output_positive.tmp -in2 ${DATA_DIR_TOPP}/RTPredict_2_output_positive.idXML )
set_tests_properties("TOPP_RTPredict_2_out1" PROPERTIES DEPENDS "TOPP_RTPredict_2")
add_test("TOPP_RTPredict_2_out2" ${DIFF} -in1 RTPredict_2_output_negative.tmp -in2 ${DATA_DIR_TOPP}/RTPredict_2_output_negative.idXML )
set_tests_properties("TOPP_RTPredict_2_out2" PROPERTIES DEPENDS "TOPP_RTPredict_2")
add_test("TOPP_RTPredict_3" ${TOPP_BIN_PATH}/RTPredict -test -in_id ${DATA_DIR_TOPP}/RTPredict_3_input.idXML -out_id:file RTPredict_3_output.tmp -total_gradient_time 1 -svm_model ${DATA_DIR_TOPP}/RTPredict_3_input.model)
add_test("TOPP_RTPredict_3_out1" ${DIFF} -whitelist "?xml-stylesheet" "predicted_RT_p_value" -in1 RTPredict_3_output.tmp -in2 ${DATA_DIR_TOPP}/RTPredict_3_output.idXML )
set_tests_properties("TOPP_RTPredict_3_out1" PROPERTIES DEPENDS "TOPP_RTPredict_3")
add_test("TOPP_RTPredict_4" ${TOPP_BIN_PATH}/RTPredict -test -in_text ${DATA_DIR_TOPP}/RTPredict_4_input.txt -out_text:file RTPredict_4_output.tmp -total_gradient_time 1 -svm_model ${DATA_DIR_TOPP}/RTPredict_4_input.model)
add_test("TOPP_RTPredict_4_out1" ${DIFF} -in1 RTPredict_4_output.tmp -in2 ${DATA_DIR_TOPP}/RTPredict_4_output.txt )
set_tests_properties("TOPP_RTPredict_4_out1" PROPERTIES DEPENDS "TOPP_RTPredict_4")
add_test("TOPP_RTPredict_5" ${TOPP_BIN_PATH}/RTPredict -test -in_id ${DATA_DIR_TOPP}/RTPredict_5_input.idXML -out_text:file RTPredict_5_output.tmp -total_gradient_time 1 -svm_model ${DATA_DIR_TOPP}/RTPredict_5_input.model)
add_test("TOPP_RTPredict_5_out1" ${DIFF} -in1 RTPredict_5_output.tmp -in2 ${DATA_DIR_TOPP}/RTPredict_5_output.txt )
set_tests_properties("TOPP_RTPredict_5_out1" PROPERTIES DEPENDS "TOPP_RTPredict_5")

#------------------------------------------------------------------------------
# InspectAdapter tests
add_test("TOPP_InspectAdapter_1" ${TOPP_BIN_PATH}/InspectAdapter -test -ini ${DATA_DIR_TOPP}/InspectAdapter_1_parameters.ini -trie_dbs ${DATA_DIR_TOPP}/Inspect_FASTAFile_test2.trie -in ${DATA_DIR_TOPP}/Inspect.mzXML -dbs ${DATA_DIR_TOPP}/Inspect_FASTAFile_test.fasta -inspect_in -out InspectAdapter_2_output.tmp)
add_test("TOPP_InspectAdapter_1_out1" ${DIFF} -whitelist "spectra," "db," -in1 InspectAdapter_2_output.tmp -in2 ${DATA_DIR_TOPP}/InspectAdapter_2_output.inspect_in )
set_tests_properties("TOPP_InspectAdapter_1_out1" PROPERTIES DEPENDS "TOPP_InspectAdapter_1")

add_test("TOPP_InspectAdapter_2" ${TOPP_BIN_PATH}/InspectAdapter -test -ini ${DATA_DIR_TOPP}/InspectAdapter_1_parameters.ini -trie_dbs ${DATA_DIR_TOPP}/Inspect_FASTAFile_test2.trie -in ${DATA_DIR_TOPP}/Inspect.mzData -dbs ${DATA_DIR_TOPP}/Inspect_FASTAFile_test.fasta -inspect_in -out InspectAdapter_3_output.tmp)
add_test("TOPP_InspectAdapter_2_out1" ${DIFF} -whitelist "spectra," "db," -in1 InspectAdapter_3_output.tmp -in2 ${DATA_DIR_TOPP}/InspectAdapter_2_output.inspect_in )
set_tests_properties("TOPP_InspectAdapter_2_out1" PROPERTIES DEPENDS "TOPP_InspectAdapter_2")

#------------------------------------------------------------------------------
# SeedListGenerator tests
add_test("TOPP_SeedListGenerator_1" ${TOPP_BIN_PATH}/SeedListGenerator -test -in ${DATA_DIR_TOPP}/PepXMLFile_test.mzML -out SeedListGenerator_1_output.tmp)
add_test("TOPP_SeedListGenerator_1_out1" ${DIFF} -in1 SeedListGenerator_1_output.tmp -in2 ${DATA_DIR_TOPP}/SeedListGenerator_1_output.featureXML )
set_tests_properties("TOPP_SeedListGenerator_1_out1" PROPERTIES DEPENDS "TOPP_SeedListGenerator_1")
add_test("TOPP_SeedListGenerator_2" ${TOPP_BIN_PATH}/SeedListGenerator -test -in ${DATA_DIR_TOPP}/IDMapper_1_output.featureXML -out SeedListGenerator_2_output.tmp)
add_test("TOPP_SeedListGenerator_2_out1" ${DIFF} -in1 SeedListGenerator_2_output.tmp -in2 ${DATA_DIR_TOPP}/SeedListGenerator_2_output.featureXML )
set_tests_properties("TOPP_SeedListGenerator_2_out1" PROPERTIES DEPENDS "TOPP_SeedListGenerator_2")
add_test("TOPP_SeedListGenerator_3" ${TOPP_BIN_PATH}/SeedListGenerator -test -in ${DATA_DIR_TOPP}/ConsensusXMLFile_1.consensusXML -out SeedListGenerator_3_output1.tmp SeedListGenerator_3_output2.tmp)
add_test("TOPP_SeedListGenerator_3_out1" ${DIFF} -in1 SeedListGenerator_3_output1.tmp -in2 ${DATA_DIR_TOPP}/SeedListGenerator_3_output1.featureXML )
set_tests_properties("TOPP_SeedListGenerator_3_out1" PROPERTIES DEPENDS "TOPP_SeedListGenerator_3")
add_test("TOPP_SeedListGenerator_3_out2" ${DIFF} -in1 SeedListGenerator_3_output2.tmp -in2 ${DATA_DIR_TOPP}/SeedListGenerator_3_output2.featureXML )
set_tests_properties("TOPP_SeedListGenerator_3_out2" PROPERTIES DEPENDS "TOPP_SeedListGenerator_3")

#------------------------------------------------------------------------------
# SpectraFilter tests
add_test("TOPP_SpectraFilterSqrtMower_1" ${TOPP_BIN_PATH}/SpectraFilterSqrtMower -test -in ${DATA_DIR_TOPP}/SpectraFilterSqrtMower_1_input.mzML -out SpectraFilterSqrtMower.tmp)
add_test("TOPP_SpectraFilterSqrtMower_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 SpectraFilterSqrtMower.tmp -in2 ${DATA_DIR_TOPP}/SpectraFilterSqrtMower_1_output.mzML )
set_tests_properties("TOPP_SpectraFilterSqrtMower_1_out1" PROPERTIES DEPENDS "TOPP_SpectraFilterSqrtMower_1")
add_test("TOPP_SpectraFilterWindowMower_1" ${TOPP_BIN_PATH}/SpectraFilterWindowMower -test -in ${DATA_DIR_TOPP}/SpectraFilterWindowMower_1_input.mzML -out SpectraFilterWindowMower_1.tmp)
add_test("TOPP_SpectraFilterWindowMower_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 SpectraFilterWindowMower_1.tmp -in2 ${DATA_DIR_TOPP}/SpectraFilterWindowMower_1_output.mzML )
set_tests_properties("TOPP_SpectraFilterWindowMower_1_out1" PROPERTIES DEPENDS "TOPP_SpectraFilterWindowMower_1")
add_test("TOPP_SpectraFilterWindowMower_2" ${TOPP_BIN_PATH}/SpectraFilterWindowMower -test -in ${DATA_DIR_TOPP}/SpectraFilterWindowMower_2_input.mzML -out SpectraFilterWindowMower_2.tmp -ini ${DATA_DIR_TOPP}/SpectraFilterWindowMower_2_parameters.ini)
add_test("TOPP_SpectraFilterWindowMower_2_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 SpectraFilterWindowMower_2.tmp -in2 ${DATA_DIR_TOPP}/SpectraFilterWindowMower_2_output.mzML )
set_tests_properties("TOPP_SpectraFilterWindowMower_2_out1" PROPERTIES DEPENDS "TOPP_SpectraFilterWindowMower_2")

#------------------------------------------------------------------------------
# InternalCalibration tests

## idXML file as calibrants
add_test("TOPP_InternalCalibration_1" ${TOPP_BIN_PATH}/InternalCalibration -test -in ${DATA_DIR_TOPP}/InternalCalibration_1_BSA1.mzML -cal:id_in ${DATA_DIR_TOPP}/InternalCalibration_1_BSA1_OMSSA.idXML -out InternalCalibration_1_BSA1_out.mzML.tmp -cal:model_type linear -RANSAC:enabled -RANSAC:iter 500 -RANSAC:threshold 1 -RT_chunking -1 -quality_control:models InternalCalibration_1_models.csv.tmp -quality_control:residuals InternalCalibration_1_residuals.csv.tmp -ms_level 1)
add_test("TOPP_InternalCalibration_1_out1" ${DIFF} -in1 InternalCalibration_1_models.csv.tmp -in2 ${DATA_DIR_TOPP}/InternalCalibration_1_models.csv)
set_tests_properties("TOPP_InternalCalibration_1_out1" PROPERTIES DEPENDS "TOPP_InternalCalibration_1")
add_test("TOPP_InternalCalibration_1_out2" ${DIFF} -in1 InternalCalibration_1_residuals.csv.tmp -in2 ${DATA_DIR_TOPP}/InternalCalibration_1_residuals.csv)
set_tests_properties("TOPP_InternalCalibration_1_out2" PROPERTIES DEPENDS "TOPP_InternalCalibration_1")

## lock masses as calibrants
add_test("TOPP_InternalCalibration_2" ${TOPP_BIN_PATH}/InternalCalibration -test -in ${DATA_DIR_TOPP}/InternalCalibration_2_lockmass.mzML.gz -cal:lock_in ${DATA_DIR_TOPP}/InternalCalibration_2_lock.csv -out InternalCalibration_2_lockmass.mzML.tmp -cal:lock_require_mono -cal:model_type linear -RT_chunking 60 -quality_control:models InternalCalibration_2_models.csv.tmp -quality_control:residuals InternalCalibration_2_residuals.csv.tmp)
add_test("TOPP_InternalCalibration_2_out1" ${DIFF} -in1 InternalCalibration_2_models.csv.tmp -in2 ${DATA_DIR_TOPP}/InternalCalibration_2_models.csv)
set_tests_properties("TOPP_InternalCalibration_2_out1" PROPERTIES DEPENDS "TOPP_InternalCalibration_2")
add_test("TOPP_InternalCalibration_2_out2" ${DIFF} -in1 InternalCalibration_2_residuals.csv.tmp -in2 ${DATA_DIR_TOPP}/InternalCalibration_2_residuals.csv)
set_tests_properties("TOPP_InternalCalibration_2_out2" PROPERTIES DEPENDS "TOPP_InternalCalibration_2")
add_test("TOPP_InternalCalibration_2_convert" ${TOPP_BIN_PATH}/FileConverter -test -in InternalCalibration_2_lockmass.mzML.tmp -no_progress -out InternalCalibration_2_out_lockmass.dta2d.tmp -out_type dta2d)
set_tests_properties("TOPP_InternalCalibration_2_convert" PROPERTIES DEPENDS "TOPP_InternalCalibration_2")
add_test("TOPP_InternalCalibration_2_out3" ${DIFF} -in1 InternalCalibration_2_out_lockmass.dta2d.tmp -in2 ${DATA_DIR_TOPP}/InternalCalibration_2_out_lockmass.dta2d)
set_tests_properties("TOPP_InternalCalibration_2_out3" PROPERTIES DEPENDS "TOPP_InternalCalibration_2_convert")

#------------------------------------------------------------------------------
## ExternalCalibration
add_test("TOPP_ExternalCalibration_1_MS1" ${TOPP_BIN_PATH}/ExternalCalibration -test -in ${DATA_DIR_TOPP}/ExternalCalibration_1_input.mzML -out ExternalCalibration_1_MS1_out.mzML.tmp -offset -5.5 -slope 0.0001 -ms_level 1)
add_test("TOPP_ExternalCalibration_1_out1" ${DIFF} -in1 ExternalCalibration_1_MS1_out.mzML.tmp -in2 ${DATA_DIR_TOPP}/ExternalCalibration_1_MS1_out.mzML)
set_tests_properties("TOPP_ExternalCalibration_1_out1" PROPERTIES DEPENDS "TOPP_ExternalCalibration_1_MS1")

add_test("TOPP_ExternalCalibration_2_MS2" ${TOPP_BIN_PATH}/ExternalCalibration -test -in ${DATA_DIR_TOPP}/ExternalCalibration_1_input.mzML -out ExternalCalibration_2_MS2_out.mzML.tmp -offset -5.5 -slope 0.0001 -ms_level 2)
add_test("TOPP_ExternalCalibration_2_out1" ${DIFF} -in1 ExternalCalibration_2_MS2_out.mzML.tmp -in2 ${DATA_DIR_TOPP}/ExternalCalibration_2_MS2_out.mzML)
set_tests_properties("TOPP_ExternalCalibration_2_out1" PROPERTIES DEPENDS "TOPP_ExternalCalibration_2_MS2")

#------------------------------------------------------------------------------
# TMTAnalyzer
add_test("TOPP_TMTAnalyzer_1" ${TOPP_BIN_PATH}/TMTAnalyzer -test -ini ${DATA_DIR_TOPP}/TMTAnalyzer.ini -in ${DATA_DIR_TOPP}/TMTAnalyzer_input_1.mzML -out TMTAnalyzer_1.tmp)
add_test("TOPP_TMTAnalyzer_1_out1" ${DIFF} -whitelist "id=" "<map" "?xml-stylesheet" -in1 TMTAnalyzer_1.tmp -in2 ${DATA_DIR_TOPP}/TMTAnalyzer_output_1.consensusXML )
set_tests_properties("TOPP_TMTAnalyzer_1_out1" PROPERTIES DEPENDS "TOPP_TMTAnalyzer_1")

#------------------------------------------------------------------------------
# TOFCalibration tests
add_test("TOPP_TOFCalibration_1" ${TOPP_BIN_PATH}/TOFCalibration -test -in ${DATA_DIR_TOPP}/TOFCalibration_1_input.mzML -out TOFCalibration.tmp -ref_masses ${DATA_DIR_TOPP}/TOFCalibration_ref_masses -ini ${DATA_DIR_TOPP}/TOFCalibration_parameters.ini -tof_const ${DATA_DIR_TOPP}/TOFCalibration_const -ext_calibrants ${DATA_DIR_TOPP}/TOFCalibration_1_calibrants.mzML)
#add_test("TOPP_TOFCalibration_1_out1" ${DIFF} -in1 TOFCalibration.tmp -in2 ${DATA_DIR_TOPP}/TOFCalibration_1_output.mzML )
#set_tests_properties("TOPP_TOFCalibration_1_out1" PROPERTIES DEPENDS "TOPP_TOFCalibration_1")
add_test("TOPP_TOFCalibration_2" ${TOPP_BIN_PATH}/TOFCalibration -test -in ${DATA_DIR_TOPP}/TOFCalibration_2_input.mzML -out TOFCalibration.tmp -ref_masses ${DATA_DIR_TOPP}/TOFCalibration_ref_masses -ini ${DATA_DIR_TOPP}/TOFCalibration_parameters.ini -tof_const ${DATA_DIR_TOPP}/TOFCalibration_const -ext_calibrants ${DATA_DIR_TOPP}/TOFCalibration_2_calibrants.mzML -peak_data)
#add_test("TOPP_TOFCalibration_2_out1" ${DIFF} -in1 TOFCalibration.tmp -in2 ${DATA_DIR_TOPP}/TOFCalibration_2_output.mzML )
#set_tests_properties("TOPP_TOFCalibration_2_out1" PROPERTIES DEPENDS "TOPP_TOFCalibration_2")

#------------------------------------------------------------------------------
# TextExporter tests
add_test("TOPP_TextExporter_1" ${TOPP_BIN_PATH}/TextExporter -test -in ${DATA_DIR_TOPP}/TextExporter_1_input.featureXML -no_progress -out TextExporter_1_output.tmp)
add_test("TOPP_TextExporter_1_out1" ${DIFF} -in1 TextExporter_1_output.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_1_output.txt )
set_tests_properties("TOPP_TextExporter_1_out1" PROPERTIES DEPENDS "TOPP_TextExporter_1")

add_test("TOPP_TextExporter_2" ${TOPP_BIN_PATH}/TextExporter -test -in ${DATA_DIR_TOPP}/TextExporter_2_input.consensusXML -ini ${DATA_DIR_TOPP}/TextExporter_2_parameters.ini)
add_test("TOPP_TextExporter_2_out1" ${DIFF} -whitelist "extracted from" -in1 TextExporter_2_consensus_centroids.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_2_consensus_centroids.wsv )
set_tests_properties("TOPP_TextExporter_2_out1" PROPERTIES DEPENDS "TOPP_TextExporter_2")
add_test("TOPP_TextExporter_2_out2" ${DIFF} -whitelist "extracted from" -in1 TextExporter_2_consensus_elements.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_2_consensus_elements.wsv )
set_tests_properties("TOPP_TextExporter_2_out2" PROPERTIES DEPENDS "TOPP_TextExporter_2")
add_test("TOPP_TextExporter_2_out3" ${DIFF} -whitelist "extracted from" -in1 TextExporter_2_consensus_features.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_2_consensus_features.wsv )
set_tests_properties("TOPP_TextExporter_2_out3" PROPERTIES DEPENDS "TOPP_TextExporter_2")
add_test("TOPP_TextExporter_2_out4" ${DIFF} -whitelist "extracted from" -in1 TextExporter_2_consensus_tsv.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_2_consensus_tsv.tsv )
set_tests_properties("TOPP_TextExporter_2_out4" PROPERTIES DEPENDS "TOPP_TextExporter_2")

add_test("TOPP_TextExporter_3" ${TOPP_BIN_PATH}/TextExporter -test -in ${DATA_DIR_TOPP}/TextExporter_3_input.idXML -no_progress -out TextExporter_3_output.tmp)
add_test("TOPP_TextExporter_3_out1" ${DIFF} -in1 TextExporter_3_output.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_3_output.txt )
set_tests_properties("TOPP_TextExporter_3_out1" PROPERTIES DEPENDS "TOPP_TextExporter_3")

add_test("TOPP_TextExporter_4" ${TOPP_BIN_PATH}/TextExporter -test -in ${DATA_DIR_TOPP}/TextExporter_3_input.idXML -no_progress -out TextExporter_4_output.tmp -id:proteins_only)
add_test("TOPP_TextExporter_4_out1" ${DIFF} -in1 TextExporter_4_output.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_4_output_proteins.txt )
set_tests_properties("TOPP_TextExporter_4_out1" PROPERTIES DEPENDS "TOPP_TextExporter_4")

add_test("TOPP_TextExporter_5" ${TOPP_BIN_PATH}/TextExporter -test -in ${DATA_DIR_TOPP}/TextExporter_5_input.idXML -no_progress -out TextExporter_5_output.tmp -id:peptides_only -id:first_dim_rt)
add_test("TOPP_TextExporter_5_out1" ${DIFF} -in1 TextExporter_5_output.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_5_output_peptides.txt )
set_tests_properties("TOPP_TextExporter_5_out1" PROPERTIES DEPENDS "TOPP_TextExporter_5")

add_test("TOPP_TextExporter_6" ${TOPP_BIN_PATH}/TextExporter -test -in ${DATA_DIR_TOPP}/TextExporter_6_input.featureXML -no_progress -out TextExporter_6_output.tmp -no_ids)
add_test("TOPP_TextExporter_6_out1" ${DIFF} -in1 TextExporter_6_output.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_6_output.txt )
set_tests_properties("TOPP_TextExporter_6_out1" PROPERTIES DEPENDS "TOPP_TextExporter_6")

add_test("TOPP_TextExporter_7" ${TOPP_BIN_PATH}/TextExporter -test -in ${DATA_DIR_TOPP}/TextExporter_7_input.consensusXML -ini ${DATA_DIR_TOPP}/TextExporter_7_parameters.ini)
add_test("TOPP_TextExporter_7_out1" ${DIFF} -whitelist "extracted from" -in1 TextExporter_7_consensus_tsv.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_7_consensus_tsv.tsv )
set_tests_properties("TOPP_TextExporter_7_out1" PROPERTIES DEPENDS "TOPP_TextExporter_7")

add_test("TOPP_TextExporter_8" ${TOPP_BIN_PATH}/TextExporter -test -in ${DATA_DIR_TOPP}/TextExporter_1_input.featureXML -no_ids -no_progress -out TextExporter_8_output.tmp)
add_test("TOPP_TextExporter_8_out1" ${DIFF} -in1 TextExporter_8_output.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_8_output.txt )
set_tests_properties("TOPP_TextExporter_8_out1" PROPERTIES DEPENDS "TOPP_TextExporter_8")

add_test("TOPP_TextExporter_9" ${TOPP_BIN_PATH}/TextExporter -test -in ${DATA_DIR_TOPP}/TextExporter_9_input.idXML -no_progress -out TextExporter_9_output.tmp -id:add_metavalues 0 -id:add_hit_metavalues 0)
add_test("TOPP_TextExporter_9_out1" ${DIFF} -in1 TextExporter_9_output.tmp -in2 ${DATA_DIR_TOPP}/TextExporter_9_output.txt )
set_tests_properties("TOPP_TextExporter_9_out1" PROPERTIES DEPENDS "TOPP_TextExporter_9")

#------------------------------------------------------------------------------
# FeatureLinker tests
# "labeled" algorithm:
add_test("TOPP_FeatureLinkerLabeled_1" ${TOPP_BIN_PATH}/FeatureLinkerLabeled -test -ini ${DATA_DIR_TOPP}/FeatureLinkerLabeled_1_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerLabeled_1_input.featureXML -out FeatureLinkerLabeled_1_output.tmp)
add_test("TOPP_FeatureLinkerLabeled_1_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerLabeled_1_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinkerLabeled_1_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerLabeled_1_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerLabeled_1")
add_test("TOPP_FeatureLinkerLabeled_2" ${TOPP_BIN_PATH}/FeatureLinkerLabeled -test -ini ${DATA_DIR_TOPP}/FeatureLinkerLabeled_2_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerLabeled_2_input.featureXML -out FeatureLinkerLabeled_2_output.tmp)
add_test("TOPP_FeatureLinkerLabeled_2_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerLabeled_2_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinkerLabeled_2_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerLabeled_2_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerLabeled_2")
# "unlabelled" algorithm:
add_test("TOPP_FeatureLinkerUnlabeled_1" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeled -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input2.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input3.featureXML -out FeatureLinkerUnlabeled_1_output.tmp)
add_test("TOPP_FeatureLinkerUnlabeled_1_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerUnlabeled_1_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerUnlabeled_1_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeled_1")
add_test("TOPP_FeatureLinkerUnlabeled_2" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeled -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_2_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_2_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_2_input2.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_2_input3.featureXML -out FeatureLinkerUnlabeled_2_output.tmp)
add_test("TOPP_FeatureLinkerUnlabeled_2_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerUnlabeled_2_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_2_output.consensusXML)
set_tests_properties("TOPP_FeatureLinkerUnlabeled_2_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeled_2")
add_test("TOPP_FeatureLinkerUnlabeled_3" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeled -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_3_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_3_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_3_input2.featureXML -out FeatureLinkerUnlabeled_3_output.tmp)
add_test("TOPP_FeatureLinkerUnlabeled_3_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerUnlabeled_3_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_3_output.consensusXML)
set_tests_properties("TOPP_FeatureLinkerUnlabeled_3_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeled_3")
# "unlabeled" algorithm with "use_identifications" option:
add_test("TOPP_FeatureLinkerUnlabeled_4" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeled -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_4_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input2.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input3.featureXML -out FeatureLinkerUnlabeled_4_output.tmp)
add_test("TOPP_FeatureLinkerUnlabeled_4_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerUnlabeled_4_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_4_output.consensusXML)
set_tests_properties("TOPP_FeatureLinkerUnlabeled_4_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeled_4")
# "unlabeled_qt" algorithm:
add_test("TOPP_FeatureLinkerUnlabeledQT_1" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeledQT -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_1_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input2.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input3.featureXML -out FeatureLinkerUnlabeledQT_1_output.tmp)
add_test("TOPP_FeatureLinkerUnlabeledQT_1_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerUnlabeledQT_1_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_1_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerUnlabeledQT_1_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeledQT_1")
# "unlabeled_qt" algorithm with consensusXML input and "keep_subelements" flag:
add_test("TOPP_FeatureLinkerUnlabeledQT_2" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeledQT -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_2_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_1_output.consensusXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_1_output.consensusXML -out FeatureLinkerUnlabeledQT_2_output.tmp)
add_test("TOPP_FeatureLinkerUnlabeledQT_2_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerUnlabeledQT_2_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_2_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerUnlabeledQT_2_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeledQT_2")
add_test("TOPP_FeatureLinkerUnlabeledQT_3" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeledQT -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_3_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_3_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_3_input2.featureXML -out FeatureLinkerUnlabeledQT_3_output.tmp)
add_test("TOPP_FeatureLinkerUnlabeledQT_3_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerUnlabeledQT_3_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_3_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerUnlabeledQT_3_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeledQT_3")
# use_identifications turned on
add_test("TOPP_FeatureLinkerUnlabeledQT_4" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeledQT -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_1_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input2.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input3.featureXML -out FeatureLinkerUnlabeledQT_4_output.tmp -algorithm:use_identifications)
add_test("TOPP_FeatureLinkerUnlabeledQT_4_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerUnlabeledQT_4_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_4_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerUnlabeledQT_4_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeledQT_4")
add_test("TOPP_FeatureLinkerUnlabeledQT_5" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeledQT -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_1_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_5_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_5_input2.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_5_input3.featureXML -out FeatureLinkerUnlabeledQT_5_output.tmp -algorithm:distance_RT:max_difference 200)
add_test("TOPP_FeatureLinkerUnlabeledQT_5_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerUnlabeledQT_5_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_5_output.consensusXML)
set_tests_properties("TOPP_FeatureLinkerUnlabeledQT_5_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeledQT_5")
add_test("TOPP_FeatureLinkerUnlabeledQT_6" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeledQT -test -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_5_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_5_input2.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_5_input3.featureXML -out FeatureLinkerUnlabeledQT_6_output.tmp -algorithm:use_identifications -algorithm:distance_RT:max_difference 200)
add_test("TOPP_FeatureLinkerUnlabeledQT_6_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerUnlabeledQT_6_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_6_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerUnlabeledQT_6_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeledQT_6")
# FeatureLinkerUnlabeledKD
add_test("TOPP_FeatureLinkerUnlabeledKD_1" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeledKD -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_1_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input2.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input3.featureXML -out FeatureLinkerUnlabeledKD_1_output.tmp)
add_test("TOPP_FeatureLinkerUnlabeledKD_1_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerUnlabeledKD_1_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_1_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerUnlabeledKD_1_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeledKD_1")
add_test("TOPP_FeatureLinkerUnlabeledKD_2" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeledKD -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_2_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_1_output.consensusXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_1_output.consensusXML -out FeatureLinkerUnlabeledKD_2_output.tmp)
add_test("TOPP_FeatureLinkerUnlabeledKD_2_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerUnlabeledKD_2_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_2_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerUnlabeledKD_2_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeledKD_2")
add_test("TOPP_FeatureLinkerUnlabeledKD_3" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeledKD -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_3_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_3_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_3_input2.featureXML -out FeatureLinkerUnlabeledKD_3_output.tmp)
add_test("TOPP_FeatureLinkerUnlabeledKD_3_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerUnlabeledKD_3_output.tmp -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_3_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerUnlabeledKD_3_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeledKD_3")

#------------------------------------------------------------------------------
# IDMapper tests
add_test("TOPP_IDMapper_1" ${TOPP_BIN_PATH}/IDMapper -test -in ${DATA_DIR_TOPP}/IDMapper_1_input.featureXML -out IDMapper_1_output.tmp -id ${DATA_DIR_TOPP}/IDMapper_1_input.idXML -mz_measure Da -mz_tolerance 1 -ignore_charge)
add_test("TOPP_IDMapper_1_out1" ${DIFF} -in1 IDMapper_1_output.tmp -in2 ${DATA_DIR_TOPP}/IDMapper_1_output.featureXML )
set_tests_properties("TOPP_IDMapper_1_out1" PROPERTIES DEPENDS "TOPP_IDMapper_1")
add_test("TOPP_IDMapper_2" ${TOPP_BIN_PATH}/IDMapper -test -in ${DATA_DIR_TOPP}/IDMapper_2_input.consensusXML -out IDMapper_2_output.tmp -id ${DATA_DIR_TOPP}/IDMapper_2_input.idXML  -mz_measure Da -mz_tolerance 1)
add_test("TOPP_IDMapper_2_out1" ${DIFF} -in1 IDMapper_2_output.tmp -in2 ${DATA_DIR_TOPP}/IDMapper_2_output.consensusXML )
set_tests_properties("TOPP_IDMapper_2_out1" PROPERTIES DEPENDS "TOPP_IDMapper_2")
add_test("TOPP_IDMapper_3" ${TOPP_BIN_PATH}/IDMapper -test -in ${DATA_DIR_TOPP}/IDMapper_3_input.featureXML -out IDMapper_3_output.tmp -id ${DATA_DIR_TOPP}/IDMapper_3_input.idXML  -mz_measure ppm -mz_tolerance 3 -rt_tolerance 4 -ignore_charge)
add_test("TOPP_IDMapper_3_out1" ${DIFF} -in1 IDMapper_3_output.tmp -in2 ${DATA_DIR_TOPP}/IDMapper_3_output.featureXML )
set_tests_properties("TOPP_IDMapper_3_out1" PROPERTIES DEPENDS "TOPP_IDMapper_3")
add_test("TOPP_IDMapper_4" ${TOPP_BIN_PATH}/IDMapper -test -in
${DATA_DIR_TOPP}/IDMapper_4_input.featureXML -out IDMapper_4_output.tmp -id
${DATA_DIR_TOPP}/IDMapper_4_input.idXML -spectra:in
${DATA_DIR_TOPP}/IDMapper_4_input.mzML -mz_measure ppm -mz_tolerance 20
-rt_tolerance 10)
add_test("TOPP_IDMapper_4_out1" ${DIFF} -in1 IDMapper_4_output.tmp -in2 ${DATA_DIR_TOPP}/IDMapper_4_output.featureXML )
set_tests_properties("TOPP_IDMapper_4_out1" PROPERTIES DEPENDS "TOPP_IDMapper_4")
add_test("TOPP_IDMapper_5" ${TOPP_BIN_PATH}/IDMapper -test -in
${DATA_DIR_TOPP}/IDMapper_5_input.featureXML -out IDMapper_5_output.tmp -id
${DATA_DIR_TOPP}/IDMapper_5_input.idXML -spectra:in
${DATA_DIR_TOPP}/IDMapper_5_input.mzML -mz_measure ppm -mz_tolerance 10
-rt_tolerance 20)
add_test("TOPP_IDMapper_5_out1" ${DIFF} -in1 IDMapper_5_output.tmp -in2 ${DATA_DIR_TOPP}/IDMapper_5_output.featureXML )
set_tests_properties("TOPP_IDMapper_5_out1" PROPERTIES DEPENDS "TOPP_IDMapper_5")

#------------------------------------------------------------------------------
# IDRipper tests
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
add_test("TOPP_IDRipper_3_prep_out1" ${DIFF} -whitelist "?xml-stylesheet" "file_origin" -in1 IDRipper_3_output.tmp -in2 ${DATA_DIR_TOPP}/IDRipper_3_output.idXML)
set_tests_properties("TOPP_IDRipper_3_prep_out1" PROPERTIES DEPENDS "TOPP_IDRipper_3_prep")
## create an output directory during configure time in the binary tree
set(TMP_RIP_PATH ${PROJECT_BINARY_DIR}/tmp_path/) ## ${PROJECT_BINARY_DIR} already includes '/src/tests/topp/'
file(MAKE_DIRECTORY ${TMP_RIP_PATH})
add_test("TOPP_IDRipper_3" ${TOPP_BIN_PATH}/IDRipper -test -in ${DATA_DIR_TOPP}/IDRipper_3_output.idXML -out_path ${TMP_RIP_PATH})
set_tests_properties("TOPP_IDRipper_3" PROPERTIES DEPENDS "TOPP_IDRipper_3_prep")
add_test("TOPP_IDRipper_3_out1" ${DIFF} -in1 ${TMP_RIP_PATH}/IDRipper_3_input1.idXML -in2 ${DATA_DIR_TOPP}/IDRipper_3_input1.idXML)
set_tests_properties("TOPP_IDRipper_3_out1" PROPERTIES DEPENDS "TOPP_IDRipper_3")
add_test("TOPP_IDRipper_3_out2" ${DIFF} -in1 ${TMP_RIP_PATH}/IDRipper_3_input2.idXML -in2 ${DATA_DIR_TOPP}/IDRipper_3_input2.idXML)
set_tests_properties("TOPP_IDRipper_3_out2" PROPERTIES DEPENDS "TOPP_IDRipper_3")

#------------------------------------------------------------------------------
# ConsensusID tests
# idXML input, "PEPMatrix" algorithm using PEP30MS matrix (example from the JPR paper):
add_test("TOPP_ConsensusID_1" ${TOPP_BIN_PATH}/ConsensusID -test -in ${DATA_DIR_TOPP}/ConsensusID_1_input.idXML -out ConsensusID_1_output.tmp -algorithm PEPMatrix -PEPMatrix:matrix PAM30MS)
add_test("TOPP_ConsensusID_1_out1" ${DIFF} -whitelist "?xml-stylesheet" "IdentificationRun date" -in1 ConsensusID_1_output.tmp -in2 ${DATA_DIR_TOPP}/ConsensusID_1_output.idXML )
set_tests_properties("TOPP_ConsensusID_1_out1" PROPERTIES DEPENDS "TOPP_ConsensusID_1")
# featureXML input, "average" algorithm:
add_test("TOPP_ConsensusID_2" ${TOPP_BIN_PATH}/ConsensusID -test -in ${DATA_DIR_TOPP}/ConsensusID_2_input.featureXML -out ConsensusID_2_output.tmp -algorithm average)
add_test("TOPP_ConsensusID_2_out1" ${DIFF} -whitelist "IdentificationRun id" -in1 ConsensusID_2_output.tmp -in2 ${DATA_DIR_TOPP}/ConsensusID_2_output.featureXML )
set_tests_properties("TOPP_ConsensusID_2_out1" PROPERTIES DEPENDS "TOPP_ConsensusID_2")
# consensusXML input, "best" algorithm:
add_test("TOPP_ConsensusID_3" ${TOPP_BIN_PATH}/ConsensusID -test -in ${DATA_DIR_TOPP}/ConsensusID_3_input.consensusXML -out ConsensusID_3_output.tmp -algorithm best)
add_test("TOPP_ConsensusID_3_out1" ${DIFF} -whitelist "?xml-stylesheet" "IdentificationRun id" -in1 ConsensusID_3_output.tmp -in2 ${DATA_DIR_TOPP}/ConsensusID_3_output.consensusXML )
set_tests_properties("TOPP_ConsensusID_3_out1" PROPERTIES DEPENDS "TOPP_ConsensusID_3")
# idXML input, "PEPMatrix" algorithm using identity matrix, limited number of considered hits:
add_test("TOPP_ConsensusID_4" ${TOPP_BIN_PATH}/ConsensusID -test -in ${DATA_DIR_TOPP}/ConsensusID_1_input.idXML -out ConsensusID_4_output.tmp -algorithm PEPMatrix -PEPMatrix:matrix identity -filter:considered_hits 6)
add_test("TOPP_ConsensusID_4_out1" ${DIFF} -whitelist "?xml-stylesheet" "IdentificationRun date" -in1 ConsensusID_4_output.tmp -in2 ${DATA_DIR_TOPP}/ConsensusID_4_output.idXML )
set_tests_properties("TOPP_ConsensusID_4_out1" PROPERTIES DEPENDS "TOPP_ConsensusID_4")
# idXML input, "PEPIons" algorithm:
add_test("TOPP_ConsensusID_5" ${TOPP_BIN_PATH}/ConsensusID -test -in ${DATA_DIR_TOPP}/ConsensusID_1_input.idXML -out ConsensusID_5_output.tmp -algorithm PEPIons)
add_test("TOPP_ConsensusID_5_out1" ${DIFF} -whitelist "?xml-stylesheet" "IdentificationRun date" -in1 ConsensusID_5_output.tmp -in2 ${DATA_DIR_TOPP}/ConsensusID_5_output.idXML )
set_tests_properties("TOPP_ConsensusID_5_out1" PROPERTIES DEPENDS "TOPP_ConsensusID_5")
# idXML input, "best" algorithm, filtering:
add_test("TOPP_ConsensusID_6" ${TOPP_BIN_PATH}/ConsensusID -test -in ${DATA_DIR_TOPP}/ConsensusID_1_input.idXML -out ConsensusID_6_output.tmp -algorithm best -filter:min_support 0.5)
add_test("TOPP_ConsensusID_6_out1" ${DIFF} -whitelist "?xml-stylesheet" "IdentificationRun date" -in1 ConsensusID_6_output.tmp -in2 ${DATA_DIR_TOPP}/ConsensusID_6_output.idXML )
set_tests_properties("TOPP_ConsensusID_6_out1" PROPERTIES DEPENDS "TOPP_ConsensusID_6")

#------------------------------------------------------------------------------
# PrecursorIonSelector tests
add_test("TOPP_PrecursorIonSelector_1" ${TOPP_BIN_PATH}/PrecursorIonSelector -test -in ${DATA_DIR_TOPP}/PrecursorIonSelector_features.featureXML -sim_results PrecursorIonSelector_1_output.tmp -ini ${DATA_DIR_TOPP}/PrecursorIonSelector_1_parameters.ini -db_path ${DATA_DIR_TOPP}/PrecursorIonSelector_db.fasta -ids ${DATA_DIR_TOPP}/PrecursorIonSelector_ids.idXML)
add_test("TOPP_PrecursorIonSelector_1_out1" ${DIFF} -in1 PrecursorIonSelector_1_output.tmp -in2 ${DATA_DIR_TOPP}/PrecursorIonSelector_1_output.txt )
set_tests_properties("TOPP_PrecursorIonSelector_1_out1" PROPERTIES DEPENDS "TOPP_PrecursorIonSelector_1")
add_test("TOPP_PrecursorIonSelector_2" ${TOPP_BIN_PATH}/PrecursorIonSelector -test -in ${DATA_DIR_TOPP}/PrecursorIonSelector_features.featureXML -out PrecursorIonSelector_2_output.tmp -ini ${DATA_DIR_TOPP}/PrecursorIonSelector_2_parameters.ini -db_path ${DATA_DIR_TOPP}/PrecursorIonSelector_db.fasta -ids ${DATA_DIR_TOPP}/PrecursorIonSelector_ids.idXML)
add_test("TOPP_PrecursorIonSelector_2_out1" ${DIFF} -in1 PrecursorIonSelector_2_output.tmp -in2 ${DATA_DIR_TOPP}/PrecursorIonSelector_2_output.featureXML )
set_tests_properties("TOPP_PrecursorIonSelector_2_out1" PROPERTIES DEPENDS "TOPP_PrecursorIonSelector_2")

#------------------------------------------------------------------------------
# CompNovo tests
add_test("TOPP_CompNovo_1" ${TOPP_BIN_PATH}/CompNovo -test -in ${DATA_DIR_TOPP}/CompNovo_1_input.mzML -ini ${DATA_DIR_TOPP}/CompNovo_1.ini -out CompNovo_1_output.tmp)
add_test("TOPP_CompNovo_1_out1" ${DIFF} -whitelist "?xml-stylesheet" "IdentificationRun date" -in1 CompNovo_1_output.tmp -in2 ${DATA_DIR_TOPP}/CompNovo_1_output.idXML )
set_tests_properties("TOPP_CompNovo_1_out1" PROPERTIES DEPENDS "TOPP_CompNovo_1")
add_test("TOPP_CompNovoCID_1" ${TOPP_BIN_PATH}/CompNovoCID -test -in ${DATA_DIR_TOPP}/CompNovoCID_1_input.mzML -ini ${DATA_DIR_TOPP}/CompNovoCID_1.ini -out CompNovoCID_1_output.tmp)
add_test("TOPP_CompNovoCID_1_out1" ${DIFF} -whitelist "?xml-stylesheet" "IdentificationRun date" -in1 CompNovoCID_1_output.tmp -in2 ${DATA_DIR_TOPP}/CompNovoCID_1_output.idXML )
set_tests_properties("TOPP_CompNovoCID_1_out1" PROPERTIES DEPENDS "TOPP_CompNovoCID_1")

#------------------------------------------------------------------------------
# PrecursorMassCorrector tests
add_test("TOPP_PrecursorMassCorrector_1" ${TOPP_BIN_PATH}/PrecursorMassCorrector -test -in ${DATA_DIR_TOPP}/PrecursorMassCorrector_1_input.mzML -ini ${DATA_DIR_TOPP}/PrecursorMassCorrector_1.ini -out PrecursorMassCorrector_1_output.tmp)
add_test("TOPP_PrecursorMassCorrector_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} "location=" -in1 PrecursorMassCorrector_1_output.tmp -in2 ${DATA_DIR_TOPP}/PrecursorMassCorrector_1_output.mzML )
set_tests_properties("TOPP_PrecursorMassCorrector_1_out1" PROPERTIES DEPENDS "TOPP_PrecursorMassCorrector_1")

#------------------------------------------------------------------------------
# ProteinInference tests
add_test("TOPP_ProteinInference_1" ${TOPP_BIN_PATH}/ProteinInference -test -in ${DATA_DIR_TOPP}/ProteinInference.idXML -out ProteinInference_output.tmp)
add_test("TOPP_ProteinInference_1_out1" ${DIFF} -in1 ProteinInference_output.tmp -in2 ${DATA_DIR_TOPP}/ProteinInference_proteins.idXML )
set_tests_properties("TOPP_ProteinInference_1_out1" PROPERTIES DEPENDS "TOPP_ProteinInference_1")

#------------------------------------------------------------------------------
# FalseDiscoveryRate_tests (1: default param, 2: with runs separated, 3: with charge variants individual, 4: with charge variants and no charge 1 decoys -> q_value 0)
add_test("TOPP_FalseDiscoveryRate_1" ${TOPP_BIN_PATH}/FalseDiscoveryRate -test -in ${DATA_DIR_TOPP}/FalseDiscoveryRate_OMSSA.idXML -out FalseDiscoveryRate_output_1.tmp -peptides_only)
add_test("TOPP_FalseDiscoveryRate_1_out1" ${DIFF} -whitelist "?xml-stylesheet" -in1 FalseDiscoveryRate_output_1.tmp -in2 ${DATA_DIR_TOPP}/FalseDiscoveryRate_output_1.idXML)
set_tests_properties("TOPP_FalseDiscoveryRate_1_out1" PROPERTIES DEPENDS "TOPP_FalseDiscoveryRate_1")
add_test("TOPP_FalseDiscoveryRate_2" ${TOPP_BIN_PATH}/FalseDiscoveryRate -test -in ${DATA_DIR_TOPP}/FalseDiscoveryRate_OMSSA.idXML -out FalseDiscoveryRate_output_2.tmp -peptides_only -algorithm:treat_runs_separately)
add_test("TOPP_FalseDiscoveryRate_2_out1" ${DIFF} -whitelist "?xml-stylesheet" -in1 FalseDiscoveryRate_output_2.tmp -in2 ${DATA_DIR_TOPP}/FalseDiscoveryRate_output_2.idXML)
set_tests_properties("TOPP_FalseDiscoveryRate_2_out1" PROPERTIES DEPENDS "TOPP_FalseDiscoveryRate_2")
add_test("TOPP_FalseDiscoveryRate_3" ${TOPP_BIN_PATH}/FalseDiscoveryRate -test -in ${DATA_DIR_TOPP}/FalseDiscoveryRate_OMSSA.idXML -out FalseDiscoveryRate_output_3.tmp -peptides_only -algorithm:split_charge_variants)
add_test("TOPP_FalseDiscoveryRate_3_out1" ${DIFF} -whitelist "?xml-stylesheet" -in1 FalseDiscoveryRate_output_3.tmp -in2 ${DATA_DIR_TOPP}/FalseDiscoveryRate_output_3.idXML)
set_tests_properties("TOPP_FalseDiscoveryRate_3_out1" PROPERTIES DEPENDS "TOPP_FalseDiscoveryRate_3")
add_test("TOPP_FalseDiscoveryRate_4" ${TOPP_BIN_PATH}/FalseDiscoveryRate -test -in ${DATA_DIR_TOPP}/FalseDiscoveryRate_OMSSA_4.idXML -out FalseDiscoveryRate_output_4.tmp -peptides_only -algorithm:split_charge_variants)
add_test("TOPP_FalseDiscoveryRate_4_out1" ${DIFF} -whitelist "?xml-stylesheet" -in1 FalseDiscoveryRate_output_4.tmp -in2 ${DATA_DIR_TOPP}/FalseDiscoveryRate_output_4.idXML)
set_tests_properties("TOPP_FalseDiscoveryRate_4_out1" PROPERTIES DEPENDS "TOPP_FalseDiscoveryRate_4")
# proteins only (input is derived from "THIRDPARTY/FidoAdapter_1_output.idXML"):
add_test("TOPP_FalseDiscoveryRate_5" ${TOPP_BIN_PATH}/FalseDiscoveryRate -test -in ${DATA_DIR_TOPP}/FalseDiscoveryRate_5_input.idXML -out FalseDiscoveryRate_5_output.tmp -proteins_only)
add_test("TOPP_FalseDiscoveryRate_5_out1" ${DIFF} -whitelist "?xml-stylesheet" -in1 FalseDiscoveryRate_5_output.tmp -in2 ${DATA_DIR_TOPP}/FalseDiscoveryRate_5_output.idXML )
set_tests_properties("TOPP_FalseDiscoveryRate_5_out1" PROPERTIES DEPENDS "TOPP_FalseDiscoveryRate_5")

#------------------------------------------------------------------------------
# IDPosteriorErrorProbability
# Mascot input:
add_test("TOPP_IDPosteriorErrorProbability_1" ${TOPP_BIN_PATH}/IDPosteriorErrorProbability -test -in ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_Mascot_input.idXML -out 	IDPosteriorErrorProbability_output_1.tmp)
add_test("TOPP_IDPosteriorErrorProbability_1_out1" ${DIFF} -in1 IDPosteriorErrorProbability_output_1.tmp -in2 ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_Mascot_output.idXML)
set_tests_properties("TOPP_IDPosteriorErrorProbability_1_out1" PROPERTIES DEPENDS "TOPP_IDPosteriorErrorProbability_1")
# X! Tandem input:
add_test("TOPP_IDPosteriorErrorProbability_2" ${TOPP_BIN_PATH}/IDPosteriorErrorProbability -test -in ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_XTandem_input.idXML -out IDPosteriorErrorProbability_output_2.tmp)
add_test("TOPP_IDPosteriorErrorProbability_2_out1" ${DIFF} -in1 IDPosteriorErrorProbability_output_2.tmp -in2 ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_XTandem_output.idXML)
set_tests_properties("TOPP_IDPosteriorErrorProbability_2_out1" PROPERTIES DEPENDS "TOPP_IDPosteriorErrorProbability_2")
# OMSSA input:
add_test("TOPP_IDPosteriorErrorProbability_3" ${TOPP_BIN_PATH}/IDPosteriorErrorProbability -test -in ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_OMSSA_input.idXML -out IDPosteriorErrorProbability_output_3.tmp)
add_test("TOPP_IDPosteriorErrorProbability_3_out1" ${DIFF} -in1 IDPosteriorErrorProbability_output_3.tmp -in2 ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_OMSSA_output.idXML)
set_tests_properties("TOPP_IDPosteriorErrorProbability_3_out1" PROPERTIES DEPENDS "TOPP_IDPosteriorErrorProbability_3")
# OMSSA input with "split_charge":
add_test("TOPP_IDPosteriorErrorProbability_4" ${TOPP_BIN_PATH}/IDPosteriorErrorProbability -test -in ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_OMSSA_input2.idXML -out IDPosteriorErrorProbability_output_4.tmp -split_charge)
add_test("TOPP_IDPosteriorErrorProbability_4_out1" ${DIFF} -in1 IDPosteriorErrorProbability_output_4.tmp -in2 ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_OMSSA_output2.idXML)
set_tests_properties("TOPP_IDPosteriorErrorProbability_4_out1" PROPERTIES DEPENDS "TOPP_IDPosteriorErrorProbability_4")
# X! Tandem input with "split_charge":
add_test("TOPP_IDPosteriorErrorProbability_5" ${TOPP_BIN_PATH}/IDPosteriorErrorProbability -test -in ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_XTandem_input2.idXML -out IDPosteriorErrorProbability_output_5.tmp -split_charge)
add_test("TOPP_IDPosteriorErrorProbability_5_out1" ${DIFF} -in1 IDPosteriorErrorProbability_output_5.tmp -in2 ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_XTandem_output2.idXML)
set_tests_properties("TOPP_IDPosteriorErrorProbability_5_out1" PROPERTIES DEPENDS "TOPP_IDPosteriorErrorProbability_5")
# Mascot input with "split_charge":
add_test("TOPP_IDPosteriorErrorProbability_6" ${TOPP_BIN_PATH}/IDPosteriorErrorProbability -test -in ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_Mascot_input2.idXML -out IDPosteriorErrorProbability_output_6.tmp -split_charge)
add_test("TOPP_IDPosteriorErrorProbability_6_out1" ${DIFF} -in1 IDPosteriorErrorProbability_output_6.tmp -in2 ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_Mascot_output2.idXML)
set_tests_properties("TOPP_IDPosteriorErrorProbability_6_out1" PROPERTIES DEPENDS "TOPP_IDPosteriorErrorProbability_6")
# with "ignore_bad_data" (no check for correctness of output, bc of unresolved platform specific errors):
add_test("TOPP_IDPosteriorErrorProbability_7" ${TOPP_BIN_PATH}/IDPosteriorErrorProbability -test -in ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_bad_data.idXML -out IDPosteriorErrorProbability_bad_data_output.tmp -ignore_bad_data)
# add_test("TOPP_IDPosteriorErrorProbability_7_out1" ${DIFF} -in1 IDPosteriorErrorProbability_bad_data_output.tmp -in2 ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_bad_data_out.idXML)
# with "prob_correct":
add_test("TOPP_IDPosteriorErrorProbability_8" ${TOPP_BIN_PATH}/IDPosteriorErrorProbability -test -in ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_OMSSA_input.idXML -out IDPosteriorErrorProbability_output_8.tmp -prob_correct)
add_test("TOPP_IDPosteriorErrorProbability_8_out1" ${DIFF} -in1 IDPosteriorErrorProbability_output_8.tmp -in2 ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_prob_correct_output.idXML)
set_tests_properties("TOPP_IDPosteriorErrorProbability_8_out1" PROPERTIES DEPENDS "TOPP_IDPosteriorErrorProbability_8")

#------------------------------------------------------------------------------
# ProteinResolver tests
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

#------------------------------------------------------------------------------
# ProteinQuantifier tests
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
#add_test("TOPP_ProteinQuantifier_10" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.consensusXML -mzTab_out ProteinQuantifier_10_output.tmp -top 0 -average sum)
#add_test("TOPP_ProteinQuantifier_10_out1" ${DIFF} -in1 ProteinQuantifier_10_output.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_10_output.mzTab)
#set_tests_properties("TOPP_ProteinQuantifier_10_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_10")
#add_test("TOPP_ProteinQuantifier_11" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.featureXML -mzTab_out ProteinQuantifier_11_output.tmp -top 0 -average sum -filter_charge)
#add_test("TOPP_ProteinQuantifier_11_out1" ${DIFF} -in1 ProteinQuantifier_11_output.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_11_output.mzTab)
#set_tests_properties("TOPP_ProteinQuantifier_11_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_11")
#ratios
add_test("TOPP_ProteinQuantifier_12" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.consensusXML -out ProteinQuantifier_12_output.tmp -top 3 -include_all -consensus:fix_peptides -average sum -ratios)
add_test("TOPP_ProteinQuantifier_12_out1" ${DIFF} -in1 ProteinQuantifier_12_output.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_12_output.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_12_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_12")
add_test("TOPP_ProteinQuantifier_13" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.consensusXML -out ProteinQuantifier_13_output.tmp -top 3 -include_all -consensus:fix_peptides -average sum -ratiosSILAC)
add_test("TOPP_ProteinQuantifier_13_out1" ${DIFF} -in1 ProteinQuantifier_13_output.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_13_output.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_13_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_13")
# idXML input (spectral counting):
add_test("TOPP_ProteinQuantifier_14" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.idXML -out ProteinQuantifier_14_output1.tmp -peptide_out ProteinQuantifier_14_output2.tmp -top 0 -average sum)
add_test("TOPP_ProteinQuantifier_14_out1" ${DIFF} -in1 ProteinQuantifier_14_output1.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_14_output1.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_14_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_14")
add_test("TOPP_ProteinQuantifier_14_out2" ${DIFF} -in1 ProteinQuantifier_14_output2.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_14_output2.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_14_out2" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_14")

#------------------------------------------------------------------------------
# InclusionExclusionListCreator tests
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
add_test("TOPP_InclusionExclusionListCreator_5" ${TOPP_BIN_PATH}/InclusionExclusionListCreator -test -exclude ${DATA_DIR_TOPP}/InclusionExclusionListCreator.idXML -out InclusionExclusionListCreator_5_output.tmp -exclusion_charges 1 2 -rt_model ${DATA_DIR_TOPP}/InclusionExclusionListCreator_rt.model)
add_test("TOPP_InclusionExclusionListCreator_5_out1" ${DIFF} -in1 InclusionExclusionListCreator_5_output.tmp -in2 ${DATA_DIR_TOPP}/InclusionExclusionListCreator_5_output.txt )
set_tests_properties("TOPP_InclusionExclusionListCreator_5_out1" PROPERTIES DEPENDS "TOPP_InclusionExclusionListCreator_5")
add_test("TOPP_InclusionExclusionListCreator_6" ${TOPP_BIN_PATH}/InclusionExclusionListCreator -test -exclude ${DATA_DIR_TOPP}/InclusionExclusionListCreator.idXML -out InclusionExclusionListCreator_6_output.tmp -exclusion_charges 1 2 -rt_model ${DATA_DIR_TOPP}/InclusionExclusionListCreator_rt.model  -ini ${DATA_DIR_TOPP}/InclusionExclusionListCreator_6.ini)
add_test("TOPP_InclusionExclusionListCreator_6_out1" ${DIFF} -in1 InclusionExclusionListCreator_6_output.tmp -in2 ${DATA_DIR_TOPP}/InclusionExclusionListCreator_6_output.txt )
set_tests_properties("TOPP_InclusionExclusionListCreator_6_out1" PROPERTIES DEPENDS "TOPP_InclusionExclusionListCreator_6")

#------------------------------------------------------------------------------
# PeptideIndexer tests
add_test("TOPP_PeptideIndexer_1" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_1.idXML -out PeptideIndexer_1_out.tmp.idXML -allow_unmatched -full_tolerant_search -enzyme:specificity none)
add_test("TOPP_PeptideIndexer_1_out" ${DIFF} -in1 PeptideIndexer_1_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/PeptideIndexer_1_out.idXML )
set_tests_properties("TOPP_PeptideIndexer_1_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_1")
add_test("TOPP_PeptideIndexer_2" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_1.idXML -out PeptideIndexer_2_out.tmp.idXML -allow_unmatched -write_protein_sequence -full_tolerant_search -enzyme:specificity none)
add_test("TOPP_PeptideIndexer_2_out" ${DIFF} -in1 PeptideIndexer_2_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/PeptideIndexer_2_out.idXML )
set_tests_properties("TOPP_PeptideIndexer_2_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_2")
add_test("TOPP_PeptideIndexer_3" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_1.idXML -out PeptideIndexer_3_out.tmp.idXML -allow_unmatched -keep_unreferenced_proteins -full_tolerant_search -enzyme:specificity none)
add_test("TOPP_PeptideIndexer_3_out" ${DIFF} -in1 PeptideIndexer_3_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/PeptideIndexer_3_out.idXML )
set_tests_properties("TOPP_PeptideIndexer_3_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_3")
add_test("TOPP_PeptideIndexer_4" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_1.idXML -out PeptideIndexer_4_out.tmp.idXML -allow_unmatched -aaa_max 0 -write_protein_sequence -enzyme:specificity none)
add_test("TOPP_PeptideIndexer_4_out" ${DIFF} -in1 PeptideIndexer_4_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/PeptideIndexer_4_out.idXML )
set_tests_properties("TOPP_PeptideIndexer_4_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_4")
add_test("TOPP_PeptideIndexer_5" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_1.idXML -out PeptideIndexer_5_out.tmp.idXML -allow_unmatched -enzyme:specificity none)
add_test("TOPP_PeptideIndexer_5_out" ${DIFF} -in1 PeptideIndexer_5_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/PeptideIndexer_5_out.idXML )
set_tests_properties("TOPP_PeptideIndexer_5_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_5")
add_test("TOPP_PeptideIndexer_6" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_2.idXML -out PeptideIndexer_6_out.tmp.idXML -allow_unmatched -enzyme:specificity none)
add_test("TOPP_PeptideIndexer_6_out" ${DIFF} -in1 PeptideIndexer_6_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/PeptideIndexer_6_out.idXML )
set_tests_properties("TOPP_PeptideIndexer_6_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_6")
# specificity tests
add_test("TOPP_PeptideIndexer_7" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_3.idXML -out PeptideIndexer_7_out.tmp.idXML -allow_unmatched)
add_test("TOPP_PeptideIndexer_7_out" ${DIFF} -in1 PeptideIndexer_7_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/PeptideIndexer_7_out.idXML )
set_tests_properties("TOPP_PeptideIndexer_7_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_7")
add_test("TOPP_PeptideIndexer_8" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_3.idXML -out PeptideIndexer_8_out.tmp.idXML -allow_unmatched -enzyme:specificity semi)
add_test("TOPP_PeptideIndexer_8_out" ${DIFF} -in1 PeptideIndexer_8_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/PeptideIndexer_8_out.idXML )
set_tests_properties("TOPP_PeptideIndexer_8_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_8")
add_test("TOPP_PeptideIndexer_9" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_3.idXML -out PeptideIndexer_9_out.tmp.idXML -allow_unmatched -enzyme:specificity none)
add_test("TOPP_PeptideIndexer_9_out" ${DIFF} -in1 PeptideIndexer_9_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/PeptideIndexer_9_out.idXML )
set_tests_properties("TOPP_PeptideIndexer_9_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_9")
# "IL_equivalent" option (and 1) make sure that I/L does not count towards ambiguous AAs; 2) that the original protein and peptide sequence (no I/L substitutions) are reported
add_test("TOPP_PeptideIndexer_10" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_10_input.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_10_input.idXML -out PeptideIndexer_10_output.tmp.idXML -IL_equivalent -aaa_max 3 -write_protein_sequence)
add_test("TOPP_PeptideIndexer_10_out" ${DIFF} -in1 PeptideIndexer_10_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/PeptideIndexer_10_output.idXML )
set_tests_properties("TOPP_PeptideIndexer_10_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_10")
# empty FASTA
add_test("TOPP_PeptideIndexer_11" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/degenerate_cases/empty.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_1.idXML -out PeptideIndexer_11_out.tmp.idXML)
set_tests_properties("TOPP_PeptideIndexer_11" PROPERTIES WILL_FAIL 1)
# with "annotate_proteins" (target/decoy):
add_test("TOPP_PeptideIndexer_12" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_1.idXML -out PeptideIndexer_12_out.tmp.idXML -allow_unmatched -enzyme:specificity none)
add_test("TOPP_PeptideIndexer_12_out" ${DIFF} -in1 PeptideIndexer_12_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/PeptideIndexer_12_out.idXML )
set_tests_properties("TOPP_PeptideIndexer_12_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_12")
# empty input
add_test("TOPP_PeptideIndexer_13" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -in ${DATA_DIR_TOPP}/degenerate_cases/empty.idXML -out PeptideIndexer_13_out.tmp.idXML)
add_test("TOPP_PeptideIndexer_13_out" ${DIFF} -in1 PeptideIndexer_13_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/degenerate_cases/empty.idXML )
set_tests_properties("TOPP_PeptideIndexer_13_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_13")
# allow mismatches (supported in full_tolerant_search mode only)
# -- first test without mismatches
add_test("TOPP_PeptideIndexer_14" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_14.idXML -out PeptideIndexer_14_out.tmp.idXML -mismatches_max 0 -full_tolerant_search -allow_unmatched)
add_test("TOPP_PeptideIndexer_14_out" ${DIFF} -in1 PeptideIndexer_14_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/PeptideIndexer_14_out.idXML )
set_tests_properties("TOPP_PeptideIndexer_14_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_14")
# -- and now with mismatches allowed
add_test("TOPP_PeptideIndexer_15" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_14.idXML -out PeptideIndexer_15_out.tmp.idXML -mismatches_max 1 -full_tolerant_search -allow_unmatched)
add_test("TOPP_PeptideIndexer_15_out" ${DIFF} -in1 PeptideIndexer_15_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/PeptideIndexer_15_out.idXML )
set_tests_properties("TOPP_PeptideIndexer_15_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_15")
add_test("TOPP_PeptideIndexer_16" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_14.idXML -out PeptideIndexer_16_out.tmp.idXML -mismatches_max 2 -full_tolerant_search -allow_unmatched)
add_test("TOPP_PeptideIndexer_16_out" ${DIFF} -in1 PeptideIndexer_16_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/PeptideIndexer_16_out.idXML )
set_tests_properties("TOPP_PeptideIndexer_16_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_16")
add_test("TOPP_PeptideIndexer_17" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_14.idXML -out PeptideIndexer_17_out.tmp.idXML -mismatches_max 3 -full_tolerant_search)
add_test("TOPP_PeptideIndexer_17_out" ${DIFF} -in1 PeptideIndexer_17_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/PeptideIndexer_17_out.idXML )
set_tests_properties("TOPP_PeptideIndexer_17_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_17")
## -- test correct DB index handling of partial results from AC vs SA
add_test("TOPP_PeptideIndexer_18" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_2.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_18.idXML -out PeptideIndexer_18_out.tmp.idXML -missing_decoy_action warn -filter_aaa_proteins)
add_test("TOPP_PeptideIndexer_18_out" ${DIFF} -in1 PeptideIndexer_18_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/PeptideIndexer_18_out.idXML )
set_tests_properties("TOPP_PeptideIndexer_18_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_18")
## -- same as 18, but using SA all the way -- results should be identical
add_test("TOPP_PeptideIndexer_19" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_2.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_18.idXML -out PeptideIndexer_19_out.tmp.idXML  -missing_decoy_action warn -filter_aaa_proteins -full_tolerant_search)
add_test("TOPP_PeptideIndexer_19_out" ${DIFF} -in1 PeptideIndexer_19_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/PeptideIndexer_19_out.idXML )
set_tests_properties("TOPP_PeptideIndexer_19_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_19")

if(WITH_GUI)
  #------------------------------------------------------------------------------
  # ExecutePipeline tests (as substitute for TOPPAS) - the ResourceFiles are in binary tree, as they have been configured from a .in file (see above)!
  add_test("TOPP_ExecutePipeline_1" ${TOPP_BIN_PATH}/ExecutePipeline -test -in ${DATA_DIR_TOPP}/ExecutePipeline_1.toppas -resource_file ${DATA_DIR_TOPP_BIN}/ExecutePipeline_1.trf -out_dir .)
  set_tests_properties("TOPP_ExecutePipeline_1" PROPERTIES WILL_FAIL 1)
endif()

#------------------------------------------------------------------------------
# search engine tests (might be optional, depending on presence of external binaries, e.g., omssacl.exe)
#------------------------------------------------------------------------------
include(THIRDPARTY/third_party_tests.cmake)

#------------------------------------------------------------------------------
# UTILS tests (if this becomes too much, we should create a separate directory)
#------------------------------------------------------------------------------
# ImageCreator:
if(WITH_GUI)
  add_test("UTILS_ImageCreator_1" ${TOPP_BIN_PATH}/ImageCreator -test -in ${DATA_DIR_TOPP}/ImageCreator_1_input.mzML -out ImageCreator1.bmp -width 20 -height 15)
  add_test("UTILS_ImageCreator_1_out1" ${DIFF} -in1 ImageCreator1.bmp -in2 ${DATA_DIR_TOPP}/ImageCreator_1_output.bmp )
  set_tests_properties("UTILS_ImageCreator_1_out1" PROPERTIES DEPENDS "UTILS_ImageCreator_1")
  # testing forced out_type
  add_test("UTILS_ImageCreator_2" ${TOPP_BIN_PATH}/ImageCreator -test -in ${DATA_DIR_TOPP}/ImageCreator_2_input.mzML -out ImageCreator2.png -out_type bmp -width 20 -height 15 -precursors -precursor_size 3 -precursor_color green -log_intensity)
  add_test("UTILS_ImageCreator_2_out1" ${DIFF} -in1 ImageCreator2.png -in2 ${DATA_DIR_TOPP}/ImageCreator_2_output.bmp )
  set_tests_properties("UTILS_ImageCreator_2_out1" PROPERTIES DEPENDS "UTILS_ImageCreator_2")
endif()

# IDScoreSwitcher - switch q-value and PEP in Percolator result:
add_test("UTILS_IDScoreSwitcher_1" ${TOPP_BIN_PATH}/IDScoreSwitcher -test -in ${DATA_DIR_TOPP}/IDFileConverter_12_output.idXML -out IDScoreSwitcher_1_output.tmp -new_score Percolator_PEP -new_score_type "Posterior Error Probability" -new_score_orientation lower_better -old_score Percolator_qvalue)
add_test("UTILS_IDScoreSwitcher_1_out1" ${DIFF} -in1 IDScoreSwitcher_1_output.tmp -in2 ${DATA_DIR_TOPP}/IDFileConverter_13_output.idXML -whitelist "IdentificationRun date")
set_tests_properties("UTILS_IDScoreSwitcher_1_out1" PROPERTIES DEPENDS "UTILS_IDScoreSwitcher_1")
# IDScoreSwitcher - switch q-value and protein probability:
add_test("UTILS_IDScoreSwitcher_2" ${TOPP_BIN_PATH}/IDScoreSwitcher -test -in ${DATA_DIR_TOPP}/FalseDiscoveryRate_5_output.idXML -out IDScoreSwitcher_2_output.tmp -new_score "Posterior Probability_score" -new_score_type "Posterior Probability" -new_score_orientation higher_better -proteins)
add_test("UTILS_IDScoreSwitcher_2_out1" ${DIFF} -in1 IDScoreSwitcher_2_output.tmp -in2 ${DATA_DIR_TOPP}/IDScoreSwitcher_2_output.idXML -whitelist "IdentificationRun date")
set_tests_properties("UTILS_IDScoreSwitcher_2_out1" PROPERTIES DEPENDS "UTILS_IDScoreSwitcher_2")

# IDSplitter:
add_test("UTILS_IDSplitter_1" ${TOPP_BIN_PATH}/IDSplitter -test -in ${DATA_DIR_TOPP}/IDMapper_1_output.featureXML -out IDSplitter_1_output1.tmp -id_out IDSplitter_1_output2.tmp)
add_test("UTILS_IDSplitter_1_out1" ${DIFF} -in1 IDSplitter_1_output1.tmp -in2 ${DATA_DIR_TOPP}/IDSplitter_1_output.featureXML )
set_tests_properties("UTILS_IDSplitter_1_out1" PROPERTIES DEPENDS "UTILS_IDSplitter_1")
add_test("UTILS_IDSplitter_1_out2" ${DIFF} -in1 IDSplitter_1_output2.tmp -in2 ${DATA_DIR_TOPP}/IDSplitter_1_output.idXML )
set_tests_properties("UTILS_IDSplitter_1_out2" PROPERTIES DEPENDS "UTILS_IDSplitter_1")

# MassCalculator:
add_test("UTILS_MassCalculator_1" ${TOPP_BIN_PATH}/MassCalculator -test -in ${DATA_DIR_TOPP}/MassCalculator_1_input.txt -out MassCalculator_1_output.tmp -charge 0 1 -format table -separator ,)
add_test("UTILS_MassCalculator_1_out1" ${DIFF} -in1 MassCalculator_1_output.tmp -in2 ${DATA_DIR_TOPP}/MassCalculator_1_output.csv )
set_tests_properties("UTILS_MassCalculator_1_out1" PROPERTIES DEPENDS "UTILS_MassCalculator_1")
add_test("UTILS_MassCalculator_2" ${TOPP_BIN_PATH}/MassCalculator -test -in_seq "LDQWLC(Carbamidomethyl)EKL" "(Glu->pyro-Glu)EAM(Oxidation)APKHK" "RANVM(Oxidation)DYR" "FGVEQDVDMVFASFIR" -out MassCalculator_2_output.tmp -charge 1 2 3)
add_test("UTILS_MassCalculator_2_out1" ${DIFF} -in1 MassCalculator_2_output.tmp -in2 ${DATA_DIR_TOPP}/MassCalculator_2_output.txt )
set_tests_properties("UTILS_MassCalculator_2_out1" PROPERTIES DEPENDS "UTILS_MassCalculator_2")

#------------------------------------------------------------------------------
# MultiplexResolver test
add_test("UTILS_MultiplexResolver_1" ${TOPP_BIN_PATH}/MultiplexResolver -test -in ${DATA_DIR_TOPP}/MultiplexResolver_1_input.consensusXML -ini ${DATA_DIR_TOPP}/MultiplexResolver_1_parameters.ini -out MultiplexResolver_1.tmp -out_conflicts MultiplexResolver_1_conflicts.tmp)
add_test("UTILS_MultiplexResolver_1_out1" ${DIFF} -whitelist "id=" "href=" -in1 MultiplexResolver_1.tmp -in2 ${DATA_DIR_TOPP}/MultiplexResolver_1_output.consensusXML )
set_tests_properties("UTILS_MultiplexResolver_1_out1" PROPERTIES DEPENDS "UTILS_MultiplexResolver_1")
add_test("UTILS_MultiplexResolver_1_out2" ${DIFF} -whitelist "id=" "href=" -in1 MultiplexResolver_1_conflicts.tmp -in2 ${DATA_DIR_TOPP}/MultiplexResolver_1_output_conflicts.consensusXML )
set_tests_properties("UTILS_MultiplexResolver_1_out2" PROPERTIES DEPENDS "UTILS_MultiplexResolver_1")
add_test("UTILS_MultiplexResolver_2" ${TOPP_BIN_PATH}/MultiplexResolver -test -in ${DATA_DIR_TOPP}/MultiplexResolver_2_input.consensusXML -ini ${DATA_DIR_TOPP}/MultiplexResolver_2_parameters.ini -out MultiplexResolver_2.tmp -out_conflicts MultiplexResolver_2_conflicts.tmp)
add_test("UTILS_MultiplexResolver_2_out1" ${DIFF} -whitelist "id=" "href=" -in1 MultiplexResolver_2.tmp -in2 ${DATA_DIR_TOPP}/MultiplexResolver_2_output.consensusXML )
set_tests_properties("UTILS_MultiplexResolver_2_out1" PROPERTIES DEPENDS "UTILS_MultiplexResolver_2")
add_test("UTILS_MultiplexResolver_2_out2" ${DIFF} -whitelist "id=" "href=" -in1 MultiplexResolver_2_conflicts.tmp -in2 ${DATA_DIR_TOPP}/MultiplexResolver_2_output_conflicts.consensusXML )
set_tests_properties("UTILS_MultiplexResolver_2_out2" PROPERTIES DEPENDS "UTILS_MultiplexResolver_2")
add_test("UTILS_MultiplexResolver_3" ${TOPP_BIN_PATH}/MultiplexResolver -test -in ${DATA_DIR_TOPP}/MultiplexResolver_3_input.consensusXML -ini ${DATA_DIR_TOPP}/MultiplexResolver_3_parameters.ini -out MultiplexResolver_3.tmp -out_conflicts MultiplexResolver_3_conflicts.tmp)
add_test("UTILS_MultiplexResolver_3_out1" ${DIFF} -whitelist "id=" "href=" -in1 MultiplexResolver_3.tmp -in2 ${DATA_DIR_TOPP}/MultiplexResolver_3_output.consensusXML )
set_tests_properties("UTILS_MultiplexResolver_3_out1" PROPERTIES DEPENDS "UTILS_MultiplexResolver_3")
add_test("UTILS_MultiplexResolver_3_out3" ${DIFF} -whitelist "id=" "href=" -in1 MultiplexResolver_3_conflicts.tmp -in2 ${DATA_DIR_TOPP}/MultiplexResolver_3_output_conflicts.consensusXML )
set_tests_properties("UTILS_MultiplexResolver_3_out3" PROPERTIES DEPENDS "UTILS_MultiplexResolver_3")

# INIUpdater:
if(WITH_GUI)
  add_test("UTILS_INIUpdater_1" ${TOPP_BIN_PATH}/INIUpdater -test -in ${DATA_DIR_TOPP}/INIUpdater_1_noupdate.toppas -out INIUpdater_1_noupdate.toppas.tmp)
  add_test("UTILS_INIUpdater_1_out" ${DIFF} -in1 INIUpdater_1_noupdate.toppas.tmp -in2 ${DATA_DIR_TOPP}/INIUpdater_1_noupdate.toppas )
  set_tests_properties("UTILS_INIUpdater_1_out" PROPERTIES DEPENDS "UTILS_INIUpdater_1")

  add_test("UTILS_INIUpdater_2" ${TOPP_BIN_PATH}/INIUpdater -test -in ${DATA_DIR_TOPP}/INIUpdater_2_broken.toppas -out INIUpdater_2_broken.toppas.tmp)
  set_tests_properties("UTILS_INIUpdater_2" PROPERTIES WILL_FAIL 1)

  add_test("UTILS_INIUpdater_3" ${TOPP_BIN_PATH}/INIUpdater -test -in ${DATA_DIR_TOPP}/INIUpdater_3_old.toppas -out INIUpdater_3_old.toppas.tmp)
  add_test("UTILS_INIUpdater_3_out" ${DIFF} -in1 INIUpdater_3_old.toppas.tmp -in2 ${DATA_DIR_TOPP}/INIUpdater_3_old_out.toppas )
  set_tests_properties("UTILS_INIUpdater_3_out" PROPERTIES DEPENDS "UTILS_INIUpdater_3")
endif()

# DatabaseFilter test:
add_test("UTILS_DatabaseFilter_1" ${TOPP_BIN_PATH}/DatabaseFilter -test -in ${DATA_DIR_TOPP}/DatabaseFilter_1.fasta -id ${DATA_DIR_TOPP}/DatabaseFilter_1.idXML -out DatabaseFilter_1_out.fasta.tmp)
add_test("UTILS_DatabaseFilter_1_out" ${DIFF} -in1 DatabaseFilter_1_out.fasta.tmp -in2 ${DATA_DIR_TOPP}/DatabaseFilter_1_out.fasta )
set_tests_properties("UTILS_DatabaseFilter_1_out" PROPERTIES DEPENDS "UTILS_DatabaseFilter_1")
add_test("UTILS_DatabaseFilter_2" ${TOPP_BIN_PATH}/DatabaseFilter -test -in ${DATA_DIR_TOPP}/DatabaseFilter_1.fasta -id ${DATA_DIR_TOPP}/DatabaseFilter_1.idXML -out DatabaseFilter_2_out.fasta.tmp -method blacklist)
add_test("UTILS_DatabaseFilter_2_out" ${DIFF} -in1 DatabaseFilter_2_out.fasta.tmp -in2 ${DATA_DIR_TOPP}/DatabaseFilter_2_out.fasta )
set_tests_properties("UTILS_DatabaseFilter_2_out" PROPERTIES DEPENDS "UTILS_DatabaseFilter_2")
add_test("UTILS_DatabaseFilter_3" ${TOPP_BIN_PATH}/DatabaseFilter -test -in ${DATA_DIR_TOPP}/DatabaseFilter_3.fasta -id ${DATA_DIR_TOPP}/DatabaseFilter_3.mzid -out DatabaseFilter_3_out.fasta.tmp)
add_test("UTILS_DatabaseFilter_3_out" ${DIFF} -in1 DatabaseFilter_3_out.fasta.tmp -in2 ${DATA_DIR_TOPP}/DatabaseFilter_3_out.fasta )
set_tests_properties("UTILS_DatabaseFilter_3_out" PROPERTIES DEPENDS "UTILS_DatabaseFilter_3")
add_test("UTILS_DatabaseFilter_4" ${TOPP_BIN_PATH}/DatabaseFilter -test -in ${DATA_DIR_TOPP}/DatabaseFilter_3.fasta -id ${DATA_DIR_TOPP}/DatabaseFilter_3.mzid -out DatabaseFilter_4_out.fasta.tmp -method blacklist)
add_test("UTILS_DatabaseFilter_4_out" ${DIFF} -in1 DatabaseFilter_4_out.fasta.tmp -in2 ${DATA_DIR_TOPP}/DatabaseFilter_4_out.fasta )
set_tests_properties("UTILS_DatabaseFilter_4_out" PROPERTIES DEPENDS "UTILS_DatabaseFilter_4")


# DecoyDatabase:
add_test("UTILS_DecoyDatabase_1" ${TOPP_BIN_PATH}/DecoyDatabase -test -in ${DATA_DIR_TOPP}/DecoyDatabase_1.fasta -out DecoyDatabase_1.fasta.tmp -only_decoy)
add_test("UTILS_DecoyDatabase_1_out" ${DIFF} -in1 DecoyDatabase_1.fasta.tmp -in2 ${DATA_DIR_TOPP}/DecoyDatabase_1_out.fasta )
set_tests_properties("UTILS_DecoyDatabase_1_out" PROPERTIES DEPENDS "UTILS_DecoyDatabase_1")
add_test("UTILS_DecoyDatabase_2" ${TOPP_BIN_PATH}/DecoyDatabase -test -in ${DATA_DIR_TOPP}/DecoyDatabase_1.fasta -out DecoyDatabase_2.fasta.tmp -decoy_string "blabla" -decoy_string_position "prefix")
add_test("UTILS_DecoyDatabase_2_out" ${DIFF} -in1 DecoyDatabase_2.fasta.tmp -in2 ${DATA_DIR_TOPP}/DecoyDatabase_2_out.fasta )
set_tests_properties("UTILS_DecoyDatabase_2_out" PROPERTIES DEPENDS "UTILS_DecoyDatabase_2")

# SimpleSearchEngine:
add_test("UTILS_SimpleSearchEngine_1" ${TOPP_BIN_PATH}/SimpleSearchEngine -test
-ini ${DATA_DIR_TOPP}/SimpleSearchEngine_1.ini -in
${DATA_DIR_TOPP}/SimpleSearchEngine_1.mzML -out SimpleSearchEngine_1_out.tmp
-database ${DATA_DIR_TOPP}/SimpleSearchEngine_1.fasta)
add_test("UTILS_SimpleSearchEngine_1_out" ${DIFF} -in1 SimpleSearchEngine_1_out.tmp -in2 ${DATA_DIR_TOPP}/SimpleSearchEngine_1_out.idXML -whitelist "IdentificationRun date" "SearchParameters id=\"SP_0\" db=")
set_tests_properties("UTILS_SimpleSearchEngine_1_out" PROPERTIES DEPENDS
"UTILS_SimpleSearchEngine_1")

# FeatureFinderSuperHirn - test on centroided data:
add_test("UTILS_FeatureFinderSuperHirn_1" ${TOPP_BIN_PATH}/FeatureFinderSuperHirn -test -in ${DATA_DIR_TOPP}/FeatureFinderSuperHirn_input_1.mzML -out FeatureFinderSuperHirn_1_output.featureXML.tmp -ini ${DATA_DIR_TOPP}/FeatureFinderSuperHirn_1_parameters.ini)
add_test("UTILS_FeatureFinderSuperHirn_1_out1" ${DIFF} -whitelist "id=" -in1 FeatureFinderSuperHirn_1_output.featureXML.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinderSuperHirn_1_output.featureXML)
set_tests_properties("UTILS_FeatureFinderSuperHirn_1_out1" PROPERTIES DEPENDS "UTILS_FeatureFinderSuperHirn_1")
# FeatureFinder of SuperHirn - test on profile data:
add_test("UTILS_FeatureFinderSuperHirn_2" ${TOPP_BIN_PATH}/FeatureFinderSuperHirn -test -in ${DATA_DIR_TOPP}/FeatureFinderSuperHirn_input_2.mzML -out FeatureFinderSuperHirn_2_output.featureXML.tmp -ini ${DATA_DIR_TOPP}/FeatureFinderSuperHirn_2_parameters.ini)
add_test("UTILS_FeatureFinderSuperHirn_2_out1" ${DIFF} -whitelist "id=" -in1 FeatureFinderSuperHirn_2_output.featureXML.tmp -in2 ${DATA_DIR_TOPP}/FeatureFinderSuperHirn_2_output.featureXML )
set_tests_properties("UTILS_FeatureFinderSuperHirn_2_out1" PROPERTIES DEPENDS "UTILS_FeatureFinderSuperHirn_2")

# MzMLSplitter test:
# specify number of parts directly:
add_test("UTILS_MzMLSplitter_1" ${TOPP_BIN_PATH}/MzMLSplitter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out MzMLSplitter_1_output -parts 2)
add_test("UTILS_MzMLSplitter_1_out1" ${DIFF} -in1 MzMLSplitter_1_output_part1of2.mzML -in2 ${DATA_DIR_TOPP}/MzMLSplitter_output_part1.mzML)
set_tests_properties("UTILS_MzMLSplitter_1_out1" PROPERTIES DEPENDS "UTILS_MzMLSplitter_1")
add_test("UTILS_MzMLSplitter_1_out2" ${DIFF} -in1 MzMLSplitter_1_output_part2of2.mzML -in2 ${DATA_DIR_TOPP}/MzMLSplitter_output_part2.mzML)
set_tests_properties("UTILS_MzMLSplitter_1_out2" PROPERTIES DEPENDS "UTILS_MzMLSplitter_1")
# specify target file size:
add_test("UTILS_MzMLSplitter_2" ${TOPP_BIN_PATH}/MzMLSplitter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out MzMLSplitter_2_output -size 40 -unit KB)
add_test("UTILS_MzMLSplitter_2_out1" ${DIFF} -in1 MzMLSplitter_2_output_part1of2.mzML -in2 ${DATA_DIR_TOPP}/MzMLSplitter_output_part1.mzML)
set_tests_properties("UTILS_MzMLSplitter_2_out1" PROPERTIES DEPENDS "UTILS_MzMLSplitter_2")
add_test("UTILS_MzMLSplitter_2_out2" ${DIFF} -in1 MzMLSplitter_2_output_part2of2.mzML -in2 ${DATA_DIR_TOPP}/MzMLSplitter_output_part2.mzML)
set_tests_properties("UTILS_MzMLSplitter_2_out2" PROPERTIES DEPENDS "UTILS_MzMLSplitter_2")

# TICCalculator test:
add_test("UTILS_TICCalculator_1" ${TOPP_BIN_PATH}/TICCalculator -test -in ${DATA_DIR_TOPP}/MapNormalizer_output.mzML -read_method regular)
add_test("UTILS_TICCalculator_2" ${TOPP_BIN_PATH}/TICCalculator -test -in ${DATA_DIR_TOPP}/MapNormalizer_output.mzML -read_method streaming -loadData true)
add_test("UTILS_TICCalculator_3" ${TOPP_BIN_PATH}/TICCalculator -test -in ${DATA_DIR_TOPP}/MapNormalizer_output.mzML -read_method streaming -loadData false)
add_test("UTILS_TICCalculator_4" ${TOPP_BIN_PATH}/TICCalculator -test -in ${DATA_DIR_TOPP}/MapNormalizer_output.mzML -read_method indexed)
add_test("UTILS_TICCalculator_5" ${TOPP_BIN_PATH}/TICCalculator -test -in ${DATA_DIR_TOPP}/MapNormalizer_output.mzML -read_method indexed_parallel)

