# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2023.
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
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

#
# README:
#
# This file contains the declarations of all TOPP tests that are executed by CTest
#
# In order to add a test do the following:
# - add specific tests for your tool (this usually involves FuzzyDiff) - see below for examples
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
set(DATA_DIR_SHARE ${OPENMS_HOST_DIRECTORY}/share/OpenMS)
set(DATA_DIR_TOPP_BIN ${PROJECT_BINARY_DIR})
set(INDEX_WHITELIST "offset" "indexListOffset")
set(DIFF ${TOPP_BIN_PATH}/FuzzyDiff -test -ini ${DATA_DIR_TOPP}/FuzzyDiff.ini)

#------------------------------------------------------------------------------
# WRITE_INI test for each TOPP tool, checking -write_ini functionality
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# TOPP
set(TOPP_executables_INI ${TOPP_TOOLS})
list(REMOVE_ITEM TOPP_executables_INI GenericWrapper) # GenericWrapper (has a -type flag)
list(REMOVE_ITEM TOPP_executables_INI OpenMSInfo)
add_test("TOPPWRITEINI_GenericWrapper" ${TOPP_BIN_PATH}/GenericWrapper -test -type RAWFileConvert -write_ini "GenericWrapper.tmp.ini")
foreach (i ${TOPP_executables_INI})
  add_test("TOPPWRITEINI_${i}" "${TOPP_BIN_PATH}/${i}" -test -write_ini "${i}.tmp.ini")
  add_test("TOPPWRITEINI_${i}_SectionName" ${CMAKE_COMMAND} -test -DNAME=${i} -DFILE=${i}.tmp.ini -P "${DATA_DIR_TOPP}/check_ini.cmake")
  set_tests_properties("TOPPWRITEINI_${i}_SectionName" PROPERTIES DEPENDS "TOPPWRITEINI_${i}")
endforeach(i)

#------------------------------------------------------------------------------
# WRITECTD test to check whether -write_ctd works
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# TOPP
set(TOPP_executables_CTD ${TOPP_TOOLS})
list(REMOVE_ITEM TOPP_executables_CTD GenericWrapper)
list(REMOVE_ITEM TOPP_executables_CTD OpenMSInfo)
foreach (i ${TOPP_executables_CTD})
  add_test("TOPPWRITECTD_${i}" "${TOPP_BIN_PATH}/${i}" -test -write_ctd)
endforeach()

#------------------------------------------------------------------------------
# test if write_ini does the correct updates when given an old INI file (using
# one arbitrary TOPP tool)
#------------------------------------------------------------------------------
add_test("TOPPWRITEINI_OVERWRITE" ${TOPP_BIN_PATH}/PeakPickerHiRes -test -write_ini "WRITE_INI.tmp.ini" -ini "${DATA_DIR_TOPP}/WRITE_INI_IN.ini")
add_test("TOPPWRITEINI_OVERWRITE_out" ${DIFF} -in1 WRITE_INI.tmp.ini -in2 "${DATA_DIR_TOPP}/WRITE_INI_OUT.ini" -whitelist "version")
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

# FuzzyDiff tests
add_test("TOPP_FuzzyDiff_1" ${TOPP_BIN_PATH}/FuzzyDiff -test -ini ${DATA_DIR_TOPP}/FuzzyDiff.ini -in1 ${DATA_DIR_TOPP}/FuzzyDiff_1_in1.featureXML -in2 ${DATA_DIR_TOPP}/FuzzyDiff_1_in1.featureXML)
add_test("TOPP_FuzzyDiff_2" ${TOPP_BIN_PATH}/FuzzyDiff -test -ini ${DATA_DIR_TOPP}/FuzzyDiff.ini -in1 ${DATA_DIR_TOPP}/FuzzyDiff_1_in1.featureXML -in2 ${DATA_DIR_TOPP}/FuzzyDiff_1_in2.featureXML)
add_test("TOPP_FuzzyDiff_3" ${TOPP_BIN_PATH}/FuzzyDiff -test -ini ${DATA_DIR_TOPP}/FuzzyDiff.ini -in1 ${DATA_DIR_TOPP}/FuzzyDiff_3_in1.featureXML -in2 ${DATA_DIR_TOPP}/FuzzyDiff_3_in2.featureXML)
add_test("TOPP_FuzzyDiff_4" ${TOPP_BIN_PATH}/FuzzyDiff -test -ini ${DATA_DIR_TOPP}/FuzzyDiff.ini -in1 ${DATA_DIR_TOPP}/lorem_ipsum.featureXML -in2 ${DATA_DIR_TOPP}/FuzzyDiff_3_in2.featureXML)
set_tests_properties("TOPP_FuzzyDiff_1" PROPERTIES WILL_FAIL 1) ## needs to fail due to different numbers
set_tests_properties("TOPP_FuzzyDiff_2" PROPERTIES WILL_FAIL 1) ## needs to fail (cheating)
set_tests_properties("TOPP_FuzzyDiff_4" PROPERTIES WILL_FAIL 1) ## needs to fail (file does not exist)

#------------------------------------------------------------------------------
# IDMerger tests
add_test("TOPP_IDMerger_1" ${TOPP_BIN_PATH}/IDMerger -test -in ${DATA_DIR_TOPP}/IDMerger_1_input1.idXML ${DATA_DIR_TOPP}/IDMerger_1_input2.idXML -out IDMerger_1_output.tmp.idXML)
add_test("TOPP_IDMerger_1_out1" ${DIFF} -in1 IDMerger_1_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDMerger_1_output.idXML)
set_tests_properties("TOPP_IDMerger_1_out1" PROPERTIES DEPENDS "TOPP_IDMerger_1")
# merge pepXML and protXML:
add_test("TOPP_IDMerger_2" ${TOPP_BIN_PATH}/IDMerger -test -pepxml_protxml -in ${DATA_DIR_TOPP}/IDMerger_2_input1.idXML ${DATA_DIR_TOPP}/IDMerger_2_input2.idXML -out IDMerger_2_output.tmp.idXML -annotate_file_origin "false")
add_test("TOPP_IDMerger_2_out1" ${DIFF} -in1 IDMerger_2_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDMerger_2_output.idXML )
set_tests_properties("TOPP_IDMerger_2_out1" PROPERTIES DEPENDS "TOPP_IDMerger_2")
# conflicting IDs:
add_test("TOPP_IDMerger_3" ${TOPP_BIN_PATH}/IDMerger -test -in ${DATA_DIR_TOPP}/IDMerger_1_input1.idXML ${DATA_DIR_TOPP}/IDMerger_1_input1.idXML -out IDMerger_3_output.tmp.idXML -annotate_file_origin "false")
add_test("TOPP_IDMerger_3_out1" ${DIFF} -in1 IDMerger_3_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDMerger_3_output.idXML )
set_tests_properties("TOPP_IDMerger_3_out1" PROPERTIES DEPENDS "TOPP_IDMerger_3")
# empty IDs:
add_test("TOPP_IDMerger_4" ${TOPP_BIN_PATH}/IDMerger -test -in ${DATA_DIR_TOPP}/degenerated_empty.idXML ${DATA_DIR_TOPP}/degenerated_empty.idXML -out IDMerger_4_output.tmp.idXML -annotate_file_origin "false")
add_test("TOPP_IDMerger_4_out1" ${DIFF} -in1 IDMerger_4_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDMerger_4_output.idXML )
set_tests_properties("TOPP_IDMerger_4_out1" PROPERTIES DEPENDS "TOPP_IDMerger_4")
# "add_to" option:
add_test("TOPP_IDMerger_5" ${TOPP_BIN_PATH}/IDMerger -test -in ${DATA_DIR_TOPP}/IDMerger_2_input1.idXML -add_to ${DATA_DIR_TOPP}/IDMerger_5_input1.idXML -out IDMerger_5_output.tmp.idXML -annotate_file_origin "false")
add_test("TOPP_IDMerger_5_out1" ${DIFF} -in1 IDMerger_5_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDMerger_5_output.idXML )
set_tests_properties("TOPP_IDMerger_5_out1" PROPERTIES DEPENDS "TOPP_IDMerger_5")

# .oms (SQLite) files - inputs based on BSA example data, generated like this:
# IDFileConverter -in share/OpenMS/examples/BSA/BSA2_OMSSA.idXML -out src/tests/topp/IDMerger_6_input2.oms
# IDFileConverter -in share/OpenMS/examples/BSA/BSA3_OMSSA.idXML -out src/tests/topp/IDMerger_6_input1.oms
add_test("TOPP_IDMerger_6" ${TOPP_BIN_PATH}/IDMerger -test -in ${DATA_DIR_TOPP}/IDMerger_6_input1.oms ${DATA_DIR_TOPP}/IDMerger_6_input2.oms -out IDMerger_6_output1.oms)
add_test("TOPP_IDMerger_6_out1" ${TOPP_BIN_PATH}/IDFileConverter -in IDMerger_6_output1.oms -out IDMerger_6_output2.tmp.idXML -out_type idXML)
set_tests_properties("TOPP_IDMerger_6_out1" PROPERTIES DEPENDS "TOPP_IDMerger_6")
add_test("TOPP_IDMerger_6_out2" ${DIFF} -in1 IDMerger_6_output2.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDMerger_6_output.idXML)
set_tests_properties("TOPP_IDMerger_6_out2" PROPERTIES DEPENDS "TOPP_IDMerger_6_out1")

#------------------------------------------------------------------------------
# BaselineFilter tests
add_test("TOPP_BaselineFilter_1" ${TOPP_BIN_PATH}/BaselineFilter -test -in ${DATA_DIR_TOPP}/BaselineFilter_input.mzML -out BaselineFilter.tmp.mzML -struc_elem_length 1.5)
add_test("TOPP_BaselineFilter_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 BaselineFilter.tmp.mzML -in2 ${DATA_DIR_TOPP}/BaselineFilter_output.mzML )
set_tests_properties("TOPP_BaselineFilter_1_out1" PROPERTIES DEPENDS "TOPP_BaselineFilter_1")

#------------------------------------------------------------------------------
# ConsensusMapNormalizer tests
add_test("TOPP_ConsensusMapNormalizer_1" ${TOPP_BIN_PATH}/ConsensusMapNormalizer -test -in ${DATA_DIR_TOPP}/ConsensusMapNormalizer_input.consensusXML -out ConsensusMapNormalizer.tmp.consensusXML)
add_test("TOPP_ConsensusMapNormalizer_1_out1" ${DIFF} -in1 ConsensusMapNormalizer.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/ConsensusMapNormalizer_output.consensusXML )
set_tests_properties("TOPP_ConsensusMapNormalizer_1_out1" PROPERTIES DEPENDS "TOPP_ConsensusMapNormalizer_1")

#------------------------------------------------------------------------------
# MapNormalizer tests
add_test("TOPP_MapNormalizer_1" ${TOPP_BIN_PATH}/MapNormalizer -test -in ${DATA_DIR_TOPP}/MapNormalizer_input.mzML -out MapNormalizer.tmp.mzML)
add_test("TOPP_MapNormalizer_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 MapNormalizer.tmp.mzML -in2 ${DATA_DIR_TOPP}/MapNormalizer_output.mzML )
set_tests_properties("TOPP_MapNormalizer_1_out1" PROPERTIES DEPENDS "TOPP_MapNormalizer_1")

#------------------------------------------------------------------------------
# DatabaseSuitability tests
# test with disabled correction (not depending on Comet)
# remaining tests are found in /THIRDPARTY and are only executed if the Comet executable was found
add_test("TOPP_DatabaseSuitability" ${TOPP_BIN_PATH}/DatabaseSuitability -test -in_id ${DATA_DIR_TOPP}/THIRDPARTY/DatabaseSuitability_in_id.idXML -in_spec ${DATA_DIR_TOPP}/THIRDPARTY/DatabaseSuitability_in_spec.mzML -in_novo ${DATA_DIR_TOPP}/THIRDPARTY/DatabaseSuitability_in_novo.idXML -database ${DATA_DIR_TOPP}/THIRDPARTY/DatabaseSuitability_database.fasta -novo_database ${DATA_DIR_TOPP}/THIRDPARTY/DatabaseSuitability_novo_database.FASTA -algorithm:disable_correction -out DatabaseSuitability.tmp.tsv)
add_test("TOPP_DatabaseSuitability_out" ${DIFF} -in1 DatabaseSuitability.tmp.tsv -in2 ${DATA_DIR_TOPP}/THIRDPARTY/DatabaseSuitability_out.tsv )
set_tests_properties("TOPP_DatabaseSuitability_out" PROPERTIES DEPENDS "TOPP_DatabaseSuitability")

#------------------------------------------------------------------------------
# Decharger tests
add_test("TOPP_Decharger_1" ${TOPP_BIN_PATH}/Decharger -test -in ${DATA_DIR_TOPP}/Decharger_input.featureXML -ini ${DATA_DIR_TOPP}/Decharger_input.ini -out_fm Decharger_output_fm.tmp.featureXML -out_cm Decharger_output.tmp.consensusXML -outpairs Decharger_p_output.tmp.consensusXML)
add_test("TOPP_Decharger_1_out1" ${DIFF} -whitelist "id=" "href=" -in1 Decharger_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/Decharger_output.consensusXML )
set_tests_properties("TOPP_Decharger_1_out1" PROPERTIES DEPENDS "TOPP_Decharger_1" )
add_test("TOPP_Decharger_1_out2" ${DIFF} -whitelist "id=" "href=" "completion_time=" "version=" -in1 Decharger_p_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/Decharger_output_pairs.consensusXML )
set_tests_properties("TOPP_Decharger_1_out2" PROPERTIES DEPENDS "TOPP_Decharger_1" )
add_test("TOPP_Decharger_1_out3" ${DIFF} -whitelist "id=" -in1 Decharger_output_fm.tmp.featureXML -in2 ${DATA_DIR_TOPP}/Decharger_output.featureXML )
set_tests_properties("TOPP_Decharger_1_out3" PROPERTIES DEPENDS "TOPP_Decharger_1" )

#------------------------------------------------------------------------------
# DTAExtractor tests
add_test("TOPP_DTAExtractor_1" ${TOPP_BIN_PATH}/DTAExtractor -test -in ${DATA_DIR_TOPP}/DTAExtractor_1_input.mzML -out DTAExtractor -rt :61)
add_test("TOPP_DTAExtractor_1_out1" ${DIFF} -in1 DTAExtractor_RT60.0.dta -in2 ${DATA_DIR_TOPP}/DTAExtractor_1_output.dta )
set_tests_properties("TOPP_DTAExtractor_1_out1" PROPERTIES DEPENDS "TOPP_DTAExtractor_1")
add_test("TOPP_DTAExtractor_2" ${TOPP_BIN_PATH}/DTAExtractor -test -in ${DATA_DIR_TOPP}/DTAExtractor_1_input.mzML -out DTAExtractor -level 1)
add_test("TOPP_DTAExtractor_2_out1" ${DIFF} -in1 DTAExtractor_RT60.0.dta -in2 ${DATA_DIR_TOPP}/DTAExtractor_2_output.dta )
set_tests_properties("TOPP_DTAExtractor_2_out1" PROPERTIES DEPENDS "TOPP_DTAExtractor_2")
add_test("TOPP_DTAExtractor_3" ${TOPP_BIN_PATH}/DTAExtractor -test -in ${DATA_DIR_TOPP}/DTAExtractor_1_input.mzML -out DTAExtractor -level 2 -mz :1000)
add_test("TOPP_DTAExtractor_3_out1" ${DIFF} -in1 DTAExtractor_RT140.0_MZ5.0.dta -in2 ${DATA_DIR_TOPP}/DTAExtractor_3_output.dta )
set_tests_properties("TOPP_DTAExtractor_3_out1" PROPERTIES DEPENDS "TOPP_DTAExtractor_3")

#------------------------------------------------------------------------------
# MassTraceExtractor test
add_test("TOPP_MassTraceExtractor_1" ${TOPP_BIN_PATH}/MassTraceExtractor -test -ini ${DATA_DIR_TOPP}/MassTraceExtractor.ini -in ${DATA_DIR_TOPP}/MassTraceExtractor_1_input.mzML -out MassTraceExtractor_1.tmp.featureXML)
add_test("TOPP_MassTraceExtractor_1_out" ${DIFF} -whitelist "id=" -in1 MassTraceExtractor_1.tmp.featureXML -in2 ${DATA_DIR_TOPP}/MassTraceExtractor_1_output.featureXML )
set_tests_properties("TOPP_MassTraceExtractor_1_out" PROPERTIES DEPENDS "TOPP_MassTraceExtractor_1")
add_test("TOPP_MassTraceExtractor_2" ${TOPP_BIN_PATH}/MassTraceExtractor -test -ini ${DATA_DIR_TOPP}/MassTraceExtractor_2.ini -in ${DATA_DIR_TOPP}/MassTraceExtractor_1_input.mzML -out MassTraceExtractor_2.tmp.featureXML)
add_test("TOPP_MassTraceExtractor_2_out" ${DIFF} -whitelist "id=" -in1 MassTraceExtractor_2.tmp.featureXML -in2 ${DATA_DIR_TOPP}/MassTraceExtractor_2_output.featureXML )
set_tests_properties("TOPP_MassTraceExtractor_2_out" PROPERTIES DEPENDS "TOPP_MassTraceExtractor_2")

#------------------------------------------------------------------------------
# FeatureFinderMetabo test
add_test("TOPP_FeatureFinderMetabo_1" ${TOPP_BIN_PATH}/FeatureFinderMetabo -test -ini ${DATA_DIR_TOPP}/FeatureFinderMetabo.ini -in ${DATA_DIR_TOPP}/FeatureFinderMetabo_1_input.mzML -out FeatureFinderMetabo_1.tmp.featureXML)
add_test("TOPP_FeatureFinderMetabo_1_out" ${DIFF} -whitelist "id=" -in1 FeatureFinderMetabo_1.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FeatureFinderMetabo_1_output.featureXML)
set_tests_properties("TOPP_FeatureFinderMetabo_1_out" PROPERTIES DEPENDS "TOPP_FeatureFinderMetabo_1")
add_test("TOPP_FeatureFinderMetabo_2" ${TOPP_BIN_PATH}/FeatureFinderMetabo -test -ini ${DATA_DIR_TOPP}/FeatureFinderMetabo_2_noEPD.ini -in ${DATA_DIR_TOPP}/FeatureFinderMetabo_2_input.mzML -out FeatureFinderMetabo_2.tmp.featureXML)
add_test("TOPP_FeatureFinderMetabo_2_out" ${DIFF} -whitelist "id=" "completion_time" -in1 FeatureFinderMetabo_2.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FeatureFinderMetabo_2_noEPD_output.featureXML)
set_tests_properties("TOPP_FeatureFinderMetabo_2_out" PROPERTIES DEPENDS "TOPP_FeatureFinderMetabo_2")
add_test("TOPP_FeatureFinderMetabo_3" ${TOPP_BIN_PATH}/FeatureFinderMetabo -test -ini ${DATA_DIR_TOPP}/FeatureFinderMetabo_3.ini -in ${DATA_DIR_TOPP}/FeatureFinderMetabo_3_input.mzML -out FeatureFinderMetabo_3.tmp.featureXML -out_chrom FeatureFinderMetabo_3_chrom.tmp.mzML)
add_test("TOPP_FeatureFinderMetabo_3_chrom_out" ${DIFF} -whitelist "<chromatogram id=" "<offset" "<userParam name=" -in1 FeatureFinderMetabo_3_chrom.tmp.mzML -in2 ${DATA_DIR_TOPP}/FeatureFinderMetabo_3_chrom_output.mzML)
set_tests_properties("TOPP_FeatureFinderMetabo_3_chrom_out" PROPERTIES DEPENDS "TOPP_FeatureFinderMetabo_3")
add_test("TOPP_FeatureFinderMetabo_4" ${TOPP_BIN_PATH}/FeatureFinderMetabo -test -ini ${DATA_DIR_TOPP}/FeatureFinderMetabo_4.ini -in ${DATA_DIR_TOPP}/FeatureFinderMetabo_3_input.mzML -out FeatureFinderMetabo_4.tmp.featureXML -out_chrom FeatureFinderMetabo_4_chrom.tmp.mzML)
add_test("TOPP_FeatureFinderMetabo_4_out" ${DIFF} -whitelist "id=" -in1 FeatureFinderMetabo_4.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FeatureFinderMetabo_4_output.featureXML)
add_test("TOPP_FeatureFinderMetabo_4_chrom_out" ${DIFF} -whitelist "<chromatogram id=" "<offset" "<userParam name=" -in1 FeatureFinderMetabo_4_chrom.tmp.mzML -in2 ${DATA_DIR_TOPP}/FeatureFinderMetabo_4_chrom_output.mzML)
set_tests_properties("TOPP_FeatureFinderMetabo_4_chrom_out" PROPERTIES DEPENDS "TOPP_FeatureFinderMetabo_4")
add_test("TOPP_FeatureFinderMetabo_5" ${TOPP_BIN_PATH}/FeatureFinderMetabo -test -ini ${DATA_DIR_TOPP}/FeatureFinderMetabo.ini -in ${DATA_DIR_TOPP}/FeatureFinderMetabo_1_input.mzML -out FeatureFinderMetabo_5.tmp.featureXML -algorithm:mtd:quant_method max_height)
add_test("TOPP_FeatureFinderMetabo_5_out" ${DIFF} -whitelist "id=" -in1 FeatureFinderMetabo_5.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FeatureFinderMetabo_5_output.featureXML)
set_tests_properties("TOPP_FeatureFinderMetabo_5_out" PROPERTIES DEPENDS "TOPP_FeatureFinderMetabo_5")


#------------------------------------------------------------------------------
# FeatureFinderCentroided test
add_test("TOPP_FeatureFinderCentroided_1" ${TOPP_BIN_PATH}/FeatureFinderCentroided -test -ini ${DATA_DIR_TOPP}/FeatureFinderCentroided_1_parameters.ini -in ${DATA_DIR_TOPP}/FeatureFinderCentroided_1_input.mzML -out FeatureFinderCentroided_1.tmp.featureXML)
add_test("TOPP_FeatureFinderCentroided_1_out1" ${DIFF} -whitelist "id=" -in1 FeatureFinderCentroided_1.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FeatureFinderCentroided_1_1_output.featureXML )
set_tests_properties("TOPP_FeatureFinderCentroided_1_out1" PROPERTIES DEPENDS "TOPP_FeatureFinderCentroided_1")

#------------------------------------------------------------------------------
# FeatureFinderIdentification test
## internal IDs only:
add_test("TOPP_FeatureFinderIdentification_1" ${TOPP_BIN_PATH}/FeatureFinderIdentification -test -in ${DATA_DIR_TOPP}/FeatureFinderIdentification_1_input.mzML -id ${DATA_DIR_TOPP}/FeatureFinderIdentification_1_input.idXML -out FeatureFinderIdentification_1.tmp.featureXML -extract:mz_window 0.1 -detect:peak_width 60 -model:type none)
add_test("TOPP_FeatureFinderIdentification_1_out1" ${DIFF} -whitelist "spectra_data" "featureMap" -in1 FeatureFinderIdentification_1.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FeatureFinderIdentification_1_output.featureXML)
set_tests_properties("TOPP_FeatureFinderIdentification_1_out1" PROPERTIES DEPENDS "TOPP_FeatureFinderIdentification_1")
## with (faked) external IDs; fix SVM parameters to avoid randomness:
## test currently produces different results on Windows, Mac, and Linux
# add_test("TOPP_FeatureFinderIdentification_2" ${TOPP_BIN_PATH}/FeatureFinderIdentification -test -in ${DATA_DIR_TOPP}/FeatureFinderIdentification_1_input.mzML -id ${DATA_DIR_TOPP}/FeatureFinderIdentification_2_input.idXML -id_ext ${DATA_DIR_TOPP}/FeatureFinderIdentification_1_input.idXML -out FeatureFinderIdentification_2.tmp.featureXML -extract:mz_window 0.1 -detect:peak_width 60 -svm:no_selection -svm:log2_C 1 -svm:log2_gamma 1 -model:type none)
## results are somewhat random, so don't be too strict when comparing:
#add_test("TOPP_FeatureFinderIdentification_2_out1" ${DIFF} -whitelist "FDR_probabilities" "FFId_category" -ratio 1.1 -absdiff 0.03 -in1 FeatureFinderIdentification_2.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FeatureFinderIdentification_2_output.featureXML)
#set_tests_properties("TOPP_FeatureFinderIdentification_2_out1" PROPERTIES DEPENDS "TOPP_FeatureFinderIdentification_2")

## with elution model fitting:
add_test("TOPP_FeatureFinderIdentification_3" ${TOPP_BIN_PATH}/FeatureFinderIdentification -test -in ${DATA_DIR_TOPP}/FeatureFinderIdentification_1_input.mzML -id ${DATA_DIR_TOPP}/FeatureFinderIdentification_1_input.idXML -out FeatureFinderIdentification_3.tmp.featureXML -extract:mz_window 0.1 -detect:peak_width 60 -model:type symmetric)
add_test("TOPP_FeatureFinderIdentification_3_out1" ${DIFF} -whitelist "spectra_data" "featureMap" -in1 FeatureFinderIdentification_3.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FeatureFinderIdentification_3_output.featureXML)
set_tests_properties("TOPP_FeatureFinderIdentification_3_out1" PROPERTIES DEPENDS "TOPP_FeatureFinderIdentification_3")
## elution model fitting for each individual mass trace:
## TODO: reenable - currently not stable on windows in CI
##add_test("TOPP_FeatureFinderIdentification_4" ${TOPP_BIN_PATH}/FeatureFinderIdentification -test -in ${DATA_DIR_TOPP}/FeatureFinderIdentification_1_input.mzML -id ${DATA_DIR_TOPP}/FeatureFinderIdentification_1_input.idXML -out FeatureFinderIdentification_4.tmp.featureXML -extract:mz_window 0.1 -detect:peak_width 60 -model:type symmetric -model:each_trace)
##add_test("TOPP_FeatureFinderIdentification_4_out1" ${DIFF} -whitelist "spectra_data" -in1 FeatureFinderIdentification_4.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FeatureFinderIdentification_4_output.featureXML)
##set_tests_properties("TOPP_FeatureFinderIdentification_4_out1" PROPERTIES DEPENDS "TOPP_FeatureFinderIdentification_4")

## with batch extraction size smaller than the nr of peptides (we need to whitelist feature ids, since they might be generated differently)
add_test("TOPP_FeatureFinderIdentification_5" ${TOPP_BIN_PATH}/FeatureFinderIdentification -test -in ${DATA_DIR_TOPP}/FeatureFinderIdentification_1_input.mzML -id ${DATA_DIR_TOPP}/FeatureFinderIdentification_1_input.idXML -out FeatureFinderIdentification_5.tmp.featureXML -candidates_out FeatureFinderIdentification_5_candidates.tmp.featureXML -extract:mz_window 0.1 -extract:batch_size 10 -detect:peak_width 60 -model:type none)
add_test("TOPP_FeatureFinderIdentification_5_out1" ${DIFF} -whitelist "feature id" "spectra_data" "featureMap" -in1 FeatureFinderIdentification_5.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FeatureFinderIdentification_1_output.featureXML)
set_tests_properties("TOPP_FeatureFinderIdentification_5_out1" PROPERTIES DEPENDS "TOPP_FeatureFinderIdentification_5")
add_test("TOPP_FeatureFinderIdentification_5_out2" ${DIFF} -whitelist "feature id" "spectra_data" "featureMap" -in1 FeatureFinderIdentification_5_candidates.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FeatureFinderIdentification_5_candidates.featureXML)
set_tests_properties("TOPP_FeatureFinderIdentification_5_out2" PROPERTIES DEPENDS "TOPP_FeatureFinderIdentification_5")



#------------------------------------------------------------------------------
# FeatureFinderMRM test
add_test("TOPP_FeatureFinderMRM_1" ${TOPP_BIN_PATH}/FeatureFinderMRM -test -in ${DATA_DIR_TOPP}/FeatureFinderMRM_1_input.mzML -ini ${DATA_DIR_TOPP}/FeatureFinderMRM_1_parameters.ini -out FeatureFinderMRM_1.tmp.featureXML)
add_test("TOPP_FeatureFinderMRM_1_out1" ${DIFF} -whitelist "id=" -in1 FeatureFinderMRM_1.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FeatureFinderMRM_1_output.featureXML )
set_tests_properties("TOPP_FeatureFinderMRM_1_out1" PROPERTIES DEPENDS "TOPP_FeatureFinderMRM_1")

#------------------------------------------------------------------------------
# FeatureFinderMultiplex test
add_test("TOPP_FeatureFinderMultiplex_1" ${TOPP_BIN_PATH}/FeatureFinderMultiplex -test -in ${DATA_DIR_TOPP}/FeatureFinderMultiplex_1_input.mzML -ini ${DATA_DIR_TOPP}/FeatureFinderMultiplex_1_parameters.ini -out FeatureFinderMultiplex_1.tmp.featureXML -out_multiplets FeatureFinderMultiplex_2.tmp.consensusXML)
add_test("TOPP_FeatureFinderMultiplex_1_out1" ${DIFF} -whitelist "id=" "href=" "spectra_data" -in1 FeatureFinderMultiplex_1.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_1_output.featureXML )
set_tests_properties("TOPP_FeatureFinderMultiplex_1_out1" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_1")
add_test("TOPP_FeatureFinderMultiplex_1_out2" ${DIFF} -whitelist "id=" "href=" -in1 FeatureFinderMultiplex_2.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_1_output.consensusXML )
set_tests_properties("TOPP_FeatureFinderMultiplex_1_out2" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_1")

add_test("TOPP_FeatureFinderMultiplex_2" ${TOPP_BIN_PATH}/FeatureFinderMultiplex -test -in ${DATA_DIR_TOPP}/FeatureFinderMultiplex_2_input.mzML -ini ${DATA_DIR_TOPP}/FeatureFinderMultiplex_2_parameters.ini -out FeatureFinderMultiplex_3.tmp.featureXML -out_multiplets FeatureFinderMultiplex_4.tmp.consensusXML)
add_test("TOPP_FeatureFinderMultiplex_2_out1" ${DIFF} -whitelist "id=" "href=" "spectra_data" -in1 FeatureFinderMultiplex_3.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_2_output.featureXML )
set_tests_properties("TOPP_FeatureFinderMultiplex_2_out1" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_2")
add_test("TOPP_FeatureFinderMultiplex_2_out2" ${DIFF} -whitelist "id=" "href=" -in1 FeatureFinderMultiplex_4.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_2_output.consensusXML )
set_tests_properties("TOPP_FeatureFinderMultiplex_2_out2" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_2")

add_test("TOPP_FeatureFinderMultiplex_3" ${TOPP_BIN_PATH}/FeatureFinderMultiplex -test -in ${DATA_DIR_TOPP}/FeatureFinderMultiplex_3_input.mzML -ini ${DATA_DIR_TOPP}/FeatureFinderMultiplex_3_parameters.ini -out FeatureFinderMultiplex_5.tmp.featureXML -out_multiplets FeatureFinderMultiplex_6.tmp.consensusXML)
add_test("TOPP_FeatureFinderMultiplex_3_out1" ${DIFF} -whitelist "id=" "href=" "spectra_data" -in1 FeatureFinderMultiplex_5.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_3_output.featureXML )
set_tests_properties("TOPP_FeatureFinderMultiplex_3_out1" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_3")
add_test("TOPP_FeatureFinderMultiplex_3_out2" ${DIFF} -whitelist "id=" "href=" -in1 FeatureFinderMultiplex_6.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_3_output.consensusXML )
set_tests_properties("TOPP_FeatureFinderMultiplex_3_out2" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_3")

add_test("TOPP_FeatureFinderMultiplex_4" ${TOPP_BIN_PATH}/FeatureFinderMultiplex -test -in ${DATA_DIR_TOPP}/FeatureFinderMultiplex_4_input.mzML -ini ${DATA_DIR_TOPP}/FeatureFinderMultiplex_4_parameters.ini -out FeatureFinderMultiplex_7.tmp.featureXML -out_multiplets FeatureFinderMultiplex_8.tmp.consensusXML)
add_test("TOPP_FeatureFinderMultiplex_4_out1" ${DIFF} -whitelist "id=" "href=" "spectra_data" -in1 FeatureFinderMultiplex_7.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_4_output.featureXML )
set_tests_properties("TOPP_FeatureFinderMultiplex_4_out1" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_4")
add_test("TOPP_FeatureFinderMultiplex_4_out2" ${DIFF} -whitelist "id=" "href=" -in1 FeatureFinderMultiplex_8.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_4_output.consensusXML )
set_tests_properties("TOPP_FeatureFinderMultiplex_4_out2" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_4")

add_test("TOPP_FeatureFinderMultiplex_5" ${TOPP_BIN_PATH}/FeatureFinderMultiplex -test -in ${DATA_DIR_TOPP}/FeatureFinderMultiplex_5_input.mzML -ini ${DATA_DIR_TOPP}/FeatureFinderMultiplex_5_parameters.ini -out FeatureFinderMultiplex_9.tmp.featureXML -out_multiplets FeatureFinderMultiplex_10.tmp.consensusXML)
add_test("TOPP_FeatureFinderMultiplex_5_out1" ${DIFF} -whitelist "id=" "href=" "spectra_data" -in1 FeatureFinderMultiplex_9.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_5_output.featureXML )
set_tests_properties("TOPP_FeatureFinderMultiplex_5_out1" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_5")
add_test("TOPP_FeatureFinderMultiplex_5_out2" ${DIFF} -whitelist "id=" "href=" -in1 FeatureFinderMultiplex_10.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_5_output.consensusXML )
set_tests_properties("TOPP_FeatureFinderMultiplex_5_out2" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_5")

add_test("TOPP_FeatureFinderMultiplex_6" ${TOPP_BIN_PATH}/FeatureFinderMultiplex -test -in ${DATA_DIR_TOPP}/FeatureFinderMultiplex_6_input.mzML -ini ${DATA_DIR_TOPP}/FeatureFinderMultiplex_6_parameters.ini -out FeatureFinderMultiplex_11.tmp.featureXML -out_multiplets FeatureFinderMultiplex_12.tmp.consensusXML)
add_test("TOPP_FeatureFinderMultiplex_6_out1" ${DIFF} -whitelist "id=" "href=" "spectra_data" -in1 FeatureFinderMultiplex_11.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_6_output.featureXML )
set_tests_properties("TOPP_FeatureFinderMultiplex_6_out1" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_6")
add_test("TOPP_FeatureFinderMultiplex_6_out2" ${DIFF} -whitelist "id=" "href=" -in1 FeatureFinderMultiplex_12.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_6_output.consensusXML )
set_tests_properties("TOPP_FeatureFinderMultiplex_6_out2" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_6")

add_test("TOPP_FeatureFinderMultiplex_7" ${TOPP_BIN_PATH}/FeatureFinderMultiplex -test -in ${DATA_DIR_TOPP}/FeatureFinderMultiplex_7_input.mzML -ini ${DATA_DIR_TOPP}/FeatureFinderMultiplex_7_parameters.ini -out FeatureFinderMultiplex_13.tmp.featureXML -out_multiplets FeatureFinderMultiplex_14.tmp.consensusXML)
add_test("TOPP_FeatureFinderMultiplex_7_out1" ${DIFF} -whitelist "id=" "href=" "spectra_data" -in1 FeatureFinderMultiplex_13.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_7_output.featureXML )
set_tests_properties("TOPP_FeatureFinderMultiplex_7_out1" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_7")
add_test("TOPP_FeatureFinderMultiplex_7_out2" ${DIFF} -whitelist "id=" "href=" -in1 FeatureFinderMultiplex_14.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_7_output.consensusXML )

add_test("TOPP_FeatureFinderMultiplex_8" ${TOPP_BIN_PATH}/FeatureFinderMultiplex -test -in ${DATA_DIR_TOPP}/FeatureFinderMultiplex_8_input.mzML -ini ${DATA_DIR_TOPP}/FeatureFinderMultiplex_8_parameters.ini -out FeatureFinderMultiplex_15.tmp.featureXML -out_multiplets FeatureFinderMultiplex_16.tmp.consensusXML)
add_test("TOPP_FeatureFinderMultiplex_8_out1" ${DIFF} -whitelist "id=" "href=" "spectra_data" -in1 FeatureFinderMultiplex_15.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_8_output.featureXML )
set_tests_properties("TOPP_FeatureFinderMultiplex_8_out1" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_8")
add_test("TOPP_FeatureFinderMultiplex_8_out2" ${DIFF} -whitelist "id=" "href=" -in1 FeatureFinderMultiplex_16.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_8_output.consensusXML )
set_tests_properties("TOPP_FeatureFinderMultiplex_8_out2" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_8")

add_test("TOPP_FeatureFinderMultiplex_9" ${TOPP_BIN_PATH}/FeatureFinderMultiplex -test -in ${DATA_DIR_TOPP}/FeatureFinderMultiplex_9_input.mzML -ini ${DATA_DIR_TOPP}/FeatureFinderMultiplex_9_parameters.ini -out FeatureFinderMultiplex_17.tmp.featureXML -out_multiplets FeatureFinderMultiplex_18.tmp.consensusXML)
add_test("TOPP_FeatureFinderMultiplex_9_out1" ${DIFF} -whitelist "id=" "href=" "spectra_data" -in1 FeatureFinderMultiplex_17.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_9_output.featureXML )
set_tests_properties("TOPP_FeatureFinderMultiplex_9_out1" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_9")
add_test("TOPP_FeatureFinderMultiplex_9_out2" ${DIFF} -whitelist "id=" "href=" -in1 FeatureFinderMultiplex_18.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_9_output.consensusXML )
set_tests_properties("TOPP_FeatureFinderMultiplex_9_out2" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_9")

add_test("TOPP_FeatureFinderMultiplex_10" ${TOPP_BIN_PATH}/FeatureFinderMultiplex -test -in ${DATA_DIR_TOPP}/FeatureFinderMultiplex_10_input.mzML -ini ${DATA_DIR_TOPP}/FeatureFinderMultiplex_10_parameters.ini -out FeatureFinderMultiplex_19.tmp.featureXML -out_multiplets FeatureFinderMultiplex_20.tmp.consensusXML)
add_test("TOPP_FeatureFinderMultiplex_10_out1" ${DIFF} -whitelist "id=" "href=" "spectra_data" -in1 FeatureFinderMultiplex_19.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_10_output.featureXML )
set_tests_properties("TOPP_FeatureFinderMultiplex_10_out1" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_10")
add_test("TOPP_FeatureFinderMultiplex_10_out2" ${DIFF} -whitelist "id=" "href=" -in1 FeatureFinderMultiplex_20.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_10_output.consensusXML )
set_tests_properties("TOPP_FeatureFinderMultiplex_10_out2" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_10")

add_test("TOPP_FeatureFinderMultiplex_11" ${TOPP_BIN_PATH}/FeatureFinderMultiplex -test -in ${DATA_DIR_TOPP}/FeatureFinderMultiplex_11_input.mzML -ini ${DATA_DIR_TOPP}/FeatureFinderMultiplex_11_parameters.ini -out FeatureFinderMultiplex_21.tmp.featureXML)
add_test("TOPP_FeatureFinderMultiplex_11_out1" ${DIFF} -whitelist "id=" "href=" "spectra_data" -in1 FeatureFinderMultiplex_21.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FeatureFinderMultiplex_11_output.featureXML )
set_tests_properties("TOPP_FeatureFinderMultiplex_11_out1" PROPERTIES DEPENDS "TOPP_FeatureFinderMultiplex_11")

#------------------------------------------------------------------------------
# FileConverter tests
add_test("TOPP_FileConverter_1" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_1_input.mzData -out FileConverter_1.tmp.mzML -out_type mzML)
add_test("TOPP_FileConverter_1_out1" ${DIFF} -whitelist "id=" "location=" "href=" -in1 FileConverter_1.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileConverter_1_output.mzML )
set_tests_properties("TOPP_FileConverter_1_out1" PROPERTIES DEPENDS "TOPP_FileConverter_1")
add_test("TOPP_FileConverter_2" ${TOPP_BIN_PATH}/FileConverter -test  -in ${DATA_DIR_TOPP}/FileConverter_2_input.someInputDTA2D -in_type dta2d -out FileConverter_2.tmp.mzML -out_type mzML)
add_test("TOPP_FileConverter_2_out1" ${DIFF} -whitelist "id=" "location=" "href=" "<offset" -in1 FileConverter_2.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileConverter_2_output.mzML )
set_tests_properties("TOPP_FileConverter_2_out1" PROPERTIES DEPENDS "TOPP_FileConverter_2")
add_test("TOPP_FileConverter_3" ${TOPP_BIN_PATH}/FileConverter -test  -in ${DATA_DIR_TOPP}/FileConverter_3_input.featureXML -out FileConverter_3.tmp.mzML -out_type mzML)
add_test("TOPP_FileConverter_3_out1" ${DIFF} -in1 FileConverter_3.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileConverter_3_output.mzML )
set_tests_properties("TOPP_FileConverter_3_out1" PROPERTIES DEPENDS "TOPP_FileConverter_3")
add_test("TOPP_FileConverter_4" ${TOPP_BIN_PATH}/FileConverter -test  -in ${DATA_DIR_TOPP}/FileConverter_4_input.mzXML -out FileConverter_4.tmp.mzML -out_type mzML)
add_test("TOPP_FileConverter_4_out1" ${DIFF} -whitelist "location=" "MS:100056" "SHA-1" -in1 FileConverter_4.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileConverter_4_output.mzML )
set_tests_properties("TOPP_FileConverter_4_out1" PROPERTIES DEPENDS "TOPP_FileConverter_4")
add_test("TOPP_FileConverter_5" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_5_input.mzML -out FileConverter_5.tmp.mzXML -out_type mzXML)
add_test("TOPP_FileConverter_5_out1" ${DIFF} -in1 FileConverter_5.tmp.mzXML -in2 ${DATA_DIR_TOPP}/FileConverter_5_output.mzXML )
set_tests_properties("TOPP_FileConverter_5_out1" PROPERTIES DEPENDS "TOPP_FileConverter_5")
add_test("TOPP_FileConverter_6" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_6_input.mzML -out FileConverter_6.tmp.mzXML -out_type mzXML)
add_test("TOPP_FileConverter_6_out1" ${DIFF} -in1 FileConverter_6.tmp.mzXML -in2 ${DATA_DIR_TOPP}/FileConverter_6_output.mzXML )
set_tests_properties("TOPP_FileConverter_6_out1" PROPERTIES DEPENDS "TOPP_FileConverter_6")
add_test("TOPP_FileConverter_7" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_7_input.consensusXML -out FileConverter_7.tmp.featureXML -out_type featureXML)
add_test("TOPP_FileConverter_7_out1" ${DIFF} -in1 FileConverter_7.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FileConverter_7_output.featureXML )
set_tests_properties("TOPP_FileConverter_7_out1" PROPERTIES DEPENDS "TOPP_FileConverter_7")
add_test("TOPP_FileConverter_8" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_8_input.mzML -out FileConverter_8.tmp.mzData -out_type mzData)
add_test("TOPP_FileConverter_8_out1" ${DIFF} -whitelist "<pathToFile>" -in1 FileConverter_8.tmp.mzData -in2 ${DATA_DIR_TOPP}/FileConverter_8_output.mzData )
set_tests_properties("TOPP_FileConverter_8_out1" PROPERTIES DEPENDS "TOPP_FileConverter_8")
add_test("TOPP_FileConverter_9" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_9_input.consensusXML -out FileConverter_9.tmp.featureXML -out_type featureXML)
add_test("TOPP_FileConverter_9_out1" ${DIFF} -in1 FileConverter_9.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FileConverter_9_output.featureXML )
set_tests_properties("TOPP_FileConverter_9_out1" PROPERTIES DEPENDS "TOPP_FileConverter_9")
add_test("TOPP_FileConverter_10" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_10_input.edta -no_progress -out FileConverter_10.tmp.featureXML -out_type featureXML)
add_test("TOPP_FileConverter_10_out1" ${DIFF} -whitelist "id=" -in1 FileConverter_10.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FileConverter_10_output.featureXML )
set_tests_properties("TOPP_FileConverter_10_out1" PROPERTIES DEPENDS "TOPP_FileConverter_10")
#SpecArray
add_test("TOPP_FileConverter_11" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_11_input.peplist -no_progress -out FileConverter_11.tmp.featureXML -out_type featureXML)
add_test("TOPP_FileConverter_11_out1" ${DIFF} -whitelist "id=" -in1 FileConverter_11.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FileConverter_11_output.featureXML )
set_tests_properties("TOPP_FileConverter_11_out1" PROPERTIES DEPENDS "TOPP_FileConverter_11")
#msInspect
add_test("TOPP_FileConverter_12" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_12_input.peptides.tsv -no_progress -out FileConverter_12.tmp.featureXML -out_type featureXML)
add_test("TOPP_FileConverter_12_out1" ${DIFF} -whitelist "id=" -in1 FileConverter_12.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FileConverter_12_output.featureXML )
set_tests_properties("TOPP_FileConverter_12_out1" PROPERTIES DEPENDS "TOPP_FileConverter_12")
#kroenik
add_test("TOPP_FileConverter_13" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_13_input.peptides.kroenik -no_progress -out FileConverter_13.tmp.featureXML -out_type featureXML)
add_test("TOPP_FileConverter_13_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileConverter_13.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FileConverter_13_output.featureXML )
set_tests_properties("TOPP_FileConverter_13_out1" PROPERTIES DEPENDS "TOPP_FileConverter_13")
# featureXML to consensusXML
add_test("TOPP_FileConverter_14" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_9_output.featureXML -no_progress -out FileConverter_14.tmp.consensusXML -out_type consensusXML)
add_test("TOPP_FileConverter_14_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileConverter_14.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FileConverter_14_output.consensusXML )
set_tests_properties("TOPP_FileConverter_14_out1" PROPERTIES DEPENDS "TOPP_FileConverter_14")
# EDTA to consensusXML
add_test("TOPP_FileConverter_15" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_10_input.edta -no_progress -out FileConverter_15.tmp.consensusXML -out_type consensusXML)
add_test("TOPP_FileConverter_15_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileConverter_15.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FileConverter_15_output.consensusXML )
set_tests_properties("TOPP_FileConverter_15_out1" PROPERTIES DEPENDS "TOPP_FileConverter_15")
# consensus EDTA to consensusXML
add_test("TOPP_FileConverter_16" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_16_input.edta -no_progress -out FileConverter_16.tmp.consensusXML -out_type consensusXML)
add_test("TOPP_FileConverter_16_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileConverter_16.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FileConverter_16_output.consensusXML )
set_tests_properties("TOPP_FileConverter_16_out1" PROPERTIES DEPENDS "TOPP_FileConverter_16")
# consensusXML to ibspectra without out_type
add_test("TOPP_FileConverter_17" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_17_input.consensusXML -no_progress -out FileConverter_17.csv)
add_test("TOPP_FileConverter_17_out1" ${DIFF} -in1 FileConverter_17.csv -in2 ${DATA_DIR_TOPP}/FileConverter_17_output.ibspectra.csv )
set_tests_properties("TOPP_FileConverter_17_out1" PROPERTIES DEPENDS "TOPP_FileConverter_17")
# consensusXML to ibspectra with out_type -> note that input and expected output are the same as in FileConverter_17
add_test("TOPP_FileConverter_18" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_17_input.consensusXML -no_progress -out FileConverter_18.tmp.ibspectra.csv -out_type csv)
add_test("TOPP_FileConverter_18_out1" ${DIFF} -in1 FileConverter_18.tmp.ibspectra.csv -in2 ${DATA_DIR_TOPP}/FileConverter_17_output.ibspectra.csv )
set_tests_properties("TOPP_FileConverter_18_out1" PROPERTIES DEPENDS "TOPP_FileConverter_18")
# Purpose: test the -process_lowmemory options (for mzML and mzXML)
add_test("TOPP_FileConverter_19" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out FileConverter_19.tmp.mzML -process_lowmemory -in_type mzML -out_type mzML)
add_test("TOPP_FileConverter_19_out" ${DIFF} -in1 FileConverter_19.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileConverter_19_output.mzML )
set_tests_properties("TOPP_FileConverter_19_out" PROPERTIES DEPENDS "TOPP_FileConverter_19")
# Purpose: test the writing and reading back the same featureXML
add_test("TOPP_FileConverter_20" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_20_input.featureXML -out FileConverter_20.tmp.featureXML -in_type featureXML -out_type featureXML)
add_test("TOPP_FileConverter_20_out" ${DIFF} -in1 FileConverter_20.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FileConverter_20_output.featureXML )
set_tests_properties("TOPP_FileConverter_20_out" PROPERTIES DEPENDS "TOPP_FileConverter_20")
add_test("TOPP_FileConverter_21" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_4_input.mzXML -out FileConverter_21.tmp.mzML -out_type mzML -process_lowmemory)
# TODO sourceFile does not seem to be handled correctly, we should fix this!
add_test("TOPP_FileConverter_21_out" ${DIFF} -whitelist "location=" "MS:100056" "SHA-1" "<offset" "<indexListOffset" -in1 FileConverter_21.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileConverter_4_output.mzML)
set_tests_properties("TOPP_FileConverter_21_out" PROPERTIES DEPENDS "TOPP_FileConverter_21")
add_test("TOPP_FileConverter_22" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_8_input.mzML -out FileConverter_22.tmp.mzData -out_type mzData -process_lowmemory)
set_tests_properties("TOPP_FileConverter_22" PROPERTIES WILL_FAIL 1) ## low memory not implemented yet
# Corrupt file with empty spectra and chromatograms, it should be handled gracefully
add_test("TOPP_FileConverter_23" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_23_input.mzML -out FileConverter_23.tmp.mzML -out_type mzML)
add_test("TOPP_FileConverter_23_out" ${DIFF} -in1 FileConverter_23.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileConverter_23_output.mzML -whitelist "location=")
set_tests_properties("TOPP_FileConverter_23_out" PROPERTIES DEPENDS "TOPP_FileConverter_23")
# Purpose: test the numpress conversion options (using full accuracy and 1ppm accuracy @ 100 m/z)
add_test("TOPP_FileConverter_24" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out FileConverter_24.tmp.mzML -process_lowmemory -in_type mzML -out_type mzML -lossy_compression -lossy_mass_accuracy 0.0001)
add_test("TOPP_FileConverter_24_out" ${DIFF} -in1 FileConverter_24.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileConverter_24_output.mzML )
set_tests_properties("TOPP_FileConverter_24_out" PROPERTIES DEPENDS "TOPP_FileConverter_24")
add_test("TOPP_FileConverter_25" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out FileConverter_25.tmp.mzML -process_lowmemory -in_type mzML -out_type mzML -lossy_compression)
add_test("TOPP_FileConverter_25_out" ${DIFF} -in1 FileConverter_25.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileConverter_25_output.mzML )
set_tests_properties("TOPP_FileConverter_25_out" PROPERTIES DEPENDS "TOPP_FileConverter_25")
# Purpose: test mzXML conversion with MaxQuant compatibility options (should add a scan index, Thermo manufacturer and other things)
add_test("TOPP_FileConverter_26" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out FileConverter_26.tmp.mzXML -force_MaxQuant_compatibility -out_type mzXML)
add_test("TOPP_FileConverter_26_out" ${DIFF} -in1 FileConverter_26.tmp.mzXML -in2 ${DATA_DIR_TOPP}/FileConverter_26_output.mzXML )
set_tests_properties("TOPP_FileConverter_26_out" PROPERTIES DEPENDS "TOPP_FileConverter_26")
#  .. and self conversion should not loose anything
add_test("TOPP_FileConverter_26-2" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_26_output.mzXML -out FileConverter_26-2.tmp.mzXML -force_MaxQuant_compatibility -out_type mzXML)
add_test("TOPP_FileConverter_26-2_out" ${DIFF} -in1 FileConverter_26-2.tmp.mzXML -in2 ${DATA_DIR_TOPP}/FileConverter_26-2_output.mzXML )
set_tests_properties("TOPP_FileConverter_26-2_out" PROPERTIES DEPENDS "TOPP_FileConverter_26-2")
# Purpose: test conversion of targeted MS1 (e.g. GC-MS data) to chromatograms
add_test("TOPP_FileConverter_27" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_27_input.mzML -out FileConverter_27.tmp.mzML -out_type mzML  -convert_to_chromatograms)
add_test("TOPP_FileConverter_27_out" ${DIFF} -in1 FileConverter_27.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileConverter_27_output.mzML -whitelist "location=" "<offset")
set_tests_properties("TOPP_FileConverter_27_out" PROPERTIES DEPENDS "TOPP_FileConverter_27")
# Purpose: test conversion of chromatogram data in minutes to data in seconds
add_test("TOPP_FileConverter_28" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_28_input.mzML -out FileConverter_28.tmp.mzML -out_type mzML  -convert_to_chromatograms)
add_test("TOPP_FileConverter_28_out" ${DIFF} -in1 FileConverter_28.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileConverter_28_output.mzML -whitelist "location=" "<offset")
set_tests_properties("TOPP_FileConverter_28_out" PROPERTIES DEPENDS "TOPP_FileConverter_28")
# add test for lossy FDA compression
add_test("TOPP_FileConverter_29" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_17_input.mzML -out FileConverter_29.tmp.mzML -out_type mzML -lossy_compression -lossy_mass_accuracy 1e-5 -process_lowmemory)
add_test("TOPP_FileConverter_29_out" ${DIFF} -in1 FileConverter_29.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileConverter_29_output.mzML -whitelist "location=" "<offset")
add_test("TOPP_FileConverter_29_back" ${TOPP_BIN_PATH}/FileConverter -test -in FileConverter_29.tmp.mzML -in_type mzML -out FileConverter_29.back.tmp.mzML -out_type mzML)
set_tests_properties("TOPP_FileConverter_29_out" PROPERTIES DEPENDS "TOPP_FileConverter_29")
set_tests_properties("TOPP_FileConverter_29_back" PROPERTIES DEPENDS "TOPP_FileConverter_29")
# Test with ion mobility
add_test("TOPP_FileConverter_30" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_30_input.mzML -change_im_format multiple_spectra -out_type mzML -out FileConverter_30.tmp.mzML)
add_test("TOPP_FileConverter_30_out" ${DIFF} -in1 FileConverter_30.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileConverter_30_output.mzML -whitelist "location=" "<offset")
set_tests_properties("TOPP_FileConverter_30_out" PROPERTIES DEPENDS "TOPP_FileConverter_30")
add_test("TOPP_FileConverter_31" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_30_output.mzML -change_im_format concatenated -out_type mzML -out FileConverter_31.tmp.mzML)
add_test("TOPP_FileConverter_31_out" ${DIFF} -in1 FileConverter_31.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileConverter_31_output.mzML -whitelist "location=" "<offset")
set_tests_properties("TOPP_FileConverter_31_out" PROPERTIES DEPENDS "TOPP_FileConverter_31")
add_test("TOPP_FileConverter_32" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/FileConverter_32_input.mzML -out_type mzML -out FileConverter_32.tmp.mzML)
add_test("TOPP_FileConverter_32_out" ${DIFF} -in1 FileConverter_32.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileConverter_32_output.mzML -whitelist "location=" "<offset")
set_tests_properties("TOPP_FileConverter_32_out" PROPERTIES DEPENDS "TOPP_FileConverter_32")

#------------------------------------------------------------------------------
# FileFilter tests
add_test("TOPP_FileFilter_1" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out FileFilter_1.tmp.mzML -rt :30 -mz :1000 -int :20000 -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_1_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_1.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_1_output.mzML )
set_tests_properties("TOPP_FileFilter_1_out1" PROPERTIES DEPENDS "TOPP_FileFilter_1")
add_test("TOPP_FileFilter_2" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out FileFilter_2.tmp.mzML -rt 30: -mz 1000: -int 100: -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_2_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_2.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_2_output.mzML )
set_tests_properties("TOPP_FileFilter_2_out1" PROPERTIES DEPENDS "TOPP_FileFilter_2")
add_test("TOPP_FileFilter_3" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out FileFilter_3.tmp.mzML -peak_options:level 2 -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_3_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_3.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_3_output.mzML )
set_tests_properties("TOPP_FileFilter_3_out1" PROPERTIES DEPENDS "TOPP_FileFilter_3")
add_test("TOPP_FileFilter_4" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_4_input.mzML -out FileFilter_4.tmp.mzML -spectra:remove_zoom -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_4_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_4.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_4_output.mzML )
set_tests_properties("TOPP_FileFilter_4_out1" PROPERTIES DEPENDS "TOPP_FileFilter_4")
add_test("TOPP_FileFilter_5" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_5_input.featureXML -out FileFilter_5.tmp.featureXML -rt :1000 -mz :480 -int :79000 -f_and_c:charge :3 -feature:q :0.6 -in_type featureXML -out_type featureXML)
add_test("TOPP_FileFilter_5_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_5.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FileFilter_5_out.featureXML )
set_tests_properties("TOPP_FileFilter_5_out1" PROPERTIES DEPENDS "TOPP_FileFilter_5")
add_test("TOPP_FileFilter_6" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_5_input.featureXML -out FileFilter_6.tmp.featureXML -rt 1000: -mz 440: -int 70000: -f_and_c:charge 3: -feature:q 0.51: -in_type featureXML -out_type featureXML)
add_test("TOPP_FileFilter_6_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_6.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FileFilter_6_out.featureXML )
set_tests_properties("TOPP_FileFilter_6_out1" PROPERTIES DEPENDS "TOPP_FileFilter_6")
add_test("TOPP_FileFilter_7" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_7_input.mzML -out FileFilter_7.tmp.mzML -int 7000: -peak_options:level 1 2 3 -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_7_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_7.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_7_output.mzML )
set_tests_properties("TOPP_FileFilter_7_out1" PROPERTIES DEPENDS "TOPP_FileFilter_7")
add_test("TOPP_FileFilter_8" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_8_input.consensusXML -out FileFilter_8.tmp.consensusXML -rt 600:1400 -mz 700:2300 -int 1100:6000 -in_type consensusXML -out_type consensusXML)
add_test("TOPP_FileFilter_8_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_8.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FileFilter_8_output.consensusXML )
set_tests_properties("TOPP_FileFilter_8_out1" PROPERTIES DEPENDS "TOPP_FileFilter_8")
add_test("TOPP_FileFilter_9" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_9_input.mzML -out FileFilter_9.tmp.mzML -spectra:remove_mode SelectedIonMonitoring -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_9_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_9.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_9_output.mzML)
set_tests_properties("TOPP_FileFilter_9_out1" PROPERTIES DEPENDS "TOPP_FileFilter_9")
add_test("TOPP_FileFilter_10" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_10_input.mzML -out FileFilter_10.tmp.mzML -spectra:remove_activation "Collision-induced dissociation" -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_10_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_10.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_10_output.mzML )
set_tests_properties("TOPP_FileFilter_10_out1" PROPERTIES DEPENDS "TOPP_FileFilter_10")
add_test("TOPP_FileFilter_11" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_11_input.mzML -out FileFilter_11.tmp.mzML -spectra:remove_activation "Plasma desorption" -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_11_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_11.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_11_output.mzML )
set_tests_properties("TOPP_FileFilter_11_out1" PROPERTIES DEPENDS "TOPP_FileFilter_11")
add_test("TOPP_FileFilter_12" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_12_input.mzML -out FileFilter_12.tmp.mzML -peak_options:remove_chromatograms -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_12_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_12.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_12_output.mzML )
set_tests_properties("TOPP_FileFilter_12_out1" PROPERTIES DEPENDS "TOPP_FileFilter_12")
add_test("TOPP_FileFilter_13" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_13_input.consensusXML -out FileFilter_13.tmp.featureXML -in_type consensusXML -out_type featureXML -consensus:map 2)
add_test("TOPP_FileFilter_13_out1" ${DIFF} -whitelist "id=" "href="  -in1 FileFilter_13.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FileFilter_13_output.featureXML )
set_tests_properties("TOPP_FileFilter_13_out1" PROPERTIES DEPENDS "TOPP_FileFilter_13")
add_test("TOPP_FileFilter_14" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_14_input.consensusXML -out FileFilter_14.tmp.consensusXML -in_type consensusXML -out_type consensusXML -consensus:map 0 2)
add_test("TOPP_FileFilter_14_out1" ${DIFF} -whitelist "id=" "href="  -in1 FileFilter_14.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FileFilter_14_output.consensusXML )
set_tests_properties("TOPP_FileFilter_14_out1" PROPERTIES DEPENDS "TOPP_FileFilter_14")
add_test("TOPP_FileFilter_15" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_15_input.featureXML -out FileFilter_15.tmp.featureXML -id:sequences_whitelist Oxidation -id:remove_unassigned_ids)
add_test("TOPP_FileFilter_15_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_15.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FileFilter_15_output.featureXML )
set_tests_properties("TOPP_FileFilter_15_out1" PROPERTIES DEPENDS "TOPP_FileFilter_15")
add_test("TOPP_FileFilter_16" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_15_input.featureXML -out FileFilter_16.tmp.featureXML -id:sequences_whitelist Oxidation -id:remove_unassigned_ids -mz 400:600 -rt 3000:4000)
add_test("TOPP_FileFilter_16_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_16.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FileFilter_16_output.featureXML )
set_tests_properties("TOPP_FileFilter_16_out1" PROPERTIES DEPENDS "TOPP_FileFilter_16")
add_test("TOPP_FileFilter_17" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_15_input.featureXML -out FileFilter_17.tmp.featureXML -id:remove_annotated_features -mz 400:600 -rt 3000:4000)
add_test("TOPP_FileFilter_17_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_17.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FileFilter_17_output.featureXML )
set_tests_properties("TOPP_FileFilter_17_out1" PROPERTIES DEPENDS "TOPP_FileFilter_17")
add_test("TOPP_FileFilter_18" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_18_input.consensusXML -out FileFilter_18.tmp.consensusXML -id:sequences_whitelist Oxidation -id:remove_unassigned_ids)
add_test("TOPP_FileFilter_18_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_18.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FileFilter_18_output.consensusXML )
set_tests_properties("TOPP_FileFilter_18_out1" PROPERTIES DEPENDS "TOPP_FileFilter_18")
add_test("TOPP_FileFilter_19" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_18_input.consensusXML -out FileFilter_19.tmp.consensusXML -id:remove_unannotated_features -mz 400:600 -rt 3000:4000)
add_test("TOPP_FileFilter_19_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_19.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FileFilter_19_output.consensusXML )
set_tests_properties("TOPP_FileFilter_19_out1" PROPERTIES DEPENDS "TOPP_FileFilter_19")
add_test("TOPP_FileFilter_20" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_15_input.featureXML -out FileFilter_20.tmp.featureXML -id:accessions_whitelist YDL217C -id:remove_unassigned_ids)
add_test("TOPP_FileFilter_20_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_20.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FileFilter_20_output.featureXML )
set_tests_properties("TOPP_FileFilter_20_out1" PROPERTIES DEPENDS "TOPP_FileFilter_20")
add_test("TOPP_FileFilter_21" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_15_input.featureXML -out FileFilter_21.tmp.featureXML -id:remove_unassigned_ids -id:remove_unannotated_features -id:keep_best_score_id)
add_test("TOPP_FileFilter_21_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_21.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FileFilter_21_output.featureXML )
set_tests_properties("TOPP_FileFilter_21_out1" PROPERTIES DEPENDS "TOPP_FileFilter_21")
add_test("TOPP_FileFilter_22" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_22_input.consensusXML -out FileFilter_22.tmp.consensusXML -f_and_c:remove_meta distinct_charges gt "1,2")
add_test("TOPP_FileFilter_22_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_22.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FileFilter_22_output.consensusXML )
set_tests_properties("TOPP_FileFilter_22_out1" PROPERTIES DEPENDS "TOPP_FileFilter_22")
add_test("TOPP_FileFilter_23" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_22_input.consensusXML -out FileFilter_23.tmp.consensusXML -f_and_c:remove_meta distinct_charges_size gt 2)
add_test("TOPP_FileFilter_23_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_23.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FileFilter_22_output.consensusXML )
set_tests_properties("TOPP_FileFilter_23_out1" PROPERTIES DEPENDS "TOPP_FileFilter_23")
add_test("TOPP_FileFilter_24" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_22_input.consensusXML -out FileFilter_24.tmp.consensusXML -f_and_c:remove_meta DOESNOTEXIST lt "whatever")
add_test("TOPP_FileFilter_24_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_24.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FileFilter_24_output.consensusXML )
set_tests_properties("TOPP_FileFilter_24_out1" PROPERTIES DEPENDS "TOPP_FileFilter_24")

add_test("TOPP_FileFilter_25" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_25_input.mzML.gz -id:blacklist ${DATA_DIR_TOPP}/FileFilter_25_input.idXML -out FileFilter_25.tmp.mzML -id:mz 0.05 -id:rt 1 -id:blacklist_imperfect)
add_test("TOPP_FileFilter_25_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_25.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_25_output.mzML )
set_tests_properties("TOPP_FileFilter_25_out1" PROPERTIES DEPENDS "TOPP_FileFilter_25")

add_test("TOPP_FileFilter_26" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_25_input.mzML.gz -id:blacklist ${DATA_DIR_TOPP}/FileFilter_25_input.idXML -out FileFilter_26.tmp.mzML -id:blacklist_imperfect)
add_test("TOPP_FileFilter_26_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_26.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_25_output.mzML ) ## with missing id:mz and id:rt -- to test the default values
set_tests_properties("TOPP_FileFilter_26_out1" PROPERTIES DEPENDS "TOPP_FileFilter_26")

add_test("TOPP_FileFilter_27" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_25_input.mzML.gz -id:blacklist ${DATA_DIR_TOPP}/FileFilter_25_input.idXML -out FileFilter_27.tmp.mzML  -id:mz 0.05 -id:rt 1 )
set_tests_properties("TOPP_FileFilter_27" PROPERTIES WILL_FAIL 1) ## has 2 imperfect matches

add_test("TOPP_FileFilter_28" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_28_input.mzML.gz -peak_options:pc_mz_range 832:836 -out FileFilter_28.tmp.mzML -peak_options:level 2) ## only one MS2 remaining
add_test("TOPP_FileFilter_28_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_28.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_28_output.mzML )
set_tests_properties("TOPP_FileFilter_28_out1" PROPERTIES DEPENDS "TOPP_FileFilter_28")

add_test("TOPP_FileFilter_29" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_28_input.mzML.gz -peak_options:pc_mz_range 832:836 -out FileFilter_29.tmp.mzML -peak_options:level 1 2) ## only one MS2 remaining, and lots of MS1
add_test("TOPP_FileFilter_29_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_29.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_29_output.mzML )
set_tests_properties("TOPP_FileFilter_29_out1" PROPERTIES DEPENDS "TOPP_FileFilter_29")

add_test("TOPP_FileFilter_30" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_28_input.mzML.gz -peak_options:pc_mz_range 832:836 -out FileFilter_30.tmp.mzML -peak_options:level 1 2) ## no MS2 remaining since -mz filter empties it and it will be deleted, and lots of MS1
add_test("TOPP_FileFilter_30_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_30.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_30_output.mzML )
set_tests_properties("TOPP_FileFilter_30_out1" PROPERTIES DEPENDS "TOPP_FileFilter_30")

add_test("TOPP_FileFilter_31" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_31_34_input.mzML -spectra:remove_isolation_window_width :2 -out FileFilter_31.tmp.mzML )
add_test("TOPP_FileFilter_31_out" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_31.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_31_remove_collision_energy.mzML )
set_tests_properties("TOPP_FileFilter_31_out" PROPERTIES DEPENDS "TOPP_FileFilter_31")

add_test("TOPP_FileFilter_32" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_31_34_input.mzML -spectra:select_isolation_window_width :2 -out FileFilter_32.tmp.mzML )
add_test("TOPP_FileFilter_32_out" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_32.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_32_select_collision_energy.mzML )
set_tests_properties("TOPP_FileFilter_32_out" PROPERTIES DEPENDS "TOPP_FileFilter_32")

add_test("TOPP_FileFilter_33" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_31_34_input.mzML -spectra:remove_collision_energy :35 -out FileFilter_33.tmp.mzML )
add_test("TOPP_FileFilter_33_out" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_33.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_33_remove_isolation_window.mzML )
set_tests_properties("TOPP_FileFilter_33_out" PROPERTIES DEPENDS "TOPP_FileFilter_33")

add_test("TOPP_FileFilter_34" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_31_34_input.mzML -spectra:select_collision_energy :35 -out FileFilter_34.tmp.mzML )
add_test("TOPP_FileFilter_34_out" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_34.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_34_select_isolation_window.mzML )
set_tests_properties("TOPP_FileFilter_34_out" PROPERTIES DEPENDS "TOPP_FileFilter_34")

# Purpose: test the -peak_options:indexed_file to write mzML files with/without indices
add_test("TOPP_FileFilter_35" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out FileFilter_35.tmp.mzML -peak_options:indexed_file true -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_35_out" ${DIFF} -in1 FileFilter_35.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_35_output.mzML )
set_tests_properties("TOPP_FileFilter_35_out" PROPERTIES DEPENDS "TOPP_FileFilter_35")
add_test("TOPP_FileFilter_36" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out FileFilter_36.tmp.mzML -peak_options:indexed_file false -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_36_out" ${DIFF} -in1 FileFilter_36.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_36_output.mzML )
set_tests_properties("TOPP_FileFilter_36_out" PROPERTIES DEPENDS "TOPP_FileFilter_36")

# Purpose: test the -numpress options for mzML files
add_test("TOPP_FileFilter_37" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out FileFilter_37.tmp.mzML -test -peak_options:numpress:intensity slof -peak_options:numpress:masstime linear -in_type mzML -peak_options:numpress:lossy_mass_accuracy 1e-4 -out_type mzML)
add_test("TOPP_FileFilter_37_out" ${DIFF} -in1 FileFilter_37.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_37_output.mzML )
set_tests_properties("TOPP_FileFilter_37_out" PROPERTIES DEPENDS "TOPP_FileFilter_37")
add_test("TOPP_FileFilter_38" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out FileFilter_38.tmp.mzML -test -peak_options:numpress:intensity pic -peak_options:numpress:masstime linear -in_type mzML -peak_options:numpress:lossy_mass_accuracy 1e-4 -out_type mzML)
add_test("TOPP_FileFilter_38_out" ${DIFF} -in1 FileFilter_38.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_38_output.mzML )
set_tests_properties("TOPP_FileFilter_38_out" PROPERTIES DEPENDS "TOPP_FileFilter_38")
add_test("TOPP_FileFilter_39" ${TOPP_BIN_PATH}/FileFilter -in ${DATA_DIR_TOPP}/FileFilter_39_input_corrupt.mzML -out FileFilter_39.tmp.mzML -test)
set_tests_properties("TOPP_FileFilter_39" PROPERTIES WILL_FAIL 1) ## this mzML has a corrupt numpress encoded string and should fail

add_test("TOPP_FileFilter_40" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_40_input.mzML -out FileFilter_40.tmp.mzML -spectra:select_polarity positive -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_40_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_40.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_40_output.mzML)
set_tests_properties("TOPP_FileFilter_40_out1" PROPERTIES DEPENDS "TOPP_FileFilter_40")
add_test("TOPP_FileFilter_41" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_40_input.mzML -out FileFilter_41.tmp.mzML -spectra:select_polarity negative -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_41_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_41.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_41_output.mzML)
set_tests_properties("TOPP_FileFilter_41_out1" PROPERTIES DEPENDS "TOPP_FileFilter_41")
add_test("TOPP_FileFilter_42" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_40_input.mzML -out FileFilter_42.tmp.mzML -spectra:select_polarity "" -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_42_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_42.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_42_output.mzML)
set_tests_properties("TOPP_FileFilter_42_out1" PROPERTIES DEPENDS "TOPP_FileFilter_42")
# Test XML escaping of our XML handlers
add_test("TOPP_FileFilter_43" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_43_input.mzML -out FileFilter_43.tmp.mzML -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_43_out1" ${DIFF} -whitelist "id=" "href=" -in1 FileFilter_43.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_43_output.mzML)
# We check whether the XML we just produced can be read again
add_test("TOPP_FileFilter_43_read_again" ${TOPP_BIN_PATH}/FileFilter -test -in FileFilter_43.tmp.mzML -out FileFilter_43.dummy.tmp.mzML  -in_type mzML -out_type mzML)
set_tests_properties("TOPP_FileFilter_43_out1" PROPERTIES DEPENDS "TOPP_FileFilter_43")
set_tests_properties("TOPP_FileFilter_43_read_again" PROPERTIES DEPENDS "TOPP_FileFilter_43")
# test corrupt numpress file, see #1641
add_test("TOPP_FileFilter_44" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_44_input.mzML -out FileFilter_44.tmp.mzML -test -in_type mzML -out_type mzML)
add_test("TOPP_FileFilter_44_out" ${DIFF} -in1 FileFilter_44.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_44_output.mzML )
set_tests_properties("TOPP_FileFilter_44_out" PROPERTIES DEPENDS "TOPP_FileFilter_44")
add_test("TOPP_FileFilter_45" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_45_input.featureXML -id:sequences_whitelist "YSFS" "STLIPPPSK(Label:13C(6)15N(2))" -out FileFilter_45.tmp.featureXML)
add_test("TOPP_FileFilter_45_out" ${DIFF} -in1 FileFilter_45.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FileFilter_45_output.featureXML)
set_tests_properties("TOPP_FileFilter_45_out" PROPERTIES DEPENDS "TOPP_FileFilter_45")
add_test("TOPP_FileFilter_46" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_46_input.featureXML -id:sequences_whitelist "YSFS" "STLIPPPSK(Label:13C(6)15N(2))" -id:sequence_comparison_method "exact" -out FileFilter_46.tmp.featureXML)
add_test("TOPP_FileFilter_46_out" ${DIFF} -in1 FileFilter_46.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FileFilter_46_output.featureXML)
set_tests_properties("TOPP_FileFilter_46_out" PROPERTIES DEPENDS "TOPP_FileFilter_46")
# test filtering by spectral similarity
# remove similar spectrum
add_test("TOPP_FileFilter_47" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_47_input.mzML -spectra:blackorwhitelist:file ${DATA_DIR_TOPP}/FileFilter_47_input_select.mzML -spectra:blackorwhitelist:similarity_threshold 0.9 -out FileFilter_47_1.tmp.mzML)
add_test("TOPP_FileFilter_47_out" ${DIFF} -in1 FileFilter_47_1.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_47_output.mzML)
set_tests_properties("TOPP_FileFilter_47_out" PROPERTIES DEPENDS "TOPP_FileFilter_47")
# keep similar spectrum (we reuse the input from 47)
add_test("TOPP_FileFilter_48" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_47_input.mzML -spectra:blackorwhitelist:file ${DATA_DIR_TOPP}/FileFilter_47_input_select.mzML -spectra:blackorwhitelist:similarity_threshold 0.9 -spectra:blackorwhitelist:blacklist false -out FileFilter_48_1.tmp.mzML)
add_test("TOPP_FileFilter_48_out" ${DIFF} -in1 FileFilter_48_1.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_48_output.mzML)
set_tests_properties("TOPP_FileFilter_48_out" PROPERTIES DEPENDS "TOPP_FileFilter_48")

# Test with ion mobility
add_test("TOPP_FileFilter_49" ${TOPP_BIN_PATH}/FileFilter -test -in ${DATA_DIR_TOPP}/FileFilter_49_input.mzML -peak_options:numpress:intensity pic -peak_options:numpress:masstime linear -peak_options:numpress:float_da slof -peak_options:zlib_compression true -out FileFilter_49_1.tmp.mzML)
add_test("TOPP_FileFilter_49_out" ${DIFF} -in1 FileFilter_49_1.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileFilter_49_output.mzML)
set_tests_properties("TOPP_FileFilter_49_out" PROPERTIES DEPENDS "TOPP_FileFilter_49")


#------------------------------------------------------------------------------
# FileInfo tests
add_test("TOPP_FileInfo_1" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_1_input.dta -in_type dta -no_progress -out FileInfo_1.tmp.txt)
add_test("TOPP_FileInfo_1_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_1.tmp.txt -in2 ${DATA_DIR_TOPP}/FileInfo_1_output.txt )
set_tests_properties("TOPP_FileInfo_1_out1" PROPERTIES DEPENDS "TOPP_FileInfo_1")
add_test("TOPP_FileInfo_2" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_2_input.dta2d -no_progress -out FileInfo_2.tmp.txt)
add_test("TOPP_FileInfo_2_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_2.tmp.txt -in2 ${DATA_DIR_TOPP}/FileInfo_2_output.txt )
set_tests_properties("TOPP_FileInfo_2_out1" PROPERTIES DEPENDS "TOPP_FileInfo_2")
add_test("TOPP_FileInfo_3" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_3_input.featureXML -m -s -p -no_progress -out FileInfo_3.tmp.txt)
add_test("TOPP_FileInfo_3_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_3.tmp.txt -in2 ${DATA_DIR_TOPP}/FileInfo_3_output.txt )
set_tests_properties("TOPP_FileInfo_3_out1" PROPERTIES DEPENDS "TOPP_FileInfo_3")
add_test("TOPP_FileInfo_4" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_4_input.mzXML -m -no_progress -out FileInfo_4.tmp.txt)
add_test("TOPP_FileInfo_4_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_4.tmp.txt -in2 ${DATA_DIR_TOPP}/FileInfo_4_output.txt )
set_tests_properties("TOPP_FileInfo_4_out1" PROPERTIES DEPENDS "TOPP_FileInfo_4")
add_test("TOPP_FileInfo_5" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_5_input.mzDat -in_type mzData -m -s -no_progress -out FileInfo_5.tmp.txt)
add_test("TOPP_FileInfo_5_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_5.tmp.txt -in2 ${DATA_DIR_TOPP}/FileInfo_5_output.txt )
set_tests_properties("TOPP_FileInfo_5_out1" PROPERTIES DEPENDS "TOPP_FileInfo_5")
add_test("TOPP_FileInfo_6" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_6_input.mzData -d -s -no_progress -out FileInfo_6.tmp.txt)
add_test("TOPP_FileInfo_6_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_6.tmp.txt -in2 ${DATA_DIR_TOPP}/FileInfo_6_output.txt )
set_tests_properties("TOPP_FileInfo_6_out1" PROPERTIES DEPENDS "TOPP_FileInfo_6")
add_test("TOPP_FileInfo_7" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_7_input.consensusXML -s -m -p -no_progress -out FileInfo_7.tmp.txt)
add_test("TOPP_FileInfo_7_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_7.tmp.txt -in2 ${DATA_DIR_TOPP}/FileInfo_7_output.txt )
set_tests_properties("TOPP_FileInfo_7_out1" PROPERTIES DEPENDS "TOPP_FileInfo_7")
add_test("TOPP_FileInfo_9" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_9_input.mzML -m -p -s -no_progress -out FileInfo_9.tmp.txt)
add_test("TOPP_FileInfo_9_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_9.tmp.txt -in2 ${DATA_DIR_TOPP}/FileInfo_9_output.txt )
set_tests_properties("TOPP_FileInfo_9_out1" PROPERTIES DEPENDS "TOPP_FileInfo_9")
add_test("TOPP_FileInfo_10" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_10_input.idXML -no_progress -out FileInfo_10.tmp.txt)
add_test("TOPP_FileInfo_10_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_10.tmp.txt -in2 ${DATA_DIR_TOPP}/FileInfo_10_output.txt )
set_tests_properties("TOPP_FileInfo_10_out1" PROPERTIES DEPENDS "TOPP_FileInfo_10")
add_test("TOPP_FileInfo_11" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_11_input.mzML -i -no_progress)
set_tests_properties("TOPP_FileInfo_11" PROPERTIES WILL_FAIL 1) ## the mzML has no index
add_test("TOPP_FileInfo_12" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_12_input.mzML -i  -no_progress)
set_tests_properties("TOPP_FileInfo_12" PROPERTIES WILL_FAIL 0) ## this mzML has an index, should succeed
add_test("TOPP_FileInfo_13" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_13_input.consensusXML -no_progress) # empty file should not crash
add_test("TOPP_FileInfo_14" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_14_input.mzid -v -no_progress -out FileInfo_14.tmp.txt)
add_test("TOPP_FileInfo_14_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_14.tmp.txt -in2 ${DATA_DIR_TOPP}/FileInfo_14_output.txt )
set_tests_properties("TOPP_FileInfo_14_out1" PROPERTIES DEPENDS "TOPP_FileInfo_14")
add_test("TOPP_FileInfo_15" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_15_input.mzid -v -no_progress -out FileInfo_15.tmp.txt)
add_test("TOPP_FileInfo_15_out1" ${DIFF} -whitelist "File name" "line 327 column 183" -in1 FileInfo_15.tmp.txt -in2 ${DATA_DIR_TOPP}/FileInfo_15_output.txt )
set_tests_properties("TOPP_FileInfo_15_out1" PROPERTIES DEPENDS "TOPP_FileInfo_15")
add_test("TOPP_FileInfo_16" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_16_input.trafoXML -out FileInfo_16_output.tmp.txt)
add_test("TOPP_FileInfo_16_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_16_output.tmp.txt -in2 ${DATA_DIR_TOPP}/FileInfo_16_output.txt )
set_tests_properties("TOPP_FileInfo_16_out1" PROPERTIES DEPENDS "TOPP_FileInfo_16")
add_test("TOPP_FileInfo_17" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_17_input.fasta -out FileInfo_17_output.tmp.txt)
add_test("TOPP_FileInfo_17_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_17_output.tmp.txt -in2 ${DATA_DIR_TOPP}/FileInfo_17_output.txt )
set_tests_properties("TOPP_FileInfo_17_out1" PROPERTIES DEPENDS "TOPP_FileInfo_17")
add_test("TOPP_FileInfo_18" ${TOPP_BIN_PATH}/FileInfo -test -in ${DATA_DIR_TOPP}/FileInfo_18_input.fasta -out FileInfo_18_output.tmp.txt)
add_test("TOPP_FileInfo_18_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_18_output.tmp.txt -in2 ${DATA_DIR_TOPP}/FileInfo_18_output.txt )
set_tests_properties("TOPP_FileInfo_18_out1" PROPERTIES DEPENDS "TOPP_FileInfo_18")
add_test("TOPP_FileInfo_19" ${TOPP_BIN_PATH}/FileInfo -test -d -in ${DATA_DIR_TOPP}/FAIMS_test_data.mzML -out FileInfo_19_output.tmp.txt)
add_test("TOPP_FileInfo_19_out1" ${DIFF} -whitelist "File name" -in1 FileInfo_19_output.tmp.txt -in2 ${DATA_DIR_TOPP}/FileInfo_19_output.txt )
set_tests_properties("TOPP_FileInfo_19_out1" PROPERTIES DEPENDS "TOPP_FileInfo_19")



#------------------------------------------------------------------------------
# FileMerger tests
add_test("TOPP_FileMerger_1" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_1_input1.dta2d ${DATA_DIR_TOPP}/FileMerger_1_input2.dta ${DATA_DIR_TOPP}/FileMerger_1_input3.dta2d -out FileMerger_1.tmp.mzML -raw:rt_custom 1 2 3)
add_test("TOPP_FileMerger_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} "id=" "location=" "href=" -in1 FileMerger_1.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileMerger_1_output.mzML )
set_tests_properties("TOPP_FileMerger_1_out1" PROPERTIES DEPENDS "TOPP_FileMerger_1")
add_test("TOPP_FileMerger_2" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_2_input1.dta ${DATA_DIR_TOPP}/FileMerger_2_input2.dta -in_type dta -out FileMerger_2.tmp.mzML -raw:rt_custom 5 10)
add_test("TOPP_FileMerger_2_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} "id=" "location=" "href=" -in1 FileMerger_2.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileMerger_2_output.mzML )
set_tests_properties("TOPP_FileMerger_2_out1" PROPERTIES DEPENDS "TOPP_FileMerger_2")
add_test("TOPP_FileMerger_3" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_3_input1.dta ${DATA_DIR_TOPP}/FileMerger_3_input2.dta -out FileMerger_3.tmp.mzML -raw:rt_auto)
add_test("TOPP_FileMerger_3_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} "id=" "location=" "href=" -in1 FileMerger_3.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileMerger_3_output.mzML )
set_tests_properties("TOPP_FileMerger_3_out1" PROPERTIES DEPENDS "TOPP_FileMerger_3")
add_test("TOPP_FileMerger_4" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_4_input1.dta2d ${DATA_DIR_TOPP}/FileMerger_4_input2.dta2d -out FileMerger_4.tmp.mzML)
add_test("TOPP_FileMerger_4_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} "id=" "location=" "href=" -in1 FileMerger_4.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileMerger_4_output.mzML )
set_tests_properties("TOPP_FileMerger_4_out1" PROPERTIES DEPENDS "TOPP_FileMerger_4")
add_test("TOPP_FileMerger_5" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_5_input_rt1023.331.dta ${DATA_DIR_TOPP}/FileMerger_5_input_rt2044.334.dta ${DATA_DIR_TOPP}/FileMerger_5_input_rt889.32.dta -raw:rt_filename -raw:ms_level 2 -out FileMerger_5.tmp.mzML)
add_test("TOPP_FileMerger_5_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} "id=" "location=" "href=" -in1 FileMerger_5.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileMerger_5_output.mzML )
set_tests_properties("TOPP_FileMerger_5_out1" PROPERTIES DEPENDS "TOPP_FileMerger_5")
add_test("TOPP_FileMerger_6" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_6_input1.mzML ${DATA_DIR_TOPP}/FileMerger_6_input2.mzML -out FileMerger_6.tmp.mzML)
add_test("TOPP_FileMerger_6_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} "id=" "location=" "href=" -in1 FileMerger_6.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileMerger_6_output.mzML )
set_tests_properties("TOPP_FileMerger_6_out1" PROPERTIES DEPENDS "TOPP_FileMerger_6")
add_test("TOPP_FileMerger_7" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_7_input1.featureXML ${DATA_DIR_TOPP}/FileMerger_7_input2.featureXML -out FileMerger_7.tmp.featureXML )
add_test("TOPP_FileMerger_7_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} "id=" "location=" "href=" -in1 FileMerger_7.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FileMerger_7_output.featureXML )
set_tests_properties("TOPP_FileMerger_7_out1" PROPERTIES DEPENDS "TOPP_FileMerger_7")
add_test("TOPP_FileMerger_8" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_8_input1.consensusXML ${DATA_DIR_TOPP}/FileMerger_8_input2.consensusXML -out FileMerger_8.tmp.consensusXML)
add_test("TOPP_FileMerger_8_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} "id=" "location=" "href=" -in1 FileMerger_8.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FileMerger_8_output.consensusXML )
set_tests_properties("TOPP_FileMerger_8_out1" PROPERTIES DEPENDS "TOPP_FileMerger_8")
add_test("TOPP_FileMerger_9" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_9_input1.traML ${DATA_DIR_TOPP}/FileMerger_9_input2.traML -out FileMerger_9.tmp.traML)
add_test("TOPP_FileMerger_9_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} "id=" "location=" "href=" -in1 FileMerger_9.tmp.traML -in2 ${DATA_DIR_TOPP}/FileMerger_9_output.traML )
set_tests_properties("TOPP_FileMerger_9_out1" PROPERTIES DEPENDS "TOPP_FileMerger_9")
# "rt_concat" options:
add_test("TOPP_FileMerger_10" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_6_input2.mzML ${DATA_DIR_TOPP}/FileMerger_6_input2.mzML -out FileMerger_10_output.tmp.mzML -rt_concat:gap 10.0 -rt_concat:trafo_out FileMerger_10_trafo1.tmp.trafoXML FileMerger_10_trafo2.tmp.trafoXML)
add_test("TOPP_FileMerger_10_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} "id=" "location=" "href=" "offset" "indexListOffset" -in1 FileMerger_10_output.tmp.mzML -in2 ${DATA_DIR_TOPP}/FileMerger_10_output.mzML )
set_tests_properties("TOPP_FileMerger_10_out1" PROPERTIES DEPENDS "TOPP_FileMerger_10")
add_test("TOPP_FileMerger_10_out2" ${DIFF} -in1 FileMerger_10_trafo1.tmp.trafoXML -in2 ${DATA_DIR_TOPP}/FileMerger_10_trafo1.trafoXML)
set_tests_properties("TOPP_FileMerger_10_out2" PROPERTIES DEPENDS "TOPP_FileMerger_10")
add_test("TOPP_FileMerger_10_out3" ${DIFF} -in1 FileMerger_10_trafo2.tmp.trafoXML -in2 ${DATA_DIR_TOPP}/FileMerger_10_trafo2.trafoXML)
set_tests_properties("TOPP_FileMerger_10_out3" PROPERTIES DEPENDS "TOPP_FileMerger_10")
# append_cols
add_test("TOPP_FileMerger_11" ${TOPP_BIN_PATH}/FileMerger -test -in ${DATA_DIR_TOPP}/FileMerger_8_input1.consensusXML ${DATA_DIR_TOPP}/FileMerger_8_input2.consensusXML -append_method append_cols -out FileMerger_11.tmp.consensusXML)
add_test("TOPP_FileMerger_11_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} "id=" "location=" "href=" -in1 FileMerger_11.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FileMerger_11_output.consensusXML )
set_tests_properties("TOPP_FileMerger_11_out1" PROPERTIES DEPENDS "TOPP_FileMerger_11")

#------------------------------------------------------------------------------
# FLASHDeconv tests
add_test("TOPP_FLASHDeconv_1" ${TOPP_BIN_PATH}/FLASHDeconv -test -in ${DATA_DIR_TOPP}/FLASHDeconv_sample_input.mzML -out FLASHDeconv_1.tmp.tsv)
# Currently just testing if the tool runs
# We cannot currently test the correctness of the tsv output
# FuzzyDiff cannot give out tolerances/whitelist for tabular format outputs (.tsv)

#------------------------------------------------------------------------------
# GNPSExport tests
# test mostintense output
add_test("TOPP_GNPSExport_1" ${TOPP_BIN_PATH}/GNPSExport -test -ini ${DATA_DIR_TOPP}/GNPSExport_1_mostint.ini -in_cm ${DATA_DIR_TOPP}/GNPSExport_cons1.consensusXML -in_mzml ${DATA_DIR_TOPP}/GNPSExport_mzml1.mzML ${DATA_DIR_TOPP}/GNPSExport_mzml2.mzML -out GNPSExport_1_out.tmp.mgf -out_quantification GNPSExport_1_out_quant.tmp.txt)
add_test("TOPP_GNPSExport_1_out" ${DIFF} -in1 GNPSExport_1_out.tmp.mgf -in2 ${DATA_DIR_TOPP}/GNPSExport_1_out.mgf)
set_tests_properties("TOPP_GNPSExport_1_out" PROPERTIES DEPENDS "TOPP_GNPSExport_1")
add_test("TOPP_GNPSExport_1_out_QuantificationTable" ${DIFF} -in1 GNPSExport_1_out_quant.tmp.txt -in2 ${DATA_DIR_TOPP}/GNPSExport_1_out_quant.txt)
set_tests_properties("TOPP_GNPSExport_1_out_QuantificationTable" PROPERTIES DEPENDS "TOPP_GNPSExport_1")
# test mergedspectra test
add_test("TOPP_GNPSExport_2" ${TOPP_BIN_PATH}/GNPSExport -test -ini ${DATA_DIR_TOPP}/GNPSExport_2_merged.ini -in_cm ${DATA_DIR_TOPP}/GNPSExport_cons1.consensusXML -in_mzml ${DATA_DIR_TOPP}/GNPSExport_mzml1.mzML ${DATA_DIR_TOPP}/GNPSExport_mzml2.mzML -out GNPSExport_2_out.tmp.mgf -out_quantification GNPSExport_2_out_quant.tmp.txt)
add_test("TOPP_GNPSExport_2_out" ${DIFF} -in1 GNPSExport_2_out.tmp.mgf -in2 ${DATA_DIR_TOPP}/GNPSExport_2_out.mgf)
set_tests_properties("TOPP_GNPSExport_2_out" PROPERTIES DEPENDS "TOPP_GNPSExport_2")
# test ms2_bin_size
add_test("TOPP_GNPSExport_3" ${TOPP_BIN_PATH}/GNPSExport -test -ini ${DATA_DIR_TOPP}/GNPSExport_3_binsize.ini -in_cm ${DATA_DIR_TOPP}/GNPSExport_cons1.consensusXML -in_mzml ${DATA_DIR_TOPP}/GNPSExport_mzml1.mzML ${DATA_DIR_TOPP}/GNPSExport_mzml2.mzML -out GNPSExport_3_out.tmp.mgf -out_quantification GNPSExport_3_out_quant.tmp.txt)
add_test("TOPP_GNPSExport_3_out" ${DIFF} -in1 GNPSExport_3_out.tmp.mgf -in2 ${DATA_DIR_TOPP}/GNPSExport_3_out.mgf)
set_tests_properties("TOPP_GNPSExport_3_out" PROPERTIES DEPENDS "TOPP_GNPSExport_3")
# test output for IonIdentityMolecularNetworking
add_test("TOPP_GNPSExport_4" ${TOPP_BIN_PATH}/GNPSExport -test -in_cm ${DATA_DIR_TOPP}/GNPSExport_cons2.consensusXML -in_mzml ${DATA_DIR_TOPP}/GNPSExport_mzml1.mzML ${DATA_DIR_TOPP}/GNPSExport_mzml2.mzML -out GNPSExport_4_out.tmp.mgf -out_quantification GNPSExport_4_out_quant.tmp.txt -out_pairs GNPSExport_4_out_pairs.tmp.csv -out_meta_values GNPSExport_4_out_meta_values.tmp.tsv)
add_test("TOPP_GNPSExport_4_out_quant" ${DIFF} -in1 GNPSExport_4_out_quant.tmp.txt -in2 ${DATA_DIR_TOPP}/GNPSExport_4_out_quant.txt)
set_tests_properties("TOPP_GNPSExport_4_out_quant" PROPERTIES DEPENDS "TOPP_GNPSExport_4")
add_test("TOPP_GNPSExport_4_out_pairs" ${DIFF} -in1 GNPSExport_4_out_pairs.tmp.csv -in2 ${DATA_DIR_TOPP}/GNPSExport_4_out_pairs.csv)
set_tests_properties("TOPP_GNPSExport_4_out_pairs" PROPERTIES DEPENDS "TOPP_GNPSExport_4")
add_test("TOPP_GNPSExport_4_out_meta_values" ${DIFF} -in1 GNPSExport_4_out_meta_values.tmp.tsv -in2 ${DATA_DIR_TOPP}/GNPSExport_4_out_meta_values.tsv)
set_tests_properties("TOPP_GNPSExport_4_out_meta_values" PROPERTIES DEPENDS "TOPP_GNPSExport_4")

# HighResPrecursorMassCorrector tests
add_test("TOPP_HighResPrecursorMassCorrector_1"
${TOPP_BIN_PATH}/HighResPrecursorMassCorrector -test -in ${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_2860_1103_3.mzML -feature:in ${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_2860_1103_3.featureXML -out HighResPrecursorMassCorrector_2860_1103_3_out.tmp.mzML )
add_test("TOPP_HighResPrecursorMassCorrector_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 HighResPrecursorMassCorrector_2860_1103_3_out.tmp.mzML  -in2 ${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_2860_1103_3_out.mzML )
set_tests_properties("TOPP_HighResPrecursorMassCorrector_1_out1" PROPERTIES DEPENDS "TOPP_HighResPrecursorMassCorrector_1")
# precursor is on the 4th isotopic trace (with 0 beeing the mono-isotopic one)
add_test("TOPP_HighResPrecursorMassCorrector_2"
	${TOPP_BIN_PATH}/HighResPrecursorMassCorrector -test -in
	${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_1035_1178_4.mzML -feature:in
	${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_1035_1178_4.featureXML
	-feature:max_trace 4
	-feature:mz_tolerance 10
	-out HighResPrecursorMassCorrector_1035_1178_4_out.tmp.mzML )
add_test("TOPP_HighResPrecursorMassCorrector_2_out1" ${DIFF} -whitelist ${INDEX_WHITELIST}  -in1 HighResPrecursorMassCorrector_1035_1178_4_out.tmp.mzML -in2 ${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_1035_1178_4_out.mzML )
set_tests_properties("TOPP_HighResPrecursorMassCorrector_2_out1" PROPERTIES DEPENDS "TOPP_HighResPrecursorMassCorrector_2")
add_test("TOPP_HighResPrecursorMassCorrector_3" ${TOPP_BIN_PATH}/HighResPrecursorMassCorrector -test -in ${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_2538_1091_2.mzML -feature:in ${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_2538_1091_2.featureXML -out HighResPrecursorMassCorrector_2538_1091_2_out.tmp.mzML )
add_test("TOPP_HighResPrecursorMassCorrector_3_out1" ${DIFF} -whitelist ${INDEX_WHITELIST}  -in1 HighResPrecursorMassCorrector_2538_1091_2_out.tmp.mzML -in2 ${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_2538_1091_2_out.mzML )
set_tests_properties("TOPP_HighResPrecursorMassCorrector_3_out1" PROPERTIES DEPENDS "TOPP_HighResPrecursorMassCorrector_3")
add_test("TOPP_HighResPrecursorMassCorrector_4" ${TOPP_BIN_PATH}/HighResPrecursorMassCorrector -test -in ${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_2810_1091_3.mzML -feature:in ${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_2810_1091_3.featureXML -out HighResPrecursorMassCorrector_2810_1091_3_out.tmp.mzML )
add_test("TOPP_HighResPrecursorMassCorrector_4_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 HighResPrecursorMassCorrector_2810_1091_3_out.tmp.mzML -in2 ${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_2810_1091_3_out.mzML )
set_tests_properties("TOPP_HighResPrecursorMassCorrector_4_out1" PROPERTIES DEPENDS "TOPP_HighResPrecursorMassCorrector_4")
add_test("TOPP_HighResPrecursorMassCorrector_5" ${TOPP_BIN_PATH}/HighResPrecursorMassCorrector -test -in ${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_3070_1191_3.mzML -feature:in ${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_3070_1191_3.featureXML -out HighResPrecursorMassCorrector_3070_1191_3_out.tmp.mzML )
add_test("TOPP_HighResPrecursorMassCorrector_5_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 HighResPrecursorMassCorrector_3070_1191_3_out.tmp.mzML -in2 ${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_3070_1191_3_out.mzML )
set_tests_properties("TOPP_HighResPrecursorMassCorrector_5_out1" PROPERTIES DEPENDS "TOPP_HighResPrecursorMassCorrector_5")
add_test("TOPP_HighResPrecursorMassCorrector_6" ${TOPP_BIN_PATH}/HighResPrecursorMassCorrector -test -in ${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_6.mzML -highest_intensity_peak:mz_tolerance 0.2 -out HighResPrecursorMassCorrector_6_out.tmp.mzML )
add_test("TOPP_HighResPrecursorMassCorrector_6_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 HighResPrecursorMassCorrector_6_out.tmp.mzML -in2 ${DATA_DIR_TOPP}/HighResPrecursorMassCorrector_6_out.mzML )
set_tests_properties("TOPP_HighResPrecursorMassCorrector_6_out1" PROPERTIES DEPENDS "TOPP_HighResPrecursorMassCorrector_6")

#------------------------------------------------------------------------------
# IDRTCalibration tests
add_test("TOPP_IDRTCalibration_1" ${TOPP_BIN_PATH}/IDRTCalibration -test -in ${DATA_DIR_TOPP}/IDRTCalibration_1_input.idXML -out IDRTCalibration_1_output.tmp.idXML -calibrant_1_input 10 -calibrant_2_input 90)
add_test("TOPP_IDRTCalibration_1_out1" ${DIFF} -in1 IDRTCalibration_1_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDRTCalibration_1_output.idXML )
set_tests_properties("TOPP_IDRTCalibration_1_out1" PROPERTIES DEPENDS "TOPP_IDRTCalibration_1")
# test empty IDs
add_test("TOPP_IDRTCalibration_2" ${TOPP_BIN_PATH}/IDRTCalibration -test -in ${DATA_DIR_TOPP}/degenerated_empty.idXML -out IDRTCalibration_2_output.tmp.idXML -calibrant_1_input 10 -calibrant_2_input 90)
add_test("TOPP_IDRTCalibration_2_out1" ${DIFF} -in1 IDRTCalibration_2_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDRTCalibration_2_output.idXML )
set_tests_properties("TOPP_IDRTCalibration_2_out1" PROPERTIES DEPENDS "TOPP_IDRTCalibration_2")


#------------------------------------------------------------------------------
# IDMassAccuracy tests
add_test("TOPP_IDMassAccuracy_1" ${TOPP_BIN_PATH}/IDMassAccuracy -test -in ${DATA_DIR_TOPP}/THIRDPARTY/spectra.mzML -id_in ${DATA_DIR_TOPP}/THIRDPARTY/MSGFPlusAdapter_1_out.idXML -number_of_bins 10 -out_fragment IDMassAccuracy_1_out_fragment.tsv -out_fragment_fit IDMassAccuracy_1_out_fragment_fit.tsv -out_precursor IDMassAccuracy_1_out_precursor.tsv -out_precursor_fit IDMassAccuracy_1_out_precursor_fit.tsv)
# Currently just testing if the tool runs

#------------------------------------------------------------------------------
# IsobaricAnalyzer tests
add_test("TOPP_IsobaricAnalyzer_1" ${TOPP_BIN_PATH}/IsobaricAnalyzer -test -in ${DATA_DIR_TOPP}/IsobaricAnalyzer_input_1.mzML -ini ${DATA_DIR_TOPP}/IsobaricAnalyzer.ini -out IsobaricAnalyzer_output_1.tmp.consensusXML)
add_test("TOPP_IsobaricAnalyzer_1_out1" ${DIFF} -whitelist "id=" "<map" "?xml-stylesheet" -in1 IsobaricAnalyzer_output_1.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/IsobaricAnalyzer_output_1.consensusXML )
set_tests_properties("TOPP_IsobaricAnalyzer_1_out1" PROPERTIES DEPENDS "TOPP_IsobaricAnalyzer_1")

add_test("TOPP_IsobaricAnalyzer_TMTTenPlexMethod_1" ${TOPP_BIN_PATH}/IsobaricAnalyzer -test -in ${DATA_DIR_TOPP}/TMTTenPlexMethod_test.mzML -ini ${DATA_DIR_TOPP}/TMTTenPlexMethod_test.ini -out TMTTenPlexMethod_output.tmp.consensusXML)
add_test("TOPP_IsobaricAnalyzer_TMTTenPlexMethod_1_out1" ${DIFF} -whitelist "id=" "<map" "?xml-stylesheet" -in1 TMTTenPlexMethod_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/TMTTenPlexMethod_test.consensusXML )
set_tests_properties("TOPP_IsobaricAnalyzer_TMTTenPlexMethod_1_out1" PROPERTIES DEPENDS "TOPP_IsobaricAnalyzer_TMTTenPlexMethod_1")

add_test("TOPP_IsobaricAnalyzer_TMTElevenPlexMethod_1" ${TOPP_BIN_PATH}/IsobaricAnalyzer -test -in ${DATA_DIR_TOPP}/TMTTenPlexMethod_test.mzML -ini ${DATA_DIR_TOPP}/TMTElevenPlexMethod_test.ini -out TMTElevenPlexMethod_output.tmp.consensusXML)
add_test("TOPP_IsobaricAnalyzer_TMTElevenPlexMethod_1_out1" ${DIFF} -whitelist "id=" "<map" "?xml-stylesheet" -in1 TMTElevenPlexMethod_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/TMTElevenPlexMethod_test.consensusXML )
set_tests_properties("TOPP_IsobaricAnalyzer_TMTElevenPlexMethod_1_out1" PROPERTIES DEPENDS "TOPP_IsobaricAnalyzer_TMTElevenPlexMethod_1")

add_test("TOPP_IsobaricAnalyzer_MS3TMT10Plex_1" ${TOPP_BIN_PATH}/IsobaricAnalyzer -test -in ${DATA_DIR_TOPP}/MS3_nonHierarchical.mzML -extraction:select_activation "Collision-induced dissociation" -type tmt10plex -out MS3TMT10Plex_output.tmp.consensusXML)
add_test("TOPP_IsobaricAnalyzer_MS3TMT10Plex_1_out1" ${DIFF} -whitelist "id=" "<map" "?xml-stylesheet" -in1 MS3TMT10Plex_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/MS3TMT10Plex_test.consensusXML )
set_tests_properties("TOPP_IsobaricAnalyzer_MS3TMT10Plex_1_out1" PROPERTIES DEPENDS "TOPP_IsobaricAnalyzer_MS3TMT10Plex_1")

#------------------------------------------------------------------------------
# IDConflictResolver tests
add_test("TOPP_IDConflictResolver_1" ${TOPP_BIN_PATH}/IDConflictResolver -test -in ${DATA_DIR_TOPP}/IDConflictResolver_1_input.featureXML -out IDConflictResolver_1_output.tmp.featureXML)
add_test("TOPP_IDConflictResolver_1_out1" ${DIFF} -in1 IDConflictResolver_1_output.tmp.featureXML -in2 ${DATA_DIR_TOPP}/IDConflictResolver_1_output.featureXML )
set_tests_properties("TOPP_IDConflictResolver_1_out1" PROPERTIES DEPENDS "TOPP_IDConflictResolver_1")
add_test("TOPP_IDConflictResolver_2" ${TOPP_BIN_PATH}/IDConflictResolver -test -in ${DATA_DIR_TOPP}/IDConflictResolver_2_input.consensusXML -out IDConflictResolver_2_output.tmp.consensusXML)
add_test("TOPP_IDConflictResolver_2_out1" ${DIFF} -in1 IDConflictResolver_2_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/IDConflictResolver_2_output.consensusXML )
set_tests_properties("TOPP_IDConflictResolver_2_out1" PROPERTIES DEPENDS "TOPP_IDConflictResolver_2")
# Tests including empty features and features with IDs but without hits.
add_test("TOPP_IDConflictResolver_3" ${TOPP_BIN_PATH}/IDConflictResolver -test -in ${DATA_DIR_TOPP}/IDConflictResolver_3_input.consensusXML -out IDConflictResolver_3_output.tmp.consensusXML)
add_test("TOPP_IDConflictResolver_3_out1" ${DIFF} -in1 IDConflictResolver_3_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/IDConflictResolver_3_output.consensusXML )
set_tests_properties("TOPP_IDConflictResolver_3_out1" PROPERTIES DEPENDS "TOPP_IDConflictResolver_3")
# Tests resolve_between_features option
add_test("TOPP_IDConflictResolver_4" ${TOPP_BIN_PATH}/IDConflictResolver -test -in ${DATA_DIR_TOPP}/IDConflictResolver_4_input.featureXML -resolve_between_features highest_intensity -out IDConflictResolver_4.tmp.featureXML)
add_test("TOPP_IDConflictResolver_4_out" ${DIFF} -in1 IDConflictResolver_4.tmp.featureXML -in2 ${DATA_DIR_TOPP}/IDConflictResolver_4_output.featureXML)
set_tests_properties("TOPP_IDConflictResolver_4_out" PROPERTIES DEPENDS "TOPP_IDConflictResolver_4")

#------------------------------------------------------------------------------
# IDFileConverter tests
# Mascot XML to idXML:
add_test("TOPP_IDFileConverter_1" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_1_input1.mascotXML -mz_file ${DATA_DIR_TOPP}/IDFileConverter_1_input2.mzML -out IDFileConverter_1_output.tmp.idXML -out_type idXML)
add_test("TOPP_IDFileConverter_1_out1" ${DIFF} -in1 IDFileConverter_1_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFileConverter_1_output.idXML )
set_tests_properties("TOPP_IDFileConverter_1_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_1")
# pepXML to idXML:
add_test("TOPP_IDFileConverter_2" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/PepXMLFile_test.pepxml -mz_file ${DATA_DIR_TOPP}/PepXMLFile_test.mzML -mz_name PepXMLFile_test -out IDFileConverter_2_output.tmp.idXML -out_type idXML)
add_test("TOPP_IDFileConverter_2_out1" ${DIFF} -in1 IDFileConverter_2_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFileConverter_2_output.idXML )
set_tests_properties("TOPP_IDFileConverter_2_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_2")
# protXML to idXML:
add_test("TOPP_IDFileConverter_3" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_3_input.protXML -out IDFileConverter_3_output.tmp.idXML -out_type idXML)
add_test("TOPP_IDFileConverter_3_out1" ${DIFF} -in1 IDFileConverter_3_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFileConverter_3_output.idXML )
set_tests_properties("TOPP_IDFileConverter_3_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_3")
# Mascot XML to idXML (with user-defined regular expression):
add_test("TOPP_IDFileConverter_5" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_1_input1.mascotXML -mz_file ${DATA_DIR_TOPP}/IDFileConverter_1_input2.mzML -out IDFileConverter_5_output.tmp.idXML -ini ${DATA_DIR_TOPP}/IDFileConverter_5_parameters.ini)
add_test("TOPP_IDFileConverter_5_out1" ${DIFF} -in1 IDFileConverter_5_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFileConverter_5_output.idXML )
set_tests_properties("TOPP_IDFileConverter_5_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_5")
# pepXML (exported from Mascot) to idXML:
add_test("TOPP_IDFileConverter_6" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_6_input1.pepXML -mz_file ${DATA_DIR_TOPP}/IDFileConverter_1_input2.mzML -mz_name F025589.dat.mzML -out IDFileConverter_6_output.tmp.idXML -out_type idXML)
add_test("TOPP_IDFileConverter_6_out1" ${DIFF} -in1 IDFileConverter_6_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFileConverter_6_output.idXML )
set_tests_properties("TOPP_IDFileConverter_6_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_6")
# X! Tandem XML to idXML:
add_test("TOPP_IDFileConverter_7" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_7_input1.xml -out IDFileConverter_7_output.tmp.idXML -out_type idXML)
add_test("TOPP_IDFileConverter_7_out1" ${DIFF} -in1 IDFileConverter_7_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFileConverter_7_output.idXML -whitelist "IdentificationRun date" )
set_tests_properties("TOPP_IDFileConverter_7_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_7")
# mzIdentML to idXML (input file is "MzIdentMLFile_whole.mzid" from the class test):
add_test("TOPP_IDFileConverter_8" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_8_input.mzid -out IDFileConverter_8_output.tmp.idXML -out_type idXML)
add_test("TOPP_IDFileConverter_8_out1" ${DIFF} -in1 IDFileConverter_8_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFileConverter_8_output.idXML -whitelist "IdentificationRun date" )
set_tests_properties("TOPP_IDFileConverter_8_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_8")
#  idXML to mzIdentML (input file is "IdXMLFile_whole.idXML" from the class test):
add_test("TOPP_IDFileConverter_9" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_9_input.idXML -out IDFileConverter_9_output.tmp.mzid -out_type mzid)
add_test("TOPP_IDFileConverter_9_out1" ${DIFF} -in1 IDFileConverter_9_output.tmp.mzid -in2 ${DATA_DIR_TOPP}/IDFileConverter_9_output.mzid -whitelist "creationDate" "id" "spectraData_ref" "searchDatabase_ref")
set_tests_properties("TOPP_IDFileConverter_9_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_9")
# pepXML with very small delta mass as modification (most probably an artefact and removed during parsing)
add_test("TOPP_IDFileConverter_10" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_10_input.pepXML -out IDFileConverter_10_output.tmp.idXML -out_type idXML)
add_test("TOPP_IDFileConverter_10_out1" ${DIFF} -in1 IDFileConverter_10_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFileConverter_10_output.idXML )
set_tests_properties("TOPP_IDFileConverter_10_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_10")

# tsv to idXML
add_test("TOPP_IDFileConverter_11" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_11_input.tsv -out IDFileConverter_11_output.tmp.idXML -out_type idXML)
add_test("TOPP_IDFileConverter_11_out1" ${DIFF} -in1 IDFileConverter_11_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFileConverter_11_output.idXML)
set_tests_properties("TOPP_IDFileConverter_11_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_11")
# Percolator .psms to idXML (score: q-value)
add_test("TOPP_IDFileConverter_12" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_12_input.psms -out IDFileConverter_12_output.tmp.idXML -out_type idXML -score_type qvalue)
add_test("TOPP_IDFileConverter_12_out1" ${DIFF} -in1 IDFileConverter_12_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFileConverter_12_output.idXML -whitelist "IdentificationRun date")
set_tests_properties("TOPP_IDFileConverter_12_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_12")
# Percolator .psms to idXML (score: PEP)
add_test("TOPP_IDFileConverter_13" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_12_input.psms -out IDFileConverter_13_output.tmp.idXML -out_type idXML -score_type PEP)
add_test("TOPP_IDFileConverter_13_out1" ${DIFF} -in1 IDFileConverter_13_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFileConverter_13_output.idXML -whitelist "IdentificationRun date")
set_tests_properties("TOPP_IDFileConverter_13_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_13")
# Percolator .psms to idXML (score: raw score)
add_test("TOPP_IDFileConverter_14" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_12_input.psms -out IDFileConverter_14_output.tmp.idXML -out_type idXML -score_type score)
add_test("TOPP_IDFileConverter_14_out1" ${DIFF} -in1 IDFileConverter_14_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFileConverter_14_output.idXML -whitelist "IdentificationRun date")
set_tests_properties("TOPP_IDFileConverter_14_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_14")
# mzIdentML to idXML with RT look-up in mzML:
add_test("TOPP_IDFileConverter_15" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/THIRDPARTY/MSGFPlusAdapter_1_out.mzid -out IDFileConverter_15_output.tmp.idXML -out_type idXML -mz_file ${DATA_DIR_TOPP}/THIRDPARTY/spectra.mzML)
add_test("TOPP_IDFileConverter_15_out1" ${DIFF} -in1 IDFileConverter_15_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFileConverter_15_output.idXML -whitelist "IdentificationRun date" "db=" "spectra_data" "IsotopeError")
set_tests_properties("TOPP_IDFileConverter_15_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_15")
# pepXML (from iProphet) to idXML:
add_test("TOPP_IDFileConverter_16" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_16_input.pepXML -out IDFileConverter_16_output.tmp.idXML -out_type idXML)
add_test("TOPP_IDFileConverter_16_out1" ${DIFF} -in1 IDFileConverter_16_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFileConverter_16_output.idXML)
set_tests_properties("TOPP_IDFileConverter_16_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_16")
# idXML to pepXML (like in LuciphorAdapter)
# Tests  E-Value, hyperscore
add_test("TOPP_IDFileConverter_17" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_17_input.idXML -out IDFileConverter_17_output.tmp.pepXML -out_type pepXML)
add_test("TOPP_IDFileConverter_17_out1" ${DIFF} -in1 IDFileConverter_17_output.tmp.pepXML -in2 ${DATA_DIR_TOPP}/IDFileConverter_17_output.pepXML)
set_tests_properties("TOPP_IDFileConverter_17_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_17")
# Tests nextscore
add_test("TOPP_IDFileConverter_18" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_18_input.idXML -out IDFileConverter_18_output.tmp.pepXML -out_type pepXML)
add_test("TOPP_IDFileConverter_18_out1" ${DIFF} -in1 IDFileConverter_18_output.tmp.pepXML -in2 ${DATA_DIR_TOPP}/IDFileConverter_18_output.pepXML)
set_tests_properties("TOPP_IDFileConverter_18_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_18")
# Tests Percolator_score
add_test("TOPP_IDFileConverter_19" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_19_input.idXML -out IDFileConverter_19_output.tmp.pepXML -out_type pepXML)
add_test("TOPP_IDFileConverter_19_out1" ${DIFF} -in1 IDFileConverter_19_output.tmp.pepXML -in2 ${DATA_DIR_TOPP}/IDFileConverter_19_output.pepXML)
set_tests_properties("TOPP_IDFileConverter_19_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_19")
# Tests hyperscore
add_test("TOPP_IDFileConverter_20" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_20_input.idXML -out IDFileConverter_20_output.tmp.pepXML -out_type pepXML)
add_test("TOPP_IDFileConverter_20_out1" ${DIFF} -in1 IDFileConverter_20_output.tmp.pepXML -in2 ${DATA_DIR_TOPP}/IDFileConverter_20_output.pepXML)
set_tests_properties("TOPP_IDFileConverter_20_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_20")
# Testing spectrum_reference addition
add_test("TOPP_IDFileConverter_21" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_21_input.idXML -out_type idXML -out IDFileConverter_21_output.tmp.idXML -mz_file ${DATA_DIR_TOPP}/IDMapper_4_input.mzML)
add_test("TOPP_IDFileConverter_21_out1" ${DIFF} -in1 IDFileConverter_21_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFileConverter_21_output.idXML -whitelist "creationDate" "id" "spectraData_ref" "searchDatabase_ref")
set_tests_properties("TOPP_IDFileConverter_21_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_21")
# Testing add_ionmatch_annotation
add_test("TOPP_IDFileConverter_22" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/THIRDPARTY/MSGFPlusAdapter_1_out.mzid -out IDFileConverter_22_output.tmp.idXML -out_type idXML -mz_file ${DATA_DIR_TOPP}/THIRDPARTY/spectra.mzML -add_ionmatch_annotation 0.01)
add_test("TOPP_IDFileConverter_22_out1" ${DIFF} -in1 IDFileConverter_22_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFileConverter_22_output.idXML -whitelist "IdentificationRun date" "IsotopeError" "id" "spectra_data")
set_tests_properties("TOPP_IDFileConverter_22_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_22")
# Crux mzid
add_test("TOPP_IDFileConverter_23" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_23_input.mzid -out IDFileConverter_23_output.tmp.idXML -out_type idXML)
add_test("TOPP_IDFileConverter_23_out1" ${DIFF} -in1 IDFileConverter_23_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFileConverter_23_output.idXML -whitelist "IdentificationRun date" )
set_tests_properties("TOPP_IDFileConverter_23_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_23")
# Comet pep.xml
add_test("TOPP_IDFileConverter_24" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_24_input.pep.xml -out IDFileConverter_24_output.tmp.idXML -out_type idXML)
add_test("TOPP_IDFileConverter_24_out1" ${DIFF} -in1 IDFileConverter_24_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFileConverter_24_output.idXML -whitelist "IdentificationRun date" )
set_tests_properties("TOPP_IDFileConverter_24_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_24")
# Comet pep.xml
add_test("TOPP_IDFileConverter_25" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_25_input.idXML -out IDFileConverter_25_output.tmp.pepXML -out_type pepXML)
add_test("TOPP_IDFileConverter_25_out1" ${DIFF} -in1 IDFileConverter_25_output.tmp.pepXML -in2 ${DATA_DIR_TOPP}/IDFileConverter_25_output.pep.xml)
set_tests_properties("TOPP_IDFileConverter_25_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_25")
# MS-GF+ mzid with non-Unimod modifications
add_test("TOPP_IDFileConverter_26" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_26_input.mzid -out IDFileConverter_26_output.tmp.idXML -out_type idXML)
add_test("TOPP_IDFileConverter_26_out1" ${DIFF} -in1 IDFileConverter_26_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFileConverter_26_output.idXML -whitelist "IdentificationRun date" )
set_tests_properties("TOPP_IDFileConverter_26_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_26")
# Test idXML to FASTA
# no concatenation
add_test("TOPP_IDFileConverter_27" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_27_input.idXML -out IDFileConverter_27_output.tmp.fasta -out_type fasta)
add_test("TOPP_IDFileConverter_27_out1" ${DIFF} -in1 IDFileConverter_27_output.tmp.fasta -in2 ${DATA_DIR_TOPP}/IDFileConverter_27_output.fasta)
set_tests_properties("TOPP_IDFileConverter_27_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_27")
add_test("TOPP_IDFileConverter_28" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_27_input.idXML -out IDFileConverter_28_output.tmp.fasta -out_type fasta -number_of_hits -1)
add_test("TOPP_IDFileConverter_28_out1" ${DIFF} -in1 IDFileConverter_28_output.tmp.fasta -in2 ${DATA_DIR_TOPP}/IDFileConverter_28_output.fasta)
set_tests_properties("TOPP_IDFileConverter_28_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_28")
# with concatenation
add_test("TOPP_IDFileConverter_29" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_27_input.idXML -out IDFileConverter_29_output.tmp.fasta -out_type fasta -concatenate_peptides)
add_test("TOPP_IDFileConverter_29_out1" ${DIFF} -in1 IDFileConverter_29_output.tmp.fasta -in2 ${DATA_DIR_TOPP}/IDFileConverter_29_output.fasta)
set_tests_properties("TOPP_IDFileConverter_29_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_29")
add_test("TOPP_IDFileConverter_30" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_27_input.idXML -out IDFileConverter_30_output.tmp.fasta -out_type fasta -concatenate_peptides -number_of_hits 2)
add_test("TOPP_IDFileConverter_30_out1" ${DIFF} -in1 IDFileConverter_30_output.tmp.fasta -in2 ${DATA_DIR_TOPP}/IDFileConverter_30_output.fasta)
set_tests_properties("TOPP_IDFileConverter_30_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_30")
# Test MzIdentML files with unofficial (non CV-term) scores
add_test("TOPP_IDFileConverter_31" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_31_input.mzid -out IDFileConverter_31_output.tmp.idXML -out_type idXML)
add_test("TOPP_IDFileConverter_31_out1" ${DIFF} -in1 IDFileConverter_31_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFileConverter_31_output.idXML)
set_tests_properties("TOPP_IDFileConverter_31_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_31")
# Test FASTA to MzML
# default parameters
add_test("TOPP_IDFileConverter_32" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_32_input.FASTA -out IDFileConverter_32_output.tmp.mzml -out_type mzML)
add_test("TOPP_IDFileConverter_32_File_Conversion" ${TOPP_BIN_PATH}/FileConverter -test -in IDFileConverter_32_output.tmp.mzml -out IDFileConverter_32_output_mgf.tmp.mgf -out_type mgf)
set_tests_properties("TOPP_IDFileConverter_32_File_Conversion" PROPERTIES DEPENDS "TOPP_IDFileConverter_32")
add_test("TOPP_IDFileConverter_32_out1" ${DIFF} -in1 IDFileConverter_32_output_mgf.tmp.mgf -in2 ${DATA_DIR_TOPP}/IDFileConverter_32_output.mgf -whitelist "TITLE")
set_tests_properties("TOPP_IDFileConverter_32_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_32_File_Conversion")
add_test("TOPP_IDFileConverter_33" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_32_input.FASTA -out IDFileConverter_33_output.tmp.mzml -out_type mzML -fasta_to_mzml:missed_cleavages 1 -fasta_to_mzml:precursor_charge 3)
add_test("TOPP_IDFileConverter_33_File_Conversion" ${TOPP_BIN_PATH}/FileConverter -test -in IDFileConverter_33_output.tmp.mzml -out IDFileConverter_33_output_mgf.tmp.mgf -out_type mgf)
set_tests_properties("TOPP_IDFileConverter_33_File_Conversion" PROPERTIES DEPENDS "TOPP_IDFileConverter_33")
add_test("TOPP_IDFileConverter_33_out1" ${DIFF} -in1 IDFileConverter_33_output_mgf.tmp.mgf -in2 ${DATA_DIR_TOPP}/IDFileConverter_33_output.mgf  -whitelist "TITLE")
set_tests_properties("TOPP_IDFileConverter_33_out1" PROPERTIES DEPENDS "TOPP_IDFileConverter_33_File_Conversion")
# OpenMS SQLite format (input file is "IdXMLFile_whole.idXML" from the class test, with one score type renamed for consistency):
add_test("TOPP_IDFileConverter_34_1" ${TOPP_BIN_PATH}/IDFileConverter -test -in ${DATA_DIR_TOPP}/IDFileConverter_34_input.idXML -out IDFileConverter_34_output1.oms)
add_test("TOPP_IDFileConverter_34_2" ${TOPP_BIN_PATH}/IDFileConverter -test -in IDFileConverter_34_output1.oms -out IDFileConverter_34_output2.tmp.idxml -out_type idXML)
add_test("TOPP_IDFileConverter_34_out" ${DIFF} -in1 IDFileConverter_34_output2.tmp.idxml -in2 ${DATA_DIR_TOPP}/IDFileConverter_34_output.idXML)
set_tests_properties("TOPP_IDFileConverter_34_2" PROPERTIES DEPENDS "TOPP_IDFileConverter_34_1")
set_tests_properties("TOPP_IDFileConverter_34_out" PROPERTIES DEPENDS "TOPP_IDFileConverter_34_2")

#------------------------------------------------------------------------------
# IDFilter tests
add_test("TOPP_IDFilter_1" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_1_input.idXML -out IDFilter_1_output.tmp.idXML -whitelist:proteins ${DATA_DIR_TOPP}/IDFilter_1_input.fas)
add_test("TOPP_IDFilter_1_out1" ${DIFF} -in1 IDFilter_1_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFilter_1_output.idXML )
set_tests_properties("TOPP_IDFilter_1_out1" PROPERTIES DEPENDS "TOPP_IDFilter_1")
add_test("TOPP_IDFilter_3" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_3_input.idXML -out IDFilter_3_output.tmp.idXML -blacklist:peptides ${DATA_DIR_TOPP}/IDFilter_3_2_input.idXML)
add_test("TOPP_IDFilter_3_out1" ${DIFF} -in1 IDFilter_3_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFilter_3_output.idXML )
set_tests_properties("TOPP_IDFilter_3_out1" PROPERTIES DEPENDS "TOPP_IDFilter_3")
add_test("TOPP_IDFilter_4" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_4_input.idXML -out IDFilter_4_output.tmp.idXML -rt:p_value 0.08)
add_test("TOPP_IDFilter_4_out1" ${DIFF} -in1 IDFilter_4_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFilter_4_output.idXML )
set_tests_properties("TOPP_IDFilter_4_out1" PROPERTIES DEPENDS "TOPP_IDFilter_4")

add_test("TOPP_IDFilter_5" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_5_input.idXML -out IDFilter_5_output.tmp.idXML -score:pep 32 -score:prot 25)
add_test("TOPP_IDFilter_5_out1" ${DIFF} -in1 IDFilter_5_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFilter_5_output.idXML )
set_tests_properties("TOPP_IDFilter_5_out1" PROPERTIES DEPENDS "TOPP_IDFilter_5")
add_test("TOPP_IDFilter_5a" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_5_input.idXML -out IDFilter_5a_output.tmp.idXML -score:pep 32 )
add_test("TOPP_IDFilter_5a_out1" ${DIFF} -in1 IDFilter_5a_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFilter_5_output.idXML )
set_tests_properties("TOPP_IDFilter_5a_out1" PROPERTIES DEPENDS "TOPP_IDFilter_5a")

add_test("TOPP_IDFilter_5b" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_5_input.idXML -out IDFilter_5b_output.tmp.idXML -score:prot 25)
add_test("TOPP_IDFilter_5b_out1" ${DIFF} -in1 IDFilter_5b_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFilter_5b_output.idXML )
set_tests_properties("TOPP_IDFilter_5b_out1" PROPERTIES DEPENDS "TOPP_IDFilter_5b")

add_test("TOPP_IDFilter_5c" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_5_input.idXML -out IDFilter_5c_output.tmp.idXML -score:prot 25 -delete_unreferenced_peptide_hits)
add_test("TOPP_IDFilter_5c_out1" ${DIFF} -in1 IDFilter_5c_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFilter_5c_output.idXML )
set_tests_properties("TOPP_IDFilter_5c_out1" PROPERTIES DEPENDS "TOPP_IDFilter_5c")

add_test("TOPP_IDFilter_6" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_6_input.idXML -out IDFilter_6_output.tmp.idXML -best:n_peptide_hits 2 -best:n_protein_hits 10)
add_test("TOPP_IDFilter_6_out1" ${DIFF} -in1 IDFilter_6_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFilter_6_output.idXML )
set_tests_properties("TOPP_IDFilter_6_out1" PROPERTIES DEPENDS "TOPP_IDFilter_6")

add_test("TOPP_IDFilter_7" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_7_input.idXML -out IDFilter_7_output.tmp.idXML -remove_duplicate_psm)
add_test("TOPP_IDFilter_7_out1" ${DIFF} -in1 IDFilter_7_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFilter_7_output.idXML )
set_tests_properties("TOPP_IDFilter_7_out1" PROPERTIES DEPENDS "TOPP_IDFilter_7")
# filter by precursor RT and m/z:
add_test("TOPP_IDFilter_8" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_8_input.idXML -out IDFilter_8_output.tmp.idXML -precursor:rt 200:350 -precursor:mz 999:1000)
add_test("TOPP_IDFilter_8_out1" ${DIFF} -in1 IDFilter_8_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFilter_8_output.idXML )
set_tests_properties("TOPP_IDFilter_8_out1" PROPERTIES DEPENDS "TOPP_IDFilter_8")
# check for issue 887 regression (multiple runs)
add_test("TOPP_IDFilter_9" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_9_input.idXML -out IDFilter_9_output.tmp.idXML -score:pep 0.05)
add_test("TOPP_IDFilter_9_out1" ${DIFF} -in1 IDFilter_9_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFilter_9_output.idXML )
set_tests_properties("TOPP_IDFilter_9_out1" PROPERTIES DEPENDS "TOPP_IDFilter_9")
# update protein groups ("IDFilter_10_input.idXML" is same as "IDFileConverter_3_output.idXML", i.e. converted protXML):
add_test("TOPP_IDFilter_10" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_10_input.idXML -out IDFilter_10_output.tmp.idXML -score:prot 0.3 -delete_unreferenced_peptide_hits)
add_test("TOPP_IDFilter_10_out1" ${DIFF} -in1 IDFilter_10_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFilter_10_output.idXML )
set_tests_properties("TOPP_IDFilter_10_out1" PROPERTIES DEPENDS "TOPP_IDFilter_10")
# remove decoys and update protein groups:
add_test("TOPP_IDFilter_11" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_11_input.idXML -out IDFilter_11_output.tmp.idXML -remove_decoys -delete_unreferenced_peptide_hits)
add_test("TOPP_IDFilter_11_out1" ${DIFF} -in1 IDFilter_11_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFilter_11_output.idXML )
set_tests_properties("TOPP_IDFilter_11_out1" PROPERTIES DEPENDS "TOPP_IDFilter_11")

# filters by digestion products
# default digestion parameters
add_test("TOPP_IDFilter_12" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_12_input.idXML -out IDFilter_12_output.tmp.idXML -in_silico_digestion:fasta ${DATA_DIR_TOPP}/IDFilter_12_input.fasta)
add_test("TOPP_IDFilter_12_out1" ${DIFF} -in1 IDFilter_12_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFilter_12_output.idXML )
set_tests_properties("TOPP_IDFilter_12_out1" PROPERTIES DEPENDS "TOPP_IDFilter_12")
# filters by digestion products considering missed cleavages
add_test("TOPP_IDFilter_13" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_13_input.idXML -out IDFilter_13_output.tmp.idXML -in_silico_digestion:fasta ${DATA_DIR_TOPP}/IDFilter_13_input.fasta -in_silico_digestion:missed_cleavages 1)
add_test("TOPP_IDFilter_13_out1" ${DIFF} -in1 IDFilter_13_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFilter_13_output.idXML)
set_tests_properties("TOPP_IDFilter_13_out1" PROPERTIES DEPENDS "TOPP_IDFilter_13")
# filters by digestion products considering missed cleavages with specific peptide
add_test("TOPP_IDFilter_14" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_14_input.idXML -out IDFilter_14_output.tmp.idXML -in_silico_digestion:fasta ${DATA_DIR_TOPP}/IDFilter_14_input.fasta -in_silico_digestion:enzyme Trypsin/P -in_silico_digestion:missed_cleavages 1)
add_test("TOPP_IDFilter_14_out1" ${DIFF} -in1 IDFilter_14_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFilter_14_output.idXML)
set_tests_properties("TOPP_IDFilter_14_out1" PROPERTIES DEPENDS "TOPP_IDFilter_14")
# filters by digestion products considering missed cleavages with semi specificity filtering
add_test("TOPP_IDFilter_15" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_15_input.idXML -out IDFilter_15_output.tmp.idXML -in_silico_digestion:fasta ${DATA_DIR_TOPP}/IDFilter_15_input.fasta -in_silico_digestion:enzyme Trypsin/P -in_silico_digestion:missed_cleavages 1 -in_silico_digestion:specificity semi)
add_test("TOPP_IDFilter_15_out1" ${DIFF} -in1 IDFilter_15_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFilter_15_output.idXML)
set_tests_properties("TOPP_IDFilter_15_out1" PROPERTIES DEPENDS "TOPP_IDFilter_15")
add_test("TOPP_IDFilter_16" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_16_input.idXML -out IDFilter_16_output.tmp.idXML -in_silico_digestion:fasta ${DATA_DIR_TOPP}/IDFilter_16_input.fasta -in_silico_digestion:enzyme Trypsin/P -in_silico_digestion:missed_cleavages 1 -in_silico_digestion:methionine_cleavage)
add_test("TOPP_IDFilter_16_out1" ${DIFF} -in1 IDFilter_16_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFilter_16_output.idXML)
set_tests_properties("TOPP_IDFilter_16_out1" PROPERTIES DEPENDS "TOPP_IDFilter_16")

#filters by peptide missed cleavages
add_test("TOPP_IDFilter_17" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_missed_cleavages_input.idXML -out IDFilter_17_output.tmp.idXML -missed_cleavages:enzyme Lys-N -missed_cleavages:number_of_missed_cleavages :2)
add_test("TOPP_IDFilter_17_out1" ${DIFF} -in1 IDFilter_17_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFilter_17_output.idXML )
set_tests_properties("TOPP_IDFilter_17_out1" PROPERTIES DEPENDS "TOPP_IDFilter_17")

add_test("TOPP_IDFilter_18" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_missed_cleavages_input.idXML -out IDFilter_18_output.tmp.idXML -missed_cleavages:enzyme Lys-N -missed_cleavages:number_of_missed_cleavages 2:)
add_test("TOPP_IDFilter_18_out1" ${DIFF} -in1 IDFilter_18_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFilter_18_output.idXML )
set_tests_properties("TOPP_IDFilter_18_out1" PROPERTIES DEPENDS "TOPP_IDFilter_18")

add_test("TOPP_IDFilter_19" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_missed_cleavages_input.idXML -out IDFilter_19_output.tmp.idXML -missed_cleavages:enzyme Lys-N -missed_cleavages:number_of_missed_cleavages 1:3)
add_test("TOPP_IDFilter_19_out1" ${DIFF} -in1 IDFilter_19_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFilter_19_output.idXML )
set_tests_properties("TOPP_IDFilter_19_out1" PROPERTIES DEPENDS "TOPP_IDFilter_19")

#this should filter all peptide hits invalid range
add_test("TOPP_IDFilter_20" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_missed_cleavages_input.idXML -out IDFilter_20_output.tmp.idXML -missed_cleavages:number_of_missed_cleavages 1:0)
add_test("TOPP_IDFilter_20_out1" ${DIFF} -in1 IDFilter_20_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFilter_20_output.idXML )
set_tests_properties("TOPP_IDFilter_20_out1" PROPERTIES DEPENDS "TOPP_IDFilter_20")

#filter peptide hits by meta value name
add_test("TOPP_IDFilter_21" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_16_input.idXML -out IDFilter_21_output.tmp.idXML -remove_peptide_hits_by_metavalue "calcMZ" "gt" "750.0")
add_test("TOPP_IDFilter_21_out1" ${DIFF} -in1 IDFilter_21_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFilter_21_output.idXML)
set_tests_properties("TOPP_IDFilter_21_out1" PROPERTIES DEPENDS "TOPP_IDFilter_21")

add_test("TOPP_IDFilter_22" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDFilter_16_input.idXML -out IDFilter_22_output.tmp.idXML -remove_peptide_hits_by_metavalue "end" "ne" "23")
add_test("TOPP_IDFilter_22_out1" ${DIFF} -in1 IDFilter_22_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFilter_22_output.idXML)
set_tests_properties("TOPP_IDFilter_22_out1" PROPERTIES DEPENDS "TOPP_IDFilter_22")

#consensusXML
add_test("TOPP_IDFilter_23" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/IDMapper_2_output.consensusXML -out IDFilter_23_output.tmp.consensusXML -whitelist:protein_accessions "Q9HP81")
add_test("TOPP_IDFilter_23_out1" ${DIFF} -in1 IDFilter_23_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/IDFilter_23_output.consensusXML)
set_tests_properties("TOPP_IDFilter_23_out1" PROPERTIES DEPENDS "TOPP_IDFilter_23")

#consenusXML and protein groups
add_test("TOPP_IDFilter_24" ${TOPP_BIN_PATH}/IDFilter -test -in ${DATA_DIR_TOPP}/Epifany_3_out.consensusXML -out IDFilter_24_output.tmp.consensusXML -score:protgroup 0.99)
add_test("TOPP_IDFilter_24_out1" ${DIFF} -in1 IDFilter_24_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/IDFilter_24_output.consensusXML)
set_tests_properties("TOPP_IDFilter_24_out1" PROPERTIES DEPENDS "TOPP_IDFilter_24")

#-----------------------------------------------------------------------------
# MapAlignerPoseClustering tests
# featureXML input:
add_test("TOPP_MapAlignerPoseClustering_1" ${TOPP_BIN_PATH}/MapAlignerPoseClustering -test -ini ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_parameters.ini -in ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_input1.featureXML ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_input2.featureXML ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_input3.featureXML -out MapAlignerPoseClustering_1_output1.tmp.featureXML MapAlignerPoseClustering_1_output2.tmp.featureXML MapAlignerPoseClustering_1_output3.tmp.featureXML -trafo_out MapAlignerPoseClustering_1_trafo1.tmp.trafoXML MapAlignerPoseClustering_1_trafo2.tmp.trafoXML MapAlignerPoseClustering_1_trafo3.tmp.trafoXML)
add_test("TOPP_MapAlignerPoseClustering_1_out1" ${DIFF} -in1 MapAlignerPoseClustering_1_output1.tmp.featureXML -in2 ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_output1.featureXML )
set_tests_properties("TOPP_MapAlignerPoseClustering_1_out1" PROPERTIES DEPENDS "TOPP_MapAlignerPoseClustering_1")
add_test("TOPP_MapAlignerPoseClustering_1_out2" ${DIFF} -in1 MapAlignerPoseClustering_1_output2.tmp.featureXML -in2 ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_output2.featureXML )
set_tests_properties("TOPP_MapAlignerPoseClustering_1_out2" PROPERTIES DEPENDS "TOPP_MapAlignerPoseClustering_1")
add_test("TOPP_MapAlignerPoseClustering_1_out3" ${DIFF} -in1 MapAlignerPoseClustering_1_output3.tmp.featureXML -in2 ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_output3.featureXML )
set_tests_properties("TOPP_MapAlignerPoseClustering_1_out3" PROPERTIES DEPENDS "TOPP_MapAlignerPoseClustering_1")
add_test("TOPP_MapAlignerPoseClustering_1_out4" ${DIFF} -in1 MapAlignerPoseClustering_1_trafo1.tmp.trafoXML -in2 ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_trafo1.trafoXML )
set_tests_properties("TOPP_MapAlignerPoseClustering_1_out4" PROPERTIES DEPENDS "TOPP_MapAlignerPoseClustering_1")
add_test("TOPP_MapAlignerPoseClustering_1_out5" ${DIFF} -in1 MapAlignerPoseClustering_1_trafo2.tmp.trafoXML -in2 ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_trafo2.trafoXML )
set_tests_properties("TOPP_MapAlignerPoseClustering_1_out5" PROPERTIES DEPENDS "TOPP_MapAlignerPoseClustering_1")
add_test("TOPP_MapAlignerPoseClustering_1_out6" ${DIFF} -in1 MapAlignerPoseClustering_1_trafo3.tmp.trafoXML -in2 ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_trafo3.trafoXML )
set_tests_properties("TOPP_MapAlignerPoseClustering_1_out6" PROPERTIES DEPENDS "TOPP_MapAlignerPoseClustering_1")
# mzML input:
add_test("TOPP_MapAlignerPoseClustering_2" ${TOPP_BIN_PATH}/MapAlignerPoseClustering -test -ini ${DATA_DIR_TOPP}/MapAlignerPoseClustering_2_parameters.ini -in ${DATA_DIR_TOPP}/MapAlignerPoseClustering_2_input1.mzML ${DATA_DIR_TOPP}/MapAlignerPoseClustering_2_input2.mzML ${DATA_DIR_TOPP}/MapAlignerPoseClustering_2_input3.mzML -out MapAlignerPoseClustering_2_output1.tmp.mzML MapAlignerPoseClustering_2_output2.tmp.mzML MapAlignerPoseClustering_2_output3.tmp.mzML)
add_test("TOPP_MapAlignerPoseClustering_2_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 MapAlignerPoseClustering_2_output1.tmp.mzML -in2 ${DATA_DIR_TOPP}/MapAlignerPoseClustering_2_output1.mzML )
set_tests_properties("TOPP_MapAlignerPoseClustering_2_out1" PROPERTIES DEPENDS "TOPP_MapAlignerPoseClustering_2")
add_test("TOPP_MapAlignerPoseClustering_2_out2" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 MapAlignerPoseClustering_2_output2.tmp.mzML -in2 ${DATA_DIR_TOPP}/MapAlignerPoseClustering_2_output2.mzML )
set_tests_properties("TOPP_MapAlignerPoseClustering_2_out2" PROPERTIES DEPENDS "TOPP_MapAlignerPoseClustering_2")
add_test("TOPP_MapAlignerPoseClustering_2_out3" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 MapAlignerPoseClustering_2_output3.tmp.mzML -in2 ${DATA_DIR_TOPP}/MapAlignerPoseClustering_2_output3.mzML )
set_tests_properties("TOPP_MapAlignerPoseClustering_2_out3" PROPERTIES DEPENDS "TOPP_MapAlignerPoseClustering_2")
# "reference:file" option:
add_test("TOPP_MapAlignerPoseClustering_3" ${TOPP_BIN_PATH}/MapAlignerPoseClustering -test -ini ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_parameters.ini -in ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_input2.featureXML ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_input3.featureXML -out MapAlignerPoseClustering_3_output1.tmp.featureXML MapAlignerPoseClustering_3_output2.tmp.featureXML -reference:file ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_input1.featureXML)
add_test("TOPP_MapAlignerPoseClustering_3_out1" ${DIFF} -in1 MapAlignerPoseClustering_3_output1.tmp.featureXML -in2 ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_output2.featureXML )
set_tests_properties("TOPP_MapAlignerPoseClustering_3_out1" PROPERTIES DEPENDS "TOPP_MapAlignerPoseClustering_3")
add_test("TOPP_MapAlignerPoseClustering_3_out2" ${DIFF} -in1 MapAlignerPoseClustering_3_output2.tmp.featureXML -in2 ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_output3.featureXML )
set_tests_properties("TOPP_MapAlignerPoseClustering_3_out2" PROPERTIES DEPENDS "TOPP_MapAlignerPoseClustering_3")
# "reference:index" option:
add_test("TOPP_MapAlignerPoseClustering_4" ${TOPP_BIN_PATH}/MapAlignerPoseClustering -test -ini ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_parameters.ini -in ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_input1.featureXML ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_input2.featureXML -trafo_out MapAlignerPoseClustering_4_trafo1.tmp.trafoXML MapAlignerPoseClustering_4_trafo2.tmp.trafoXML -reference:index 2)
add_test("TOPP_MapAlignerPoseClustering_4_out1" ${DIFF} -in1 MapAlignerPoseClustering_4_trafo1.tmp.trafoXML -in2 ${DATA_DIR_TOPP}/MapAlignerPoseClustering_4_trafo1.trafoXML )
set_tests_properties("TOPP_MapAlignerPoseClustering_4_out1" PROPERTIES DEPENDS "TOPP_MapAlignerPoseClustering_4")
add_test("TOPP_MapAlignerPoseClustering_4_out2" ${DIFF} -in1 MapAlignerPoseClustering_4_trafo2.tmp.trafoXML -in2 ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_trafo1.trafoXML )
set_tests_properties("TOPP_MapAlignerPoseClustering_4_out2" PROPERTIES DEPENDS "TOPP_MapAlignerPoseClustering_4")

#------------------------------------------------------------------------------
# MapAlignerIdentification tests:
add_test("TOPP_MapAlignerIdentification_1" ${TOPP_BIN_PATH}/MapAlignerIdentification -test -ini ${DATA_DIR_TOPP}/MapAlignerIdentification_parameters.ini -in ${DATA_DIR_TOPP}/MapAlignerIdentification_1_input1.featureXML ${DATA_DIR_TOPP}/MapAlignerIdentification_1_input2.featureXML -out MapAlignerIdentification_1_output1.tmp.featureXML MapAlignerIdentification_1_output2.tmp.featureXML)
add_test("TOPP_MapAlignerIdentification_1_out1" ${DIFF} -in1 MapAlignerIdentification_1_output1.tmp.featureXML -in2 ${DATA_DIR_TOPP}/MapAlignerIdentification_1_output1.featureXML )
set_tests_properties("TOPP_MapAlignerIdentification_1_out1" PROPERTIES DEPENDS "TOPP_MapAlignerIdentification_1")
add_test("TOPP_MapAlignerIdentification_1_out2" ${DIFF} -in1 MapAlignerIdentification_1_output2.tmp.featureXML -in2 ${DATA_DIR_TOPP}/MapAlignerIdentification_1_output2.featureXML )
set_tests_properties("TOPP_MapAlignerIdentification_1_out2" PROPERTIES DEPENDS "TOPP_MapAlignerIdentification_1")
# "reference:file" option:
add_test("TOPP_MapAlignerIdentification_2" ${TOPP_BIN_PATH}/MapAlignerIdentification -test -ini ${DATA_DIR_TOPP}/MapAlignerIdentification_parameters.ini -in ${DATA_DIR_TOPP}/MapAlignerIdentification_1_input1.featureXML -out MapAlignerIdentification_2_output1.tmp.featureXML -reference:file ${DATA_DIR_TOPP}/MapAlignerIdentification_1_input2.featureXML)
add_test("TOPP_MapAlignerIdentification_2_out1" ${DIFF} -in1 MapAlignerIdentification_2_output1.tmp.featureXML -in2 ${DATA_DIR_TOPP}/MapAlignerIdentification_2_output1.featureXML )
set_tests_properties("TOPP_MapAlignerIdentification_2_out1" PROPERTIES DEPENDS "TOPP_MapAlignerIdentification_2")
# "reference:index" option:
add_test("TOPP_MapAlignerIdentification_3" ${TOPP_BIN_PATH}/MapAlignerIdentification -test -ini ${DATA_DIR_TOPP}/MapAlignerIdentification_parameters.ini -in ${DATA_DIR_TOPP}/MapAlignerIdentification_1_input2.featureXML ${DATA_DIR_TOPP}/MapAlignerIdentification_1_input1.featureXML -out MapAlignerIdentification_3_output1.tmp.featureXML MapAlignerIdentification_3_output2.tmp.featureXML -reference:index 1)
add_test("TOPP_MapAlignerIdentification_3_out1" ${DIFF} -in1 MapAlignerIdentification_3_output1.tmp.featureXML -in2 ${DATA_DIR_TOPP}/MapAlignerIdentification_3_output1.featureXML )
set_tests_properties("TOPP_MapAlignerIdentification_3_out1" PROPERTIES DEPENDS "TOPP_MapAlignerIdentification_3")
add_test("TOPP_MapAlignerIdentification_3_out2" ${DIFF} -in1 MapAlignerIdentification_3_output2.tmp.featureXML -in2 ${DATA_DIR_TOPP}/MapAlignerIdentification_2_output1.featureXML )
set_tests_properties("TOPP_MapAlignerIdentification_3_out2" PROPERTIES DEPENDS "TOPP_MapAlignerIdentification_3")
add_test("TOPP_MapAlignerIdentification_4" ${TOPP_BIN_PATH}/MapAlignerIdentification -test -ini ${DATA_DIR_TOPP}/MapAlignerIdentification_parameters.ini -in  ${DATA_DIR_TOPP}/MapAlignerIdentification_1_input1.featureXML ${DATA_DIR_TOPP}/MapAlignerIdentification_1_input2.featureXML -out MapAlignerIdentification_4_output1.tmp.featureXML MapAlignerIdentification_4_output2.tmp.featureXML -reference:index 2)
add_test("TOPP_MapAlignerIdentification_4_out1" ${DIFF} -in1 MapAlignerIdentification_4_output1.tmp.featureXML -in2 ${DATA_DIR_TOPP}/MapAlignerIdentification_2_output1.featureXML )
set_tests_properties("TOPP_MapAlignerIdentification_4_out1" PROPERTIES DEPENDS "TOPP_MapAlignerIdentification_4")
add_test("TOPP_MapAlignerIdentification_4_out2" ${DIFF} -in1 MapAlignerIdentification_4_output2.tmp.featureXML -in2 ${DATA_DIR_TOPP}/MapAlignerIdentification_3_output1.featureXML )
set_tests_properties("TOPP_MapAlignerIdentification_4_out2" PROPERTIES DEPENDS "TOPP_MapAlignerIdentification_4")
# consensusXML input:
add_test("TOPP_MapAlignerIdentification_5" ${TOPP_BIN_PATH}/MapAlignerIdentification -test -ini ${DATA_DIR_TOPP}/MapAlignerIdentification_parameters.ini -in ${DATA_DIR_TOPP}/MapAlignerIdentification_5_input1.consensusXML ${DATA_DIR_TOPP}/MapAlignerIdentification_5_input2.consensusXML -out MapAlignerIdentification_5_output1.tmp.consensusXML MapAlignerIdentification_5_output2.tmp.consensusXML)
add_test("TOPP_MapAlignerIdentification_5_out1" ${DIFF} -in1 MapAlignerIdentification_5_output1.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/MapAlignerIdentification_5_output1.consensusXML )
set_tests_properties("TOPP_MapAlignerIdentification_5_out1" PROPERTIES DEPENDS "TOPP_MapAlignerIdentification_5")
add_test("TOPP_MapAlignerIdentification_5_out2" ${DIFF} -in1 MapAlignerIdentification_5_output2.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/MapAlignerIdentification_5_output2.consensusXML )
set_tests_properties("TOPP_MapAlignerIdentification_5_out2" PROPERTIES DEPENDS "TOPP_MapAlignerIdentification_5")
# trafoXML output:
add_test("TOPP_MapAlignerIdentification_6" ${TOPP_BIN_PATH}/MapAlignerIdentification -test -ini ${DATA_DIR_TOPP}/MapAlignerIdentification_parameters.ini -in ${DATA_DIR_TOPP}/MapAlignerIdentification_1_input1.featureXML -trafo_out MapAlignerIdentification_6_output1.tmp.trafoXML -reference:file ${DATA_DIR_TOPP}/MapAlignerIdentification_1_input2.featureXML)
add_test("TOPP_MapAlignerIdentification_6_out1" ${DIFF} -in1 MapAlignerIdentification_6_output1.tmp.trafoXML -in2 ${DATA_DIR_TOPP}/MapAlignerIdentification_6_output1.trafoXML )
set_tests_properties("TOPP_MapAlignerIdentification_6_out1" PROPERTIES DEPENDS "TOPP_MapAlignerIdentification_6")
# "force" option (not enough RT data):
add_test("TOPP_MapAlignerIdentification_7" ${TOPP_BIN_PATH}/MapAlignerIdentification -test -in ${DATA_DIR_TOPP}/MapAlignerIdentification_7_input1.idXML -out MapAlignerIdentification_7_output1.tmp.idXML -trafo_out MapAlignerIdentification_7_output2.tmp.trafoXML -reference:file ${DATA_DIR_TOPP}/MapAlignerIdentification_7_input2.idXML -force)
add_test("TOPP_MapAlignerIdentification_7_out1" ${DIFF} -in1 MapAlignerIdentification_7_output1.tmp.idXML -in2 ${DATA_DIR_TOPP}/MapAlignerIdentification_7_input1.idXML)
add_test("TOPP_MapAlignerIdentification_7_out2" ${DIFF} -in1 MapAlignerIdentification_7_output2.tmp.trafoXML -in2 ${DATA_DIR_TOPP}/MapAlignerIdentification_7_output2.trafoXML -whitelist "TrafoXML version=")
set_tests_properties("TOPP_MapAlignerIdentification_7_out1" PROPERTIES DEPENDS "TOPP_MapAlignerIdentification_7")
set_tests_properties("TOPP_MapAlignerIdentification_7_out2" PROPERTIES DEPENDS "TOPP_MapAlignerIdentification_7")
# oms (database file) input - files based on BSA example data, generated like this:
# IDFileConverter -in share/OpenMS/examples/BSA/BSA2_OMSSA.idXML -out src/tests/topp/MapAlignerIdentification_8_input2.oms
# IDFileConverter -in share/OpenMS/examples/BSA/BSA3_OMSSA.idXML -out src/tests/topp/MapAlignerIdentification_8_input1.oms
add_test("TOPP_MapAlignerIdentification_8" ${TOPP_BIN_PATH}/MapAlignerIdentification -test -in ${DATA_DIR_TOPP}/MapAlignerIdentification_8_input1.oms -trafo_out MapAlignerIdentification_8_output1.tmp.trafoXML -out MapAlignerIdentification_8_output2.tmp.oms -reference:file ${DATA_DIR_TOPP}/MapAlignerIdentification_8_input2.oms -store_original_rt)
# don't compare binary SQLite files
add_test("TOPP_MapAlignerIdentification_8_out1" ${DIFF} -in1 MapAlignerIdentification_8_output1.tmp.trafoXML -in2 ${DATA_DIR_TOPP}/MapAlignerIdentification_8_output1.trafoXML)
set_tests_properties("TOPP_MapAlignerIdentification_8_out1" PROPERTIES DEPENDS "TOPP_MapAlignerIdentification_8")

#------------------------------------------------------------------------------
# MapAlignerTreeGuided tests:
add_test("TOPP_MapAlignerTreeGuided_1" ${TOPP_BIN_PATH}/MapAlignerTreeGuided -test -ini ${DATA_DIR_TOPP}/MapAlignerTreeGuided_parameters.ini -in ${DATA_DIR_TOPP}/MapAlignerTreeGuided_1_input1.featureXML ${DATA_DIR_TOPP}/MapAlignerTreeGuided_1_input2.featureXML ${DATA_DIR_TOPP}/MapAlignerTreeGuided_1_input3.featureXML -out MapAlignerTreeGuided_1_output1.tmp.featureXML MapAlignerTreeGuided_1_output2.tmp.featureXML MapAlignerTreeGuided_1_output3.tmp.featureXML)
add_test("TOPP_MapAlignerTreeGuided_1_out1" ${DIFF} -in1 MapAlignerTreeGuided_1_output1.tmp.featureXML -in2 ${DATA_DIR_TOPP}/MapAlignerTreeGuided_1_output1.featureXML )
set_tests_properties("TOPP_MapAlignerTreeGuided_1_out1" PROPERTIES DEPENDS "TOPP_MapAlignerTreeGuided_1")
add_test("TOPP_MapAlignerTreeGuided_1_out2" ${DIFF} -in1 MapAlignerTreeGuided_1_output2.tmp.featureXML -in2 ${DATA_DIR_TOPP}/MapAlignerTreeGuided_1_output2.featureXML )
set_tests_properties("TOPP_MapAlignerTreeGuided_1_out2" PROPERTIES DEPENDS "TOPP_MapAlignerTreeGuided_1")
add_test("TOPP_MapAlignerTreeGuided_1_out3" ${DIFF} -in1 MapAlignerTreeGuided_1_output3.tmp.featureXML -in2 ${DATA_DIR_TOPP}/MapAlignerTreeGuided_1_output3.featureXML )
set_tests_properties("TOPP_MapAlignerTreeGuided_1_out3" PROPERTIES DEPENDS "TOPP_MapAlignerTreeGuided_1")
# trafoXML output:
add_test("TOPP_MapAlignerTreeGuided_2" ${TOPP_BIN_PATH}/MapAlignerTreeGuided -test -ini ${DATA_DIR_TOPP}/MapAlignerTreeGuided_parameters.ini -in ${DATA_DIR_TOPP}/MapAlignerTreeGuided_1_input1.featureXML ${DATA_DIR_TOPP}/MapAlignerTreeGuided_1_input2.featureXML ${DATA_DIR_TOPP}/MapAlignerTreeGuided_1_input3.featureXML -trafo_out MapAlignerTreeGuided_2_output1.tmp.trafoXML MapAlignerTreeGuided_2_output2.tmp.trafoXML MapAlignerTreeGuided_2_output3.tmp.trafoXML)
add_test("TOPP_MapAlignerTreeGuided_2_out1" ${DIFF} -in1 MapAlignerTreeGuided_2_output1.tmp.trafoXML -in2 ${DATA_DIR_TOPP}/MapAlignerTreeGuided_2_output1.trafoXML )
set_tests_properties("TOPP_MapAlignerTreeGuided_2_out1" PROPERTIES DEPENDS "TOPP_MapAlignerTreeGuided_2")
add_test("TOPP_MapAlignerTreeGuided_2_out2" ${DIFF} -in1 MapAlignerTreeGuided_2_output2.tmp.trafoXML -in2 ${DATA_DIR_TOPP}/MapAlignerTreeGuided_2_output2.trafoXML )
set_tests_properties("TOPP_MapAlignerTreeGuided_2_out2" PROPERTIES DEPENDS "TOPP_MapAlignerTreeGuided_2")
add_test("TOPP_MapAlignerTreeGuided_2_out3" ${DIFF} -in1 MapAlignerTreeGuided_2_output3.tmp.trafoXML -in2 ${DATA_DIR_TOPP}/MapAlignerTreeGuided_2_output3.trafoXML )
set_tests_properties("TOPP_MapAlignerTreeGuided_2_out3" PROPERTIES DEPENDS "TOPP_MapAlignerTreeGuided_2")
# loading files twice option
add_test("TOPP_MapAlignerTreeGuided_3" ${TOPP_BIN_PATH}/MapAlignerTreeGuided -test -ini ${DATA_DIR_TOPP}/MapAlignerTreeGuided_parameters2.ini -in ${DATA_DIR_TOPP}/MapAlignerTreeGuided_1_input1.featureXML ${DATA_DIR_TOPP}/MapAlignerTreeGuided_1_input2.featureXML ${DATA_DIR_TOPP}/MapAlignerTreeGuided_1_input3.featureXML -out MapAlignerTreeGuided_3_output1.tmp.featureXML MapAlignerTreeGuided_3_output2.tmp.featureXML MapAlignerTreeGuided_3_output3.tmp.featureXML)
add_test("TOPP_MapAlignerTreeGuided_3_out1" ${DIFF} -in1 MapAlignerTreeGuided_3_output1.tmp.featureXML -in2 ${DATA_DIR_TOPP}/MapAlignerTreeGuided_1_output1.featureXML )
set_tests_properties("TOPP_MapAlignerTreeGuided_3_out1" PROPERTIES DEPENDS "TOPP_MapAlignerTreeGuided_3")
add_test("TOPP_MapAlignerTreeGuided_3_out2" ${DIFF} -in1 MapAlignerTreeGuided_3_output2.tmp.featureXML -in2 ${DATA_DIR_TOPP}/MapAlignerTreeGuided_1_output2.featureXML )
set_tests_properties("TOPP_MapAlignerTreeGuided_3_out2" PROPERTIES DEPENDS "TOPP_MapAlignerTreeGuided_3")
add_test("TOPP_MapAlignerTreeGuided_3_out3" ${DIFF} -in1 MapAlignerTreeGuided_3_output3.tmp.featureXML -in2 ${DATA_DIR_TOPP}/MapAlignerTreeGuided_1_output3.featureXML )
set_tests_properties("TOPP_MapAlignerTreeGuided_3_out3" PROPERTIES DEPENDS "TOPP_MapAlignerTreeGuided_3")

#------------------------------------------------------------------------------
# MapRTTransformer tests
# featureXML input:
add_test("TOPP_MapRTTransformer_1" ${TOPP_BIN_PATH}/MapRTTransformer -test -in ${DATA_DIR_TOPP}/MapRTTransformer_1_input.featureXML -trafo_in ${DATA_DIR_TOPP}/MapRTTransformer_trafo_linear.trafoXML -out MapRTTransformer_1_output.tmp.featureXML)
add_test("TOPP_MapRTTransformer_1_out1" ${DIFF} -in1 MapRTTransformer_1_output.tmp.featureXML -in2 ${DATA_DIR_TOPP}/MapRTTransformer_1_output.featureXML )
set_tests_properties("TOPP_MapRTTransformer_1_out1" PROPERTIES DEPENDS "TOPP_MapRTTransformer_1")
# mzML input:
add_test("TOPP_MapRTTransformer_2" ${TOPP_BIN_PATH}/MapRTTransformer -test -in ${DATA_DIR_TOPP}/MapRTTransformer_2_input.mzML -trafo_in ${DATA_DIR_TOPP}/MapRTTransformer_trafo_linear.trafoXML -out MapRTTransformer_2_output.tmp.mzML)
add_test("TOPP_MapRTTransformer_2_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 MapRTTransformer_2_output.tmp.mzML -in2 ${DATA_DIR_TOPP}/MapRTTransformer_2_output.mzML )
set_tests_properties("TOPP_MapRTTransformer_2_out1" PROPERTIES DEPENDS "TOPP_MapRTTransformer_2")
# "invert" option:
add_test("TOPP_MapRTTransformer_3" ${TOPP_BIN_PATH}/MapRTTransformer -test -invert -trafo_in ${DATA_DIR_TOPP}/MapAlignerPoseClustering_1_trafo2.trafoXML -trafo_out MapRTTransformer_3_output.tmp.trafoXML)
add_test("TOPP_MapRTTransformer_3_out1" ${DIFF} -in1 MapRTTransformer_3_output.tmp.trafoXML -in2 ${DATA_DIR_TOPP}/MapRTTransformer_3_trafo.trafoXML )
set_tests_properties("TOPP_MapRTTransformer_3_out1" PROPERTIES DEPENDS "TOPP_MapRTTransformer_3")
# chromatograms
add_test("TOPP_MapRTTransformer_4" ${TOPP_BIN_PATH}/MapRTTransformer -test -in ${DATA_DIR_TOPP}/MapRTTransformer_4_input.chrom.mzML -trafo_in ${DATA_DIR_TOPP}/MapRTTransformer_trafo_linear.trafoXML -out MapRTTransformer_4_output.tmp.mzML)
add_test("TOPP_MapRTTransformer_4_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 MapRTTransformer_4_output.tmp.mzML -in2 ${DATA_DIR_TOPP}/MapRTTransformer_4_output.chrom.mzML )
set_tests_properties("TOPP_MapRTTransformer_4_out1" PROPERTIES DEPENDS "TOPP_MapRTTransformer_4")
# no transformation:
add_test("TOPP_MapRTTransformer_5" ${TOPP_BIN_PATH}/MapRTTransformer -test -in ${DATA_DIR_TOPP}/MapRTTransformer_1_input.featureXML -trafo_in ${DATA_DIR_TOPP}/MapRTTransformer_trafo_none.trafoXML -out MapRTTransformer_5_output.tmp.featureXML)
add_test("TOPP_MapRTTransformer_5_out1" ${DIFF} -in1 MapRTTransformer_5_output.tmp.featureXML -in2 ${DATA_DIR_TOPP}/MapRTTransformer_5_output.featureXML )
set_tests_properties("TOPP_MapRTTransformer_5_out1" PROPERTIES DEPENDS "TOPP_MapRTTransformer_5")
# "store_original_rt" option:
add_test("TOPP_MapRTTransformer_6" ${TOPP_BIN_PATH}/MapRTTransformer -test -in ${DATA_DIR_TOPP}/MapRTTransformer_1_input.featureXML -trafo_in ${DATA_DIR_TOPP}/MapRTTransformer_trafo_linear.trafoXML -out MapRTTransformer_6_output.tmp.featureXML -store_original_rt)
add_test("TOPP_MapRTTransformer_6_out1" ${DIFF} -in1 MapRTTransformer_6_output.tmp.featureXML -in2 ${DATA_DIR_TOPP}/MapRTTransformer_6_output.featureXML )
set_tests_properties("TOPP_MapRTTransformer_6_out1" PROPERTIES DEPENDS "TOPP_MapRTTransformer_6")

#------------------------------------------------------------------------------
# MetaProSIP tests (single SIP feature with approx. RIA of 1% and 37%)
add_test("TOPP_MetaProSIP_1" ${TOPP_BIN_PATH}/MetaProSIP -test -in_mzML ${DATA_DIR_TOPP}/MetaProSIP_1_input.mzML -in_fasta ${DATA_DIR_TOPP}/MetaProSIP_1_input.fasta -in_featureXML ${DATA_DIR_TOPP}/MetaProSIP_1_input.featureXML -out_csv MetaProSIP_1_output_1.tmp.csv -out_peptide_centric_csv MetaProSIP_1_output_2.tmp.csv)
add_test("TOPP_MetaProSIP_1_out1" ${DIFF} -in1 MetaProSIP_1_output_1.tmp.csv -in2 ${DATA_DIR_TOPP}/MetaProSIP_1_output_1.csv )
set_tests_properties("TOPP_MetaProSIP_1_out1" PROPERTIES DEPENDS "TOPP_MetaProSIP_1")
add_test("TOPP_MetaProSIP_1_out2" ${DIFF} -in1 MetaProSIP_1_output_2.tmp.csv -in2 ${DATA_DIR_TOPP}/MetaProSIP_1_output_2.csv )
set_tests_properties("TOPP_MetaProSIP_1_out2" PROPERTIES DEPENDS "TOPP_MetaProSIP_1")

#------------------------------------------------------------------------------
# MascotAdapter tests
add_test("TOPP_MascotAdapter_1" ${TOPP_BIN_PATH}/MascotAdapter -test -ini ${DATA_DIR_TOPP}/MascotAdapter_1_parameters.ini -in ${DATA_DIR_TOPP}/MascotAdapter_1_input.mzData -out MascotAdapter_1_output.tmp.mgf -out_type mgf)
add_test("TOPP_MascotAdapter_1_out1" ${DIFF} -in1 MascotAdapter_1_output.tmp.mgf -in2 ${DATA_DIR_TOPP}/MascotAdapter_1_output.mascot_in )
set_tests_properties("TOPP_MascotAdapter_1_out1" PROPERTIES DEPENDS "TOPP_MascotAdapter_1")
add_test("TOPP_MascotAdapter_2" ${TOPP_BIN_PATH}/MascotAdapter -test -in ${DATA_DIR_TOPP}/MascotAdapter_2_input.mascotXML -out MascotAdapter_2_output.tmp.idXML -out_type idXML)
add_test("TOPP_MascotAdapter_2_out1" ${DIFF} -in1 MascotAdapter_2_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/MascotAdapter_2_output.idXML -whitelist "UserParam type=\"string\" name=\"MascotAdapter:1:in\" value=")
set_tests_properties("TOPP_MascotAdapter_2_out1" PROPERTIES DEPENDS "TOPP_MascotAdapter_2")

#------------------------------------------------------------------------------
# MetaboliteAdductDecharger tests
add_test("TOPP_MetaboliteAdductDecharger_1" ${TOPP_BIN_PATH}/MetaboliteAdductDecharger -test -in ${DATA_DIR_TOPP}/MetaboliteAdductDecharger_input.featureXML -out_fm MetaboliteAdductDecharger_1_output.tmp.featureXML)
add_test("TOPP_MetaboliteAdductDecharger_1_out1" ${DIFF} -in1 MetaboliteAdductDecharger_1_output.tmp.featureXML -in2 ${DATA_DIR_TOPP}/MetaboliteAdductDecharger_1_output.featureXML -whitelist "xml-stylesheet" "consensusXML version=" "consensusElement id=" "<UserParam type=" "<dataProcessing completion_time=" "<software name=")
set_tests_properties("TOPP_MetaboliteAdductDecharger_1_out1" PROPERTIES DEPENDS "TOPP_MetaboliteAdductDecharger_1")
add_test("TOPP_MetaboliteAdductDecharger_2" ${TOPP_BIN_PATH}/MetaboliteAdductDecharger -test -in ${DATA_DIR_TOPP}/MetaboliteAdductDecharger_input.featureXML -out_cm MetaboliteAdductDecharger_2_output_1.tmp.consensusXML -outpairs MetaboliteAdductDecharger_2_output_2.tmp.consensusXML)
add_test("TOPP_MetaboliteAdductDecharger_2_out1" ${DIFF} -in1 MetaboliteAdductDecharger_2_output_1.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/MetaboliteAdductDecharger_2_output_1.consensusXML  -whitelist "xml-stylesheet" "consensusXML version=" "consensusElement id=" "<UserParam type=" "<dataProcessing completion_time=" "<software name=")
add_test("TOPP_MetaboliteAdductDecharger_2_out2" ${DIFF} -in1 MetaboliteAdductDecharger_2_output_2.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/MetaboliteAdductDecharger_2_output_2.consensusXML  -whitelist "xml-stylesheet" "consensusXML version=" "consensusElement id=" "<UserParam type=" "<dataProcessing completion_time=" "<software name=")
set_tests_properties("TOPP_MetaboliteAdductDecharger_2_out1" PROPERTIES DEPENDS "TOPP_MetaboliteAdductDecharger_2")
set_tests_properties("TOPP_MetaboliteAdductDecharger_2_out2" PROPERTIES DEPENDS "TOPP_MetaboliteAdductDecharger_2")

#------------------------------------------------------------------------------
# MSstatsConverter tests
add_test("TOPP_MSstatsConverter_1" ${TOPP_BIN_PATH}/MSstatsConverter -test -in ${DATA_DIR_TOPP}/MSstatsConverter_1_in.consensusXML -retention_time_summarization_method max -in_design ${DATA_DIR_TOPP}/MSstatsConverter_1_design.tsv -out MSstatsConverter_1_out.tmp.csv)
add_test("TOPP_MSstatsConverter_1_out1" ${DIFF} -in1 MSstatsConverter_1_out.tmp.csv -in2 ${DATA_DIR_TOPP}/MSstatsConverter_1_out.csv )
set_tests_properties("TOPP_MSstatsConverter_1_out1" PROPERTIES DEPENDS "TOPP_MSstatsConverter_1")

add_test("TOPP_MSstatsConverter_2" ${TOPP_BIN_PATH}/MSstatsConverter -test -in ${DATA_DIR_TOPP}/MSstatsConverter_2_in.consensusXML -method "ISO" -retention_time_summarization_method manual -in_design ${DATA_DIR_TOPP}/MSstatsConverter_2_design.tsv -out MSstatsConverter_2_out.tmp.csv -msstats_bioreplicate MSstats_BioReplicate -msstats_condition MSstats_Condition -msstats_mixture MSstats_Mixture)
add_test("TOPP_MSstatsConverter_2_out1" ${DIFF} -in1 MSstatsConverter_2_out.tmp.csv -in2 ${DATA_DIR_TOPP}/MSstatsConverter_2_out.csv )
set_tests_properties("TOPP_MSstatsConverter_2_out1" PROPERTIES DEPENDS "TOPP_MSstatsConverter_2")

add_test("TOPP_MSstatsConverter_3" ${TOPP_BIN_PATH}/MSstatsConverter -test -in ${DATA_DIR_TOPP}/MSstatsConverter_3_in.consensusXML -method "ISO" -retention_time_summarization_method manual -in_design ${DATA_DIR_TOPP}/MSstatsConverter_3_design.tsv -out MSstatsConverter_3_out.tmp.csv -msstats_bioreplicate MSstats_BioReplicate -msstats_condition MSstats_Condition -msstats_mixture MSstats_Mixture)
add_test("TOPP_MSstatsConverter_3_out1" ${DIFF} -in1 MSstatsConverter_3_out.tmp.csv -in2 ${DATA_DIR_TOPP}/MSstatsConverter_3_out.csv )
set_tests_properties("TOPP_MSstatsConverter_3_out1" PROPERTIES DEPENDS "TOPP_MSstatsConverter_3")

### SpecLibSearcher tests
add_test("TOPP_SpecLibSearcher_1" ${TOPP_BIN_PATH}/SpecLibSearcher -test -ini ${DATA_DIR_TOPP}/SpecLibSearcher_1_parameters.ini -in ${DATA_DIR_TOPP}/SpecLibSearcher_1.mzML -lib ${DATA_DIR_TOPP}/SpecLibSearcher_1.MSP -out SpecLibSearcher_1.tmp.idXML)
add_test("TOPP_SpecLibSearcher_1_out1" ${DIFF} -in1 SpecLibSearcher_1.tmp.idXML  -in2 ${DATA_DIR_TOPP}/SpecLibSearcher_1.idXML -whitelist "?xml-stylesheet" "IdentificationRun date" "db=")
set_tests_properties("TOPP_SpecLibSearcher_1_out1" PROPERTIES DEPENDS "TOPP_SpecLibSearcher_1")

if(NOT DISABLE_OPENSWATH)
  #------------------------------------------------------------------------------
  # MRM / SWATH tests (from OpenSWATH)
  add_test("TOPP_OpenSwathAssayGenerator_test_1"
${TOPP_BIN_PATH}/OpenSwathAssayGenerator -test -in ${DATA_DIR_TOPP}/OpenSwathAssayGenerator_input.TraML -out OpenSwathAssayGenerator_output.tmp.TraML -out_type TraML -min_transitions 6 -max_transitions 6 -allowed_fragment_types b,y -allowed_fragment_charges 2,3 -enable_detection_specific_losses)
  add_test("TOPP_OpenSwathAssayGenerator_test_1_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathAssayGenerator_output.tmp.TraML -in2 ${DATA_DIR_TOPP}/OpenSwathAssayGenerator_output.TraML)
  set_tests_properties("TOPP_OpenSwathAssayGenerator_test_1_out1" PROPERTIES DEPENDS "TOPP_OpenSwathAssayGenerator_test_1")
  add_test("TOPP_OpenSwathAssayGenerator_test_2"
${TOPP_BIN_PATH}/OpenSwathAssayGenerator -test -in ${DATA_DIR_TOPP}/OpenSwathAssayGenerator_input_2.TraML -out OpenSwathAssayGenerator_output_2.tmp.TraML -out_type TraML -allowed_fragment_charges 1,2,3,4 -enable_ipf -unimod_file ${DATA_DIR_TOPP}/OpenSwathAssayGenerator_input_2_unimod.xml)
  add_test("TOPP_OpenSwathAssayGenerator_test_2_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathAssayGenerator_output_2.tmp.TraML -in2 ${DATA_DIR_TOPP}/OpenSwathAssayGenerator_output_2.TraML)
  set_tests_properties("TOPP_OpenSwathAssayGenerator_test_2_out1" PROPERTIES DEPENDS "TOPP_OpenSwathAssayGenerator_test_2")
  add_test("TOPP_OpenSwathAssayGenerator_test_3"
${TOPP_BIN_PATH}/OpenSwathAssayGenerator -test -in ${DATA_DIR_TOPP}/OpenSwathAssayGenerator_input_2.TraML -out OpenSwathAssayGenerator_output_3.tmp.TraML -out_type TraML -allowed_fragment_charges 1,2,3,4 -enable_ipf -unimod_file ${DATA_DIR_TOPP}/OpenSwathAssayGenerator_input_3_unimod.xml)
  add_test("TOPP_OpenSwathAssayGenerator_test_3_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathAssayGenerator_output_3.tmp.TraML -in2 ${DATA_DIR_TOPP}/OpenSwathAssayGenerator_output_3.TraML)
  set_tests_properties("TOPP_OpenSwathAssayGenerator_test_3_out1" PROPERTIES DEPENDS "TOPP_OpenSwathAssayGenerator_test_3")

  add_test("TOPP_OpenSwathDecoyGenerator_test_1"
${TOPP_BIN_PATH}/OpenSwathDecoyGenerator -test -in ${DATA_DIR_TOPP}/OpenSwathDecoyGenerator_input.TraML -out OpenSwathDecoyGenerator.tmp.TraML -out_type TraML -method pseudo-reverse -separate -switchKR false)
  add_test("TOPP_OpenSwathDecoyGenerator_test_1_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathDecoyGenerator.tmp.TraML -in2 ${DATA_DIR_TOPP}/OpenSwathDecoyGenerator_output.TraML)
  set_tests_properties("TOPP_OpenSwathDecoyGenerator_test_1_out1" PROPERTIES DEPENDS "TOPP_OpenSwathDecoyGenerator_test_1")
  add_test("TOPP_OpenSwathDecoyGenerator_test_2"
${TOPP_BIN_PATH}/OpenSwathDecoyGenerator -test -in ${DATA_DIR_TOPP}/OpenSwathDecoyGenerator_input_2.TraML -out OpenSwathDecoyGenerator_2.tmp.TraML -out_type TraML -method pseudo-reverse -product_mz_threshold 0.8 -switchKR false)
  add_test("TOPP_OpenSwathDecoyGenerator_test_2_out1" ${DIFF} -in1 OpenSwathDecoyGenerator_2.tmp.TraML -in2 ${DATA_DIR_TOPP}/OpenSwathDecoyGenerator_output_2.TraML)
  set_tests_properties("TOPP_OpenSwathDecoyGenerator_test_2_out1" PROPERTIES DEPENDS "TOPP_OpenSwathDecoyGenerator_test_2")
  add_test("TOPP_OpenSwathDecoyGenerator_test_3"
${TOPP_BIN_PATH}/OpenSwathDecoyGenerator -test -in ${DATA_DIR_TOPP}/OpenSwathDecoyGenerator_input_3.TraML -out OpenSwathDecoyGenerator_3.tmp.TraML -out_type TraML -method pseudo-reverse -separate -switchKR false)
  add_test("TOPP_OpenSwathDecoyGenerator_test_3_out1" ${DIFF} -in1 OpenSwathDecoyGenerator_3.tmp.TraML -in2 ${DATA_DIR_TOPP}/OpenSwathDecoyGenerator_output_3.TraML)
  set_tests_properties("TOPP_OpenSwathDecoyGenerator_test_3_out1" PROPERTIES DEPENDS "TOPP_OpenSwathDecoyGenerator_test_3")
  add_test("TOPP_OpenSwathDecoyGenerator_test_4"
    ${TOPP_BIN_PATH}/OpenSwathDecoyGenerator -test -in ${DATA_DIR_TOPP}/OpenSwathDecoyGenerator_input_4.tsv -out OpenSwathDecoyGenerator_4.tmp.TraML -out_type TraML -method pseudo-reverse -separate -switchKR true -enable_detection_unspecific_losses  -enable_detection_specific_losses -min_decoy_fraction 0.4)
  add_test("TOPP_OpenSwathDecoyGenerator_test_4_out1" ${DIFF} -in1 OpenSwathDecoyGenerator_4.tmp.TraML -in2 ${DATA_DIR_TOPP}/OpenSwathDecoyGenerator_output_4.TraML)
  set_tests_properties("TOPP_OpenSwathDecoyGenerator_test_4_out1" PROPERTIES DEPENDS "TOPP_OpenSwathDecoyGenerator_test_4")
  add_test("TOPP_OpenSwathDecoyGenerator_test_5"
    ${TOPP_BIN_PATH}/OpenSwathDecoyGenerator -test -in ${DATA_DIR_TOPP}/OpenSwathDecoyGenerator_input_5.tsv -out OpenSwathDecoyGenerator_5.tmp.TraML -out_type TraML -method pseudo-reverse -switchKR false -min_decoy_fraction 0.4)
  add_test("TOPP_OpenSwathDecoyGenerator_test_5_out1" ${DIFF} -in1 OpenSwathDecoyGenerator_5.tmp.TraML -in2 ${DATA_DIR_TOPP}/OpenSwathDecoyGenerator_output_5.TraML)
  set_tests_properties("TOPP_OpenSwathDecoyGenerator_test_5_out1" PROPERTIES DEPENDS "TOPP_OpenSwathDecoyGenerator_test_5")

  # TargetedFileConverter : convert from and to different targeted files
  # Contains heavy and light transitions in 2 different notations
  add_test("TOPP_ConvertTSVToTraML_test_1" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ${DATA_DIR_TOPP}/ConvertTSVToTraML_1_input.tsv -out ConvertTSVToTraML_output.tmp.TraML -out_type TraML)
  add_test("TOPP_ConvertTSVToTraML_test_1_out1" ${DIFF} -in1 ConvertTSVToTraML_output.tmp.TraML -in2 ${DATA_DIR_TOPP}/ConvertTSVToTraML_output.TraML)
  set_tests_properties("TOPP_ConvertTSVToTraML_test_1_out1" PROPERTIES DEPENDS "TOPP_ConvertTSVToTraML_test_1")
  add_test("TOPP_ConvertTSVToTraML_test_2" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ${DATA_DIR_TOPP}/ConvertTSVToTraML_2_input.tsv -out ConvertTSVToTraML_2_output.tmp.TraML -algorithm:retentionTimeInterpretation minutes -out_type TraML)
  add_test("TOPP_ConvertTSVToTraML_test_2_out1" ${DIFF} -in1 ConvertTSVToTraML_2_output.tmp.TraML -in2 ${DATA_DIR_TOPP}/ConvertTSVToTraML_2_output.TraML)
  set_tests_properties("TOPP_ConvertTSVToTraML_test_2_out1" PROPERTIES DEPENDS "TOPP_ConvertTSVToTraML_test_2")

  add_test("TOPP_ConvertTSVToTraML_test_3" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ${DATA_DIR_TOPP}/ConvertTSVToTraML_2_input.tsv -out ConvertTSVToTraML_3_output.tmp.TraML -algorithm:retentionTimeInterpretation seconds -out_type TraML)
  add_test("TOPP_ConvertTSVToTraML_test_3_out1" ${DIFF} -in1 ConvertTSVToTraML_3_output.tmp.TraML -in2 ${DATA_DIR_TOPP}/ConvertTSVToTraML_3_output.TraML)
  set_tests_properties("TOPP_ConvertTSVToTraML_test_3_out1" PROPERTIES DEPENDS "TOPP_ConvertTSVToTraML_test_3")
  add_test("TOPP_ConvertTSVToTraML_test_4" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ${DATA_DIR_TOPP}/ConvertTSVToTraML_4_input.mrm -out ConvertTSVToTraML_4_output.tmp.TraML -out_type TraML)
  add_test("TOPP_ConvertTSVToTraML_test_4_out1" ${DIFF} -in1 ConvertTSVToTraML_4_output.tmp.TraML -in2 ${DATA_DIR_TOPP}/ConvertTSVToTraML_4_output.TraML)
  set_tests_properties("TOPP_ConvertTSVToTraML_test_4_out1" PROPERTIES DEPENDS "TOPP_ConvertTSVToTraML_test_4")
  # Metabolomics
  add_test("TOPP_ConvertTSVToTraML_test_5" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ${DATA_DIR_TOPP}/ConvertTSVToTraML_5_input.tsv -out ConvertTSVToTraML_5_output.tmp.TraML -out_type TraML)
  add_test("TOPP_ConvertTSVToTraML_test_5_out1" ${DIFF} -in1 ConvertTSVToTraML_5_output.tmp.TraML -in2 ${DATA_DIR_TOPP}/ConvertTSVToTraML_5_output.TraML)
  set_tests_properties("TOPP_ConvertTSVToTraML_test_5_out1" PROPERTIES DEPENDS "TOPP_ConvertTSVToTraML_test_5")
  add_test("TOPP_ConvertTSVToTraML_test_6" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ${DATA_DIR_TOPP}/ConvertTSVToTraML_6_input.tsv -out ConvertTSVToTraML_6_output.tmp.TraML -out_type TraML)
  add_test("TOPP_ConvertTSVToTraML_test_6_out1" ${DIFF} -in1 ConvertTSVToTraML_6_output.tmp.TraML -in2 ${DATA_DIR_TOPP}/ConvertTSVToTraML_6_output.TraML)
  set_tests_properties("TOPP_ConvertTSVToTraML_test_6_out1" PROPERTIES DEPENDS "TOPP_ConvertTSVToTraML_test_6")
  add_test("TOPP_ConvertTSVToTraML_test_7" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ${DATA_DIR_TOPP}/ConvertTSVToTraML_7_input_Skyline.tsv -out ConvertTSVToTraML_7_output.tmp.TraML -out_type TraML)
  add_test("TOPP_ConvertTSVToTraML_test_7_out1" ${DIFF} -in1 ConvertTSVToTraML_7_output.tmp.TraML -in2 ${DATA_DIR_TOPP}/ConvertTSVToTraML_7_output.TraML)
  set_tests_properties("TOPP_ConvertTSVToTraML_test_7_out1" PROPERTIES DEPENDS "TOPP_ConvertTSVToTraML_test_7")

  # Proteomics PQP - TraML Conversion
  add_test("TOPP_TargetedFileConverter_test_1_prepare" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ${DATA_DIR_TOPP}/TargetedFileConverter_1_input.TraML -out TargetedFileConverter_1_output.tmp.pqp -out_type pqp)
  add_test("TOPP_TargetedFileConverter_test_1" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in TargetedFileConverter_1_output.tmp.pqp -in_type pqp -out TargetedFileConverter_1_output.tmp.TraML -out_type TraML)
  add_test("TOPP_TargetedFileConverter_test_1_out1" ${DIFF} -in1 TargetedFileConverter_1_output.tmp.TraML -in2 ${DATA_DIR_TOPP}/TargetedFileConverter_1_output.TraML)
  set_tests_properties("TOPP_TargetedFileConverter_test_1" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_1_prepare")
  set_tests_properties("TOPP_TargetedFileConverter_test_1_out1" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_1")
  add_test("TOPP_TargetedFileConverter_test_2" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in TargetedFileConverter_1_output.tmp.pqp -in_type pqp -out TargetedFileConverter_2_output.tmp.TraML -out_type TraML -legacy_traml_id)
  add_test("TOPP_TargetedFileConverter_test_2_out1" ${DIFF} -in1 TargetedFileConverter_2_output.tmp.TraML -in2 ${DATA_DIR_TOPP}/TargetedFileConverter_2_output.TraML)
  set_tests_properties("TOPP_TargetedFileConverter_test_2" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_1_prepare")
  set_tests_properties("TOPP_TargetedFileConverter_test_2_out1" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_2")

  # Metabolomics PQP - TraML Conversion
  add_test("TOPP_TargetedFileConverter_test_3_prepare" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ${DATA_DIR_TOPP}/TargetedFileConverter_3_input.TraML -out TargetedFileConverter_3_output.tmp.pqp -out_type pqp)
  add_test("TOPP_TargetedFileConverter_test_3" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in TargetedFileConverter_3_output.tmp.pqp -in_type pqp -out TargetedFileConverter_3_output.tmp.TraML -out_type TraML)
  add_test("TOPP_TargetedFileConverter_test_3_out1" ${DIFF} -in1 TargetedFileConverter_3_output.tmp.TraML -in2 ${DATA_DIR_TOPP}/TargetedFileConverter_3_output.TraML)
  set_tests_properties("TOPP_TargetedFileConverter_test_3" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_3_prepare")
  set_tests_properties("TOPP_TargetedFileConverter_test_3_out1" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_3")
  add_test("TOPP_TargetedFileConverter_test_4" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in TargetedFileConverter_3_output.tmp.pqp -in_type pqp -out TargetedFileConverter_4_output.tmp.TraML -out_type TraML -legacy_traml_id)
  add_test("TOPP_TargetedFileConverter_test_4_out1" ${DIFF} -in1 TargetedFileConverter_4_output.tmp.TraML -in2 ${DATA_DIR_TOPP}/TargetedFileConverter_4_output.TraML)
  set_tests_properties("TOPP_TargetedFileConverter_test_4_out1" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_4")
  set_tests_properties("TOPP_TargetedFileConverter_test_4" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_3")

  # Metabolomics with Adducts PQP/TSV/TraML conversion
  add_test("TOPP_TargetedFileConverter_test_10_prepare" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ${DATA_DIR_TOPP}/TargetedFileConverter_10_input.pqp  -out TargetedFileConverter_10_output.tmp.TraML -out_type TraML -legacy_traml_id)
  add_test("TOPP_TargetedFileConverter_test_10" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in TargetedFileConverter_10_output.tmp.TraML -in_type TraML -out TargetedFileConverter_10_output.tmp.tsv -out_type tsv)
  add_test("TOPP_TargetedFileConverter_test_10_out1" ${DIFF} -in1 TargetedFileConverter_10_output.tmp.tsv -in2 ${DATA_DIR_TOPP}/TargetedFileConverter_10_output.tsv)
  set_tests_properties("TOPP_TargetedFileConverter_test_10_out1" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_10")
  set_tests_properties("TOPP_TargetedFileConverter_test_10" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_10_prepare")

  add_test("TOPP_TargetedFileConverter_test_11_prepare" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ${DATA_DIR_TOPP}/TargetedFileConverter_11_input.TraML -out TargetedFileConverter_11_output.tmp.pqp -out_type pqp)
  add_test("TOPP_TargetedFileConverter_test_11" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in TargetedFileConverter_11_output.tmp.pqp -in_type pqp -out TargetedFileConverter_11_output.tmp.TraML -out_type TraML)
  add_test("TOPP_TargetedFileConverter_test_11_out1" ${DIFF} -whitelist "id=" -in1 TargetedFileConverter_11_output.tmp.TraML -in2 ${DATA_DIR_TOPP}/TargetedFileConverter_11_input.TraML)
  set_tests_properties("TOPP_TargetedFileConverter_test_11_out1" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_11")
  set_tests_properties("TOPP_TargetedFileConverter_test_11" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_11_prepare")


  # Proteomics PQP/TSV/TraML peptidoforms conversion
  add_test("TOPP_TargetedFileConverter_test_8_prepare" ${TOPP_BIN_PATH}/TargetedFileConverter -algorithm:override_group_label_check -test -in ${DATA_DIR_TOPP}/TargetedFileConverter_8_input.tsv -out TargetedFileConverter_8_output.tmp.TraML -out_type TraML)
  add_test("TOPP_TargetedFileConverter_test_8" ${TOPP_BIN_PATH}/TargetedFileConverter -algorithm:override_group_label_check -test -in TargetedFileConverter_8_output.tmp.TraML -in_type TraML -out TargetedFileConverter_8_output.tmp.tsv -out_type tsv)
  add_test("TOPP_TargetedFileConverter_test_8_out1" ${DIFF} -in1 TargetedFileConverter_8_output.tmp.tsv -in2 ${DATA_DIR_TOPP}/TargetedFileConverter_8_input.tsv)
  set_tests_properties("TOPP_TargetedFileConverter_test_8_out1" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_8")
  set_tests_properties("TOPP_TargetedFileConverter_test_8" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_8_prepare")

  add_test("TOPP_TargetedFileConverter_test_9_prepare" ${TOPP_BIN_PATH}/TargetedFileConverter -algorithm:override_group_label_check -test -in ${DATA_DIR_TOPP}/TargetedFileConverter_9_input.TraML -out TargetedFileConverter_9_output.tmp.pqp -out_type pqp)
  add_test("TOPP_TargetedFileConverter_test_9" ${TOPP_BIN_PATH}/TargetedFileConverter -algorithm:override_group_label_check -test -in TargetedFileConverter_9_output.tmp.pqp -in_type pqp -out TargetedFileConverter_9_output.tmp.TraML -out_type TraML)
  add_test("TOPP_TargetedFileConverter_test_9_out1" ${DIFF} -whitelist "id=" -in1 TargetedFileConverter_9_output.tmp.TraML -in2 ${DATA_DIR_TOPP}/TargetedFileConverter_9_input.TraML)
  set_tests_properties("TOPP_TargetedFileConverter_test_9_out1" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_9")
  set_tests_properties("TOPP_TargetedFileConverter_test_9" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_9_prepare")

  # additional tests (previously ConvertTSVToTraML tests)
  add_test("TOPP_TargetedFileConverter_test_5" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ${DATA_DIR_TOPP}/ConvertTSVToTraML_output.TraML -out ConvertTraMLToTSV_output.tmp.tsv -out_type tsv)
  add_test("TOPP_TargetedFileConverter_test_5_out1" ${DIFF} -in1 ConvertTraMLToTSV_output.tmp.tsv -in2 ${DATA_DIR_TOPP}/ConvertTraMLToTSV_output.tsv)
  set_tests_properties("TOPP_TargetedFileConverter_test_5_out1" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_5")

  # protein group conversion test PQP
  add_test("TOPP_TargetedFileConverter_test_12_prepare" ${TOPP_BIN_PATH}/TargetedFileConverter -algorithm:override_group_label_check -test -in ${DATA_DIR_TOPP}/TargetedFileConverter_12_input.TraML -out TargetedFileConverter_12_output.tmp.pqp -out_type pqp)
  add_test("TOPP_TargetedFileConverter_test_12" ${TOPP_BIN_PATH}/TargetedFileConverter -algorithm:override_group_label_check -test -in TargetedFileConverter_12_output.tmp.pqp -in_type pqp -out TargetedFileConverter_12_output.tmp.TraML -out_type TraML)
  add_test("TOPP_TargetedFileConverter_test_12_out1" ${DIFF} -whitelist "id=" -in1 TargetedFileConverter_12_output.tmp.TraML -in2 ${DATA_DIR_TOPP}/TargetedFileConverter_12_input.TraML)
  set_tests_properties("TOPP_TargetedFileConverter_test_12_out1" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_12")
  set_tests_properties("TOPP_TargetedFileConverter_test_12" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_12_prepare")

  # protein group conversion test TraML
  add_test("TOPP_TargetedFileConverter_test_13" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ${DATA_DIR_TOPP}/TargetedFileConverter_12_input.tsv -out TargetedFileConverter_13_output.tmp.TraML -out_type TraML)
  add_test("TOPP_TargetedFileConverter_test_13_out1" ${DIFF} -in1 TargetedFileConverter_13_output.tmp.TraML -in2 ${DATA_DIR_TOPP}/TargetedFileConverter_13_output.TraML)
  set_tests_properties("TOPP_TargetedFileConverter_test_13_out1" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_13")

  # Metabolomics
  add_test("TOPP_TargetedFileConverter_test_6" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ${DATA_DIR_TOPP}/ConvertTSVToTraML_5_output.TraML -out ConvertTraMLToTSV_output_2.tmp.tsv -out_type tsv)
  add_test("TOPP_TargetedFileConverter_test_6_out1" ${DIFF} -in1 ConvertTraMLToTSV_output_2.tmp.tsv -in2 ${DATA_DIR_TOPP}/ConvertTraMLToTSV_output_2.tsv)
  set_tests_properties("TOPP_TargetedFileConverter_test_6_out1" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_6")

  add_test("TOPP_TargetedFileConverter_test_6_convert_back" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ConvertTraMLToTSV_output_2.tmp.tsv -in_type tsv -out ConvertTraMLToTSV_output_2.tsv.back.TraML -out_type TraML)
  add_test("TOPP_TargetedFileConverter_test_6_out2" ${DIFF} -in1 ConvertTraMLToTSV_output_2.tmp.tsv -in2 ${DATA_DIR_TOPP}/ConvertTraMLToTSV_output_2.tsv)
  set_tests_properties("TOPP_TargetedFileConverter_test_6_convert_back" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_6")
  set_tests_properties("TOPP_TargetedFileConverter_test_6_out2" PROPERTIES DEPENDS "TOPP_TargetedFileConverter_test_6_convert_back")

  add_test("TOPP_MRMMapper_test_1" ${TOPP_BIN_PATH}/MRMMapper -in ${DATA_DIR_TOPP}/MRMMapping_input.chrom.mzML -tr ${DATA_DIR_TOPP}/MRMMapping_input.TraML -out MRMMapping_output_1.chrom.tmp.mzML -test -algorithm:precursor_tolerance 0.3 -algorithm:product_tolerance 0.3)
  add_test("TOPP_MRMMapper_test_1_out_1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 MRMMapping_output_1.chrom.tmp.mzML -in2 ${DATA_DIR_TOPP}/MRMMapping_output.chrom.mzML)
  set_tests_properties("TOPP_MRMMapper_test_1_out_1" PROPERTIES DEPENDS "TOPP_MRMMapper_test_1")

  add_test("TOPP_MRMMapper_test_2" ${TOPP_BIN_PATH}/MRMMapper -in ${DATA_DIR_TOPP}/MRMMapping_input.chrom.mzML -tr ${DATA_DIR_TOPP}/MRMMapping_input_2.TraML -out MRMMapping_output_2.chrom.tmp.mzML -test -algorithm:precursor_tolerance 0.01 -algorithm:product_tolerance 0.01)
  add_test("TOPP_MRMMapper_test_2_out_1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 MRMMapping_output_2.chrom.tmp.mzML -in2 ${DATA_DIR_TOPP}/MRMMapping_output_2.chrom.mzML)
  set_tests_properties("TOPP_MRMMapper_test_2_out_1" PROPERTIES DEPENDS "TOPP_MRMMapper_test_2")

  add_test("TOPP_OpenSwathFeatureXMLToTSV_test_1" ${TOPP_BIN_PATH}/OpenSwathFeatureXMLToTSV -in ${DATA_DIR_TOPP}/OpenSwathFeatureXMLToTSV_input.featureXML -tr ${DATA_DIR_TOPP}/OpenSwathFeatureXMLToTSV_input.TraML -out OpenSwathFeatureXMLToTSV_output.short.tmp.csv -short_format -test)
  add_test("TOPP_OpenSwathFeatureXMLToTSV_test_1_out1" ${DIFF} -in1 OpenSwathFeatureXMLToTSV_output.short.tmp.csv -in2 ${DATA_DIR_TOPP}/OpenSwathFeatureXMLToTSV_output.short.csv)
  set_tests_properties("TOPP_OpenSwathFeatureXMLToTSV_test_1_out1" PROPERTIES DEPENDS "TOPP_OpenSwathFeatureXMLToTSV_test_1")
  add_test("TOPP_OpenSwathFeatureXMLToTSV_test_2" ${TOPP_BIN_PATH}/OpenSwathFeatureXMLToTSV -in ${DATA_DIR_TOPP}/OpenSwathFeatureXMLToTSV_input.featureXML -tr ${DATA_DIR_TOPP}/OpenSwathFeatureXMLToTSV_input.TraML -out OpenSwathFeatureXMLToTSV_output.long.tmp.csv -test)
  add_test("TOPP_OpenSwathFeatureXMLToTSV_test_2_out1" ${DIFF} -in1 OpenSwathFeatureXMLToTSV_output.long.tmp.csv -in2 ${DATA_DIR_TOPP}/OpenSwathFeatureXMLToTSV_output.long.csv)
  set_tests_properties("TOPP_OpenSwathFeatureXMLToTSV_test_2_out1" PROPERTIES DEPENDS "TOPP_OpenSwathFeatureXMLToTSV_test_2")
  add_test("TOPP_OpenSwathFeatureXMLToTSV_test_3" ${TOPP_BIN_PATH}/OpenSwathFeatureXMLToTSV -in ${DATA_DIR_TOPP}/OpenSwathFeatureXMLToTSV_input.featureXML -tr ${DATA_DIR_TOPP}/OpenSwathFeatureXMLToTSV_input.TraML -out OpenSwathFeatureXMLToTSV_3_output.short.tmp.csv -short_format -best_scoring_peptide main_var_xx_lda_prelim_score -test)
  add_test("TOPP_OpenSwathFeatureXMLToTSV_test_3_out1" ${DIFF} -in1 OpenSwathFeatureXMLToTSV_3_output.short.tmp.csv -in2 ${DATA_DIR_TOPP}/OpenSwathFeatureXMLToTSV_3_output.short.csv)
  set_tests_properties("TOPP_OpenSwathFeatureXMLToTSV_test_3_out1" PROPERTIES DEPENDS "TOPP_OpenSwathFeatureXMLToTSV_test_3")

  add_test("TOPP_OpenSwathChromatogramExtractor_test_1" ${TOPP_BIN_PATH}/OpenSwathChromatogramExtractor -in ${DATA_DIR_TOPP}/OpenSwathChromatogramExtractor_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathChromatogramExtractor_input.TraML -out OpenSwathChromatogramExtractor_output.tmp.mzML -test)
  add_test("TOPP_OpenSwathChromatogramExtractor_test_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 OpenSwathChromatogramExtractor_output.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathChromatogramExtractor_output.mzML)
  set_tests_properties("TOPP_OpenSwathChromatogramExtractor_test_1_out1" PROPERTIES DEPENDS "TOPP_OpenSwathChromatogramExtractor_test_1")
  add_test("TOPP_OpenSwathChromatogramExtractor_test_2" ${TOPP_BIN_PATH}/OpenSwathChromatogramExtractor -in ${DATA_DIR_TOPP}/OpenSwathChromatogramExtractor_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathChromatogramExtractor_input.TraML -rt_norm ${DATA_DIR_TOPP}/OpenSwathChromatogramExtractor_input.trafoXML -out OpenSwathChromatogramExtractor_output_2.tmp.mzML -test -rt_window 50)
  add_test("TOPP_OpenSwathChromatogramExtractor_test_2_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 OpenSwathChromatogramExtractor_output_2.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathChromatogramExtractor_output_2.mzML)
  set_tests_properties("TOPP_OpenSwathChromatogramExtractor_test_2_out1" PROPERTIES DEPENDS "TOPP_OpenSwathChromatogramExtractor_test_2")
  # test for MS1 extraction
  add_test("TOPP_OpenSwathChromatogramExtractor_test_3" ${TOPP_BIN_PATH}/OpenSwathChromatogramExtractor -in ${DATA_DIR_TOPP}/OpenSwathChromatogramExtractor_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathChromatogramExtractor_input.TraML -out OpenSwathChromatogramExtractor_output_3.tmp.mzML -test -extract_MS1)
  add_test("TOPP_OpenSwathChromatogramExtractor_test_3_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 OpenSwathChromatogramExtractor_output_3.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathChromatogramExtractor_output_3.mzML)
  set_tests_properties("TOPP_OpenSwathChromatogramExtractor_test_3_out1" PROPERTIES DEPENDS "TOPP_OpenSwathChromatogramExtractor_test_3")
  # test ion mobility extraction
  add_test("TOPP_OpenSwathChromatogramExtractor_test_4" ${TOPP_BIN_PATH}/OpenSwathChromatogramExtractor -in ${DATA_DIR_TOPP}/OpenSwathChromatogramExtractor_4_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathChromatogramExtractor_4_input.TraML -out OpenSwathChromatogramExtractor_output_4.tmp.mzML -ion_mobility_window 0.05 -is_swath -test)
  add_test("TOPP_OpenSwathChromatogramExtractor_test_4_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 OpenSwathChromatogramExtractor_output_4.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathChromatogramExtractor_4_output.mzML)
  set_tests_properties("TOPP_OpenSwathChromatogramExtractor_test_4_out1" PROPERTIES DEPENDS "TOPP_OpenSwathChromatogramExtractor_test_4")
  # test metabolomic extraction for OpenSwathChromatogramExtractor
  add_test("TOPP_OpenSwathChromatogramExtractor_test_5" ${TOPP_BIN_PATH}/OpenSwathChromatogramExtractor -in ${DATA_DIR_TOPP}/OpenSwathChromatogramExtractor_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathChromatogramExtractor_5_input.TraML -out OpenSwathChromatogramExtractor_5_output.tmp.mzML -test -extract_MS1)
  add_test("TOPP_OpenSwathChromatogramExtractor_test_5_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 OpenSwathChromatogramExtractor_5_output.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathChromatogramExtractor_5_output.mzML)
  set_tests_properties("TOPP_OpenSwathChromatogramExtractor_test_5_out1" PROPERTIES DEPENDS "TOPP_OpenSwathChromatogramExtractor_test_5")

  add_test("TOPP_OpenSwathAnalyzer_test_1" ${TOPP_BIN_PATH}/OpenSwathAnalyzer -in ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input_chrom.mzML -tr ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input.TraML -out OpenSwathAnalyzer_1_output.featureXML.tmp.featureXML -rt_norm ${DATA_DIR_TOPP}/OpenSwathAnalyzer_input.trafoXML -algorithm:TransitionGroupPicker:PeakPickerChromatogram:peak_width 40.0 -algorithm:TransitionGroupPicker:PeakPickerChromatogram:method legacy -test)
  add_test("TOPP_OpenSwathAnalyzer_test_1_out1" ${DIFF} -whitelist "id="  -in1 OpenSwathAnalyzer_1_output.featureXML.tmp.featureXML -in2 ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_output.featureXML)
  set_tests_properties("TOPP_OpenSwathAnalyzer_test_1_out1" PROPERTIES DEPENDS "TOPP_OpenSwathAnalyzer_test_1")
  add_test("TOPP_OpenSwathAnalyzer_test_2" ${TOPP_BIN_PATH}/OpenSwathAnalyzer -in ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input_chrom.mzML -tr ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input.TraML -out OpenSwathAnalyzer_2_output.featureXML.tmp.featureXML -swath_files ${DATA_DIR_TOPP}/OpenSwathAnalyzer_2_swathfile.mzML -algorithm:TransitionGroupPicker:PeakPickerChromatogram:peak_width 40.0 -algorithm:TransitionGroupPicker:PeakPickerChromatogram:method legacy -test)
  add_test("TOPP_OpenSwathAnalyzer_test_2_out1" ${DIFF} -whitelist "id="  -in1 OpenSwathAnalyzer_2_output.featureXML.tmp.featureXML -in2 ${DATA_DIR_TOPP}/OpenSwathAnalyzer_2_output.featureXML)
  set_tests_properties("TOPP_OpenSwathAnalyzer_test_2_out1" PROPERTIES DEPENDS "TOPP_OpenSwathAnalyzer_test_2")
  add_test("TOPP_OpenSwathAnalyzer_test_5_mod" ${TOPP_BIN_PATH}/OpenSwathAnalyzer -in ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input_chrom.mzML -tr ${DATA_DIR_TOPP}/OpenSwathAnalyzer_mod_input.TraML -out OpenSwathAnalyzer_5_output.featureXML.tmp.featureXML -swath_files ${DATA_DIR_TOPP}/OpenSwathAnalyzer_2_swathfile.mzML -ini ${DATA_DIR_TOPP}/OpenSwathAnalyzer_5.ini -test)
  add_test("TOPP_OpenSwathAnalyzer_test_5_mod_out1" ${DIFF} -whitelist "id="  -in1 OpenSwathAnalyzer_5_output.featureXML.tmp.featureXML -in2 ${DATA_DIR_TOPP}/OpenSwathAnalyzer_5_output.featureXML)
  set_tests_properties("TOPP_OpenSwathAnalyzer_test_5_mod_out1" PROPERTIES DEPENDS "TOPP_OpenSwathAnalyzer_test_5_mod")
  add_test("TOPP_OpenSwathAnalyzer_test_6_nomod" ${TOPP_BIN_PATH}/OpenSwathAnalyzer -in ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input_chrom.mzML -tr ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input.TraML -out OpenSwathAnalyzer_6_output.featureXML.tmp.featureXML -swath_files ${DATA_DIR_TOPP}/OpenSwathAnalyzer_2_swathfile.mzML -ini ${DATA_DIR_TOPP}/OpenSwathAnalyzer_5.ini -test)
  add_test("TOPP_OpenSwathAnalyzer_test_6_nomod_out1" ${DIFF} -whitelist "id="  -in1 OpenSwathAnalyzer_6_output.featureXML.tmp.featureXML -in2 ${DATA_DIR_TOPP}/OpenSwathAnalyzer_6_output.featureXML)
  set_tests_properties("TOPP_OpenSwathAnalyzer_test_6_nomod_out1" PROPERTIES DEPENDS "TOPP_OpenSwathAnalyzer_test_6_nomod")
  #// TODO write test with TrafoXML


  add_test("TOPP_OpenSwathAnalyzer_test_7_backgroundSubtraction" ${TOPP_BIN_PATH}/OpenSwathAnalyzer -in ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input_chrom.mzML -tr ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input.TraML -out OpenSwathAnalyzer_7_output.featureXML.tmp.featureXML -ini ${DATA_DIR_TOPP}/OpenSwathAnalyzer_7_backgroundSubtraction.ini -test)
  add_test("TOPP_OpenSwathAnalyzer_test_7_backgroundSubtraction_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathAnalyzer_7_output.featureXML.tmp.featureXML -in2 ${DATA_DIR_TOPP}/OpenSwathAnalyzer_7_output.featureXML)
  set_tests_properties("TOPP_OpenSwathAnalyzer_test_7_backgroundSubtraction_out1" PROPERTIES DEPENDS "TOPP_OpenSwathAnalyzer_test_7_backgroundSubtraction")
  add_test("TOPP_OpenSwathAnalyzer_test_8" ${TOPP_BIN_PATH}/OpenSwathAnalyzer -in ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input_chrom.mzML -tr ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input.TraML -out OpenSwathAnalyzer_8_output.featureXML.tmp.featureXML -rt_norm ${DATA_DIR_TOPP}/OpenSwathAnalyzer_input.trafoXML -algorithm:TransitionGroupPicker:PeakPickerChromatogram:peak_width 40.0 -algorithm:TransitionGroupPicker:PeakPickerChromatogram:method legacy -algorithm:Scores:use_ms1_mi -algorithm:Scores:use_mi_score -algorithm:Scores:use_total_mi_score -test)
  add_test("TOPP_OpenSwathAnalyzer_test_8_out1" ${DIFF} -whitelist "id="  -in1 OpenSwathAnalyzer_8_output.featureXML.tmp.featureXML -in2 ${DATA_DIR_TOPP}/OpenSwathAnalyzer_8_output.featureXML)
  set_tests_properties("TOPP_OpenSwathAnalyzer_test_8_out1" PROPERTIES DEPENDS "TOPP_OpenSwathAnalyzer_test_8")

  add_test("TOPP_OpenSwathAnalyzer_test_9" ${TOPP_BIN_PATH}/OpenSwathAnalyzer -in ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input_chrom.mzML -tr ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input.TraML -out OpenSwathAnalyzer_9_output.featureXML.tmp.featureXML -rt_norm ${DATA_DIR_TOPP}/OpenSwathAnalyzer_input.trafoXML -algorithm:TransitionGroupPicker:PeakPickerChromatogram:peak_width 40.0 -algorithm:TransitionGroupPicker:PeakPickerChromatogram:method legacy -algorithm:Scores:use_mi_score -test)
  add_test("TOPP_OpenSwathAnalyzer_test_9_out1" ${DIFF} -whitelist "id="  -in1 OpenSwathAnalyzer_9_output.featureXML.tmp.featureXML -in2 ${DATA_DIR_TOPP}/OpenSwathAnalyzer_9_output.featureXML)
  set_tests_properties("TOPP_OpenSwathAnalyzer_test_9_out1" PROPERTIES DEPENDS "TOPP_OpenSwathAnalyzer_test_9")
  add_test("TOPP_OpenSwathAnalyzer_test_10" ${TOPP_BIN_PATH}/OpenSwathAnalyzer -in ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input_chrom.mzML -tr ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input.TraML -out OpenSwathAnalyzer_10_output.featureXML.tmp.featureXML -rt_norm ${DATA_DIR_TOPP}/OpenSwathAnalyzer_input.trafoXML -algorithm:TransitionGroupPicker:PeakPickerChromatogram:peak_width 40.0 -algorithm:TransitionGroupPicker:PeakPickerChromatogram:method legacy -algorithm:Scores:use_mi_score -algorithm:Scores:use_total_mi_score -test)
  add_test("TOPP_OpenSwathAnalyzer_test_10_out1" ${DIFF} -whitelist "id="  -in1 OpenSwathAnalyzer_10_output.featureXML.tmp.featureXML -in2 ${DATA_DIR_TOPP}/OpenSwathAnalyzer_10_output.featureXML)
  set_tests_properties("TOPP_OpenSwathAnalyzer_test_10_out1" PROPERTIES DEPENDS "TOPP_OpenSwathAnalyzer_test_10")
  add_test("TOPP_OpenSwathAnalyzer_test_11" ${TOPP_BIN_PATH}/OpenSwathAnalyzer -in ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input_chrom.mzML -tr ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input.TraML -out OpenSwathAnalyzer_11_output.featureXML.tmp.featureXML -rt_norm ${DATA_DIR_TOPP}/OpenSwathAnalyzer_input.trafoXML -algorithm:TransitionGroupPicker:PeakPickerChromatogram:peak_width 40.0 -algorithm:TransitionGroupPicker:PeakPickerChromatogram:method legacy -algorithm:Scores:use_total_mi_score -test)
  add_test("TOPP_OpenSwathAnalyzer_test_11_out1" ${DIFF} -whitelist "id="  -in1 OpenSwathAnalyzer_11_output.featureXML.tmp.featureXML -in2 ${DATA_DIR_TOPP}/OpenSwathAnalyzer_11_output.featureXML)
  set_tests_properties("TOPP_OpenSwathAnalyzer_test_11_out1" PROPERTIES DEPENDS "TOPP_OpenSwathAnalyzer_test_11")

  add_test("TOPP_OpenSwathRTNormalizer_test_1" ${TOPP_BIN_PATH}/OpenSwathRTNormalizer -in ${DATA_DIR_TOPP}/OpenSwathRTNormalizer_1_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathRTNormalizer_1_input.TraML -out OpenSwathRTNormalizer_1_output.tmp.trafoXML -test)
  add_test("TOPP_OpenSwathRTNormalizer_test_1_out1" ${DIFF} -in1 OpenSwathRTNormalizer_1_output.tmp.trafoXML -in2 ${DATA_DIR_TOPP}/OpenSwathRTNormalizer_1_output.trafoXML)
  set_tests_properties("TOPP_OpenSwathRTNormalizer_test_1_out1" PROPERTIES DEPENDS "TOPP_OpenSwathRTNormalizer_test_1")
  add_test("TOPP_OpenSwathRTNormalizer_test_2" ${TOPP_BIN_PATH}/OpenSwathRTNormalizer -in ${DATA_DIR_TOPP}/OpenSwathRTNormalizer_1_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathRTNormalizer_1_input.TraML -out OpenSwathRTNormalizer_1_output.tmp.trafoXML -test -estimateBestPeptides -peptideEstimation:NrRTBins 30 -peptideEstimation:MinBinsFilled 30)
  set_tests_properties("TOPP_OpenSwathRTNormalizer_test_2" PROPERTIES WILL_FAIL 1) ## will not reach bin coverage
  add_test("TOPP_OpenSwathRTNormalizer_test_3" ${TOPP_BIN_PATH}/OpenSwathRTNormalizer -in ${DATA_DIR_TOPP}/OpenSwathRTNormalizer_1_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathRTNormalizer_1_input.TraML -out OpenSwathRTNormalizer_3_output.tmp.trafoXML -test -estimateBestPeptides -peptideEstimation:NrRTBins 3 -peptideEstimation:MinBinsFilled 3)
  add_test("TOPP_OpenSwathRTNormalizer_test_3_out1" ${DIFF} -in1 OpenSwathRTNormalizer_3_output.tmp.trafoXML -in2 ${DATA_DIR_TOPP}/OpenSwathRTNormalizer_3_output.trafoXML)
  set_tests_properties("TOPP_OpenSwathRTNormalizer_test_3_out1" PROPERTIES DEPENDS "TOPP_OpenSwathRTNormalizer_test_3")
  add_test("TOPP_OpenSwathRTNormalizer_test_4" ${TOPP_BIN_PATH}/OpenSwathRTNormalizer -in ${DATA_DIR_TOPP}/OpenSwathRTNormalizer_1_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathRTNormalizer_1_input.TraML -out OpenSwathRTNormalizer_4_output.tmp.trafoXML -test -RTNormalization:outlierMethod iter_residual)
  add_test("TOPP_OpenSwathRTNormalizer_test_4_out1" ${DIFF} -in1 OpenSwathRTNormalizer_4_output.tmp.trafoXML -in2 ${DATA_DIR_TOPP}/OpenSwathRTNormalizer_4_output.trafoXML)
  set_tests_properties("TOPP_OpenSwathRTNormalizer_test_4_out1" PROPERTIES DEPENDS "TOPP_OpenSwathRTNormalizer_test_4")

  add_test("TOPP_OpenSwathConfidenceScoring_1" ${TOPP_BIN_PATH}/OpenSwathConfidenceScoring -test -in ${DATA_DIR_TOPP}/OpenSwathFeatureXMLToTSV_input.featureXML -lib ${DATA_DIR_TOPP}/OpenSwathFeatureXMLToTSV_input.TraML -trafo ${DATA_DIR_TOPP}/OpenSwathConfidenceScoring_1_input.trafoXML -transitions 2 -decoys 1 -out OpenSwathConfidenceScoring_1_output.tmp.featureXML)
  add_test("TOPP_OpenSwathConfidenceScoring_1_out1" ${DIFF} -in1 OpenSwathConfidenceScoring_1_output.tmp.featureXML -in2 ${DATA_DIR_TOPP}/OpenSwathConfidenceScoring_1_output.featureXML)
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

  # Test conversion to sqMass with lossy conversion
  add_test("TOPP_OpenSwathMzMLFileCacher_test_3_step1" ${TOPP_BIN_PATH}/OpenSwathMzMLFileCacher -in ${DATA_DIR_TOPP}/OpenSwathMzMLFileCacher_1_input.mzML -out OpenSwathMzMLFileCacher_3_input.tmp.sqMass -test -lossy_compression true -lossy_mass_accuracy 1e-4)
  add_test("TOPP_OpenSwathMzMLFileCacher_test_3_step2" ${TOPP_BIN_PATH}/OpenSwathMzMLFileCacher -in OpenSwathMzMLFileCacher_3_input.tmp.sqMass -out OpenSwathMzMLFileCacher_3_output.tmp.mzML -test)
  add_test("TOPP_OpenSwathMzMLFileCacher_test_3_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 OpenSwathMzMLFileCacher_3_output.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathMzMLFileCacher_3_output.mzML)
  set_tests_properties("TOPP_OpenSwathMzMLFileCacher_test_3_step2" PROPERTIES DEPENDS "TOPP_OpenSwathMzMLFileCacher_test_3_step1")
  set_tests_properties("TOPP_OpenSwathMzMLFileCacher_test_3_out1" PROPERTIES DEPENDS "TOPP_OpenSwathMzMLFileCacher_test_3_step2")

  # Test conversion to sqMass without lossy conversion
  add_test("TOPP_OpenSwathMzMLFileCacher_test_4_step1" ${TOPP_BIN_PATH}/OpenSwathMzMLFileCacher -in ${DATA_DIR_TOPP}/OpenSwathMzMLFileCacher_2_input.chrom.mzML -out OpenSwathMzMLFileCacher_4_input.tmp.sqMass -test)
  add_test("TOPP_OpenSwathMzMLFileCacher_test_4_step2" ${TOPP_BIN_PATH}/OpenSwathMzMLFileCacher -in OpenSwathMzMLFileCacher_4_input.tmp.sqMass -out OpenSwathMzMLFileCacher_4_output.tmp.mzML -test)
  add_test("TOPP_OpenSwathMzMLFileCacher_test_4_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 OpenSwathMzMLFileCacher_4_output.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathMzMLFileCacher_4_output.mzML)
  set_tests_properties("TOPP_OpenSwathMzMLFileCacher_test_4_step2" PROPERTIES DEPENDS "TOPP_OpenSwathMzMLFileCacher_test_4_step1")
  set_tests_properties("TOPP_OpenSwathMzMLFileCacher_test_4_out1" PROPERTIES DEPENDS "TOPP_OpenSwathMzMLFileCacher_test_4_step2")

  # Test with low-memory option
  add_test("TOPP_OpenSwathMzMLFileCacher_test_5_step1" ${TOPP_BIN_PATH}/OpenSwathMzMLFileCacher -in ${DATA_DIR_TOPP}/OpenSwathMzMLFileCacher_1_input.mzML -out OpenSwathMzMLFileCacher_1_input.cached.tmp.mzML -test -process_lowmemory)
  add_test("TOPP_OpenSwathMzMLFileCacher_test_5_step2" ${TOPP_BIN_PATH}/OpenSwathMzMLFileCacher -in OpenSwathMzMLFileCacher_1_input.cached.tmp.mzML -out OpenSwathMzMLFileCacher_1_output.tmp.mzML -convert_back -test)
  add_test("TOPP_OpenSwathMzMLFileCacher_test_5_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 OpenSwathMzMLFileCacher_1_input.cached.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathMzMLFileCacher_1_output_1.mzML)
  add_test("TOPP_OpenSwathMzMLFileCacher_test_5_out2" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 OpenSwathMzMLFileCacher_1_output.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathMzMLFileCacher_1_output_2.mzML)
  set_tests_properties("TOPP_OpenSwathMzMLFileCacher_test_5_step2" PROPERTIES DEPENDS "TOPP_OpenSwathMzMLFileCacher_test_5_step1")
  set_tests_properties("TOPP_OpenSwathMzMLFileCacher_test_5_out1" PROPERTIES DEPENDS "TOPP_OpenSwathMzMLFileCacher_test_5_step1")
  set_tests_properties("TOPP_OpenSwathMzMLFileCacher_test_5_out2" PROPERTIES DEPENDS "TOPP_OpenSwathMzMLFileCacher_test_5_step2")

  # Testing the OpenSwathAnalyzer together with the FileCacher (for spectra and chromatograms)
  add_test("TOPP_OpenSwathAnalyzer_test_3_prepare" ${TOPP_BIN_PATH}/OpenSwathMzMLFileCacher -in ${DATA_DIR_TOPP}/OpenSwathAnalyzer_2_swathfile.mzML -out OpenSwathAnalyzer_3_swathfile.mzML.cached.tmp.mzML -out_type mzML -test)
  add_test("TOPP_OpenSwathAnalyzer_test_3" ${TOPP_BIN_PATH}/OpenSwathAnalyzer -in ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input_chrom.mzML -tr ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input.TraML -out MRMFeatureFinderScore_output_3.featureXML.tmp.featureXML -swath_files OpenSwathAnalyzer_3_swathfile.mzML.cached.tmp.mzML -algorithm:TransitionGroupPicker:PeakPickerChromatogram:peak_width 40.0 -algorithm:TransitionGroupPicker:PeakPickerChromatogram:method legacy -test)
  add_test("TOPP_OpenSwathAnalyzer_test_3_out1" ${DIFF} -whitelist "id=" -in1 MRMFeatureFinderScore_output_3.featureXML.tmp.featureXML -in2 ${DATA_DIR_TOPP}/OpenSwathAnalyzer_2_output.featureXML)
  set_tests_properties("TOPP_OpenSwathAnalyzer_test_3" PROPERTIES DEPENDS "TOPP_OpenSwathAnalyzer_test_3_prepare")
  set_tests_properties("TOPP_OpenSwathAnalyzer_test_3_out1" PROPERTIES DEPENDS "TOPP_OpenSwathAnalyzer_test_3")
  add_test("TOPP_OpenSwathAnalyzer_test_4_prepare1" ${TOPP_BIN_PATH}/OpenSwathMzMLFileCacher -in ${DATA_DIR_TOPP}/OpenSwathAnalyzer_2_swathfile.mzML -out OpenSwathAnalyzer_4_swathfile.mzML.cached.tmp.mzML -out_type mzML -test)
  add_test("TOPP_OpenSwathAnalyzer_test_4_prepare2" ${TOPP_BIN_PATH}/OpenSwathMzMLFileCacher -in ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input_chrom.mzML -out OpenSwathAnalyzer_4_input_chrom.mzML.cached.tmp.mzML -out_type mzML -test)
  add_test("TOPP_OpenSwathAnalyzer_test_4" ${TOPP_BIN_PATH}/OpenSwathAnalyzer -in OpenSwathAnalyzer_4_input_chrom.mzML.cached.tmp.mzML -tr ${DATA_DIR_TOPP}/OpenSwathAnalyzer_1_input.TraML -out MRMFeatureFinderScore_output_4.featureXML.tmp.featureXML -swath_files OpenSwathAnalyzer_4_swathfile.mzML.cached.tmp.mzML -algorithm:TransitionGroupPicker:PeakPickerChromatogram:peak_width 40.0 -algorithm:TransitionGroupPicker:PeakPickerChromatogram:method legacy -test)
  add_test("TOPP_OpenSwathAnalyzer_test_4_out1" ${DIFF} -whitelist "id=" -in1 MRMFeatureFinderScore_output_4.featureXML.tmp.featureXML -in2 ${DATA_DIR_TOPP}/OpenSwathAnalyzer_2_output.featureXML)
  set_tests_properties("TOPP_OpenSwathAnalyzer_test_4" PROPERTIES DEPENDS "TOPP_OpenSwathAnalyzer_test_4_prepare1;TOPP_OpenSwathAnalyzer_test_4_prepare2")
  set_tests_properties("TOPP_OpenSwathAnalyzer_test_4_out1" PROPERTIES DEPENDS "TOPP_OpenSwathAnalyzer_test_4")

  add_test("TOPP_MRMTransitionGroupPicker_test_1" ${TOPP_BIN_PATH}/MRMTransitionGroupPicker -in ${DATA_DIR_TOPP}/MRMTransitionGroupPicker_1_input.mzML -tr ${DATA_DIR_TOPP}/MRMTransitionGroupPicker_1_input.TraML -out MRMTransitionGroupPicker_test_1.featureXML.tmp.featureXML -test -algorithm:PeakPickerChromatogram:remove_overlapping_peaks true -algorithm:PeakPickerChromatogram:method legacy -algorithm:PeakPickerChromatogram:peak_width 40.0)
  add_test("TOPP_MRMTransitionGroupPicker_test_1_out" ${DIFF} -whitelist "id=" -in1 MRMTransitionGroupPicker_test_1.featureXML.tmp.featureXML -in2 ${DATA_DIR_TOPP}/MRMTransitionGroupPicker_1_output.featureXML)
  set_tests_properties("TOPP_MRMTransitionGroupPicker_test_1_out" PROPERTIES DEPENDS "TOPP_MRMTransitionGroupPicker_test_1")

  # test MI scoring
  add_test("TOPP_MRMTransitionGroupPicker_test_2" ${TOPP_BIN_PATH}/MRMTransitionGroupPicker -in ${DATA_DIR_TOPP}/MRMTransitionGroupPicker_1_input.mzML -tr ${DATA_DIR_TOPP}/MRMTransitionGroupPicker_1_input.TraML -out MRMTransitionGroupPicker_test_2.featureXML.tmp.featureXML -test -algorithm:PeakPickerChromatogram:remove_overlapping_peaks true -algorithm:PeakPickerChromatogram:method legacy -algorithm:PeakPickerChromatogram:peak_width 40.0 -algorithm:compute_total_mi)
  add_test("TOPP_MRMTransitionGroupPicker_test_2_out" ${DIFF} -whitelist "id=" -in1 MRMTransitionGroupPicker_test_2.featureXML.tmp.featureXML -in2 ${DATA_DIR_TOPP}/MRMTransitionGroupPicker_2_output.featureXML)
  set_tests_properties("TOPP_MRMTransitionGroupPicker_test_2_out" PROPERTIES DEPENDS "TOPP_MRMTransitionGroupPicker_test_2")

  #------------------------------------------------------------------------------
  # OpenSwath Workflow tests
  # NOTE: OpenMS will use the *first* occurrence of any parameter
  set(OLD_OSW_PARAM "-test" "-mz_extraction_window" "0.05" "-mz_extraction_window_unit" "Th" "-ms1_isotopes" "0"
     "-Scoring:TransitionGroupPicker:compute_peak_quality"
     "-Scoring:Scores:use_ms1_mi" "false"
     "-Scoring:Scores:use_mi_score" "false")
  add_test("TOPP_OpenSwathWorkflow_1" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.TraML -rt_norm ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.trafoXML -out_chrom OpenSwathWorkflow_1.chrom.tmp.mzML -out_features OpenSwathWorkflow_1.tmp.featureXML -out_qc OpenSwathWorkflow_1.tmp.json
  -enable_ms1 "false" ${OLD_OSW_PARAM} )
  add_test("TOPP_OpenSwathWorkflow_1_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_1.tmp.featureXML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_output.featureXML)
  add_test("TOPP_OpenSwathWorkflow_1_out2" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_1.chrom.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_output.chrom.mzML)
  add_test("TOPP_OpenSwathWorkflow_1_out3" ${DIFF} -in1 OpenSwathWorkflow_1.tmp.json -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_output.json)
  set_tests_properties("TOPP_OpenSwathWorkflow_1_out1" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_1")
  set_tests_properties("TOPP_OpenSwathWorkflow_1_out2" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_1")
  set_tests_properties("TOPP_OpenSwathWorkflow_1_out3" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_1")

  add_test("TOPP_OpenSwathWorkflow_2" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_2_input.mzXML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_2_input.TraML -rt_norm ${DATA_DIR_TOPP}/OpenSwathWorkflow_2_input.trafoXML -out_chrom OpenSwathWorkflow_2.chrom.tmp.mzML -out_features OpenSwathWorkflow_2.tmp.featureXML
  -enable_ms1 "false" ${OLD_OSW_PARAM} )
  add_test("TOPP_OpenSwathWorkflow_2_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_2.tmp.featureXML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_2_output.featureXML)
  add_test("TOPP_OpenSwathWorkflow_2_out2" ${DIFF} -in1 OpenSwathWorkflow_2.chrom.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_2_output.chrom.mzML)
  set_tests_properties("TOPP_OpenSwathWorkflow_2_out1" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_2")
  set_tests_properties("TOPP_OpenSwathWorkflow_2_out2" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_2")

  # Also test with ms1 traces (precursor ion trace)
  add_test("TOPP_OpenSwathWorkflow_3" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.TraML -rt_norm ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.trafoXML -out_chrom OpenSwathWorkflow_3.chrom.tmp.mzML -out_features OpenSwathWorkflow_3.tmp.featureXML
    ${OLD_OSW_PARAM})
  add_test("TOPP_OpenSwathWorkflow_3_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_3.tmp.featureXML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_3_output.featureXML)
  add_test("TOPP_OpenSwathWorkflow_3_out2" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_3.chrom.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_3_output.chrom.mzML)
  set_tests_properties("TOPP_OpenSwathWorkflow_3_out1" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_3")
  set_tests_properties("TOPP_OpenSwathWorkflow_3_out2" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_3")

  # Also test whether it is able to write csv output
  add_test("TOPP_OpenSwathWorkflow_4" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.TraML -rt_norm ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.trafoXML -out_chrom OpenSwathWorkflow_4.chrom.tmp.mzML -out_tsv OpenSwathWorkflow_4.tmp.tsv
    ${OLD_OSW_PARAM})
  add_test("TOPP_OpenSwathWorkflow_4_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_4.chrom.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_3_output.chrom.mzML)
  set_tests_properties("TOPP_OpenSwathWorkflow_4_out1" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_4")
  # We cannot currently test the correctness of the csv output

  # Also test with readoptions cache
  add_test("TOPP_OpenSwathWorkflow_5" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.TraML -rt_norm ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.trafoXML -out_chrom OpenSwathWorkflow_5.chrom.tmp.mzML -out_features OpenSwathWorkflow_5.tmp.featureXML
  -readOptions cache -tempDirectory "." ${OLD_OSW_PARAM})
  add_test("TOPP_OpenSwathWorkflow_5_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_5.tmp.featureXML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_3_output.featureXML)
  add_test("TOPP_OpenSwathWorkflow_5_out2" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_5.chrom.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_3_output.chrom.mzML)
  set_tests_properties("TOPP_OpenSwathWorkflow_5_out1" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_5")
  set_tests_properties("TOPP_OpenSwathWorkflow_5_out2" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_5")

  # Also test with readoptions cacheWorkingInMemory
  add_test("TOPP_OpenSwathWorkflow_6" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.TraML -rt_norm ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.trafoXML -out_chrom OpenSwathWorkflow_6.chrom.tmp.mzML -out_features OpenSwathWorkflow_6.tmp.featureXML -readOptions cacheWorkingInMemory -tempDirectory "."
  ${OLD_OSW_PARAM} )
  add_test("TOPP_OpenSwathWorkflow_6_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_6.tmp.featureXML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_3_output.featureXML)
  add_test("TOPP_OpenSwathWorkflow_6_out2" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_6.chrom.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_3_output.chrom.mzML)
  set_tests_properties("TOPP_OpenSwathWorkflow_6_out1" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_6")
  set_tests_properties("TOPP_OpenSwathWorkflow_6_out2" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_6")

  # Test with input swath windows file
  add_test("TOPP_OpenSwathWorkflow_7" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.TraML -rt_norm ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.trafoXML -out_chrom OpenSwathWorkflow_7.chrom.tmp.mzML -out_features OpenSwathWorkflow_7.tmp.featureXML -swath_windows_file ${DATA_DIR_TOPP}/swath_windows.txt
  ${OLD_OSW_PARAM} )
  add_test("TOPP_OpenSwathWorkflow_7_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_7.tmp.featureXML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_3_output.featureXML)
  add_test("TOPP_OpenSwathWorkflow_7_out2" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_7.chrom.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_3_output.chrom.mzML)
  set_tests_properties("TOPP_OpenSwathWorkflow_7_out1" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_7")
  set_tests_properties("TOPP_OpenSwathWorkflow_7_out2" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_7")

  # Test with faulty swath windows file
  add_test("TOPP_OpenSwathWorkflow_8" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.TraML -rt_norm ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.trafoXML -out_chrom OpenSwathWorkflow_8.chrom.tmp.mzML -out_features OpenSwathWorkflow_8.tmp.featureXML -swath_windows_file ${DATA_DIR_TOPP}/swath_windows_overlap.txt
  ${OLD_OSW_PARAM} )
  set_tests_properties("TOPP_OpenSwathWorkflow_8" PROPERTIES WILL_FAIL 1) # file contains overlaps
  add_test("TOPP_OpenSwathWorkflow_9" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.TraML -rt_norm ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.trafoXML -out_chrom OpenSwathWorkflow_9.chrom.tmp.mzML -out_features OpenSwathWorkflow_9.tmp.featureXML -swath_windows_file ${DATA_DIR_TOPP}/swath_windows_tooshort.txt
  ${OLD_OSW_PARAM} )
  set_tests_properties("TOPP_OpenSwathWorkflow_9" PROPERTIES WILL_FAIL 1) # file does not have enough entries
  # works with force (does not check for overlaps and missing swath windows any more)
  add_test("TOPP_OpenSwathWorkflow_10" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.TraML -rt_norm ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.trafoXML -out_chrom OpenSwathWorkflow_10.chrom.tmp.mzML -out_features OpenSwathWorkflow_10.tmp.featureXML -swath_windows_file ${DATA_DIR_TOPP}/swath_windows_overlap.txt -force
  ${OLD_OSW_PARAM} )

  # test SONAR (also test iRT correction ... )
  # - use iRT correction (but no outlier removal since only 2 datapoints)
  # - use quadratic fitting for m/z correction
  # - use 550 ppm for extraction (0.05 Da @ 100 m/z)
  add_test("TOPP_OpenSwathWorkflow_11" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_11_input.mzML -tr_irt ${DATA_DIR_TOPP}/OpenSwathWorkflow_11_input.TraML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_11_input.TraML -out_chrom OpenSwathWorkflow_11.chrom.tmp.mzML -out_features OpenSwathWorkflow_11.tmp.featureXML
    -mz_extraction_window 0.2 -rt_extraction_window -1 -Scoring:Scores:use_sonar_scores -sonar -RTNormalization:outlierMethod none -mz_correction_function quadratic_regression_delta_ppm -irt_mz_extraction_window 550 -irt_mz_extraction_window_unit ppm
  -enable_ms1 "false" ${OLD_OSW_PARAM} )
  add_test("TOPP_OpenSwathWorkflow_11_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_11.chrom.tmp.mzML  -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_11_output.chrom.mzML)
  add_test("TOPP_OpenSwathWorkflow_11_out2" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_11.tmp.featureXML  -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_11_output.featureXML)
  set_tests_properties("TOPP_OpenSwathWorkflow_11_out1" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_11")
  set_tests_properties("TOPP_OpenSwathWorkflow_11_out2" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_11")

  # 20 ppm is too small for extraction and no iRT peptides are found (thus the tool fails)
  add_test("TOPP_OpenSwathWorkflow_12" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_11_input.mzML -tr_irt ${DATA_DIR_TOPP}/OpenSwathWorkflow_11_input.TraML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_11_input.TraML -mz_extraction_window 0.2 -rt_extraction_window -1 -Scoring:Scores:use_sonar_scores -sonar -out_chrom OpenSwathWorkflow_11.chrom.tmp.mzML -out_features OpenSwathWorkflow_11.tmp.featureXML -RTNormalization:outlierMethod none -mz_correction_function quadratic_regression_delta_ppm -irt_mz_extraction_window 20 -irt_mz_extraction_window_unit ppm -test)
  set_tests_properties("TOPP_OpenSwathWorkflow_12" PROPERTIES WILL_FAIL 1) # there are no points for iRT regression

  # Test whether PQP files can be read and OSW output can be written
  add_test("TOPP_OpenSwathWorkflow_13_prepare" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.TraML -out OpenSwathWorkflow_13_input.tmp.pqp -out_type pqp)
  add_test("TOPP_OpenSwathWorkflow_13" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.mzML -tr OpenSwathWorkflow_13_input.tmp.pqp -tr_type pqp -rt_norm ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.trafoXML -out_chrom OpenSwathWorkflow_13.chrom.tmp.mzML -out_osw OpenSwathWorkflow_13.osw
  ${OLD_OSW_PARAM} )
  add_test("TOPP_OpenSwathWorkflow_13_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_13.chrom.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_13_output.chrom.mzML)
  set_tests_properties("TOPP_OpenSwathWorkflow_13" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_13_prepare")
  set_tests_properties("TOPP_OpenSwathWorkflow_13_out1" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_13")
  # We cannot currently test the correctness of the PQP and OSW files

  # Test whether we can write out sqMass files
  add_test("TOPP_OpenSwathWorkflow_14_prepare" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.TraML -out OpenSwathWorkflow_14_input.tmp.pqp -out_type pqp)
  add_test("TOPP_OpenSwathWorkflow_14" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.mzML -tr OpenSwathWorkflow_14_input.tmp.pqp -tr_type pqp -rt_norm ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.trafoXML -out_chrom OpenSwathWorkflow_14.chrom.tmp.sqMass -out_osw OpenSwathWorkflow_14.osw
  ${OLD_OSW_PARAM} )
  # add_test("TOPP_OpenSwathWorkflow_14_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_4.chrom.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_14_output.chrom.mzML)
  add_test("TOPP_OpenSwathWorkflow_14_step2" ${TOPP_BIN_PATH}/OpenSwathMzMLFileCacher -in OpenSwathWorkflow_14.chrom.tmp.sqMass -out OpenSwathWorkflow_14.chrom.tmp.mzML -test -lossy_compression false -lossy_mass_accuracy 1e-4 -full_meta false)
  add_test("TOPP_OpenSwathWorkflow_14_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_14.chrom.tmp.mzML  -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_14_output.chrom.mzML)
  set_tests_properties("TOPP_OpenSwathWorkflow_14" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_14_prepare")
  set_tests_properties("TOPP_OpenSwathWorkflow_14_step2" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_14")
  set_tests_properties("TOPP_OpenSwathWorkflow_14_out1" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_14_step2")
  # We cannot currently test the correctness of the PQP and OSW files

  # Test with ms1 only (precursor ion trace)
  add_test("TOPP_OpenSwathWorkflow_15" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_15_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.TraML -rt_norm ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.trafoXML -out_chrom OpenSwathWorkflow_15.chrom.tmp.mzML -out_features OpenSwathWorkflow_15.tmp.featureXML
    -Scoring:TransitionGroupPicker:use_precursors -ms1_isotopes 2
  "-test" "-mz_extraction_window" "0.05" "-mz_extraction_window_unit" "Th" # "-ms1_isotopes" "0"
  #   "-Scoring:TransitionGroupPicker:compute_peak_quality"
     "-Scoring:Scores:use_ms1_mi" "false"
     "-Scoring:Scores:use_mi_score" "false")
  add_test("TOPP_OpenSwathWorkflow_15_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_15.tmp.featureXML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_15_output.featureXML)
  add_test("TOPP_OpenSwathWorkflow_15_out2" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_15.chrom.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_15_output.chrom.mzML)
  set_tests_properties("TOPP_OpenSwathWorkflow_15_out1" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_15")
  set_tests_properties("TOPP_OpenSwathWorkflow_15_out2" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_15")

  # Also test with ms1 traces (precursor ion trace)
  add_test("TOPP_OpenSwathWorkflow_16_prepare" ${TOPP_BIN_PATH}/OpenSwathMzMLFileCacher -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.mzML -out OpenSwathWorkflow_16_input.sqMass -test -lossy_compression true -lossy_mass_accuracy 1e-4 -full_meta true)
  add_test("TOPP_OpenSwathWorkflow_16" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in OpenSwathWorkflow_16_input.sqMass -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.TraML -rt_norm ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.trafoXML -out_chrom OpenSwathWorkflow_16.chrom.tmp.mzML -out_features OpenSwathWorkflow_16.tmp.featureXML -readOptions workingInMemory -ms1_isotopes 2
  ${OLD_OSW_PARAM} )
  add_test("TOPP_OpenSwathWorkflow_16_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_16.tmp.featureXML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_16_output.featureXML)
  add_test("TOPP_OpenSwathWorkflow_16_out2" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_16.chrom.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_16_output.chrom.mzML)
  set_tests_properties("TOPP_OpenSwathWorkflow_16" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_16_prepare")
  set_tests_properties("TOPP_OpenSwathWorkflow_16_out1" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_16")
  set_tests_properties("TOPP_OpenSwathWorkflow_16_out2" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_16")

  # Test with ion mobility
  add_test("TOPP_OpenSwathWorkflow_17" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_17_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_17_input.tsv -out_chrom OpenSwathWorkflow_17.chrom.tmp.mzML -out_features OpenSwathWorkflow_17.tmp.featureXML
    -readOptions workingInMemory -ion_mobility_window 0.05 -use_ms1_ion_mobility "false" -Scoring:Scores:use_ion_mobility_scores
    ${OLD_OSW_PARAM} )
  add_test("TOPP_OpenSwathWorkflow_17_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_17.tmp.featureXML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_17_output.featureXML)
  add_test("TOPP_OpenSwathWorkflow_17_out2" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_17.chrom.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_17_output.chrom.mzML)
  set_tests_properties("TOPP_OpenSwathWorkflow_17_out1" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_17")
  set_tests_properties("TOPP_OpenSwathWorkflow_17_out2" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_17")
  # Test with ion mobility
  add_test("TOPP_OpenSwathWorkflow_17_cache" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_17_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_17_input.tsv -out_chrom OpenSwathWorkflow_17.chrom.tmp.mzML -out_features OpenSwathWorkflow_17.tmp.featureXML -readOptions cache -ion_mobility_window 0.05 -use_ms1_ion_mobility "false" -Scoring:Scores:use_ion_mobility_scores
    ${OLD_OSW_PARAM} )
  add_test("TOPP_OpenSwathWorkflow_17_cache_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_17.tmp.featureXML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_17_output.featureXML)
  add_test("TOPP_OpenSwathWorkflow_17_cache_out2" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_17.chrom.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_17_output.chrom.mzML)
  set_tests_properties("TOPP_OpenSwathWorkflow_17_cache_out1" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_17_cache")
  set_tests_properties("TOPP_OpenSwathWorkflow_17_cache_out2" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_17_cache")

  # Test with ion mobility and the SQLite workflow
  add_test("TOPP_OpenSwathWorkflow_17_b_prepare" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_17_input.tsv -out OpenSwathWorkflow_17_input.tmp.pqp -out_type pqp)
  add_test("TOPP_OpenSwathWorkflow_17_b" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_17_input.mzML -tr_type pqp -tr OpenSwathWorkflow_17_input.tmp.pqp -out_chrom OpenSwathWorkflow_17_b.chrom.tmp.mzML -out_features OpenSwathWorkflow_17_b.tmp.featureXML -readOptions workingInMemory -ion_mobility_window 0.05 -use_ms1_ion_mobility "false"
    ${OLD_OSW_PARAM} )
  add_test("TOPP_OpenSwathWorkflow_17_b_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_17_b.tmp.featureXML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_17_b_output.featureXML)
  add_test("TOPP_OpenSwathWorkflow_17_b_out2" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_17_b.chrom.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_17_b_output.chrom.mzML)
  set_tests_properties("TOPP_OpenSwathWorkflow_17_b" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_17_b_prepare")
  set_tests_properties("TOPP_OpenSwathWorkflow_17_b_out1" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_17_b")
  set_tests_properties("TOPP_OpenSwathWorkflow_17_b_out2" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_17_b")

  add_test("TOPP_OpenSwathFileSplitter_1" ${TOPP_BIN_PATH}/OpenSwathFileSplitter -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.mzML -outputDirectory ""  -out_qc OpenSwathFileSplitter_1.tmp.json -test)
  add_test("TOPP_OpenSwathFileSplitter_1_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_1_input_ms1.mzML -in2 ${DATA_DIR_TOPP}/openswath_tmpfile_ms1.mzML)
  add_test("TOPP_OpenSwathFileSplitter_1_out2" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_1_input_4.mzML -in2 ${DATA_DIR_TOPP}/openswath_tmpfile_4.mzML)
  add_test("TOPP_OpenSwathFileSplitter_1_out3" ${DIFF} -in1 OpenSwathFileSplitter_1.tmp.json -in2 ${DATA_DIR_TOPP}/OpenSwathFileSplitter_1.json)
  set_tests_properties("TOPP_OpenSwathFileSplitter_1_out1" PROPERTIES DEPENDS "TOPP_OpenSwathFileSplitter_1")
  set_tests_properties("TOPP_OpenSwathFileSplitter_1_out2" PROPERTIES DEPENDS "TOPP_OpenSwathFileSplitter_1")
  set_tests_properties("TOPP_OpenSwathFileSplitter_1_out3" PROPERTIES DEPENDS "TOPP_OpenSwathFileSplitter_1")

  # Test MS1 MI scoring
  add_test("TOPP_OpenSwathWorkflow_18" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.TraML -rt_norm ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.trafoXML -out_features OpenSwathWorkflow_18.tmp.featureXML
  "-test" "-mz_extraction_window" "0.05" "-mz_extraction_window_unit" "Th" "-ms1_isotopes" "0"
     "-Scoring:TransitionGroupPicker:compute_peak_quality"
     # "-Scoring:Scores:use_ms1_mi" "false"
     "-Scoring:Scores:use_mi_score" "false")
  add_test("TOPP_OpenSwathWorkflow_18_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_18.tmp.featureXML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_18_output.featureXML)
  set_tests_properties("TOPP_OpenSwathWorkflow_18_out1" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_18")

  # Test MS2 MI scoring
  add_test("TOPP_OpenSwathWorkflow_19" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.TraML -rt_norm ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.trafoXML -out_features OpenSwathWorkflow_19.tmp.featureXML
  "-test" "-mz_extraction_window" "0.05" "-mz_extraction_window_unit" "Th" "-ms1_isotopes" "0"
     "-Scoring:TransitionGroupPicker:compute_peak_quality"
     "-Scoring:Scores:use_ms1_mi" "false")
     #"-Scoring:Scores:use_mi_score" "false")
  add_test("TOPP_OpenSwathWorkflow_19_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_19.tmp.featureXML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_19_output.featureXML)
  set_tests_properties("TOPP_OpenSwathWorkflow_19_out1" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_19")

  # Test MS2 MI scoring and total MI scoring
  add_test("TOPP_OpenSwathWorkflow_20" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.TraML -rt_norm ${DATA_DIR_TOPP}/OpenSwathWorkflow_1_input.trafoXML -out_features OpenSwathWorkflow_20.tmp.featureXML
  "-test" "-mz_extraction_window" "0.05" "-mz_extraction_window_unit" "Th" "-ms1_isotopes" "0"
     "-Scoring:TransitionGroupPicker:compute_peak_quality"
     "-Scoring:Scores:use_total_mi_score"
     "-Scoring:Scores:use_ms1_mi" "false"
  )
  add_test("TOPP_OpenSwathWorkflow_20_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_20.tmp.featureXML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_20_output.featureXML)
  set_tests_properties("TOPP_OpenSwathWorkflow_20_out1" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_20")

  # Test iRT alignment, full workflow
  add_test("TOPP_OpenSwathWorkflow_21" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_21_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_21_input.tsv -tr_irt ${DATA_DIR_TOPP}/OpenSwathWorkflow_21_input.irt.TraML -out_features OpenSwathWorkflow_21.tmp.featureXML -Debugging:irt_trafo OpenSwathWorkflow_21.tmp.trafoXML -out_chrom OpenSwathWorkflow_21.tmp.mzML
    -RTNormalization:lowess:span 0.666666666666666666666666666
    "-test" "-mz_extraction_window" "0.05" "-mz_extraction_window_unit" "Th" "-ms1_isotopes" "0"
    "-Scoring:Scores:use_total_mi_score"
    "-Scoring:Scores:use_ms1_mi" "false"
    "-Scoring:TransitionGroupPicker:compute_peak_quality")
  add_test("TOPP_OpenSwathWorkflow_21_out1" ${DIFF} -in1 OpenSwathWorkflow_21.tmp.featureXML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_21_output.featureXML)
  add_test("TOPP_OpenSwathWorkflow_21_out2" ${DIFF} -in1 OpenSwathWorkflow_21.tmp.trafoXML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_21_output.trafoXML)
  # add_test("TOPP_OpenSwathWorkflow_21_out3" ${DIFF} -in1 OpenSwathWorkflow_21.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_21_output.mzML)
  set_tests_properties("TOPP_OpenSwathWorkflow_21_out1" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_21")
  set_tests_properties("TOPP_OpenSwathWorkflow_21_out2" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_21")
  # set_tests_properties("TOPP_OpenSwathWorkflow_21_out3" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_21")

  # Test PRM workflow
  add_test("TOPP_OpenSwathWorkflow_22" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_22_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_22_input.tsv -out_chrom OpenSwathWorkflow_22.chrom.tmp.mzML -out_features OpenSwathWorkflow_22.tmp.featureXML -readOptions workingInMemory -matching_window_only "true" -ms1_isotopes 3
  ${OLD_OSW_PARAM} )
  add_test("TOPP_OpenSwathWorkflow_22_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_22.tmp.featureXML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_22_output.featureXML)
  add_test("TOPP_OpenSwathWorkflow_22_out2" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_22.chrom.tmp.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_22_output.chrom.mzML)
  set_tests_properties("TOPP_OpenSwathWorkflow_22_out1" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_22")
  set_tests_properties("TOPP_OpenSwathWorkflow_22_out2" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_22")

  # Test with PASEF ion mobility data, overlapping m/z SWATHs with distinct IM
  add_test("TOPP_OpenSwathWorkflow_23" ${TOPP_BIN_PATH}/OpenSwathWorkflow -in ${DATA_DIR_TOPP}/OpenSwathWorkflow_23_input.mzML -tr ${DATA_DIR_TOPP}/OpenSwathWorkflow_23_input.tsv -out_chrom OpenSwathWorkflow_23.tmp.chrom.mzML -out_features OpenSwathWorkflow_23.tmp.featureXML -readOptions workingInMemory -ion_mobility_window 0.1 -Scoring:Scores:use_ion_mobility_scores -pasef ${OLD_OSW_PARAM})
  add_test("TOPP_OpenSwathWorkflow_23_out1" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_23.tmp.featureXML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_23_output.featureXML)
  add_test("TOPP_OpenSwathWorkflow_23_out2" ${DIFF} -whitelist "id=" -in1 OpenSwathWorkflow_23.tmp.chrom.mzML -in2 ${DATA_DIR_TOPP}/OpenSwathWorkflow_23_output.chrom.mzML)
  set_tests_properties("TOPP_OpenSwathWorkflow_23_out1" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_23")
  set_tests_properties("TOPP_OpenSwathWorkflow_23_out2" PROPERTIES DEPENDS "TOPP_OpenSwathWorkflow_23")



endif(NOT DISABLE_OPENSWATH)

#------------------------------------------------------------------------------
# NoiseFilter tests
add_test("TOPP_NoiseFilterGaussian_1" ${TOPP_BIN_PATH}/NoiseFilterGaussian -test -ini ${DATA_DIR_TOPP}/NoiseFilterGaussian_1_parameters.ini -in ${DATA_DIR_TOPP}/NoiseFilterGaussian_1_input.mzML -out NoiseFilterGaussian_1.tmp.mzML)
add_test("TOPP_NoiseFilterGaussian_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 NoiseFilterGaussian_1.tmp.mzML -in2 ${DATA_DIR_TOPP}/NoiseFilterGaussian_1_output.mzML )
set_tests_properties("TOPP_NoiseFilterGaussian_1_out1" PROPERTIES DEPENDS "TOPP_NoiseFilterGaussian_1")
add_test("TOPP_NoiseFilterGaussian_2" ${TOPP_BIN_PATH}/NoiseFilterGaussian -test -ini ${DATA_DIR_TOPP}/NoiseFilterGaussian_2_parameters.ini -in ${DATA_DIR_TOPP}/NoiseFilterGaussian_2_input.chrom.mzML -out NoiseFilterGaussian_2.tmp.mzML)
add_test("TOPP_NoiseFilterGaussian_2_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 NoiseFilterGaussian_2.tmp.mzML -in2 ${DATA_DIR_TOPP}/NoiseFilterGaussian_2_output.chrom.mzML )
set_tests_properties("TOPP_NoiseFilterGaussian_2_out1" PROPERTIES DEPENDS "TOPP_NoiseFilterGaussian_2")

# lowmem variation of the algorithm
add_test("TOPP_NoiseFilterGaussian_3" ${TOPP_BIN_PATH}/NoiseFilterGaussian -test -ini ${DATA_DIR_TOPP}/NoiseFilterGaussian_1_parameters.ini -in ${DATA_DIR_TOPP}/NoiseFilterGaussian_1_input.mzML -out NoiseFilterGaussian_3.tmp.mzML -processOption lowmemory)
add_test("TOPP_NoiseFilterGaussian_3_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 NoiseFilterGaussian_3.tmp.mzML -in2 ${DATA_DIR_TOPP}/NoiseFilterGaussian_1_output.mzML )
set_tests_properties("TOPP_NoiseFilterGaussian_3_out1" PROPERTIES DEPENDS "TOPP_NoiseFilterGaussian_3")
add_test("TOPP_NoiseFilterGaussian_4" ${TOPP_BIN_PATH}/NoiseFilterGaussian -test -ini ${DATA_DIR_TOPP}/NoiseFilterGaussian_2_parameters.ini -in ${DATA_DIR_TOPP}/NoiseFilterGaussian_2_input.chrom.mzML -out NoiseFilterGaussian_4.tmp.mzML -processOption lowmemory)
add_test("TOPP_NoiseFilterGaussian_4_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 NoiseFilterGaussian_4.tmp.mzML -in2 ${DATA_DIR_TOPP}/NoiseFilterGaussian_2_output.chrom.mzML )
set_tests_properties("TOPP_NoiseFilterGaussian_4_out1" PROPERTIES DEPENDS "TOPP_NoiseFilterGaussian_4")

add_test("TOPP_NoiseFilterSGolay_1" ${TOPP_BIN_PATH}/NoiseFilterSGolay -test -ini ${DATA_DIR_TOPP}/NoiseFilterSGolay_1_parameters.ini -in ${DATA_DIR_TOPP}/NoiseFilterSGolay_1_input.mzML -out NoiseFilterSGolay_1.tmp.mzML)
add_test("TOPP_NoiseFilterSGolay_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 NoiseFilterSGolay_1.tmp.mzML -in2 ${DATA_DIR_TOPP}/NoiseFilterSGolay_1_output.mzML )
set_tests_properties("TOPP_NoiseFilterSGolay_1_out1" PROPERTIES DEPENDS "TOPP_NoiseFilterSGolay_1")
add_test("TOPP_NoiseFilterSGolay_2" ${TOPP_BIN_PATH}/NoiseFilterSGolay -test -ini ${DATA_DIR_TOPP}/NoiseFilterSGolay_2_parameters.ini -in ${DATA_DIR_TOPP}/NoiseFilterSGolay_2_input.chrom.mzML -out NoiseFilterSGolay_2.tmp.mzML)
add_test("TOPP_NoiseFilterSGolay_2_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 NoiseFilterSGolay_2.tmp.mzML -in2 ${DATA_DIR_TOPP}/NoiseFilterSGolay_2_output.chrom.mzML )
set_tests_properties("TOPP_NoiseFilterSGolay_2_out1" PROPERTIES DEPENDS "TOPP_NoiseFilterSGolay_2")

# lowmem variation of the algorithm
add_test("TOPP_NoiseFilterSGolay_3" ${TOPP_BIN_PATH}/NoiseFilterSGolay -test -ini ${DATA_DIR_TOPP}/NoiseFilterSGolay_1_parameters.ini -in ${DATA_DIR_TOPP}/NoiseFilterSGolay_1_input.mzML -out NoiseFilterSGolay_3.tmp.mzML -processOption lowmemory)
add_test("TOPP_NoiseFilterSGolay_3_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 NoiseFilterSGolay_3.tmp.mzML -in2 ${DATA_DIR_TOPP}/NoiseFilterSGolay_1_output.mzML )
set_tests_properties("TOPP_NoiseFilterSGolay_3_out1" PROPERTIES DEPENDS "TOPP_NoiseFilterSGolay_3")
add_test("TOPP_NoiseFilterSGolay_4" ${TOPP_BIN_PATH}/NoiseFilterSGolay -test -ini ${DATA_DIR_TOPP}/NoiseFilterSGolay_2_parameters.ini -in ${DATA_DIR_TOPP}/NoiseFilterSGolay_2_input.chrom.mzML -out NoiseFilterSGolay_4.tmp.mzML -processOption lowmemory)
add_test("TOPP_NoiseFilterSGolay_4_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 NoiseFilterSGolay_4.tmp.mzML -in2 ${DATA_DIR_TOPP}/NoiseFilterSGolay_2_output.chrom.mzML )
set_tests_properties("TOPP_NoiseFilterSGolay_4_out1" PROPERTIES DEPENDS "TOPP_NoiseFilterSGolay_4")

#------------------------------------------------------------------------------
# PeakPicker tests

# "high_res" algorithm with "ms_levels" option set to only 1:
# .. spectra
add_test("TOPP_PeakPickerHiRes_1" ${TOPP_BIN_PATH}/PeakPickerHiRes -test -ini ${DATA_DIR_TOPP}/PeakPickerHiRes_parameters.ini -in ${DATA_DIR_TOPP}/PeakPickerHiRes_input.mzML -out PeakPickerHiRes_1.tmp.mzML)
add_test("TOPP_PeakPickerHiRes_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 PeakPickerHiRes_1.tmp.mzML -in2 ${DATA_DIR_TOPP}/PeakPickerHiRes_output.mzML)
set_tests_properties("TOPP_PeakPickerHiRes_1_out1" PROPERTIES DEPENDS "TOPP_PeakPickerHiRes_1")
# .. chromatogram
add_test("TOPP_PeakPickerHiRes_2" ${TOPP_BIN_PATH}/PeakPickerHiRes -test -ini ${DATA_DIR_TOPP}/PeakPickerHiRes_parameters.ini -in ${DATA_DIR_TOPP}/PeakPickerHiRes_2_input.mzML -out PeakPickerHiRes_2.tmp.mzML)
add_test("TOPP_PeakPickerHiRes_2_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 PeakPickerHiRes_2.tmp.mzML -in2 ${DATA_DIR_TOPP}/PeakPickerHiRes_2_output.mzML)
set_tests_properties("TOPP_PeakPickerHiRes_2_out1" PROPERTIES DEPENDS "TOPP_PeakPickerHiRes_2")

# "high_res" peakpicker with the lowmemory option
# .. spectra
add_test("TOPP_PeakPickerHiRes_3" ${TOPP_BIN_PATH}/PeakPickerHiRes -test -ini ${DATA_DIR_TOPP}/PeakPickerHiRes_parameters.ini -in ${DATA_DIR_TOPP}/PeakPickerHiRes_input.mzML -out PeakPickerHiRes_3.tmp.mzML -processOption lowmemory)
## Note: the output SHOULD be identical to 'PeakPickerHiRes_output.mzML', but due to a missing 'Dataprocessing' entry (which is not known when writing the mzML header), we need an extra output file
add_test("TOPP_PeakPickerHiRes_3_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 PeakPickerHiRes_3.tmp.mzML -in2 ${DATA_DIR_TOPP}/PeakPickerHiRes_output_lowMem.mzML)
set_tests_properties("TOPP_PeakPickerHiRes_3_out1" PROPERTIES DEPENDS "TOPP_PeakPickerHiRes_3")
# .. chromatogram
add_test("TOPP_PeakPickerHiRes_4" ${TOPP_BIN_PATH}/PeakPickerHiRes -test -ini ${DATA_DIR_TOPP}/PeakPickerHiRes_parameters.ini -in ${DATA_DIR_TOPP}/PeakPickerHiRes_2_input.mzML -out PeakPickerHiRes_4.tmp.mzML -processOption lowmemory)
add_test("TOPP_PeakPickerHiRes_4_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 PeakPickerHiRes_4.tmp.mzML -in2 ${DATA_DIR_TOPP}/PeakPickerHiRes_2_output.mzML)
set_tests_properties("TOPP_PeakPickerHiRes_4_out1" PROPERTIES DEPENDS "TOPP_PeakPickerHiRes_4")

# peak picking of a smoothed chromatogram (like in PeakPickerChromatogram; requires disabling of spacing constraints):
add_test("TOPP_PeakPickerHiRes_5" ${TOPP_BIN_PATH}/PeakPickerHiRes -test -in ${DATA_DIR_TOPP}/PeakPickerHiRes_5_input.mzML -out PeakPickerHiRes_5.tmp.mzML)
add_test("TOPP_PeakPickerHiRes_5_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 PeakPickerHiRes_5.tmp.mzML -in2 ${DATA_DIR_TOPP}/PeakPickerHiRes_5_output.mzML)
set_tests_properties("TOPP_PeakPickerHiRes_5_out1" PROPERTIES DEPENDS "TOPP_PeakPickerHiRes_5")

#------------------------------------------------------------------------------
# TOPP_PeakPickerIterative
add_test("TOPP_PeakPickerIterative_1" ${TOPP_BIN_PATH}/PeakPickerIterative -in ${DATA_DIR_TOPP}/PeakPickerIterative_1_input.mzML -ini ${DATA_DIR_TOPP}/PeakPickerIterative_1.ini -out PeakPickerIterative.tmp.mzML -test)
add_test("TOPP_PeakPickerIterative_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 PeakPickerIterative.tmp.mzML -in2 ${DATA_DIR_TOPP}/PeakPickerIterative_1_output.mzML)
set_tests_properties("TOPP_PeakPickerIterative_1_out1" PROPERTIES DEPENDS "TOPP_PeakPickerIterative_1")
add_test("TOPP_PeakPickerIterative_2" ${TOPP_BIN_PATH}/PeakPickerIterative -in ${DATA_DIR_TOPP}/PeakPickerIterative_2_input.mzML -ini ${DATA_DIR_TOPP}/PeakPickerIterative_2.ini -out PeakPickerIterative_2.tmp.mzML -test)
add_test("TOPP_PeakPickerIterative_2_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 PeakPickerIterative_2.tmp.mzML -in2 ${DATA_DIR_TOPP}/PeakPickerIterative_2_output.mzML)
set_tests_properties("TOPP_PeakPickerIterative_2_out1" PROPERTIES DEPENDS "TOPP_PeakPickerIterative_2")


#------------------------------------------------------------------------------
# Basic ProteinInference
add_test("TOPP_ProteinInference_1" ${TOPP_BIN_PATH}/ProteinInference -test -in ${DATA_DIR_TOPP}/ProteinInference_1_input.idXML -out ProteinInference_1_output.tmp.idXML -Algorithm:use_shared_peptides false -merge_runs all -Merging:annotate_origin false)
add_test("TOPP_ProteinInference_1_out1" ${DIFF} -whitelist "InferenceEngineVersion" -in1 ProteinInference_1_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/ProteinInference_1_output.idXML)
set_tests_properties("TOPP_ProteinInference_1_out1" PROPERTIES DEPENDS "TOPP_ProteinInference_1")

add_test("TOPP_ProteinInference_2" ${TOPP_BIN_PATH}/ProteinInference -test -in ${DATA_DIR_TOPP}/Epifany_2_input.consensusXML -out ProteinInference_2_output.tmp.consensusXML -Algorithm:use_shared_peptides false -out_type consensusXML)
add_test("TOPP_ProteinInference_2_out1" ${DIFF} -whitelist "InferenceEngineVersion" -in1 ProteinInference_2_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/ProteinInference_2_output.consensusXML)
set_tests_properties("TOPP_ProteinInference_2_out1" PROPERTIES DEPENDS "TOPP_ProteinInference_2")

add_test("TOPP_ProteinInference_3" ${TOPP_BIN_PATH}/ProteinInference -test -in ${DATA_DIR_TOPP}/MSstatsConverter_3_in.consensusXML -out ProteinInference_3_output.tmp.consensusXML -Algorithm:use_shared_peptides true -picked_decoy_string "_rev" -picked_decoy_prefix "suffix" -protein_fdr true -out_type consensusXML)
add_test("TOPP_ProteinInference_3_out1" ${DIFF} -whitelist "InferenceEngineVersion" -in1 ProteinInference_3_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/ProteinInference_3_output.consensusXML)
set_tests_properties("TOPP_ProteinInference_3_out1" PROPERTIES DEPENDS "TOPP_ProteinInference_3")

#------------------------------------------------------------------------------
# Resampler tests
if(WITH_GUI)
  add_test("TOPP_Resampler_1" ${TOPP_BIN_PATH}/Resampler -test -in ${DATA_DIR_TOPP}/Resampler_1_input.mzML -out Resampler.mzML -sampling_rate 0.3)
  add_test("TOPP_Resampler_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 Resampler.mzML -in2 ${DATA_DIR_TOPP}/Resampler_1_output.mzML )
  set_tests_properties("TOPP_Resampler_1_out1" PROPERTIES DEPENDS "TOPP_Resampler_1")
  
  add_test("TOPP_Resampler_2" ${TOPP_BIN_PATH}/Resampler -test -in ${DATA_DIR_TOPP}/Resampler_1_input.mzML -out Resampler2_ppm.mzML -ppm -sampling_rate 2000 -min_int_cutoff 3900)
  add_test("TOPP_Resampler_2_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 Resampler2_ppm.mzML -in2 ${DATA_DIR_TOPP}/Resampler_2_output.mzML )
  set_tests_properties("TOPP_Resampler_2_out1" PROPERTIES DEPENDS "TOPP_Resampler_2")
endif()

#------------------------------------------------------------------------------
# SeedListGenerator tests
add_test("TOPP_SeedListGenerator_1" ${TOPP_BIN_PATH}/SeedListGenerator -test -in ${DATA_DIR_TOPP}/PepXMLFile_test.mzML -out_prefix SeedListGenerator_1_output_tmp)
add_test("TOPP_SeedListGenerator_1_out1" ${DIFF} -in1 SeedListGenerator_1_output_tmp_0.featureXML -in2 ${DATA_DIR_TOPP}/SeedListGenerator_1_output.featureXML )
set_tests_properties("TOPP_SeedListGenerator_1_out1" PROPERTIES DEPENDS "TOPP_SeedListGenerator_1")
add_test("TOPP_SeedListGenerator_2" ${TOPP_BIN_PATH}/SeedListGenerator -test -in ${DATA_DIR_TOPP}/IDMapper_1_output.featureXML -out_prefix SeedListGenerator_2_output_tmp)
add_test("TOPP_SeedListGenerator_2_out1" ${DIFF} -in1 SeedListGenerator_2_output_tmp_0.featureXML -in2 ${DATA_DIR_TOPP}/SeedListGenerator_2_output.featureXML )
set_tests_properties("TOPP_SeedListGenerator_2_out1" PROPERTIES DEPENDS "TOPP_SeedListGenerator_2")
add_test("TOPP_SeedListGenerator_3" ${TOPP_BIN_PATH}/SeedListGenerator -test -in ${DATA_DIR_TOPP}/ConsensusXMLFile_1.consensusXML -out_prefix SeedListGenerator_3_output_tmp)
add_test("TOPP_SeedListGenerator_3_out1" ${DIFF} -in1 SeedListGenerator_3_output_tmp_0.featureXML -in2 ${DATA_DIR_TOPP}/SeedListGenerator_3_output1.featureXML )
set_tests_properties("TOPP_SeedListGenerator_3_out1" PROPERTIES DEPENDS "TOPP_SeedListGenerator_3")
add_test("TOPP_SeedListGenerator_3_out2" ${DIFF} -in1 SeedListGenerator_3_output_tmp_1.featureXML -in2 ${DATA_DIR_TOPP}/SeedListGenerator_3_output2.featureXML )
set_tests_properties("TOPP_SeedListGenerator_3_out2" PROPERTIES DEPENDS "TOPP_SeedListGenerator_3")

#------------------------------------------------------------------------------
# SpectraFilter tests
add_test("TOPP_SpectraFilterSqrtMower_1" ${TOPP_BIN_PATH}/SpectraFilterSqrtMower -test -in ${DATA_DIR_TOPP}/SpectraFilterSqrtMower_1_input.mzML -out SpectraFilterSqrtMower.tmp.mzML)
add_test("TOPP_SpectraFilterSqrtMower_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 SpectraFilterSqrtMower.tmp.mzML -in2 ${DATA_DIR_TOPP}/SpectraFilterSqrtMower_1_output.mzML )
set_tests_properties("TOPP_SpectraFilterSqrtMower_1_out1" PROPERTIES DEPENDS "TOPP_SpectraFilterSqrtMower_1")
add_test("TOPP_SpectraFilterWindowMower_1" ${TOPP_BIN_PATH}/SpectraFilterWindowMower -test -in ${DATA_DIR_TOPP}/SpectraFilterWindowMower_1_input.mzML -out SpectraFilterWindowMower_1.tmp.mzML)
add_test("TOPP_SpectraFilterWindowMower_1_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 SpectraFilterWindowMower_1.tmp.mzML -in2 ${DATA_DIR_TOPP}/SpectraFilterWindowMower_1_output.mzML )
set_tests_properties("TOPP_SpectraFilterWindowMower_1_out1" PROPERTIES DEPENDS "TOPP_SpectraFilterWindowMower_1")
add_test("TOPP_SpectraFilterWindowMower_2" ${TOPP_BIN_PATH}/SpectraFilterWindowMower -test -in ${DATA_DIR_TOPP}/SpectraFilterWindowMower_2_input.mzML -out SpectraFilterWindowMower_2.tmp.mzML -ini ${DATA_DIR_TOPP}/SpectraFilterWindowMower_2_parameters.ini)
add_test("TOPP_SpectraFilterWindowMower_2_out1" ${DIFF} -whitelist ${INDEX_WHITELIST} -in1 SpectraFilterWindowMower_2.tmp.mzML -in2 ${DATA_DIR_TOPP}/SpectraFilterWindowMower_2_output.mzML )
set_tests_properties("TOPP_SpectraFilterWindowMower_2_out1" PROPERTIES DEPENDS "TOPP_SpectraFilterWindowMower_2")

#------------------------------------------------------------------------------
# InternalCalibration tests

## idXML file as calibrants
add_test("TOPP_InternalCalibration_1" ${TOPP_BIN_PATH}/InternalCalibration -test -in ${DATA_DIR_TOPP}/InternalCalibration_1_BSA1.mzML -cal:id_in ${DATA_DIR_TOPP}/InternalCalibration_1_BSA1_OMSSA.idXML -out InternalCalibration_1_BSA1_out.tmp.mzML -cal:model_type linear -RANSAC:enabled -RANSAC:iter 500 -RANSAC:threshold 1 -RT_chunking -1 -quality_control:models InternalCalibration_1_models.tmp.csv -quality_control:residuals InternalCalibration_1_residuals.tmp.csv -ms_level 1)
add_test("TOPP_InternalCalibration_1_out1" ${DIFF} -in1 InternalCalibration_1_models.tmp.csv -in2 ${DATA_DIR_TOPP}/InternalCalibration_1_models.csv)
set_tests_properties("TOPP_InternalCalibration_1_out1" PROPERTIES DEPENDS "TOPP_InternalCalibration_1")
add_test("TOPP_InternalCalibration_1_out2" ${DIFF} -in1 InternalCalibration_1_residuals.tmp.csv -in2 ${DATA_DIR_TOPP}/InternalCalibration_1_residuals.csv)
set_tests_properties("TOPP_InternalCalibration_1_out2" PROPERTIES DEPENDS "TOPP_InternalCalibration_1")

## lock masses as calibrants
add_test("TOPP_InternalCalibration_2" ${TOPP_BIN_PATH}/InternalCalibration -test -in ${DATA_DIR_TOPP}/InternalCalibration_2_lockmass.mzML.gz -cal:lock_in ${DATA_DIR_TOPP}/InternalCalibration_2_lock.csv -out InternalCalibration_2_lockmass.tmp.mzML -cal:lock_require_mono -cal:model_type linear -RT_chunking 60 -quality_control:models InternalCalibration_2_models.tmp.csv -quality_control:residuals InternalCalibration_2_residuals.tmp.csv)
add_test("TOPP_InternalCalibration_2_out1" ${DIFF} -in1 InternalCalibration_2_models.tmp.csv -in2 ${DATA_DIR_TOPP}/InternalCalibration_2_models.csv)
set_tests_properties("TOPP_InternalCalibration_2_out1" PROPERTIES DEPENDS "TOPP_InternalCalibration_2")
add_test("TOPP_InternalCalibration_2_out2" ${DIFF} -in1 InternalCalibration_2_residuals.tmp.csv -in2 ${DATA_DIR_TOPP}/InternalCalibration_2_residuals.csv)
set_tests_properties("TOPP_InternalCalibration_2_out2" PROPERTIES DEPENDS "TOPP_InternalCalibration_2")
add_test("TOPP_InternalCalibration_2_convert" ${TOPP_BIN_PATH}/FileConverter -test -in InternalCalibration_2_lockmass.tmp.mzML -no_progress -out InternalCalibration_2_out_lockmass.tmp.dta2d -out_type dta2d)
set_tests_properties("TOPP_InternalCalibration_2_convert" PROPERTIES DEPENDS "TOPP_InternalCalibration_2")
add_test("TOPP_InternalCalibration_2_out3" ${DIFF} -in1 InternalCalibration_2_out_lockmass.tmp.dta2d -in2 ${DATA_DIR_TOPP}/InternalCalibration_2_out_lockmass.dta2d)
set_tests_properties("TOPP_InternalCalibration_2_out3" PROPERTIES DEPENDS "TOPP_InternalCalibration_2_convert")

#------------------------------------------------------------------------------
## ExternalCalibration
add_test("TOPP_ExternalCalibration_1_MS1" ${TOPP_BIN_PATH}/ExternalCalibration -test -in ${DATA_DIR_TOPP}/ExternalCalibration_1_input.mzML -out ExternalCalibration_1_MS1_out.tmp.mzML -offset -5.5 -slope 0.0001 -ms_level 1)
add_test("TOPP_ExternalCalibration_1_out1" ${DIFF} -in1 ExternalCalibration_1_MS1_out.tmp.mzML -in2 ${DATA_DIR_TOPP}/ExternalCalibration_1_MS1_out.mzML)
set_tests_properties("TOPP_ExternalCalibration_1_out1" PROPERTIES DEPENDS "TOPP_ExternalCalibration_1_MS1")

add_test("TOPP_ExternalCalibration_2_MS2" ${TOPP_BIN_PATH}/ExternalCalibration -test -in ${DATA_DIR_TOPP}/ExternalCalibration_1_input.mzML -out ExternalCalibration_2_MS2_out.tmp.mzML -offset -5.5 -slope 0.0001 -ms_level 2)
add_test("TOPP_ExternalCalibration_2_out1" ${DIFF} -in1 ExternalCalibration_2_MS2_out.tmp.mzML -in2 ${DATA_DIR_TOPP}/ExternalCalibration_2_MS2_out.mzML)
set_tests_properties("TOPP_ExternalCalibration_2_out1" PROPERTIES DEPENDS "TOPP_ExternalCalibration_2_MS2")

#------------------------------------------------------------------------------
# TextExporter tests
add_test("TOPP_TextExporter_1" ${TOPP_BIN_PATH}/TextExporter -test -in ${DATA_DIR_TOPP}/TextExporter_1_input.featureXML -no_progress -out TextExporter_1_output.tmp.txt)
add_test("TOPP_TextExporter_1_out1" ${DIFF} -in1 TextExporter_1_output.tmp.txt -in2 ${DATA_DIR_TOPP}/TextExporter_1_output.txt )
set_tests_properties("TOPP_TextExporter_1_out1" PROPERTIES DEPENDS "TOPP_TextExporter_1")

add_test("TOPP_TextExporter_2" ${TOPP_BIN_PATH}/TextExporter -test -in ${DATA_DIR_TOPP}/TextExporter_2_input.consensusXML -ini ${DATA_DIR_TOPP}/TextExporter_2_parameters.ini)
add_test("TOPP_TextExporter_2_out1" ${DIFF} -whitelist "extracted from" -in1 TextExporter_2_consensus_centroids.tmp.csv -in2 ${DATA_DIR_TOPP}/TextExporter_2_consensus_centroids.wsv )
set_tests_properties("TOPP_TextExporter_2_out1" PROPERTIES DEPENDS "TOPP_TextExporter_2")
add_test("TOPP_TextExporter_2_out2" ${DIFF} -whitelist "extracted from" -in1 TextExporter_2_consensus_elements.tmp.csv -in2 ${DATA_DIR_TOPP}/TextExporter_2_consensus_elements.wsv )
set_tests_properties("TOPP_TextExporter_2_out2" PROPERTIES DEPENDS "TOPP_TextExporter_2")
add_test("TOPP_TextExporter_2_out3" ${DIFF} -whitelist "extracted from" -in1 TextExporter_2_consensus_features.tmp.csv -in2 ${DATA_DIR_TOPP}/TextExporter_2_consensus_features.wsv )
set_tests_properties("TOPP_TextExporter_2_out3" PROPERTIES DEPENDS "TOPP_TextExporter_2")
add_test("TOPP_TextExporter_2_out4" ${DIFF} -whitelist "extracted from" -in1 TextExporter_2_consensus_tsv.tmp.tsv -in2 ${DATA_DIR_TOPP}/TextExporter_2_consensus_tsv.tsv )
set_tests_properties("TOPP_TextExporter_2_out4" PROPERTIES DEPENDS "TOPP_TextExporter_2")

add_test("TOPP_TextExporter_3" ${TOPP_BIN_PATH}/TextExporter -test -in ${DATA_DIR_TOPP}/TextExporter_3_input.idXML -no_progress -out TextExporter_3_output.tmp.txt)
add_test("TOPP_TextExporter_3_out1" ${DIFF} -in1 TextExporter_3_output.tmp.txt -in2 ${DATA_DIR_TOPP}/TextExporter_3_output.txt )
set_tests_properties("TOPP_TextExporter_3_out1" PROPERTIES DEPENDS "TOPP_TextExporter_3")

add_test("TOPP_TextExporter_4" ${TOPP_BIN_PATH}/TextExporter -test -in ${DATA_DIR_TOPP}/TextExporter_3_input.idXML -no_progress -out TextExporter_4_output.tmp.txt -id:proteins_only)
add_test("TOPP_TextExporter_4_out1" ${DIFF} -in1 TextExporter_4_output.tmp.txt -in2 ${DATA_DIR_TOPP}/TextExporter_4_output_proteins.txt )
set_tests_properties("TOPP_TextExporter_4_out1" PROPERTIES DEPENDS "TOPP_TextExporter_4")

add_test("TOPP_TextExporter_5" ${TOPP_BIN_PATH}/TextExporter -test -in ${DATA_DIR_TOPP}/TextExporter_5_input.idXML -no_progress -out TextExporter_5_output.tmp.txt -id:peptides_only -id:first_dim_rt)
add_test("TOPP_TextExporter_5_out1" ${DIFF} -in1 TextExporter_5_output.tmp.txt -in2 ${DATA_DIR_TOPP}/TextExporter_5_output_peptides.txt )
set_tests_properties("TOPP_TextExporter_5_out1" PROPERTIES DEPENDS "TOPP_TextExporter_5")

add_test("TOPP_TextExporter_6" ${TOPP_BIN_PATH}/TextExporter -test -in ${DATA_DIR_TOPP}/TextExporter_6_input.featureXML -no_progress -out TextExporter_6_output.tmp.txt -no_ids)
add_test("TOPP_TextExporter_6_out1" ${DIFF} -in1 TextExporter_6_output.tmp.txt -in2 ${DATA_DIR_TOPP}/TextExporter_6_output.txt )
set_tests_properties("TOPP_TextExporter_6_out1" PROPERTIES DEPENDS "TOPP_TextExporter_6")

add_test("TOPP_TextExporter_7" ${TOPP_BIN_PATH}/TextExporter -test -in ${DATA_DIR_TOPP}/TextExporter_7_input.consensusXML -ini ${DATA_DIR_TOPP}/TextExporter_7_parameters.ini)
add_test("TOPP_TextExporter_7_out1" ${DIFF} -whitelist "extracted from" -in1 TextExporter_7_consensus_tsv.tmp.tsv -in2 ${DATA_DIR_TOPP}/TextExporter_7_consensus_tsv.tsv )
set_tests_properties("TOPP_TextExporter_7_out1" PROPERTIES DEPENDS "TOPP_TextExporter_7")

add_test("TOPP_TextExporter_8" ${TOPP_BIN_PATH}/TextExporter -test -in ${DATA_DIR_TOPP}/TextExporter_1_input.featureXML -no_ids -no_progress -out TextExporter_8_output.tmp.txt)
add_test("TOPP_TextExporter_8_out1" ${DIFF} -in1 TextExporter_8_output.tmp.txt -in2 ${DATA_DIR_TOPP}/TextExporter_8_output.txt )
set_tests_properties("TOPP_TextExporter_8_out1" PROPERTIES DEPENDS "TOPP_TextExporter_8")

add_test("TOPP_TextExporter_9" ${TOPP_BIN_PATH}/TextExporter -test -in ${DATA_DIR_TOPP}/TextExporter_9_input.idXML -no_progress -out TextExporter_9_output.tmp.txt -id:add_metavalues 0 -id:add_hit_metavalues 0)
add_test("TOPP_TextExporter_9_out1" ${DIFF} -in1 TextExporter_9_output.tmp.txt -in2 ${DATA_DIR_TOPP}/TextExporter_9_output.txt )
set_tests_properties("TOPP_TextExporter_9_out1" PROPERTIES DEPENDS "TOPP_TextExporter_9")

#------------------------------------------------------------------------------
# FeatureLinker tests
# "labeled" algorithm:
add_test("TOPP_FeatureLinkerLabeled_1" ${TOPP_BIN_PATH}/FeatureLinkerLabeled -test -ini ${DATA_DIR_TOPP}/FeatureLinkerLabeled_1_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerLabeled_1_input.featureXML -out FeatureLinkerLabeled_1_output.tmp.consensusXML)
add_test("TOPP_FeatureLinkerLabeled_1_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerLabeled_1_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FeatureLinkerLabeled_1_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerLabeled_1_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerLabeled_1")
add_test("TOPP_FeatureLinkerLabeled_2" ${TOPP_BIN_PATH}/FeatureLinkerLabeled -test -ini ${DATA_DIR_TOPP}/FeatureLinkerLabeled_2_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerLabeled_2_input.featureXML -out FeatureLinkerLabeled_2_output.tmp.consensusXML)
add_test("TOPP_FeatureLinkerLabeled_2_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerLabeled_2_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FeatureLinkerLabeled_2_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerLabeled_2_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerLabeled_2")
# "unlabelled" algorithm:
add_test("TOPP_FeatureLinkerUnlabeled_1" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeled -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input2.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input3.featureXML -out FeatureLinkerUnlabeled_1_output.tmp.consensusXML)
add_test("TOPP_FeatureLinkerUnlabeled_1_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerUnlabeled_1_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerUnlabeled_1_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeled_1")
add_test("TOPP_FeatureLinkerUnlabeled_2" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeled -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_2_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_2_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_2_input2.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_2_input3.featureXML -out FeatureLinkerUnlabeled_2_output.tmp.consensusXML)
add_test("TOPP_FeatureLinkerUnlabeled_2_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerUnlabeled_2_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_2_output.consensusXML)
set_tests_properties("TOPP_FeatureLinkerUnlabeled_2_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeled_2")
add_test("TOPP_FeatureLinkerUnlabeled_3" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeled -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_3_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_3_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_3_input2.featureXML -out FeatureLinkerUnlabeled_3_output.tmp.consensusXML)
add_test("TOPP_FeatureLinkerUnlabeled_3_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerUnlabeled_3_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_3_output.consensusXML)
set_tests_properties("TOPP_FeatureLinkerUnlabeled_3_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeled_3")
# "unlabeled" algorithm with "use_identifications" option:
add_test("TOPP_FeatureLinkerUnlabeled_4" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeled -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_4_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input2.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input3.featureXML -out FeatureLinkerUnlabeled_4_output.tmp.consensusXML)
add_test("TOPP_FeatureLinkerUnlabeled_4_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerUnlabeled_4_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_4_output.consensusXML)
set_tests_properties("TOPP_FeatureLinkerUnlabeled_4_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeled_4")
# "unlabeled_qt" algorithm:
add_test("TOPP_FeatureLinkerUnlabeledQT_1" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeledQT -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_1_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input2.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input3.featureXML -out FeatureLinkerUnlabeledQT_1_output.tmp.consensusXML)
add_test("TOPP_FeatureLinkerUnlabeledQT_1_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerUnlabeledQT_1_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_1_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerUnlabeledQT_1_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeledQT_1")
# "unlabeled_qt" algorithm with consensusXML input and "keep_subelements" flag:
add_test("TOPP_FeatureLinkerUnlabeledQT_2" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeledQT -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_2_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_1_output.consensusXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_1_output.consensusXML -out FeatureLinkerUnlabeledQT_2_output.tmp.consensusXML)
add_test("TOPP_FeatureLinkerUnlabeledQT_2_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerUnlabeledQT_2_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_2_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerUnlabeledQT_2_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeledQT_2")
add_test("TOPP_FeatureLinkerUnlabeledQT_3" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeledQT -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_3_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_3_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_3_input2.featureXML -out FeatureLinkerUnlabeledQT_3_output.tmp.consensusXML)
add_test("TOPP_FeatureLinkerUnlabeledQT_3_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerUnlabeledQT_3_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_3_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerUnlabeledQT_3_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeledQT_3")
# use_identifications turned on
add_test("TOPP_FeatureLinkerUnlabeledQT_4" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeledQT -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_1_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input2.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input3.featureXML -out FeatureLinkerUnlabeledQT_4_output.tmp.consensusXML -algorithm:use_identifications)
add_test("TOPP_FeatureLinkerUnlabeledQT_4_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerUnlabeledQT_4_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_4_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerUnlabeledQT_4_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeledQT_4")
add_test("TOPP_FeatureLinkerUnlabeledQT_5" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeledQT -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_1_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_5_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_5_input2.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_5_input3.featureXML -out FeatureLinkerUnlabeledQT_5_output.tmp.consensusXML -algorithm:distance_RT:max_difference 200)
add_test("TOPP_FeatureLinkerUnlabeledQT_5_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerUnlabeledQT_5_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_5_output.consensusXML)
set_tests_properties("TOPP_FeatureLinkerUnlabeledQT_5_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeledQT_5")
add_test("TOPP_FeatureLinkerUnlabeledQT_6" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeledQT -test -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_5_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_5_input2.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_5_input3.featureXML -out FeatureLinkerUnlabeledQT_6_output.tmp.consensusXML -algorithm:use_identifications -algorithm:distance_RT:max_difference 200)
add_test("TOPP_FeatureLinkerUnlabeledQT_6_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerUnlabeledQT_6_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_6_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerUnlabeledQT_6_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeledQT_6")
# FeatureLinkerUnlabeledKD
add_test("TOPP_FeatureLinkerUnlabeledKD_1" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeledKD -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_1_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input2.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeled_1_input3.featureXML -out FeatureLinkerUnlabeledKD_1_output.tmp.consensusXML)
add_test("TOPP_FeatureLinkerUnlabeledKD_1_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerUnlabeledKD_1_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_1_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerUnlabeledKD_1_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeledKD_1")
add_test("TOPP_FeatureLinkerUnlabeledKD_2" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeledKD -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_2_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_1_output.consensusXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_1_output.consensusXML -out FeatureLinkerUnlabeledKD_2_output.tmp.consensusXML)
add_test("TOPP_FeatureLinkerUnlabeledKD_2_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerUnlabeledKD_2_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_2_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerUnlabeledKD_2_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeledKD_2")
add_test("TOPP_FeatureLinkerUnlabeledKD_3" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeledKD -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_3_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_3_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledQT_3_input2.featureXML -out FeatureLinkerUnlabeledKD_3_output.tmp.consensusXML)
add_test("TOPP_FeatureLinkerUnlabeledKD_3_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerUnlabeledKD_3_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_3_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerUnlabeledKD_3_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeledKD_3")
# use charge/adduct information different ways
add_test("TOPP_FeatureLinkerUnlabeledKD_4" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeledKD -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_4_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_dc_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_dc_input2.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_dc_input3.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_dc_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_dc_input2.featureXML -out FeatureLinkerUnlabeledKD_4_output.tmp.consensusXML -algorithm:link:charge_merging Identical -algorithm:link:adduct_merging Any)
 add_test("TOPP_FeatureLinkerUnlabeledKD_4_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerUnlabeledKD_4_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_4_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerUnlabeledKD_4_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeledKD_4")
add_test("TOPP_FeatureLinkerUnlabeledKD_5" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeledKD -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_4_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_dc_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_dc_input2.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_dc_input3.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_dc_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_dc_input2.featureXML -out FeatureLinkerUnlabeledKD_5_output.tmp.consensusXML -algorithm:link:charge_merging With_charge_zero -algorithm:link:adduct_merging Any)
 add_test("TOPP_FeatureLinkerUnlabeledKD_5_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerUnlabeledKD_5_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_5_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerUnlabeledKD_5_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeledKD_5")
add_test("TOPP_FeatureLinkerUnlabeledKD_6" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeledKD -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_4_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_dc_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_dc_input2.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_dc_input3.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_dc_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_dc_input2.featureXML -out FeatureLinkerUnlabeledKD_6_output.tmp.consensusXML -algorithm:link:charge_merging Any -algorithm:link:adduct_merging With_unknown_adducts)
 add_test("TOPP_FeatureLinkerUnlabeledKD_6_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerUnlabeledKD_6_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_6_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerUnlabeledKD_6_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeledKD_6")
add_test("TOPP_FeatureLinkerUnlabeledKD_7" ${TOPP_BIN_PATH}/FeatureLinkerUnlabeledKD -test -ini ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_4_parameters.ini -in ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_dc_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_dc_input2.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_dc_input3.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_dc_input1.featureXML ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_dc_input2.featureXML -out FeatureLinkerUnlabeledKD_7_output.tmp.consensusXML -algorithm:link:charge_merging Any -algorithm:link:adduct_merging Identical)
 add_test("TOPP_FeatureLinkerUnlabeledKD_7_out1" ${DIFF} -whitelist "id=" "href=" -in1 FeatureLinkerUnlabeledKD_7_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/FeatureLinkerUnlabeledKD_7_output.consensusXML )
set_tests_properties("TOPP_FeatureLinkerUnlabeledKD_7_out1" PROPERTIES DEPENDS "TOPP_FeatureLinkerUnlabeledKD_7")



#------------------------------------------------------------------------------
# IDMapper tests
add_test("TOPP_IDMapper_1" ${TOPP_BIN_PATH}/IDMapper -test -in ${DATA_DIR_TOPP}/IDMapper_1_input.featureXML -out IDMapper_1_output.tmp.featureXML -id ${DATA_DIR_TOPP}/IDMapper_1_input.idXML -mz_measure Da -mz_tolerance 1 -ignore_charge -mz_reference precursor -feature:use_centroid_mz false)
add_test("TOPP_IDMapper_1_out1" ${DIFF} -in1 IDMapper_1_output.tmp.featureXML -in2 ${DATA_DIR_TOPP}/IDMapper_1_output.featureXML )
set_tests_properties("TOPP_IDMapper_1_out1" PROPERTIES DEPENDS "TOPP_IDMapper_1")
add_test("TOPP_IDMapper_2" ${TOPP_BIN_PATH}/IDMapper -test -in ${DATA_DIR_TOPP}/IDMapper_2_input.consensusXML -out IDMapper_2_output.tmp.consensusXML -id ${DATA_DIR_TOPP}/IDMapper_2_input.idXML  -mz_measure Da -mz_tolerance 1 -mz_reference precursor)
add_test("TOPP_IDMapper_2_out1" ${DIFF} -in1 IDMapper_2_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/IDMapper_2_output.consensusXML )
set_tests_properties("TOPP_IDMapper_2_out1" PROPERTIES DEPENDS "TOPP_IDMapper_2")
add_test("TOPP_IDMapper_3" ${TOPP_BIN_PATH}/IDMapper -test -in ${DATA_DIR_TOPP}/IDMapper_3_input.featureXML -out IDMapper_3_output.tmp.featureXML -id ${DATA_DIR_TOPP}/IDMapper_3_input.idXML  -mz_measure ppm -mz_tolerance 3 -rt_tolerance 4 -ignore_charge -mz_reference precursor)
add_test("TOPP_IDMapper_3_out1" ${DIFF} -in1 IDMapper_3_output.tmp.featureXML -in2 ${DATA_DIR_TOPP}/IDMapper_3_output.featureXML )
set_tests_properties("TOPP_IDMapper_3_out1" PROPERTIES DEPENDS "TOPP_IDMapper_3")
add_test("TOPP_IDMapper_4" ${TOPP_BIN_PATH}/IDMapper -test -in
${DATA_DIR_TOPP}/IDMapper_4_input.featureXML -out IDMapper_4_output.tmp.featureXML -id
${DATA_DIR_TOPP}/IDMapper_4_input.idXML -spectra:in
${DATA_DIR_TOPP}/IDMapper_4_input.mzML -mz_measure ppm -mz_tolerance 20
-rt_tolerance 10)
add_test("TOPP_IDMapper_4_out1" ${DIFF} -in1 IDMapper_4_output.tmp.featureXML -in2 ${DATA_DIR_TOPP}/IDMapper_4_output.featureXML )
set_tests_properties("TOPP_IDMapper_4_out1" PROPERTIES DEPENDS "TOPP_IDMapper_4")
add_test("TOPP_IDMapper_5" ${TOPP_BIN_PATH}/IDMapper -test -in
${DATA_DIR_TOPP}/IDMapper_5_input.featureXML -out IDMapper_5_output.tmp.featureXML -id
${DATA_DIR_TOPP}/IDMapper_5_input.idXML -spectra:in
${DATA_DIR_TOPP}/IDMapper_5_input.mzML -mz_measure ppm -mz_tolerance 10
-rt_tolerance 20 -feature:use_centroid_mz false)
add_test("TOPP_IDMapper_5_out1" ${DIFF} -in1 IDMapper_5_output.tmp.featureXML -in2 ${DATA_DIR_TOPP}/IDMapper_5_output.featureXML )
set_tests_properties("TOPP_IDMapper_5_out1" PROPERTIES DEPENDS "TOPP_IDMapper_5")

add_test("TOPP_IDMapper_6" ${TOPP_BIN_PATH}/IDMapper -test -in
${DATA_DIR_TOPP}/IDMapper_6_input.consensusXML -out IDMapper_6_output.tmp.consensusXML -id
${DATA_DIR_TOPP}/IDMapper_6_input.idXML)
add_test("TOPP_IDMapper_6_out1" ${DIFF} -in1 IDMapper_6_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/IDMapper_6_output.consensusXML )
set_tests_properties("TOPP_IDMapper_6_out1" PROPERTIES DEPENDS "TOPP_IDMapper_6")

#------------------------------------------------------------------------------
# IDRipper tests
add_test("TOPP_IDRipper_1" ${TOPP_BIN_PATH}/IDRipper -test -in ${DATA_DIR_TOPP}/IDRipper_1_input.idXML -out ./ )
add_test("TOPP_IDRipper_1_out1" ${DIFF} -in1 IDRipper_1_output_1.idXML  -in2 ${DATA_DIR_TOPP}/IDRipper_1_output1.idXML )
set_tests_properties("TOPP_IDRipper_1_out1" PROPERTIES DEPENDS "TOPP_IDRipper_1")
add_test("TOPP_IDRipper_1_out2" ${DIFF} -in1 IDRipper_1_output_2.idXML -in2 ${DATA_DIR_TOPP}/IDRipper_1_output2.idXML )
set_tests_properties("TOPP_IDRipper_1_out2" PROPERTIES DEPENDS "TOPP_IDRipper_1")
add_test("TOPP_IDRipper_2" ${TOPP_BIN_PATH}/IDRipper -test -in ${DATA_DIR_TOPP}/IDRipper_2_input.idXML -out ./ -numeric_filenames -split_ident_runs )
add_test("TOPP_IDRipper_2_out1" ${DIFF} -in1 IDRipper_2_input_0_0.idXML -in2 ${DATA_DIR_TOPP}/IDRipper_2_input_0_0.idXML )
set_tests_properties("TOPP_IDRipper_2_out1" PROPERTIES DEPENDS "TOPP_IDRipper_2")
add_test("TOPP_IDRipper_2_out2" ${DIFF} -in1 IDRipper_2_input_1_0.idXML -in2 ${DATA_DIR_TOPP}/IDRipper_2_input_1_0.idXML )
set_tests_properties("TOPP_IDRipper_2_out2" PROPERTIES DEPENDS "TOPP_IDRipper_2")
add_test("TOPP_IDRipper_2_out3" ${DIFF} -in1 IDRipper_2_input_2_1.idXML -in2 ${DATA_DIR_TOPP}/IDRipper_2_input_2_1.idXML )
set_tests_properties("TOPP_IDRipper_2_out3" PROPERTIES DEPENDS "TOPP_IDRipper_2")
add_test("TOPP_IDRipper_2a" ${TOPP_BIN_PATH}/IDRipper -test -in ${DATA_DIR_TOPP}/IDRipper_2_input.idXML -out ./ -numeric_filenames )
add_test("TOPP_IDRipper_2a_out1" ${DIFF} -in1 IDRipper_2_input_0.idXML -in2 ${DATA_DIR_TOPP}/IDRipper_2a_output_0.idXML )
set_tests_properties("TOPP_IDRipper_2a_out1" PROPERTIES DEPENDS "TOPP_IDRipper_2a")
add_test("TOPP_IDRipper_2a_out2" ${DIFF} -in1 IDRipper_2_input_1.idXML -in2 ${DATA_DIR_TOPP}/IDRipper_2a_output_1.idXML )
set_tests_properties("TOPP_IDRipper_2a_out2" PROPERTIES DEPENDS "TOPP_IDRipper_2a")

# use IDMerger results to test. Input -> IDMerger -> IDRipper -> Output == Input
add_test("TOPP_IDRipper_3_prep" ${TOPP_BIN_PATH}/IDMerger -test -in ${DATA_DIR_TOPP}/IDRipper_3_input1.idXML ${DATA_DIR_TOPP}/IDRipper_3_input2.idXML -out IDRipper_3_output.tmp.idXML)
add_test("TOPP_IDRipper_3_prep_out1" ${DIFF} -whitelist "?xml-stylesheet" "file_origin" -in1 IDRipper_3_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDRipper_3_output.idXML)
set_tests_properties("TOPP_IDRipper_3_prep_out1" PROPERTIES DEPENDS "TOPP_IDRipper_3_prep")
## create an output directory during configure time in the binary tree
set(TMP_RIP_PATH ${PROJECT_BINARY_DIR}/tmp_path/) ## ${PROJECT_BINARY_DIR} already includes '/src/tests/topp/'
file(MAKE_DIRECTORY ${TMP_RIP_PATH})
add_test("TOPP_IDRipper_3" ${TOPP_BIN_PATH}/IDRipper -test -in ${DATA_DIR_TOPP}/IDRipper_3_output.idXML -out ${TMP_RIP_PATH})
set_tests_properties("TOPP_IDRipper_3" PROPERTIES DEPENDS "TOPP_IDRipper_3_prep")
add_test("TOPP_IDRipper_3_out1" ${DIFF} -in1 ${TMP_RIP_PATH}/IDRipper_3_input1.idXML -in2 ${DATA_DIR_TOPP}/IDRipper_3_input1.idXML)
set_tests_properties("TOPP_IDRipper_3_out1" PROPERTIES DEPENDS "TOPP_IDRipper_3")
add_test("TOPP_IDRipper_3_out2" ${DIFF} -in1 ${TMP_RIP_PATH}/IDRipper_3_input2.idXML -in2 ${DATA_DIR_TOPP}/IDRipper_3_input2.idXML)
set_tests_properties("TOPP_IDRipper_3_out2" PROPERTIES DEPENDS "TOPP_IDRipper_3")


#------------------------------------------------------------------------------
# ConsensusID tests
# idXML input, "PEPMatrix" algorithm using PEP30MS matrix (example from the JPR paper):
add_test("TOPP_ConsensusID_1" ${TOPP_BIN_PATH}/ConsensusID -test -in ${DATA_DIR_TOPP}/ConsensusID_1_input.idXML -out ConsensusID_1_output.tmp.idXML -algorithm PEPMatrix -PEPMatrix:matrix PAM30MS)
add_test("TOPP_ConsensusID_1_out1" ${DIFF} -whitelist "?xml-stylesheet" "IdentificationRun date" -in1 ConsensusID_1_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/ConsensusID_1_output.idXML )
set_tests_properties("TOPP_ConsensusID_1_out1" PROPERTIES DEPENDS "TOPP_ConsensusID_1")
# featureXML input, "average" algorithm:
add_test("TOPP_ConsensusID_2" ${TOPP_BIN_PATH}/ConsensusID -test -in ${DATA_DIR_TOPP}/ConsensusID_2_input.featureXML -out ConsensusID_2_output.tmp.featureXML -algorithm average)
add_test("TOPP_ConsensusID_2_out1" ${DIFF} -whitelist "IdentificationRun id" -in1 ConsensusID_2_output.tmp.featureXML -in2 ${DATA_DIR_TOPP}/ConsensusID_2_output.featureXML )
set_tests_properties("TOPP_ConsensusID_2_out1" PROPERTIES DEPENDS "TOPP_ConsensusID_2")
# consensusXML input, "best" algorithm:
add_test("TOPP_ConsensusID_3" ${TOPP_BIN_PATH}/ConsensusID -test -in ${DATA_DIR_TOPP}/ConsensusID_3_input.consensusXML -out ConsensusID_3_output.tmp.consensusXML -algorithm best)
add_test("TOPP_ConsensusID_3_out1" ${DIFF} -whitelist "?xml-stylesheet" "IdentificationRun id" -in1 ConsensusID_3_output.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/ConsensusID_3_output.consensusXML )
set_tests_properties("TOPP_ConsensusID_3_out1" PROPERTIES DEPENDS "TOPP_ConsensusID_3")
# idXML input, "PEPMatrix" algorithm using identity matrix, limited number of considered hits:
add_test("TOPP_ConsensusID_4" ${TOPP_BIN_PATH}/ConsensusID -test -in ${DATA_DIR_TOPP}/ConsensusID_1_input.idXML -out ConsensusID_4_output.tmp.idXML -algorithm PEPMatrix -PEPMatrix:matrix identity -filter:considered_hits 6)
add_test("TOPP_ConsensusID_4_out1" ${DIFF} -whitelist "?xml-stylesheet" "IdentificationRun date" -in1 ConsensusID_4_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/ConsensusID_4_output.idXML )
set_tests_properties("TOPP_ConsensusID_4_out1" PROPERTIES DEPENDS "TOPP_ConsensusID_4")
# idXML input, "PEPIons" algorithm:
add_test("TOPP_ConsensusID_5" ${TOPP_BIN_PATH}/ConsensusID -test -in ${DATA_DIR_TOPP}/ConsensusID_1_input.idXML -out ConsensusID_5_output.tmp.idXML -algorithm PEPIons)
add_test("TOPP_ConsensusID_5_out1" ${DIFF} -whitelist "?xml-stylesheet" "IdentificationRun date" -in1 ConsensusID_5_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/ConsensusID_5_output.idXML )
set_tests_properties("TOPP_ConsensusID_5_out1" PROPERTIES DEPENDS "TOPP_ConsensusID_5")
# idXML input, "best" algorithm, filtering:
add_test("TOPP_ConsensusID_6" ${TOPP_BIN_PATH}/ConsensusID -test -in ${DATA_DIR_TOPP}/ConsensusID_1_input.idXML -out ConsensusID_6_output.tmp.idXML -algorithm best -filter:min_support 0.5)
add_test("TOPP_ConsensusID_6_out1" ${DIFF} -whitelist "?xml-stylesheet" "IdentificationRun date" -in1 ConsensusID_6_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/ConsensusID_6_output.idXML )
set_tests_properties("TOPP_ConsensusID_6_out1" PROPERTIES DEPENDS "TOPP_ConsensusID_6")
# idXML input, "best" algorithm, per_spectrum:
add_test("TOPP_ConsensusID_7" ${TOPP_BIN_PATH}/ConsensusID -test -in ${DATA_DIR_TOPP}/ConsensusID_6_input.idXML -out ConsensusID_7_output.tmp.idXML -algorithm best -per_spectrum -filter:keep_old_scores)
add_test("TOPP_ConsensusID_7_out1" ${DIFF} -whitelist "?xml-stylesheet" "IdentificationRun date" -in1 ConsensusID_7_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/ConsensusID_7_output.idXML )
set_tests_properties("TOPP_ConsensusID_7_out1" PROPERTIES DEPENDS "TOPP_ConsensusID_7")
# idXML input, with Percolator, "best" algorithm, per_spectrum:
add_test("TOPP_ConsensusID_8" ${TOPP_BIN_PATH}/ConsensusID -test -in ${DATA_DIR_TOPP}/ConsensusID_8_input.idXML -out ConsensusID_8_output.tmp.idXML -algorithm best -per_spectrum -filter:keep_old_scores)
add_test("TOPP_ConsensusID_8_out1" ${DIFF} -whitelist "?xml-stylesheet" "IdentificationRun date" -in1 ConsensusID_8_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/ConsensusID_8_output.idXML )
set_tests_properties("TOPP_ConsensusID_8_out1" PROPERTIES DEPENDS "TOPP_ConsensusID_8")

#------------------------------------------------------------------------------
# FalseDiscoveryRate_tests (1: default param, 2: with runs separated, 3: with charge variants individual, 4: with charge variants and no charge 1 decoys -> q_value 0)
add_test("TOPP_FalseDiscoveryRate_1" ${TOPP_BIN_PATH}/FalseDiscoveryRate -test -in ${DATA_DIR_TOPP}/FalseDiscoveryRate_OMSSA.idXML -out FalseDiscoveryRate_output_1.tmp.idXML -PSM true -protein false)
add_test("TOPP_FalseDiscoveryRate_1_out1" ${DIFF} -whitelist "?xml-stylesheet" -in1 FalseDiscoveryRate_output_1.tmp.idXML -in2 ${DATA_DIR_TOPP}/FalseDiscoveryRate_output_1.idXML)
set_tests_properties("TOPP_FalseDiscoveryRate_1_out1" PROPERTIES DEPENDS "TOPP_FalseDiscoveryRate_1")
add_test("TOPP_FalseDiscoveryRate_2" ${TOPP_BIN_PATH}/FalseDiscoveryRate -test -in ${DATA_DIR_TOPP}/FalseDiscoveryRate_OMSSA.idXML -out FalseDiscoveryRate_output_2.tmp.idXML -algorithm:treat_runs_separately -PSM true -protein false)
add_test("TOPP_FalseDiscoveryRate_2_out1" ${DIFF} -whitelist "?xml-stylesheet" -in1 FalseDiscoveryRate_output_2.tmp.idXML -in2 ${DATA_DIR_TOPP}/FalseDiscoveryRate_output_2.idXML)
set_tests_properties("TOPP_FalseDiscoveryRate_2_out1" PROPERTIES DEPENDS "TOPP_FalseDiscoveryRate_2")
add_test("TOPP_FalseDiscoveryRate_3" ${TOPP_BIN_PATH}/FalseDiscoveryRate -test -in ${DATA_DIR_TOPP}/FalseDiscoveryRate_OMSSA.idXML -out FalseDiscoveryRate_output_3.tmp.idXML -algorithm:split_charge_variants -PSM true -protein false)
add_test("TOPP_FalseDiscoveryRate_3_out1" ${DIFF} -whitelist "?xml-stylesheet" -in1 FalseDiscoveryRate_output_3.tmp.idXML -in2 ${DATA_DIR_TOPP}/FalseDiscoveryRate_output_3.idXML)
set_tests_properties("TOPP_FalseDiscoveryRate_3_out1" PROPERTIES DEPENDS "TOPP_FalseDiscoveryRate_3")
add_test("TOPP_FalseDiscoveryRate_4" ${TOPP_BIN_PATH}/FalseDiscoveryRate -test -in ${DATA_DIR_TOPP}/FalseDiscoveryRate_OMSSA_4.idXML -out FalseDiscoveryRate_output_4.tmp.idXML -algorithm:split_charge_variants -PSM true -protein false)
add_test("TOPP_FalseDiscoveryRate_4_out1" ${DIFF} -whitelist "?xml-stylesheet" -in1 FalseDiscoveryRate_output_4.tmp.idXML -in2 ${DATA_DIR_TOPP}/FalseDiscoveryRate_output_4.idXML)
set_tests_properties("TOPP_FalseDiscoveryRate_4_out1" PROPERTIES DEPENDS "TOPP_FalseDiscoveryRate_4")
# proteins only, keep decoys (input is derived from "THIRDPARTY/FidoAdapter_1_output.idXML"):
# use force flag since there is no inference score in it and FDRs are basically meaningless
add_test("TOPP_FalseDiscoveryRate_5" ${TOPP_BIN_PATH}/FalseDiscoveryRate -test -in ${DATA_DIR_TOPP}/FalseDiscoveryRate_5_input.idXML -out FalseDiscoveryRate_5_output.tmp.idXML -PSM false -protein true -algorithm:add_decoy_proteins -force)
add_test("TOPP_FalseDiscoveryRate_5_out1" ${DIFF} -whitelist "?xml-stylesheet" -in1 FalseDiscoveryRate_5_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/FalseDiscoveryRate_5_output.idXML)
set_tests_properties("TOPP_FalseDiscoveryRate_5_out1" PROPERTIES DEPENDS "TOPP_FalseDiscoveryRate_5")
# FDR filtering
add_test("TOPP_FalseDiscoveryRate_6" ${TOPP_BIN_PATH}/FalseDiscoveryRate -test -in ${DATA_DIR_TOPP}/FalseDiscoveryRate_6_input.idXML -out FalseDiscoveryRate_6_output.tmp.idXML -PSM true -protein false -FDR:PSM 0.05)
add_test("TOPP_FalseDiscoveryRate_6_out1" ${DIFF} -whitelist "?xml-stylesheet" -in1 FalseDiscoveryRate_6_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/FalseDiscoveryRate_6_output.idXML)
set_tests_properties("TOPP_FalseDiscoveryRate_6_out1" PROPERTIES DEPENDS "TOPP_FalseDiscoveryRate_6")
# FDR filtering, proteins only, decoys removed (default)
add_test("TOPP_FalseDiscoveryRate_7" ${TOPP_BIN_PATH}/FalseDiscoveryRate -test -in ${DATA_DIR_TOPP}/FalseDiscoveryRate_7_input.idXML -out FalseDiscoveryRate_7_output.tmp.idXML -PSM false -protein true -FDR:protein 0.30 -force)
add_test("TOPP_FalseDiscoveryRate_7_out1" ${DIFF} -whitelist "?xml-stylesheet" -in1 FalseDiscoveryRate_7_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/FalseDiscoveryRate_7_output.idXML)
set_tests_properties("TOPP_FalseDiscoveryRate_7_out1" PROPERTIES DEPENDS "TOPP_FalseDiscoveryRate_7")


#------------------------------------------------------------------------------
# IDPosteriorErrorProbability
# Mascot input:
add_test("TOPP_IDPosteriorErrorProbability_1" ${TOPP_BIN_PATH}/IDPosteriorErrorProbability -test -in ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_Mascot_input.idXML -out 	IDPosteriorErrorProbability_output_1.tmp.idXML)
add_test("TOPP_IDPosteriorErrorProbability_1_out1" ${DIFF} -in1 IDPosteriorErrorProbability_output_1.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_Mascot_output.idXML)
set_tests_properties("TOPP_IDPosteriorErrorProbability_1_out1" PROPERTIES DEPENDS "TOPP_IDPosteriorErrorProbability_1")
# X! Tandem input:
add_test("TOPP_IDPosteriorErrorProbability_2" ${TOPP_BIN_PATH}/IDPosteriorErrorProbability -test -in ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_XTandem_input.idXML -out IDPosteriorErrorProbability_output_2.tmp.idXML)
add_test("TOPP_IDPosteriorErrorProbability_2_out1" ${DIFF} -in1 IDPosteriorErrorProbability_output_2.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_XTandem_output.idXML)
set_tests_properties("TOPP_IDPosteriorErrorProbability_2_out1" PROPERTIES DEPENDS "TOPP_IDPosteriorErrorProbability_2")
# OMSSA input:
add_test("TOPP_IDPosteriorErrorProbability_3" ${TOPP_BIN_PATH}/IDPosteriorErrorProbability -test -in ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_OMSSA_input.idXML -out IDPosteriorErrorProbability_output_3.tmp.idXML)
add_test("TOPP_IDPosteriorErrorProbability_3_out1" ${DIFF} -in1 IDPosteriorErrorProbability_output_3.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_OMSSA_output.idXML)
set_tests_properties("TOPP_IDPosteriorErrorProbability_3_out1" PROPERTIES DEPENDS "TOPP_IDPosteriorErrorProbability_3")
# OMSSA input with "split_charge":
add_test("TOPP_IDPosteriorErrorProbability_4" ${TOPP_BIN_PATH}/IDPosteriorErrorProbability -test -in ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_OMSSA_input2.idXML -out IDPosteriorErrorProbability_output_4.tmp.idXML -split_charge)
add_test("TOPP_IDPosteriorErrorProbability_4_out1" ${DIFF} -in1 IDPosteriorErrorProbability_output_4.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_OMSSA_output2.idXML)
set_tests_properties("TOPP_IDPosteriorErrorProbability_4_out1" PROPERTIES DEPENDS "TOPP_IDPosteriorErrorProbability_4")
# X! Tandem input with "split_charge":
add_test("TOPP_IDPosteriorErrorProbability_5" ${TOPP_BIN_PATH}/IDPosteriorErrorProbability -test -in ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_XTandem_input2.idXML -out IDPosteriorErrorProbability_output_5.tmp.idXML -split_charge)
add_test("TOPP_IDPosteriorErrorProbability_5_out1" ${DIFF} -in1 IDPosteriorErrorProbability_output_5.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_XTandem_output2.idXML)
set_tests_properties("TOPP_IDPosteriorErrorProbability_5_out1" PROPERTIES DEPENDS "TOPP_IDPosteriorErrorProbability_5")
# Mascot input with "split_charge":
add_test("TOPP_IDPosteriorErrorProbability_6" ${TOPP_BIN_PATH}/IDPosteriorErrorProbability -test -in ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_Mascot_input2.idXML -out IDPosteriorErrorProbability_output_6.tmp.idXML -split_charge)
add_test("TOPP_IDPosteriorErrorProbability_6_out1" ${DIFF} -in1 IDPosteriorErrorProbability_output_6.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_Mascot_output2.idXML)
set_tests_properties("TOPP_IDPosteriorErrorProbability_6_out1" PROPERTIES DEPENDS "TOPP_IDPosteriorErrorProbability_6")
# with "ignore_bad_data" (no check for correctness of output, bc of unresolved platform specific errors):
add_test("TOPP_IDPosteriorErrorProbability_7" ${TOPP_BIN_PATH}/IDPosteriorErrorProbability -test -in ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_bad_data.idXML -out IDPosteriorErrorProbability_bad_data_output.tmp.idXML -ignore_bad_data)
# add_test("TOPP_IDPosteriorErrorProbability_7_out1" ${DIFF} -in1 IDPosteriorErrorProbability_bad_data_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_bad_data_out.idXML)
# with "prob_correct":
add_test("TOPP_IDPosteriorErrorProbability_8" ${TOPP_BIN_PATH}/IDPosteriorErrorProbability -test -in ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_OMSSA_input.idXML -out IDPosteriorErrorProbability_output_8.tmp.idXML -prob_correct)
add_test("TOPP_IDPosteriorErrorProbability_8_out1" ${DIFF} -in1 IDPosteriorErrorProbability_output_8.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDPosteriorErrorProbability_prob_correct_output.idXML)
set_tests_properties("TOPP_IDPosteriorErrorProbability_8_out1" PROPERTIES DEPENDS "TOPP_IDPosteriorErrorProbability_8")

#------------------------------------------------------------------------------
# ProteinResolver tests
### TODO test for idXML
add_test("TOPP_ProteinResolver_1" ${TOPP_BIN_PATH}/ProteinResolver -test -in ${DATA_DIR_TOPP}/ProteinResolver_1_input.consensusXML -fasta ${DATA_DIR_TOPP}/ProteinResolver_1_input.fasta -protein_groups ProteinResolver_1_output1.tmp.csv -peptide_table ProteinResolver_1_output2.tmp.csv -protein_table ProteinResolver_1_output3.tmp.csv)
add_test("TOPP_ProteinResolver_1_out1" ${DIFF} -in1 ProteinResolver_1_output1.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteinResolver_1_output1.txt)
set_tests_properties("TOPP_ProteinResolver_1_out1" PROPERTIES DEPENDS "TOPP_ProteinResolver_1")
add_test("TOPP_ProteinResolver_1_out2" ${DIFF} -in1 ProteinResolver_1_output2.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteinResolver_1_output2.txt)
set_tests_properties("TOPP_ProteinResolver_1_out2" PROPERTIES DEPENDS "TOPP_ProteinResolver_1")
add_test("TOPP_ProteinResolver_1_out3" ${DIFF} -in1 ProteinResolver_1_output3.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteinResolver_1_output3.txt)
set_tests_properties("TOPP_ProteinResolver_1_out3" PROPERTIES DEPENDS "TOPP_ProteinResolver_1")
#add_test("TOPP_ProteinResolver_2" ${TOPP_BIN_PATH}/ProteinResolver -test -in ${DATA_DIR_TOPP}/ProteinResolver_2_input.consensusXML -fasta ${DATA_DIR_TOPP}/ProteinResolver_2_input.fasta -design ${DATA_DIR_TOPP}/ProteinResolver_2_input.txt -protein_groups ProteinResolver_2_output1.tmp.csv -peptide_table ProteinResolver_2_output2.tmp.csv -protein_table ProteinResolver_2_output3.tmp.csv)
#add_test("TOPP_ProteinResolver_2_out1" ${DIFF} -in1 ProteinResolver_2_output1.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteinResolver_2_output1.txt)
#set_tests_properties("TOPP_ProteinResolver_1_out1" PROPERTIES DEPENDS "TOPP_ProteinResolver_1")
#add_test("TOPP_ProteinResolver_2_out2" ${DIFF} -in1 ProteinResolver_2_output2.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteinResolver_2_output2.txt)
#set_tests_properties("TOPP_ProteinResolver_1_out1" PROPERTIES DEPENDS "TOPP_ProteinResolver_1")
#add_test("TOPP_ProteinResolver_2_out3" ${DIFF} -in1 ProteinResolver_2_output3.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteinResolver_2_output3.txt)
#set_tests_properties("TOPP_ProteinResolver_1_out1" PROPERTIES DEPENDS "TOPP_ProteinResolver_1")

#------------------------------------------------------------------------------
# ProteinQuantifier tests
add_test("TOPP_ProteinQuantifier_1" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.featureXML -out ProteinQuantifier_1_output1.tmp.csv -peptide_out ProteinQuantifier_1_output2.tmp.csv)
add_test("TOPP_ProteinQuantifier_1_out1" ${DIFF} -in1 ProteinQuantifier_1_output1.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_1_output1.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_1_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_1")
add_test("TOPP_ProteinQuantifier_1_out2" ${DIFF} -in1 ProteinQuantifier_1_output2.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_1_output2.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_1_out2" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_1")
add_test("TOPP_ProteinQuantifier_2" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.featureXML -out ProteinQuantifier_2_output1.tmp.csv -peptide_out ProteinQuantifier_2_output2.tmp.csv -top:N 2 -top:include_all -best_charge_and_fraction -top:aggregate sum)
add_test("TOPP_ProteinQuantifier_2_out1" ${DIFF} -in1 ProteinQuantifier_2_output1.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_2_output1.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_2_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_2")
add_test("TOPP_ProteinQuantifier_2_out2" ${DIFF} -in1 ProteinQuantifier_2_output2.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_2_output2.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_2_out2" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_2")
# handling of modified peptides:
add_test("TOPP_ProteinQuantifier_3" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_3_input.featureXML -out ProteinQuantifier_3_output1.tmp.csv -peptide_out ProteinQuantifier_3_output2.tmp.csv -top:N 2 -top:include_all -top:aggregate mean)
add_test("TOPP_ProteinQuantifier_3_out1" ${DIFF} -in1 ProteinQuantifier_3_output1.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_3_output1.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_3_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_3")
add_test("TOPP_ProteinQuantifier_3_out2" ${DIFF} -in1 ProteinQuantifier_3_output2.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_3_output2.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_3_out2" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_3")
# consensusXML input, different options for protein quantification:
add_test("TOPP_ProteinQuantifier_4" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.consensusXML -out ProteinQuantifier_4_output.tmp.csv -top:N 0 -top:aggregate sum)
add_test("TOPP_ProteinQuantifier_4_out1" ${DIFF} -in1 ProteinQuantifier_4_output.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_4_output.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_4_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_4")
add_test("TOPP_ProteinQuantifier_5" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.consensusXML -out ProteinQuantifier_5_output.tmp.csv -top:N 3 -top:aggregate sum)
add_test("TOPP_ProteinQuantifier_5_out1" ${DIFF} -in1 ProteinQuantifier_5_output.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_5_output.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_5_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_5")
add_test("TOPP_ProteinQuantifier_6" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.consensusXML -out ProteinQuantifier_6_output.tmp.csv -top:N 3 -top:include_all -top:aggregate sum)
add_test("TOPP_ProteinQuantifier_6_out1" ${DIFF} -in1 ProteinQuantifier_6_output.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_6_output.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_6_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_6")
add_test("TOPP_ProteinQuantifier_7" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.consensusXML -out ProteinQuantifier_7_output.tmp.csv -top:N 0 -consensus:fix_peptides -top:aggregate sum)
add_test("TOPP_ProteinQuantifier_7_out1" ${DIFF} -in1 ProteinQuantifier_7_output.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_7_output.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_7_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_7")
add_test("TOPP_ProteinQuantifier_8" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.consensusXML -out ProteinQuantifier_8_output.tmp.csv -top:N 3 -consensus:fix_peptides -top:aggregate sum)
add_test("TOPP_ProteinQuantifier_8_out1" ${DIFF} -in1 ProteinQuantifier_8_output.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_8_output.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_8_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_8")
add_test("TOPP_ProteinQuantifier_9" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.consensusXML -out ProteinQuantifier_9_output.tmp.csv -mztab ProteinQuantifier_9_output_mztab.tmp.mzTab -top:N 3 -top:include_all -consensus:fix_peptides -top:aggregate sum)
add_test("TOPP_ProteinQuantifier_9_out1" ${DIFF} -in1 ProteinQuantifier_9_output.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_9_output.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_9_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_9")
add_test("TOPP_ProteinQuantifier_9_out2" ${DIFF} -in1 ProteinQuantifier_9_output_mztab.tmp.mzTab -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_9_output_mztab.mzTab -whitelist "software")
set_tests_properties("TOPP_ProteinQuantifier_9_out2" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_9")

#ratios
add_test("TOPP_ProteinQuantifier_12" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.consensusXML -out ProteinQuantifier_12_output.tmp.csv -top:N 3 -top:include_all -consensus:fix_peptides -top:aggregate sum -ratios)
add_test("TOPP_ProteinQuantifier_12_out1" ${DIFF} -in1 ProteinQuantifier_12_output.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_12_output.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_12_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_12")
add_test("TOPP_ProteinQuantifier_13" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.consensusXML -out ProteinQuantifier_13_output.tmp.csv -top:N 3 -top:include_all -consensus:fix_peptides -top:aggregate sum -ratiosSILAC)
add_test("TOPP_ProteinQuantifier_13_out1" ${DIFF} -in1 ProteinQuantifier_13_output.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_13_output.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_13_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_13")
# idXML input (spectral counting):
add_test("TOPP_ProteinQuantifier_14" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.idXML -out ProteinQuantifier_14_output1.tmp.csv -peptide_out ProteinQuantifier_14_output2.tmp.csv -top:N 0 -top:aggregate sum)
add_test("TOPP_ProteinQuantifier_14_out1" ${DIFF} -in1 ProteinQuantifier_14_output1.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_14_output1.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_14_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_14")
add_test("TOPP_ProteinQuantifier_14_out2" ${DIFF} -in1 ProteinQuantifier_14_output2.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_14_output2.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_14_out2" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_14")
# iBAQ
add_test("TOPP_ProteinQuantifier_15" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_input.consensusXML -out ProteinQuantifier_15_output.tmp.csv -method iBAQ)
add_test("TOPP_ProteinQuantifier_15_out1" ${DIFF} -in1 ProteinQuantifier_15_output.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_15_output.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_15_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_15")
# fractionated, TMT example
add_test("TOPP_ProteinQuantifier_16" ${TOPP_BIN_PATH}/ProteinQuantifier -test -in ${DATA_DIR_TOPP}/ProteinQuantifier_16_input.consensusXML -design ${DATA_DIR_TOPP}/ProteinQuantifier_16_input.tsv -out ProteinQuantifier_16_output.tmp -top:include_all)
add_test("TOPP_ProteinQuantifier_16_out1" ${DIFF} -in1 ProteinQuantifier_16_output.tmp -in2 ${DATA_DIR_TOPP}/ProteinQuantifier_16_output.txt -whitelist "computed from file" )
set_tests_properties("TOPP_ProteinQuantifier_16_out1" PROPERTIES DEPENDS "TOPP_ProteinQuantifier_16")

#------------------------------------------------------------------------------
# PeptideIndexer tests
add_test("TOPP_PeptideIndexer_1" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_1.idXML -out PeptideIndexer_1_out.tmp.idXML -unmatched_action warn -enzyme:specificity none -aaa_max 4)
add_test("TOPP_PeptideIndexer_1_out" ${DIFF} -in1 PeptideIndexer_1_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/PeptideIndexer_1_out.idXML )
set_tests_properties("TOPP_PeptideIndexer_1_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_1")
add_test("TOPP_PeptideIndexer_2" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_1.idXML -out PeptideIndexer_2_out.tmp.idXML -unmatched_action warn -write_protein_sequence -enzyme:specificity none -aaa_max 4)
add_test("TOPP_PeptideIndexer_2_out" ${DIFF} -in1 PeptideIndexer_2_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/PeptideIndexer_2_out.idXML )
set_tests_properties("TOPP_PeptideIndexer_2_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_2")
add_test("TOPP_PeptideIndexer_3" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_1.idXML -out PeptideIndexer_3_out.tmp.idXML -unmatched_action warn -keep_unreferenced_proteins -enzyme:specificity none -aaa_max 4)
add_test("TOPP_PeptideIndexer_3_out" ${DIFF} -in1 PeptideIndexer_3_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/PeptideIndexer_3_out.idXML )
set_tests_properties("TOPP_PeptideIndexer_3_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_3")
add_test("TOPP_PeptideIndexer_4" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_1.idXML -out PeptideIndexer_4_out.tmp.idXML -unmatched_action warn -aaa_max 0 -write_protein_sequence -enzyme:specificity none -aaa_max 4)
add_test("TOPP_PeptideIndexer_4_out" ${DIFF} -in1 PeptideIndexer_4_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/PeptideIndexer_4_out.idXML )
set_tests_properties("TOPP_PeptideIndexer_4_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_4")
add_test("TOPP_PeptideIndexer_5" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_1.idXML -out PeptideIndexer_5_out.tmp.idXML -unmatched_action warn -enzyme:specificity none -aaa_max 4)
add_test("TOPP_PeptideIndexer_5_out" ${DIFF} -in1 PeptideIndexer_5_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/PeptideIndexer_5_out.idXML )
set_tests_properties("TOPP_PeptideIndexer_5_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_5")
add_test("TOPP_PeptideIndexer_6" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_2.idXML -out PeptideIndexer_6_out.tmp.idXML -unmatched_action warn -enzyme:specificity none -aaa_max 3)
add_test("TOPP_PeptideIndexer_6_out" ${DIFF} -in1 PeptideIndexer_6_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/PeptideIndexer_6_out.idXML )
set_tests_properties("TOPP_PeptideIndexer_6_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_6")
# specificity tests
add_test("TOPP_PeptideIndexer_7" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_3.idXML -out PeptideIndexer_7_out.tmp.idXML -unmatched_action warn -aaa_max 4 -enzyme:specificity full -decoy_string DECOY_)
add_test("TOPP_PeptideIndexer_7_out" ${DIFF} -in1 PeptideIndexer_7_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/PeptideIndexer_7_out.idXML )
set_tests_properties("TOPP_PeptideIndexer_7_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_7")
add_test("TOPP_PeptideIndexer_8" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_3.idXML -out PeptideIndexer_8_out.tmp.idXML -unmatched_action warn -enzyme:specificity semi -aaa_max 4)
add_test("TOPP_PeptideIndexer_8_out" ${DIFF} -in1 PeptideIndexer_8_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/PeptideIndexer_8_out.idXML )
set_tests_properties("TOPP_PeptideIndexer_8_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_8")
add_test("TOPP_PeptideIndexer_9" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_3.idXML -out PeptideIndexer_9_out.tmp.idXML -unmatched_action remove -enzyme:specificity none -aaa_max 4)
add_test("TOPP_PeptideIndexer_9_out" ${DIFF} -in1 PeptideIndexer_9_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/PeptideIndexer_9_out.idXML)
set_tests_properties("TOPP_PeptideIndexer_9_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_9")
# "IL_equivalent" option (and 1) make sure that I/L does not count towards ambiguous AAs; 2) that the original protein and peptide sequence (no I/L substitutions) are reported
add_test("TOPP_PeptideIndexer_10" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_10_input.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_10_input.idXML -out PeptideIndexer_10_output.tmp.idXML -IL_equivalent -aaa_max 3 -write_protein_sequence)
add_test("TOPP_PeptideIndexer_10_out" ${DIFF} -in1 PeptideIndexer_10_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/PeptideIndexer_10_output.idXML )
set_tests_properties("TOPP_PeptideIndexer_10_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_10")
# empty FASTA
add_test("TOPP_PeptideIndexer_11" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/degenerate_cases/empty.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_1.idXML -out PeptideIndexer_11_out.tmp.idXML)
set_tests_properties("TOPP_PeptideIndexer_11" PROPERTIES WILL_FAIL 1)
# with "annotate_proteins" (target/decoy):
add_test("TOPP_PeptideIndexer_12" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_1.idXML -out PeptideIndexer_12_out.tmp.idXML -unmatched_action warn -enzyme:specificity none -aaa_max 4)
add_test("TOPP_PeptideIndexer_12_out" ${DIFF} -in1 PeptideIndexer_12_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/PeptideIndexer_12_out.idXML)
set_tests_properties("TOPP_PeptideIndexer_12_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_12")
# empty input
add_test("TOPP_PeptideIndexer_13" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -in ${DATA_DIR_TOPP}/degenerate_cases/empty.idXML -out PeptideIndexer_13_out.tmp.idXML -aaa_max 4)
add_test("TOPP_PeptideIndexer_13_out" ${DIFF} -in1 PeptideIndexer_13_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/degenerate_cases/empty.idXML )
set_tests_properties("TOPP_PeptideIndexer_13_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_13")
# FASTA with X, overlapping peptides and some * in protein sequences (which should be removed internally)
add_test("TOPP_PeptideIndexer_14" ${TOPP_BIN_PATH}/PeptideIndexer -test -fasta ${DATA_DIR_TOPP}/PeptideIndexer_2.fasta -in ${DATA_DIR_TOPP}/PeptideIndexer_14.idXML -out PeptideIndexer_14_out.tmp.idXML -enzyme:specificity none -aaa_max 4 -write_protein_sequence)
add_test("TOPP_PeptideIndexer_14_out" ${DIFF} -in1 PeptideIndexer_14_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/PeptideIndexer_14_out.idXML )
set_tests_properties("TOPP_PeptideIndexer_14_out" PROPERTIES DEPENDS "TOPP_PeptideIndexer_14")

#------------------------------------------------------------------------------
# MzTabExporter tests
add_test("TOPP_MzTabExporter_1" ${TOPP_BIN_PATH}/MzTabExporter -test -in ${DATA_DIR_TOPP}/MzTabExporter_1_input.consensusXML -out MzTabExporter_1_output.tmp.mzTab)
add_test("TOPP_MzTabExporter_1_out" ${DIFF} -in1 MzTabExporter_1_output.tmp.mzTab -in2 ${DATA_DIR_TOPP}/MzTabExporter_1_output.mzTab -whitelist "MTD\turi[1]" "MS:1000752" "MTD\tsoftware")
set_tests_properties("TOPP_MzTabExporter_1_out" PROPERTIES DEPENDS "TOPP_MzTabExporter_1")

add_test("TOPP_MzTabExporter_2" ${TOPP_BIN_PATH}/MzTabExporter -test -in ${DATA_DIR_TOPP}/MzTabExporter_2_input.idXML -out MzTabExporter_2_output.tmp.mzTab)
add_test("TOPP_MzTabExporter_2_out" ${DIFF} -in1 MzTabExporter_2_output.tmp.mzTab -in2 ${DATA_DIR_TOPP}/MzTabExporter_2_output.mzTab -whitelist "MTD\turi[1]" "MS:1000752" "MTD\tsoftware")
set_tests_properties("TOPP_MzTabExporter_2_out" PROPERTIES DEPENDS "TOPP_MzTabExporter_2")

add_test("TOPP_MzTabExporter_3" ${TOPP_BIN_PATH}/MzTabExporter -test -in ${DATA_DIR_TOPP}/MzTabExporter_3_input.featureXML -out MzTabExporter_3_output.tmp.mzTab)
add_test("TOPP_MzTabExporter_3_out" ${DIFF} -in1 MzTabExporter_3_output.tmp.mzTab -in2 ${DATA_DIR_TOPP}/MzTabExporter_3_output.mzTab -whitelist "MTD\turi[1]" "MS:1000752" "MTD\tsoftware")
set_tests_properties("TOPP_MzTabExporter_3_out" PROPERTIES DEPENDS "TOPP_MzTabExporter_3")

add_test("TOPP_MzTabExporter_4" ${TOPP_BIN_PATH}/MzTabExporter -test -in ${DATA_DIR_TOPP}/Epifany_2_out.consensusXML -out MzTabExporter_4_output.tmp.mzTab)
add_test("TOPP_MzTabExporter_4_out" ${DIFF} -in1 MzTabExporter_4_output.tmp.mzTab -in2 ${DATA_DIR_TOPP}/MzTabExporter_4_output.mzTab -whitelist "MTD\turi[1]" "MS:1000752" "MTD\tsoftware" "protein_search_engine_score")
set_tests_properties("TOPP_MzTabExporter_4_out" PROPERTIES DEPENDS "TOPP_MzTabExporter_4")

add_test("TOPP_MzTabExporter_5" ${TOPP_BIN_PATH}/MzTabExporter -test -in ${DATA_DIR_TOPP}/MzTabExporter_5_in.consensusXML -out MzTabExporter_5_output.tmp.mzTab -first_run_inference_only)
add_test("TOPP_MzTabExporter_5_out" ${DIFF} -in1 MzTabExporter_5_output.tmp.mzTab -in2 ${DATA_DIR_TOPP}/MzTabExporter_5_output.mzTab -whitelist "MTD\turi[1]" "MS:1000752" "MTD\tsoftware" "protein_search_engine_score")
set_tests_properties("TOPP_MzTabExporter_5_out" PROPERTIES DEPENDS "TOPP_MzTabExporter_5")

add_test("TOPP_MzTabExporter_6" ${TOPP_BIN_PATH}/MzTabExporter -test -in ${DATA_DIR_TOPP}/MzTabExporter_6_input.idXML -out MzTabExporter_6_output.tmp.mzTab)
add_test("TOPP_MzTabExporter_6_out" ${DIFF} -in1 MzTabExporter_6_output.tmp.mzTab -in2 ${DATA_DIR_TOPP}/MzTabExporter_6_output.mzTab -whitelist "MTD\turi[1]" "MS:1000752" "MTD\tsoftware" "protein_search_engine_score")
set_tests_properties("TOPP_MzTabExporter_6_out" PROPERTIES DEPENDS "TOPP_MzTabExporter_6")

add_test("TOPP_MzTabExporter_7" ${TOPP_BIN_PATH}/MzTabExporter -test -in ${DATA_DIR_TOPP}/MzTabExporter_7_input.consensusXML -out MzTabExporter_7_output.tmp.mzTab)
add_test("TOPP_MzTabExporter_7_out" ${DIFF} -in1 MzTabExporter_7_output.tmp.mzTab -in2 ${DATA_DIR_TOPP}/MzTabExporter_7_output.mzTab -whitelist "MTD\turi[1]" "MS:1000752" "MTD\tsoftware" "protein_search_engine_score")
set_tests_properties("TOPP_MzTabExporter_7_out" PROPERTIES DEPENDS "TOPP_MzTabExporter_7")

add_test("TOPP_MzTabExporter_8" ${TOPP_BIN_PATH}/MzTabExporter -test -in ${DATA_DIR_TOPP}/MzTabExporter_6_input.idXML -export_all_psms -out MzTabExporter_8_output.tmp.mzTab)
add_test("TOPP_MzTabExporter_8_out" ${DIFF} -in1 MzTabExporter_8_output.tmp.mzTab -in2 ${DATA_DIR_TOPP}/MzTabExporter_8_output.mzTab -whitelist "MTD\turi[1]" "MS:1000752" "MTD\tsoftware" "protein_search_engine_score")
set_tests_properties("TOPP_MzTabExporter_8_out" PROPERTIES DEPENDS "TOPP_MzTabExporter_8")
#------------------------------------------------------------------------------
# OpenPepXL test
add_test("TOPP_OpenPepXL_1" ${TOPP_BIN_PATH}/OpenPepXL -test -ini ${DATA_DIR_TOPP}/OpenPepXL_input.ini -in ${DATA_DIR_TOPP}/OpenPepXL_input.mzML -consensus ${DATA_DIR_TOPP}/OpenPepXL_input.consensusXML -database ${DATA_DIR_TOPP}/OpenPepXL_input.fasta -out_xquestxml OpenPepXL_output.tmp.xquest.xml -out_xquest_specxml OpenPepXL_output.tmp.spec.xml -out_mzIdentML OpenPepXL_output.tmp.mzid -out_idXML OpenPepXL_output.tmp.idXML)
add_test("TOPP_OpenPepXL_1_out_1" ${DIFF} -whitelist "date=" -in1 OpenPepXL_output.tmp.xquest.xml -in2 ${DATA_DIR_TOPP}/OpenPepXL_output.xquest.xml )
# unstable due to differences in double to string conversion, caused by different boost versions (between 1.71.0 and 1.78.0)
# this file format is only important for compatibility with the xQuest pipeline's result visualization and not essential
#add_test("TOPP_OpenPepXL_1_out_2" ${DIFF} -in1 OpenPepXL_output.tmp.spec.xml -in2 ${DATA_DIR_TOPP}/OpenPepXL_output.spec.xml )
add_test("TOPP_OpenPepXL_1_out_3" ${DIFF} -whitelist "creationDate=" "id=" "spectraData_ref=" "searchDatabase_ref=" "OpenPepXL_input" -in1 OpenPepXL_output.tmp.mzid -in2 ${DATA_DIR_TOPP}/OpenPepXL_output.mzid )
add_test("TOPP_OpenPepXL_1_out_4" ${DIFF} -whitelist "db=" "input_consensusXML" "input_mzML" "date" "OpenPepXL_input" -in1 OpenPepXL_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/OpenPepXL_output.idXML )
set_tests_properties("TOPP_OpenPepXL_1_out_1" PROPERTIES DEPENDS "TOPP_OpenPepXL_1")
#set_tests_properties("TOPP_OpenPepXL_1_out_2" PROPERTIES DEPENDS "TOPP_OpenPepXL_1")
set_tests_properties("TOPP_OpenPepXL_1_out_3" PROPERTIES DEPENDS "TOPP_OpenPepXL_1")
set_tests_properties("TOPP_OpenPepXL_1_out_4" PROPERTIES DEPENDS "TOPP_OpenPepXL_1")

#------------------------------------------------------------------------------
# OpenPepXLLF test
add_test("TOPP_OpenPepXLLF_1" ${TOPP_BIN_PATH}/OpenPepXLLF -test -decoy_string "decoy" -in ${DATA_DIR_TOPP}/OpenPepXLLF_input.mzML -database ${DATA_DIR_TOPP}/OpenPepXLLF_input.fasta -out_xquestxml OpenPepXLLF_output.tmp.xquest.xml -out_xquest_specxml OpenPepXLLF_output.tmp.spec.xml -out_mzIdentML OpenPepXLLF_output.tmp.mzid -out_idXML OpenPepXLLF_output.tmp.idXML)
add_test("TOPP_OpenPepXLLF_1_out_1" ${DIFF} -whitelist "date=" -in1 OpenPepXLLF_output.tmp.xquest.xml -in2 ${DATA_DIR_TOPP}/OpenPepXLLF_output.xquest.xml )
# unstable due to differences in double to string conversion, caused by different boost versions (between 1.71.0 and 1.78.0)
# this file format is only important for compatibility with the xQuest pipeline's result visualization and not essential
#add_test("TOPP_OpenPepXLLF_1_out_2" ${DIFF} -in1 OpenPepXLLF_output.tmp.spec.xml -in2 ${DATA_DIR_TOPP}/OpenPepXLLF_output.spec.xml )
add_test("TOPP_OpenPepXLLF_1_out_3" ${DIFF} -whitelist "creationDate=" "id=" "spectraData_ref=" "searchDatabase_ref=" "input_mzML" "OpenPepXLLF_input" -in1 OpenPepXLLF_output.tmp.mzid -in2 ${DATA_DIR_TOPP}/OpenPepXLLF_output.mzid )
add_test("TOPP_OpenPepXLLF_1_out_4" ${DIFF} -whitelist "db=" "input_consensusXML" "input_mzML" "date" "OpenPepXLLF_input" -in1 OpenPepXLLF_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/OpenPepXLLF_output.idXML )
set_tests_properties("TOPP_OpenPepXLLF_1_out_1" PROPERTIES DEPENDS "TOPP_OpenPepXLLF_1")
#set_tests_properties("TOPP_OpenPepXLLF_1_out_2" PROPERTIES DEPENDS "TOPP_OpenPepXLLF_1")
set_tests_properties("TOPP_OpenPepXLLF_1_out_3" PROPERTIES DEPENDS "TOPP_OpenPepXLLF_1")
set_tests_properties("TOPP_OpenPepXLLF_1_out_4" PROPERTIES DEPENDS "TOPP_OpenPepXLLF_1")
add_test("TOPP_OpenPepXLLF_2" ${TOPP_BIN_PATH}/OpenPepXLLF -test -ini ${DATA_DIR_TOPP}/OpenPepXLLF_input2.ini -in ${DATA_DIR_TOPP}/OpenPepXLLF_input2.mzML -database ${DATA_DIR_TOPP}/OpenPepXLLF_input2.fasta -out_idXML OpenPepXLLF_output2.tmp.idXML)
add_test("TOPP_OpenPepXLLF_2_out_1" ${DIFF} -whitelist "db=" "input_consensusXML" "input_mzML" "date" "OpenPepXLLF_input" -in1 OpenPepXLLF_output2.tmp.idXML -in2 ${DATA_DIR_TOPP}/OpenPepXLLF_output2.idXML )
set_tests_properties("TOPP_OpenPepXLLF_2_out_1" PROPERTIES DEPENDS "TOPP_OpenPepXLLF_2")

#------------------------------------------------------------------------------
# XFDR tests
add_test("TOPP_XFDR_1" ${TOPP_BIN_PATH}/XFDR -test -binsize 0.1 -in ${DATA_DIR_TOPP}/XFDR_test_in1.idXML -out_idXML XFDR_test_out1_temp.idXML -out_mzIdentML XFDR_test_out1_temp.mzid -out_xquest XFDR_test_out1_temp.xquest.xml)
add_test("TOPP_XFDR_1_out11" ${DIFF} -whitelist "date=" -in1 XFDR_test_out1_temp.idXML  -in2 ${DATA_DIR_TOPP}/XFDR_test_out1.idXML)
add_test("TOPP_XFDR_1_out12" ${DIFF} -whitelist "creationDate=" "version=" -in1 XFDR_test_out1_temp.mzid -in2 ${DATA_DIR_TOPP}/XFDR_test_out1.mzid)
add_test("TOPP_XFDR_1_out13" ${DIFF} -whitelist "date=" -in1 XFDR_test_out1_temp.xquest.xml -in2 ${DATA_DIR_TOPP}/XFDR_test_out1.xquest.xml)
set_tests_properties("TOPP_XFDR_1_out11" PROPERTIES DEPENDS "TOPP_XFDR_1")
set_tests_properties("TOPP_XFDR_1_out12" PROPERTIES DEPENDS "TOPP_XFDR_1")
set_tests_properties("TOPP_XFDR_1_out13" PROPERTIES DEPENDS "TOPP_XFDR_1")
add_test("TOPP_XFDR_2" ${TOPP_BIN_PATH}/XFDR -test -binsize 0.1 -in ${DATA_DIR_TOPP}/XFDR_test_in1.idXML -uniquexl -no_qvalues -out_idXML XFDR_test_out2_temp.idXML -out_mzIdentML XFDR_test_out2_temp.mzid -out_xquest XFDR_test_out2_temp.xquest.xml)
add_test("TOPP_XFDR_2_out21" ${DIFF} -whitelist "date=" -in1 XFDR_test_out2_temp.idXML -in2 ${DATA_DIR_TOPP}/XFDR_test_out2.idXML)
add_test("TOPP_XFDR_2_out22" ${DIFF} -whitelist "creationDate=" "version=" "id=\"OpenMS_" -in1 XFDR_test_out2_temp.mzid -in2 ${DATA_DIR_TOPP}/XFDR_test_out2.mzid)
add_test("TOPP_XFDR_2_out23" ${DIFF} -whitelist "date=" -in1 XFDR_test_out2_temp.xquest.xml -in2 ${DATA_DIR_TOPP}/XFDR_test_out2.xquest.xml)
set_tests_properties("TOPP_XFDR_2_out21" PROPERTIES DEPENDS "TOPP_XFDR_2")
set_tests_properties("TOPP_XFDR_2_out22" PROPERTIES DEPENDS "TOPP_XFDR_2")
set_tests_properties("TOPP_XFDR_2_out23" PROPERTIES DEPENDS "TOPP_XFDR_2")
add_test("TOPP_XFDR_3" ${TOPP_BIN_PATH}/XFDR -test -binsize 0.1 -minscore 0.0 -in ${DATA_DIR_TOPP}/XFDR_test_in2.xquest.xml -out_idXML XFDR_test_out3_temp.idXML -out_mzIdentML XFDR_test_out3_temp.mzid -out_xquest XFDR_test_out3_temp.xquest.xml)
add_test("TOPP_XFDR_3_out31" ${DIFF} -whitelist "date=" -in1 XFDR_test_out3_temp.idXML -in2 ${DATA_DIR_TOPP}/XFDR_test_out3.idXML)
add_test("TOPP_XFDR_3_out32" ${DIFF} -whitelist "creationDate=" "version=" "id=\"OpenMS_" -in1 XFDR_test_out3_temp.mzid -in2 ${DATA_DIR_TOPP}/XFDR_test_out3.mzid)
add_test("TOPP_XFDR_3_out33" ${DIFF} -whitelist "date=" -in1 XFDR_test_out3_temp.xquest.xml -in2 ${DATA_DIR_TOPP}/XFDR_test_out3.xquest.xml)
set_tests_properties("TOPP_XFDR_3_out31" PROPERTIES DEPENDS "TOPP_XFDR_3")
set_tests_properties("TOPP_XFDR_3_out32" PROPERTIES DEPENDS "TOPP_XFDR_3")
set_tests_properties("TOPP_XFDR_3_out33" PROPERTIES DEPENDS "TOPP_XFDR_3")
add_test("TOPP_XFDR_4" ${TOPP_BIN_PATH}/XFDR -test -binsize 0.1 -minborder -3 -maxborder 3 -in ${DATA_DIR_TOPP}/XFDR_test_in2.xquest.xml -out_idXML XFDR_test_out4_temp.idXML -out_mzIdentML XFDR_test_out4_temp.mzid -out_xquest XFDR_test_out4_temp.xquest.xml)
add_test("TOPP_XFDR_4_out41" ${DIFF} -whitelist "date=" -in1 XFDR_test_out4_temp.idXML -in2 ${DATA_DIR_TOPP}/XFDR_test_out4.idXML)
add_test("TOPP_XFDR_4_out42" ${DIFF} -whitelist "creationDate=" "version=" "id=\"OpenMS_" -in1 XFDR_test_out4_temp.mzid -in2 ${DATA_DIR_TOPP}/XFDR_test_out4.mzid)
add_test("TOPP_XFDR_4_out43" ${DIFF} -whitelist "date=" -in1 XFDR_test_out4_temp.xquest.xml -in2 ${DATA_DIR_TOPP}/XFDR_test_out4.xquest.xml)
set_tests_properties("TOPP_XFDR_4_out41" PROPERTIES DEPENDS "TOPP_XFDR_4")
set_tests_properties("TOPP_XFDR_4_out42" PROPERTIES DEPENDS "TOPP_XFDR_4")
set_tests_properties("TOPP_XFDR_4_out43" PROPERTIES DEPENDS "TOPP_XFDR_4")
add_test("TOPP_XFDR_5" ${TOPP_BIN_PATH}/XFDR -test -binsize 0.1 -minborder -3 -maxborder 3 -in ${DATA_DIR_TOPP}/XFDR_test_in3.xquest.xml -out_idXML XFDR_test_out5_temp.idXML -out_mzIdentML XFDR_test_out5_temp.mzid -out_xquest XFDR_test_out5_temp.xquest.xml)
add_test("TOPP_XFDR_5_out51" ${DIFF} -whitelist "date=" -in1 XFDR_test_out5_temp.idXML -in2 ${DATA_DIR_TOPP}/XFDR_test_out5.idXML)
add_test("TOPP_XFDR_5_out52" ${DIFF} -whitelist "creationDate=" "version=" "id=\"OpenMS_" -in1 XFDR_test_out5_temp.mzid -in2 ${DATA_DIR_TOPP}/XFDR_test_out5.mzid)
add_test("TOPP_XFDR_5_out53" ${DIFF} -whitelist "date=" -in1 XFDR_test_out5_temp.xquest.xml -in2 ${DATA_DIR_TOPP}/XFDR_test_out5.xquest.xml)
set_tests_properties("TOPP_XFDR_5_out51" PROPERTIES DEPENDS "TOPP_XFDR_5")
set_tests_properties("TOPP_XFDR_5_out52" PROPERTIES DEPENDS "TOPP_XFDR_5")
set_tests_properties("TOPP_XFDR_5_out53" PROPERTIES DEPENDS "TOPP_XFDR_5")
add_test("TOPP_XFDR_6" ${TOPP_BIN_PATH}/XFDR -test -binsize 0.1 -in ${DATA_DIR_TOPP}/XFDR_test_in4.idXML -out_idXML XFDR_test_out6_temp.idXML)
add_test("TOPP_XFDR_6_out6" ${DIFF} -whitelist "date=" -in1 XFDR_test_out6_temp.idXML -in2 ${DATA_DIR_TOPP}/XFDR_test_out6.idXML)
set_tests_properties("TOPP_XFDR_6_out6" PROPERTIES DEPENDS "TOPP_XFDR_6")
add_test("TOPP_XFDR_7" ${TOPP_BIN_PATH}/XFDR -test -uniquexl -binsize 0.1 -minscore 0.0 -in ${DATA_DIR_TOPP}/XFDR_test_in2.xquest.xml -out_idXML XFDR_test_out7_temp.idXML -out_mzIdentML XFDR_test_out7_temp.mzid -out_xquest XFDR_test_out7_temp.xquest.xml)
add_test("TOPP_XFDR_7_out31" ${DIFF} -whitelist "date=" -in1 XFDR_test_out7_temp.idXML -in2 ${DATA_DIR_TOPP}/XFDR_test_out7.idXML)
add_test("TOPP_XFDR_7_out32" ${DIFF} -whitelist "creationDate=" "version=" "id=\"OpenMS_" -in1 XFDR_test_out7_temp.mzid -in2 ${DATA_DIR_TOPP}/XFDR_test_out7.mzid)
add_test("TOPP_XFDR_7_out33" ${DIFF} -whitelist "date=" -in1 XFDR_test_out7_temp.xquest.xml -in2 ${DATA_DIR_TOPP}/XFDR_test_out7.xquest.xml)
set_tests_properties("TOPP_XFDR_7_out31" PROPERTIES DEPENDS "TOPP_XFDR_7")
set_tests_properties("TOPP_XFDR_7_out32" PROPERTIES DEPENDS "TOPP_XFDR_7")
set_tests_properties("TOPP_XFDR_7_out33" PROPERTIES DEPENDS "TOPP_XFDR_7")

#------------------------------------------------------------------------------
# QualityControl tests
add_test("TOPP_QualityControl_1" ${TOPP_BIN_PATH}/QualityControl -test -in_raw ${DATA_DIR_TOPP}/QualityControl_1_in1.mzML.gz ${DATA_DIR_TOPP}/QualityControl_1_in2.mzML.gz ${DATA_DIR_TOPP}/QualityControl_1_in3.mzML.gz -in_postFDR ${DATA_DIR_TOPP}/QualityControl_1_in1.featureXML ${DATA_DIR_TOPP}/QualityControl_1_in2.featureXML ${DATA_DIR_TOPP}/QualityControl_1_in3.featureXML -in_trafo ${DATA_DIR_TOPP}/QualityControl_1_in1.trafoXML ${DATA_DIR_TOPP}/QualityControl_1_in2.trafoXML ${DATA_DIR_TOPP}/QualityControl_1_in3.trafoXML  -in_contaminants ${DATA_DIR_TOPP}/QualityControl_1.fasta -in_cm ${DATA_DIR_TOPP}/QualityControl_1_in.consensusXML -in_fasta ${DATA_DIR_TOPP}/PeptideIndexer_1.fasta -out_cm QualityControl_1_out.tmp.consensusXML -out QualityControl_1_out.tmp.mzTab)
add_test("TOPP_QualityControl_1_out_consXML" ${DIFF} -in1 QualityControl_1_out.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/QualityControl_1_out.consensusXML -whitelist "href=")
add_test("TOPP_QualityControl_1_out_mzTab" ${DIFF} -in1 QualityControl_1_out.tmp.mzTab -in2 ${DATA_DIR_TOPP}/QualityControl_1_out.mzTab -whitelist "MS:1000752")
set_tests_properties("TOPP_QualityControl_1_out_consXML" PROPERTIES DEPENDS "TOPP_QualityControl_1")
set_tests_properties("TOPP_QualityControl_1_out_mzTab" PROPERTIES DEPENDS "TOPP_QualityControl_1")

## Tests for TOPPAS/ExecutePipeline are in ./src/tests/TOPPAS/...

#------------------------------------------------------------------------------
# search engine tests (might be optional, depending on presence of external binaries, e.g., omssacl.exe)
#------------------------------------------------------------------------------
include(THIRDPARTY/third_party_tests.cmake)

# AccurateMassSearch (with default params, which forces lookup of HMDB DB in share/)
add_test("TOPP_AccurateMassSearch_1" ${TOPP_BIN_PATH}/AccurateMassSearch -test -in ${DATA_DIR_TOPP}/ConsensusMapNormalizer_input.consensusXML -out AccurateMassSearch_1_output.tmp.mzTab)
add_test("TOPP_AccurateMassSearch_1_out1" ${DIFF} -in1 AccurateMassSearch_1_output.tmp.mzTab -in2 ${DATA_DIR_TOPP}/AccurateMassSearch_1_output.csv -whitelist "")
set_tests_properties("TOPP_AccurateMassSearch_1_out1" PROPERTIES DEPENDS "TOPP_AccurateMassSearch_1")
# use FFM and MAD featurexml with custom adducts and database
# mztab, featureXML (with isotope intensities)
add_test("TOPP_AccurateMassSearch_2" ${TOPP_BIN_PATH}/AccurateMassSearch -test -in ${DATA_DIR_TOPP}/AccurateMassSearch_2_input.featureXML -out AccurateMassSearch_2_output.tmp.mzTab -out_annotation AccurateMassSearch_2_output.tmp.featureXML -db:mapping ${DATA_DIR_TOPP}/AMS_test_Mapping.tsv -db:struct ${DATA_DIR_TOPP}/AMS_test_Struct.tsv -positive_adducts ${DATA_DIR_TOPP}/AMS_PositiveAdducts.tsv -negative_adducts ${DATA_DIR_TOPP}/AMS_NegativeAdducts.tsv -algorithm:mzTab:exportIsotopeIntensities true -algorithm:keep_unidentified_masses false)
add_test("TOPP_AccurateMassSearch_2_out1" ${DIFF} -in1 AccurateMassSearch_2_output.tmp.mzTab -in2 ${DATA_DIR_TOPP}/AccurateMassSearch_2_output.mzTab -whitelist "")
add_test("TOPP_AccurateMassSearch_2_out2" ${DIFF} -in1 AccurateMassSearch_2_output.tmp.featureXML -in2 ${DATA_DIR_TOPP}/AccurateMassSearch_2_output.featureXML -whitelist "IdentificationRun id=\"PI_0\" date=" "database_location")
set_tests_properties("TOPP_AccurateMassSearch_2_out1" PROPERTIES DEPENDS "TOPP_AccurateMassSearch_2")
set_tests_properties("TOPP_AccurateMassSearch_2_out2" PROPERTIES DEPENDS "TOPP_AccurateMassSearch_2")
# mztab (default)
add_test("TOPP_AccurateMassSearch_3" ${TOPP_BIN_PATH}/AccurateMassSearch -test -in ${DATA_DIR_TOPP}/AccurateMassSearch_2_input.featureXML -out AccurateMassSearch_3_output.tmp.mzTab -out_annotation AccurateMassSearch_2_output.tmp.featureXML -db:mapping ${DATA_DIR_TOPP}/AMS_test_Mapping.tsv -db:struct ${DATA_DIR_TOPP}/AMS_test_Struct.tsv -positive_adducts ${DATA_DIR_TOPP}/AMS_PositiveAdducts.tsv -negative_adducts ${DATA_DIR_TOPP}/AMS_NegativeAdducts.tsv)
add_test("TOPP_AccurateMassSearch_3_out1" ${DIFF} -in1 AccurateMassSearch_3_output.tmp.mzTab -in2 ${DATA_DIR_TOPP}/AccurateMassSearch_3_output.mzTab -whitelist "")
set_tests_properties("TOPP_AccurateMassSearch_3_out1" PROPERTIES DEPENDS "TOPP_AccurateMassSearch_3")
# TODO: How to test the oms file?
# oms (default)
#add_test("TOPP_AccurateMassSearch_4" ${TOPP_BIN_PATH}/AccurateMassSearch -test -in ${DATA_DIR_TOPP}/AccurateMassSearch_2_input.featureXML -out AccurateMassSearch_4_output.tmp.mzTab -out_annotation AccurateMassSearch_2_output.tmp.oms -db:mapping ${DATA_DIR_TOPP}/AMS_test_Mapping.tsv -db:struct ${DATA_DIR_TOPP}/AMS_test_Struct.tsv -positive_adducts ${DATA_DIR_TOPP}/AMS_PositiveAdducts.tsv -negative_adducts ${DATA_DIR_TOPP}/AMS_NegativeAdducts.tsv)
#add_test("TOPP_AccurateMassSearch_4_out1" ${DIFF} -in1 AccurateMassSearch_4_output.tmp.oms -in2 ${DATA_DIR_TOPP}/AccurateMassSearch_3_output.mzTab -whitelist "")
#set_tests_properties("TOPP_AccurateMassSearch_4_out1" PROPERTIES DEPENDS "TOPP_AccurateMassSearch_4")
# mztab_m + oms (annotation)
add_test("TOPP_AccurateMassSearch_5" ${TOPP_BIN_PATH}/AccurateMassSearch -test -in ${DATA_DIR_TOPP}/AccurateMassSearch_2_input.featureXML -out AccurateMassSearch_5_output.tmp.mzTab -out_annotation AccurateMassSearch_5_output.tmp.oms -db:mapping ${DATA_DIR_TOPP}/AMS_test_Mapping.tsv -db:struct ${DATA_DIR_TOPP}/AMS_test_Struct.tsv -positive_adducts ${DATA_DIR_TOPP}/AMS_PositiveAdducts.tsv -negative_adducts ${DATA_DIR_TOPP}/AMS_NegativeAdducts.tsv -algorithm:id_format ID)
add_test("TOPP_AccurateMassSearch_5_out1" ${DIFF} -in1 AccurateMassSearch_5_output.tmp.mzTab -in2 ${DATA_DIR_TOPP}/AccurateMassSearch_5_output.mzTab -whitelist "")
set_tests_properties("TOPP_AccurateMassSearch_5_out1" PROPERTIES DEPENDS "TOPP_AccurateMassSearch_5")
# mztab_m + featureXML (annotation)
add_test("TOPP_AccurateMassSearch_6" ${TOPP_BIN_PATH}/AccurateMassSearch -test -in ${DATA_DIR_TOPP}/AccurateMassSearch_2_input.featureXML -out AccurateMassSearch_6_output.tmp.mzTab -out_annotation AccurateMassSearch_6_output.tmp.featureXML -db:mapping ${DATA_DIR_TOPP}/AMS_test_Mapping.tsv -db:struct ${DATA_DIR_TOPP}/AMS_test_Struct.tsv -positive_adducts ${DATA_DIR_TOPP}/AMS_PositiveAdducts.tsv -negative_adducts ${DATA_DIR_TOPP}/AMS_NegativeAdducts.tsv -algorithm:id_format ID)
add_test("TOPP_AccurateMassSearch_6_out1" ${DIFF} -in1 AccurateMassSearch_6_output.tmp.mzTab -in2 ${DATA_DIR_TOPP}/AccurateMassSearch_6_output.mzTab -whitelist "")
add_test("TOPP_AccurateMassSearch_6_out2" ${DIFF} -in1 AccurateMassSearch_6_output.tmp.featureXML -in2 ${DATA_DIR_TOPP}/AccurateMassSearch_6_output.featureXML -whitelist "IdentificationRun id=\"PI_0\" date=" "database_location")
set_tests_properties("TOPP_AccurateMassSearch_6_out1" PROPERTIES DEPENDS "TOPP_AccurateMassSearch_6")
set_tests_properties("TOPP_AccurateMassSearch_6_out2" PROPERTIES DEPENDS "TOPP_AccurateMassSearch_6")

# SiriusExport
# without featureinfo
add_test("TOPP_SiriusExport_1" ${TOPP_BIN_PATH}/SiriusExport -test -in ${DATA_DIR_TOPP}/SiriusExport_1_input.mzML -out SiriusExport_1_output.tmp.ms)
add_test("TOPP_SiriusExport_1_out" ${DIFF} -in1 SiriusExport_1_output.tmp.ms -in2 ${DATA_DIR_TOPP}/SiriusExport_1_output.ms)
set_tests_properties("TOPP_SiriusExport_1_out" PROPERTIES DEPENDS "TOPP_SiriusExport_1")

# with featureinfo, feature_only and compound info table (input for AssayGeneratorMetaboSirius tests)
add_test("TOPP_SiriusExport_2" ${TOPP_BIN_PATH}/SiriusExport -test -in ${DATA_DIR_TOPP}/AssayGeneratorMetabo_input.mzML -in_featureinfo ${DATA_DIR_TOPP}/AssayGeneratorMetabo_ams_input.featureXML -out SiriusExport_2_output.tmp.ms -out_compoundinfo SiriusExport_2_output_compoundinfo.tmp.tsv -feature_only true)
add_test("TOPP_SiriusExport_2_out1" ${DIFF} -in1 SiriusExport_2_output.tmp.ms -in2 ${DATA_DIR_TOPP}/SiriusExport_2_output.ms)
add_test("TOPP_SiriusExport_2_out2" ${DIFF} -in1 SiriusExport_2_output_compoundinfo.tmp.tsv -in2 ${DATA_DIR_TOPP}/SiriusExport_2_output_compoundinfo.tsv)
set_tests_properties("TOPP_SiriusExport_2_out1" PROPERTIES DEPENDS "TOPP_SiriusExport_2")
set_tests_properties("TOPP_SiriusExport_2_out2" PROPERTIES DEPENDS "TOPP_SiriusExport_2")

# with multiple files
add_test("TOPP_SiriusExport_3" ${TOPP_BIN_PATH}/SiriusExport -test -in ${DATA_DIR_TOPP}/SiriusExport_1_input.mzML ${DATA_DIR_TOPP}/AssayGeneratorMetabo_input.mzML -in_featureinfo ${DATA_DIR_TOPP}/SiriusExport_1_input.featureXML ${DATA_DIR_TOPP}/AssayGeneratorMetabo_ams_input.featureXML -out SiriusExport_3_output.tmp.ms -out_compoundinfo SiriusExport_3_output_compoundinfo.tmp.tsv)
add_test("TOPP_SiriusExport_3_out1" ${DIFF} -in1 SiriusExport_3_output.tmp.ms -in2 ${DATA_DIR_TOPP}/SiriusExport_3_output.ms)
add_test("TOPP_SiriusExport_3_out2" ${DIFF} -in1 SiriusExport_3_output_compoundinfo.tmp.tsv -in2 ${DATA_DIR_TOPP}/SiriusExport_3_output_compoundinfo.tsv)
set_tests_properties("TOPP_SiriusExport_3_out1" PROPERTIES DEPENDS "TOPP_SiriusExport_3")
set_tests_properties("TOPP_SiriusExport_3_out2" PROPERTIES DEPENDS "TOPP_SiriusExport_3")

# AssayGeneratorMetaboSirius
# SIRIUS project directory generated with SIRIUS 5.8.5 from SiriusExport_2 test .ms file, minimal example with one UNKNOWN 2 two named compounds from which only Myclobutanil will pass filtering
add_test("TOPP_AssayGeneratorMetaboSirius_1" ${TOPP_BIN_PATH}/AssayGeneratorMetaboSirius -test -in ${DATA_DIR_TOPP}/AssayGeneratorMetabo-in-SIRIUS-project -in_compoundinfo ${DATA_DIR_TOPP}/SiriusExport_2_output_compoundinfo.tsv -out AssayGeneratorMetaboSirius_output_1.tmp.tsv -exclude_ms2_precursor -max_transitions 3)
add_test("TOPP_AssayGeneratorMetaboSirius_1_out" ${DIFF} -in1 AssayGeneratorMetaboSirius_output_1.tmp.tsv -in2 ${DATA_DIR_TOPP}/AssayGeneratorMetaboSirius_output_1.tsv)
set_tests_properties("TOPP_AssayGeneratorMetaboSirius_1_out" PROPERTIES DEPENDS "TOPP_AssayGeneratorMetaboSirius_1")

# AssayGeneratorMetaboSirius: use_known_unknowns
add_test("TOPP_AssayGeneratorMetaboSirius_2" ${TOPP_BIN_PATH}/AssayGeneratorMetaboSirius -test -in ${DATA_DIR_TOPP}/AssayGeneratorMetabo-in-SIRIUS-project -in_compoundinfo ${DATA_DIR_TOPP}/SiriusExport_2_output_compoundinfo.tsv -out AssayGeneratorMetaboSirius_output_2.tmp.tsv -exclude_ms2_precursor -max_transitions 3 -use_known_unknowns)
add_test("TOPP_AssayGeneratorMetaboSirius_2_out" ${DIFF} -in1 AssayGeneratorMetaboSirius_output_2.tmp.tsv -in2 ${DATA_DIR_TOPP}/AssayGeneratorMetaboSirius_output_2.tsv)
set_tests_properties("TOPP_AssayGeneratorMetaboSirius_2_out" PROPERTIES DEPENDS "TOPP_AssayGeneratorMetaboSirius_2")

# AssayGeneratorMetaboSirius: decoys resolving overlaps
add_test("TOPP_AssayGeneratorMetaboSirius_3" ${TOPP_BIN_PATH}/AssayGeneratorMetaboSirius -test -in ${DATA_DIR_TOPP}/AssayGeneratorMetabo-in-SIRIUS-project -in_compoundinfo ${DATA_DIR_TOPP}/SiriusExport_2_output_compoundinfo.tsv -out AssayGeneratorMetaboSirius_output_3.tmp.tsv -exclude_ms2_precursor -max_transitions 3 -decoy_generation -decoy_generation_method resolve_overlap)
add_test("TOPP_AssayGeneratorMetaboSirius_3_out" ${DIFF} -in1 AssayGeneratorMetaboSirius_output_3.tmp.tsv -in2 ${DATA_DIR_TOPP}/AssayGeneratorMetaboSirius_output_3.tsv)
set_tests_properties("TOPP_AssayGeneratorMetaboSirius_3_out" PROPERTIES DEPENDS "TOPP_AssayGeneratorMetaboSirius_3")

# AssayGeneratorMetabo
# use FeatureFinderMetabo data
add_test("TOPP_AssayGeneratorMetabo_1" ${TOPP_BIN_PATH}/AssayGeneratorMetabo -test -in ${DATA_DIR_TOPP}/AssayGeneratorMetabo_input.mzML -in_featureinfo ${DATA_DIR_TOPP}/AssayGeneratorMetabo_ffm_input.featureXML -out AssayGeneratorMetabo_ffm_output.tmp.tsv -min_transitions 1 -max_transitions 3)
add_test("TOPP_AssayGeneratorMetabo_1_out1" ${DIFF} -in1 AssayGeneratorMetabo_ffm_output.tmp.tsv -in2 ${DATA_DIR_TOPP}/AssayGeneratorMetabo_ffm_output.tsv)
set_tests_properties("TOPP_AssayGeneratorMetabo_1_out1" PROPERTIES DEPENDS "TOPP_AssayGeneratorMetabo_1")

# use AccurateMassSearch data
add_test("TOPP_AssayGeneratorMetabo_2" ${TOPP_BIN_PATH}/AssayGeneratorMetabo -test -in ${DATA_DIR_TOPP}/AssayGeneratorMetabo_input.mzML -in_featureinfo ${DATA_DIR_TOPP}/AssayGeneratorMetabo_ams_input.featureXML -out AssayGeneratorMetabo_ams_output.tmp.tsv -min_transitions 1 -max_transitions 3)
add_test("TOPP_AssayGeneratorMetabo_2_out1" ${DIFF} -in1 AssayGeneratorMetabo_ams_output.tmp.tsv -in2 ${DATA_DIR_TOPP}/AssayGeneratorMetabo_ams_output.tsv)
set_tests_properties("TOPP_AssayGeneratorMetabo_2_out1" PROPERTIES DEPENDS "TOPP_AssayGeneratorMetabo_2")

# use FeatureFinderMetabo data with method consenesus_spectrum
add_test("TOPP_AssayGeneratorMetabo_3" ${TOPP_BIN_PATH}/AssayGeneratorMetabo -test -in ${DATA_DIR_TOPP}/AssayGeneratorMetabo_input.mzML -in_featureinfo ${DATA_DIR_TOPP}/AssayGeneratorMetabo_ffm_input.featureXML -out AssayGeneratorMetabo_ffm_output_consensus.tmp.tsv -method consensus_spectrum -min_transitions 1 -max_transitions 3)
add_test("TOPP_AssayGeneratorMetabo_3_out1" ${DIFF} -in1 AssayGeneratorMetabo_ffm_output_consensus.tmp.tsv -in2 ${DATA_DIR_TOPP}/AssayGeneratorMetabo_ffm_output_consensus.tsv)
set_tests_properties("TOPP_AssayGeneratorMetabo_3_out1" PROPERTIES DEPENDS "TOPP_AssayGeneratorMetabo_3")

# use AccurateMassSearch data with method consensus_spectrum
add_test("TOPP_AssayGeneratorMetabo_4" ${TOPP_BIN_PATH}/AssayGeneratorMetabo -test -in ${DATA_DIR_TOPP}/AssayGeneratorMetabo_input.mzML -in_featureinfo ${DATA_DIR_TOPP}/AssayGeneratorMetabo_ams_input.featureXML -out AssayGeneratorMetabo_ams_output_consensus.tmp.tsv -method consensus_spectrum -min_transitions 1 -max_transitions 3)
add_test("TOPP_AssayGeneratorMetabo_4_out1" ${DIFF} -in1 AssayGeneratorMetabo_ams_output_consensus.tmp.tsv -in2 ${DATA_DIR_TOPP}/AssayGeneratorMetabo_ams_output_consensus.tsv)
set_tests_properties("TOPP_AssayGeneratorMetabo_4_out1" PROPERTIES DEPENDS "TOPP_AssayGeneratorMetabo_4")

# use AccurateMassSearch data with method consensus_spectrum and use_known_unkowns
add_test("TOPP_AssayGeneratorMetabo_5" ${TOPP_BIN_PATH}/AssayGeneratorMetabo -test -in ${DATA_DIR_TOPP}/AssayGeneratorMetabo_input.mzML -in_featureinfo ${DATA_DIR_TOPP}/AssayGeneratorMetabo_ams_input.featureXML -out AssayGeneratorMetabo_ams_uku_output_consensus.tmp.tsv -method consensus_spectrum -use_known_unknowns -min_transitions 1 -max_transitions 3)
add_test("TOPP_AssayGeneratorMetabo_5_out1" ${DIFF} -in1 AssayGeneratorMetabo_ams_uku_output_consensus.tmp.tsv -in2 ${DATA_DIR_TOPP}/AssayGeneratorMetabo_ams_uku_output_consensus.tsv)
set_tests_properties("TOPP_AssayGeneratorMetabo_5_out1" PROPERTIES DEPENDS "TOPP_AssayGeneratorMetabo_5")

# export TOPP_AssayGeneratorMetabo_5 output to TraML
add_test("TOPP_AssayGeneratorMetabo_6" ${TOPP_BIN_PATH}/TargetedFileConverter -test -in ${DATA_DIR_TOPP}/AssayGeneratorMetabo_ams_uku_output_consensus.tsv  -out AssayGeneratorMetabo_ams_uku_output_consensus_traml.tmp.TraML)
add_test("TOPP_AssayGeneratorMetabo_6_out1" ${DIFF} -in1 AssayGeneratorMetabo_ams_uku_output_consensus_traml.tmp.TraML -in2 ${DATA_DIR_TOPP}/AssayGeneratorMetabo_ams_uku_output_consensus_traml.TraML)
set_tests_properties("TOPP_AssayGeneratorMetabo_6_out1" PROPERTIES DEPENDS "TOPP_AssayGeneratorMetabo_6")

# ImageCreator:
if(WITH_GUI)
  add_test("TOPP_ImageCreator_1" ${TOPP_BIN_PATH}/ImageCreator -test -in ${DATA_DIR_TOPP}/ImageCreator_1_input.mzML -out ImageCreator1.bmp -width 20 -height 15)
  add_test("TOPP_ImageCreator_1_out1" ${DIFF} -in1 ImageCreator1.bmp -in2 ${DATA_DIR_TOPP}/ImageCreator_1_output.bmp )
  set_tests_properties("TOPP_ImageCreator_1_out1" PROPERTIES DEPENDS "TOPP_ImageCreator_1")
  # testing forced out_type
  add_test("TOPP_ImageCreator_2" ${TOPP_BIN_PATH}/ImageCreator -test -in ${DATA_DIR_TOPP}/ImageCreator_2_input.mzML -out ImageCreator2.png -out_type bmp -width 20 -height 15 -precursors -precursor_size 3 -precursor_color green -log_intensity)
  add_test("TOPP_ImageCreator_2_out1" ${DIFF} -in1 ImageCreator2.png -in2 ${DATA_DIR_TOPP}/ImageCreator_2_output.bmp )
  set_tests_properties("TOPP_ImageCreator_2_out1" PROPERTIES DEPENDS "TOPP_ImageCreator_2")
endif()

# IDScoreSwitcher - switch q-value and PEP in Percolator result:
add_test("TOPP_IDScoreSwitcher_1" ${TOPP_BIN_PATH}/IDScoreSwitcher -test -in ${DATA_DIR_TOPP}/IDFileConverter_12_output.idXML -out IDScoreSwitcher_1_output.tmp.idXML -new_score Percolator_PEP -new_score_type "Posterior Error Probability" -new_score_orientation lower_better -old_score Percolator_qvalue)
add_test("TOPP_IDScoreSwitcher_1_out1" ${DIFF} -in1 IDScoreSwitcher_1_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDFileConverter_13_output.idXML -whitelist "IdentificationRun date")
set_tests_properties("TOPP_IDScoreSwitcher_1_out1" PROPERTIES DEPENDS "TOPP_IDScoreSwitcher_1")
# IDScoreSwitcher - switch q-value and protein probability:
add_test("TOPP_IDScoreSwitcher_2" ${TOPP_BIN_PATH}/IDScoreSwitcher -test -in ${DATA_DIR_TOPP}/IDScoreSwitcher_2_input.idXML -out IDScoreSwitcher_2_output.tmp.idXML -new_score "Posterior Probability_score" -new_score_type "Posterior Probability" -new_score_orientation higher_better -proteins)
add_test("TOPP_IDScoreSwitcher_2_out1" ${DIFF} -in1 IDScoreSwitcher_2_output.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDScoreSwitcher_2_output.idXML -whitelist "IdentificationRun date")
set_tests_properties("TOPP_IDScoreSwitcher_2_out1" PROPERTIES DEPENDS "TOPP_IDScoreSwitcher_2")

# IDSplitter:
add_test("TOPP_IDSplitter_1" ${TOPP_BIN_PATH}/IDSplitter -test -in ${DATA_DIR_TOPP}/IDMapper_1_output.featureXML -out IDSplitter_1_output1.tmp.featureXML -id_out IDSplitter_1_output2.tmp.idXML)
add_test("TOPP_IDSplitter_1_out1" ${DIFF} -in1 IDSplitter_1_output1.tmp.featureXML -in2 ${DATA_DIR_TOPP}/IDSplitter_1_output.featureXML )
set_tests_properties("TOPP_IDSplitter_1_out1" PROPERTIES DEPENDS "TOPP_IDSplitter_1")
add_test("TOPP_IDSplitter_1_out2" ${DIFF} -in1 IDSplitter_1_output2.tmp.idXML -in2 ${DATA_DIR_TOPP}/IDSplitter_1_output.idXML )
set_tests_properties("TOPP_IDSplitter_1_out2" PROPERTIES DEPENDS "TOPP_IDSplitter_1")

# JSONExporter:
# proteomics test file (= IDFileConverter_34_output1.oms):
add_test("TOPP_JSONExporter_1" ${TOPP_BIN_PATH}/JSONExporter -test -in ${DATA_DIR_TOPP}/JSONExporter_protein.oms -out JSONExporter_1_output.tmp.json)
add_test("TOPP_JSONExporter_1_out" ${DIFF} -in1 JSONExporter_1_output.tmp.json -in2 ${DATA_DIR_TOPP}/JSONExporter_1_output.json)
set_tests_properties("TOPP_JSONExporter_1_out" PROPERTIES DEPENDS "TOPP_JSONExporter_1")
# RNA test file (= NucleicAcidSearchEngine_13_out.tmp.mzml):
add_test("TOPP_JSONExporter_2" ${TOPP_BIN_PATH}/JSONExporter -test -in ${DATA_DIR_TOPP}/JSONExporter_RNA.oms -out JSONExporter_2_output.tmp.json)
add_test("TOPP_JSONExporter_2_out" ${DIFF} -in1 JSONExporter_2_output.tmp.json -in2 ${DATA_DIR_TOPP}/JSONExporter_2_output.json)
set_tests_properties("TOPP_JSONExporter_2_out" PROPERTIES DEPENDS "TOPP_JSONExporter_2")

# MassCalculator:
add_test("TOPP_MassCalculator_1" ${TOPP_BIN_PATH}/MassCalculator -test -in ${DATA_DIR_TOPP}/MassCalculator_1_input.tsv -out MassCalculator_1_output.tmp.csv -charge 0 1 -format table -separator ,)
add_test("TOPP_MassCalculator_1_out1" ${DIFF} -in1 MassCalculator_1_output.tmp.csv -in2 ${DATA_DIR_TOPP}/MassCalculator_1_output.csv )
set_tests_properties("TOPP_MassCalculator_1_out1" PROPERTIES DEPENDS "TOPP_MassCalculator_1")
add_test("TOPP_MassCalculator_2" ${TOPP_BIN_PATH}/MassCalculator -test -in_seq "LDQWLC(Carbamidomethyl)EKL" "(Glu->pyro-Glu)EAM(Oxidation)APKHK" "RANVM(Oxidation)DYR" "FGVEQDVDMVFASFIR" -out MassCalculator_2_output.tmp.csv -charge 1 2 3)
add_test("TOPP_MassCalculator_2_out1" ${DIFF} -in1 MassCalculator_2_output.tmp.csv -in2 ${DATA_DIR_TOPP}/MassCalculator_2_output.txt )
set_tests_properties("TOPP_MassCalculator_2_out1" PROPERTIES DEPENDS "TOPP_MassCalculator_2")

#------------------------------------------------------------------------------
# MultiplexResolver test
add_test("TOPP_MultiplexResolver_1" ${TOPP_BIN_PATH}/MultiplexResolver -test -in ${DATA_DIR_TOPP}/MultiplexResolver_1_input.consensusXML -ini ${DATA_DIR_TOPP}/MultiplexResolver_1_parameters.ini -out MultiplexResolver_1.tmp.consensusXML -out_conflicts MultiplexResolver_1_conflicts.tmp.consensusXML)
add_test("TOPP_MultiplexResolver_1_out1" ${DIFF} -whitelist "id=" "href=" -in1 MultiplexResolver_1.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/MultiplexResolver_1_output.consensusXML )
set_tests_properties("TOPP_MultiplexResolver_1_out1" PROPERTIES DEPENDS "TOPP_MultiplexResolver_1")
add_test("TOPP_MultiplexResolver_1_out2" ${DIFF} -whitelist "id=" "href=" -in1 MultiplexResolver_1_conflicts.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/MultiplexResolver_1_output_conflicts.consensusXML )
set_tests_properties("TOPP_MultiplexResolver_1_out2" PROPERTIES DEPENDS "TOPP_MultiplexResolver_1")

add_test("TOPP_MultiplexResolver_2" ${TOPP_BIN_PATH}/MultiplexResolver -test -in ${DATA_DIR_TOPP}/MultiplexResolver_2_input.consensusXML -ini ${DATA_DIR_TOPP}/MultiplexResolver_2_parameters.ini -out MultiplexResolver_2.tmp.consensusXML -out_conflicts MultiplexResolver_2_conflicts.tmp.consensusXML)
add_test("TOPP_MultiplexResolver_2_out1" ${DIFF} -whitelist "id=" "href=" -in1 MultiplexResolver_2.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/MultiplexResolver_2_output.consensusXML )
set_tests_properties("TOPP_MultiplexResolver_2_out1" PROPERTIES DEPENDS "TOPP_MultiplexResolver_2")
add_test("TOPP_MultiplexResolver_2_out2" ${DIFF} -whitelist "id=" "href=" -in1 MultiplexResolver_2_conflicts.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/MultiplexResolver_2_output_conflicts.consensusXML )
set_tests_properties("TOPP_MultiplexResolver_2_out2" PROPERTIES DEPENDS "TOPP_MultiplexResolver_2")

add_test("TOPP_MultiplexResolver_3" ${TOPP_BIN_PATH}/MultiplexResolver -test -in ${DATA_DIR_TOPP}/MultiplexResolver_3_input.consensusXML -ini ${DATA_DIR_TOPP}/MultiplexResolver_3_parameters.ini -out MultiplexResolver_3.tmp.consensusXML -out_conflicts MultiplexResolver_3_conflicts.tmp.consensusXML)
add_test("TOPP_MultiplexResolver_3_out1" ${DIFF} -whitelist "id=" "href=" -in1 MultiplexResolver_3.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/MultiplexResolver_3_output.consensusXML )
set_tests_properties("TOPP_MultiplexResolver_3_out1" PROPERTIES DEPENDS "TOPP_MultiplexResolver_3")
add_test("TOPP_MultiplexResolver_3_out3" ${DIFF} -whitelist "id=" "href=" -in1 MultiplexResolver_3_conflicts.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/MultiplexResolver_3_output_conflicts.consensusXML )
set_tests_properties("TOPP_MultiplexResolver_3_out3" PROPERTIES DEPENDS "TOPP_MultiplexResolver_3")

add_test("TOPP_MultiplexResolver_4" ${TOPP_BIN_PATH}/MultiplexResolver -test -in ${DATA_DIR_TOPP}/MultiplexResolver_4_input.consensusXML -in_blacklist ${DATA_DIR_TOPP}/MultiplexResolver_4_input.mzML -ini ${DATA_DIR_TOPP}/MultiplexResolver_4_parameters.ini -out MultiplexResolver_4.tmp.consensusXML -out_conflicts MultiplexResolver_4_conflicts.tmp.consensusXML)
add_test("TOPP_MultiplexResolver_4_out1" ${DIFF} -whitelist "id=" "href=" -in1 MultiplexResolver_4.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/MultiplexResolver_4_output.consensusXML )
set_tests_properties("TOPP_MultiplexResolver_4_out1" PROPERTIES DEPENDS "TOPP_MultiplexResolver_4")
add_test("TOPP_MultiplexResolver_4_out3" ${DIFF} -whitelist "id=" "href=" -in1 MultiplexResolver_4_conflicts.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/MultiplexResolver_4_output_conflicts.consensusXML )
set_tests_properties("TOPP_MultiplexResolver_4_out3" PROPERTIES DEPENDS "TOPP_MultiplexResolver_4")

# INIUpdater:
if(WITH_GUI AND HAS_XSERVER)
  add_test("TOPP_INIUpdater_1" ${TOPP_BIN_PATH}/INIUpdater -test -in ${DATA_DIR_TOPP}/INIUpdater_1_noupdate.toppas -out INIUpdater_1_noupdate.tmp.toppas)
  add_test("TOPP_INIUpdater_1_out" ${DIFF} -in1 INIUpdater_1_noupdate.tmp.toppas -in2 ${DATA_DIR_TOPP}/INIUpdater_1_noupdate.toppas -whitelist "version")
  set_tests_properties("TOPP_INIUpdater_1_out" PROPERTIES DEPENDS "TOPP_INIUpdater_1")

  add_test("TOPP_INIUpdater_2" ${TOPP_BIN_PATH}/INIUpdater -test -in ${DATA_DIR_TOPP}/INIUpdater_2_broken.toppas -out INIUpdater_2_broken.tmp.toppas)
  set_tests_properties("TOPP_INIUpdater_2" PROPERTIES WILL_FAIL 1)

  add_test("TOPP_INIUpdater_3" ${TOPP_BIN_PATH}/INIUpdater -test -in ${DATA_DIR_TOPP}/INIUpdater_3_old.toppas -out INIUpdater_3_old.tmp.toppas)
  add_test("TOPP_INIUpdater_3_out" ${DIFF} -in1 INIUpdater_3_old.tmp.toppas -in2 ${DATA_DIR_TOPP}/INIUpdater_3_old_out.toppas -whitelist "version")
  set_tests_properties("TOPP_INIUpdater_3_out" PROPERTIES DEPENDS "TOPP_INIUpdater_3")
endif()

# DatabaseFilter test:
add_test("TOPP_DatabaseFilter_1" ${TOPP_BIN_PATH}/DatabaseFilter -test -in ${DATA_DIR_TOPP}/DatabaseFilter_1.fasta -id ${DATA_DIR_TOPP}/DatabaseFilter_1.idXML -out DatabaseFilter_1_out.tmp.fasta)
add_test("TOPP_DatabaseFilter_1_out" ${DIFF} -in1 DatabaseFilter_1_out.tmp.fasta -in2 ${DATA_DIR_TOPP}/DatabaseFilter_1_out.fasta )
set_tests_properties("TOPP_DatabaseFilter_1_out" PROPERTIES DEPENDS "TOPP_DatabaseFilter_1")
add_test("TOPP_DatabaseFilter_2" ${TOPP_BIN_PATH}/DatabaseFilter -test -in ${DATA_DIR_TOPP}/DatabaseFilter_1.fasta -id ${DATA_DIR_TOPP}/DatabaseFilter_1.idXML -out DatabaseFilter_2_out.tmp.fasta -method blacklist)
add_test("TOPP_DatabaseFilter_2_out" ${DIFF} -in1 DatabaseFilter_2_out.tmp.fasta -in2 ${DATA_DIR_TOPP}/DatabaseFilter_2_out.fasta )
set_tests_properties("TOPP_DatabaseFilter_2_out" PROPERTIES DEPENDS "TOPP_DatabaseFilter_2")
add_test("TOPP_DatabaseFilter_3" ${TOPP_BIN_PATH}/DatabaseFilter -test -in ${DATA_DIR_TOPP}/DatabaseFilter_3.fasta -id ${DATA_DIR_TOPP}/DatabaseFilter_3.mzid -out DatabaseFilter_3_out.tmp.fasta)
add_test("TOPP_DatabaseFilter_3_out" ${DIFF} -in1 DatabaseFilter_3_out.tmp.fasta -in2 ${DATA_DIR_TOPP}/DatabaseFilter_3_out.fasta )
set_tests_properties("TOPP_DatabaseFilter_3_out" PROPERTIES DEPENDS "TOPP_DatabaseFilter_3")
add_test("TOPP_DatabaseFilter_4" ${TOPP_BIN_PATH}/DatabaseFilter -test -in ${DATA_DIR_TOPP}/DatabaseFilter_3.fasta -id ${DATA_DIR_TOPP}/DatabaseFilter_3.mzid -out DatabaseFilter_4_out.tmp.fasta -method blacklist)
add_test("TOPP_DatabaseFilter_4_out" ${DIFF} -in1 DatabaseFilter_4_out.tmp.fasta -in2 ${DATA_DIR_TOPP}/DatabaseFilter_4_out.fasta )
set_tests_properties("TOPP_DatabaseFilter_4_out" PROPERTIES DEPENDS "TOPP_DatabaseFilter_4")

# DecoyDatabase:
add_test("TOPP_DecoyDatabase_1" ${TOPP_BIN_PATH}/DecoyDatabase -test -in ${DATA_DIR_TOPP}/DecoyDatabase_1.fasta -out DecoyDatabase_1.tmp.fasta -only_decoy)
add_test("TOPP_DecoyDatabase_1_out" ${DIFF} -in1 DecoyDatabase_1.tmp.fasta -in2 ${DATA_DIR_TOPP}/DecoyDatabase_1_out.fasta )
set_tests_properties("TOPP_DecoyDatabase_1_out" PROPERTIES DEPENDS "TOPP_DecoyDatabase_1")
add_test("TOPP_DecoyDatabase_2" ${TOPP_BIN_PATH}/DecoyDatabase -test -in ${DATA_DIR_TOPP}/DecoyDatabase_1.fasta -out DecoyDatabase_2.tmp.fasta -decoy_string "DECOY_SEQ_" -decoy_string_position "prefix" -method shuffle -Decoy:non_shuffle_pattern "KRP" -seed 42)
add_test("TOPP_DecoyDatabase_2_out" ${DIFF} -in1 DecoyDatabase_2.tmp.fasta -in2 ${DATA_DIR_TOPP}/DecoyDatabase_2_out.fasta )
set_tests_properties("TOPP_DecoyDatabase_2_out" PROPERTIES DEPENDS "TOPP_DecoyDatabase_2")
add_test("TOPP_DecoyDatabase_3" ${TOPP_BIN_PATH}/DecoyDatabase -test -in ${DATA_DIR_TOPP}/DecoyDatabase_1.fasta -out DecoyDatabase_3.tmp.fasta -decoy_string "DECOY_SEQ_" -decoy_string_position "prefix" -method shuffle -Decoy:non_shuffle_pattern "KR" -seed 42 -enzyme "Chymotrypsin")
add_test("TOPP_DecoyDatabase_3_out" ${DIFF} -in1 DecoyDatabase_3.tmp.fasta -in2 ${DATA_DIR_TOPP}/DecoyDatabase_3_out.fasta )
set_tests_properties("TOPP_DecoyDatabase_3_out" PROPERTIES DEPENDS "TOPP_DecoyDatabase_3")
add_test("TOPP_DecoyDatabase_4" ${TOPP_BIN_PATH}/DecoyDatabase -test -type RNA -in ${DATA_DIR_TOPP}/DecoyDatabase_4.fasta -out DecoyDatabase_4.tmp.fasta -decoy_string "DECOY_SEQ_" -decoy_string_position "prefix" -method reverse -seed 42 )
add_test("TOPP_DecoyDatabase_4_out" ${DIFF} -in1 DecoyDatabase_4.tmp.fasta -in2 ${DATA_DIR_TOPP}/DecoyDatabase_4_out.fasta )
set_tests_properties("TOPP_DecoyDatabase_4_out" PROPERTIES DEPENDS "TOPP_DecoyDatabase_4")
# tests with already decoyed input data (should fail)
add_test("TOPP_DecoyDatabase_5" ${TOPP_BIN_PATH}/DecoyDatabase -test -in ${DATA_DIR_TOPP}/DecoyDatabase_5.fasta -out DecoyDatabase_5.tmp.fasta)
set_tests_properties("TOPP_DecoyDatabase_5" PROPERTIES WILL_FAIL 1)
add_test("TOPP_DecoyDatabase_6" ${TOPP_BIN_PATH}/DecoyDatabase -test -in ${DATA_DIR_TOPP}/DecoyDatabase_6.fasta -out DecoyDatabase_6.tmp.fasta)
set_tests_properties("TOPP_DecoyDatabase_6" PROPERTIES WILL_FAIL 1)

# SimpleSearchEngine:
add_test("TOPP_SimpleSearchEngine_1" ${TOPP_BIN_PATH}/SimpleSearchEngine -test
-ini ${DATA_DIR_TOPP}/SimpleSearchEngine_1.ini -in
${DATA_DIR_TOPP}/SimpleSearchEngine_1.mzML -out SimpleSearchEngine_1_out.tmp.idXML 
-database ${DATA_DIR_TOPP}/SimpleSearchEngine_1.fasta)
add_test("TOPP_SimpleSearchEngine_1_out" ${DIFF} -in1 SimpleSearchEngine_1_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/SimpleSearchEngine_1_out.idXML -whitelist "IdentificationRun date" "SearchParameters id=\"SP_0\" db=")
set_tests_properties("TOPP_SimpleSearchEngine_1_out" PROPERTIES DEPENDS
"TOPP_SimpleSearchEngine_1")

# extra percolator features
add_test("TOPP_SimpleSearchEngine_2" ${TOPP_BIN_PATH}/SimpleSearchEngine -test
-ini ${DATA_DIR_TOPP}/SimpleSearchEngine_2.ini -in
${DATA_DIR_TOPP}/SimpleSearchEngine_1.mzML -out SimpleSearchEngine_2_out.tmp.idXML
-database ${DATA_DIR_TOPP}/SimpleSearchEngine_1.fasta)
add_test("TOPP_SimpleSearchEngine_2_out" ${DIFF} -in1 SimpleSearchEngine_2_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/SimpleSearchEngine_2_out.idXML -whitelist "IdentificationRun date" "SearchParameters id=\"SP_0\" db=")
set_tests_properties("TOPP_SimpleSearchEngine_2_out" PROPERTIES DEPENDS
"TOPP_SimpleSearchEngine_2")


# FeatureFinderMetaboIdent:
add_test("TOPP_FeatureFinderMetaboIdent_1" ${TOPP_BIN_PATH}/FeatureFinderMetaboIdent -test -in ${DATA_DIR_TOPP}/FeatureFinderMetaboIdent_1_input.mzML -id ${DATA_DIR_TOPP}/FeatureFinderMetaboIdent_1_input.tsv -out FeatureFinderMetaboIdent_1_output.tmp.featureXML -extract:mz_window 5 -extract:rt_window 20 -detect:peak_width 3)
add_test("TOPP_FeatureFinderMetaboIdent_1_out1" ${DIFF} -whitelist "id=" -in1 FeatureFinderMetaboIdent_1_output.tmp.featureXML -in2 ${DATA_DIR_TOPP}/FeatureFinderMetaboIdent_1_output.featureXML)
set_tests_properties("TOPP_FeatureFinderMetaboIdent_1_out1" PROPERTIES DEPENDS "TOPP_FeatureFinderMetaboIdent_1")

# MzMLSplitter test:
# specify number of parts directly:
add_test("TOPP_MzMLSplitter_1" ${TOPP_BIN_PATH}/MzMLSplitter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out MzMLSplitter_1_output -parts 2)
add_test("TOPP_MzMLSplitter_1_out1" ${DIFF} -in1 MzMLSplitter_1_output_part1of2.mzML -in2 ${DATA_DIR_TOPP}/MzMLSplitter_output_part1.mzML)
set_tests_properties("TOPP_MzMLSplitter_1_out1" PROPERTIES DEPENDS "TOPP_MzMLSplitter_1")
add_test("TOPP_MzMLSplitter_1_out2" ${DIFF} -in1 MzMLSplitter_1_output_part2of2.mzML -in2 ${DATA_DIR_TOPP}/MzMLSplitter_output_part2.mzML)
set_tests_properties("TOPP_MzMLSplitter_1_out2" PROPERTIES DEPENDS "TOPP_MzMLSplitter_1")
# specify target file size:
add_test("TOPP_MzMLSplitter_2" ${TOPP_BIN_PATH}/MzMLSplitter -test -in ${DATA_DIR_TOPP}/FileFilter_1_input.mzML -out MzMLSplitter_2_output -size 40 -unit KB)
add_test("TOPP_MzMLSplitter_2_out1" ${DIFF} -in1 MzMLSplitter_2_output_part1of2.mzML -in2 ${DATA_DIR_TOPP}/MzMLSplitter_output_part1.mzML)
set_tests_properties("TOPP_MzMLSplitter_2_out1" PROPERTIES DEPENDS "TOPP_MzMLSplitter_2")
add_test("TOPP_MzMLSplitter_2_out2" ${DIFF} -in1 MzMLSplitter_2_output_part2of2.mzML -in2 ${DATA_DIR_TOPP}/MzMLSplitter_output_part2.mzML)
set_tests_properties("TOPP_MzMLSplitter_2_out2" PROPERTIES DEPENDS "TOPP_MzMLSplitter_2")

# TICCalculator test:
add_test("TOPP_TICCalculator_1" ${TOPP_BIN_PATH}/TICCalculator -test -in ${DATA_DIR_TOPP}/MapNormalizer_output.mzML -read_method regular)
add_test("TOPP_TICCalculator_2" ${TOPP_BIN_PATH}/TICCalculator -test -in ${DATA_DIR_TOPP}/MapNormalizer_output.mzML -read_method streaming -loadData true)
add_test("TOPP_TICCalculator_3" ${TOPP_BIN_PATH}/TICCalculator -test -in ${DATA_DIR_TOPP}/MapNormalizer_output.mzML -read_method streaming -loadData false)
add_test("TOPP_TICCalculator_4" ${TOPP_BIN_PATH}/TICCalculator -test -in ${DATA_DIR_TOPP}/MapNormalizer_output.mzML -read_method indexed)
add_test("TOPP_TICCalculator_5" ${TOPP_BIN_PATH}/TICCalculator -test -in ${DATA_DIR_TOPP}/MapNormalizer_output.mzML -read_method indexed_parallel)

# ProteomicsLFQ test with a subset of the experimental design. Also tests the nonconsecutiveness of the design:
add_test("TOPP_ProteomicsLFQ_1_subset" ${TOPP_BIN_PATH}/ProteomicsLFQ
         -in
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F1.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F2.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F1.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F2.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA3_F1.mzML
         -ids
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F1.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F2.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F1.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F2.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA3_F1.idXML
         -design
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA_design_onetable_nonconsec.tsv
         -Alignment:align_algorithm:max_rt_shift 0
         -fasta
         ${DATA_DIR_SHARE}/examples/TOPPAS/data/BSA_Identification/18Protein_SoCe_Tr_detergents_trace_target_decoy.fasta
         -targeted_only true
         -mass_recalibration false
         -out_cxml BSA_subset.tmp.consensusXML
         -out_msstats BSA_subset.tmp.csv
         -out BSA_subset.tmp.mzTab
         -out_triqler BSA_subset.tmp.tsv
         -threads 1
         -proteinFDR 0.3
         -test
         )
add_test("TOPP_ProteomicsLFQ_1_subset_out_1" ${DIFF} -whitelist "spectra_data" "map id=" "InferenceEngineVersion" -in1 BSA_subset.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_1_subset_out.consensusXML )
set_tests_properties("TOPP_ProteomicsLFQ_1_subset_out_1" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_1_subset")
add_test("TOPP_ProteomicsLFQ_1_subset_out_2" ${DIFF}  -whitelist "software" "location" -in1 BSA_subset.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_1_subset_out.csv )
set_tests_properties("TOPP_ProteomicsLFQ_1_subset_out_2" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_1_subset")
add_test("TOPP_ProteomicsLFQ_1_subset_out_3" ${DIFF} -whitelist "software" "location" "InferenceEngineVersion" "TOPPProteinInference q-value" -in1 BSA_subset.tmp.mzTab -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_1_subset_out.mzTab )
set_tests_properties("TOPP_ProteomicsLFQ_1_subset_out_3" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_1_subset")
add_test("TOPP_ProteomicsLFQ_1_subset_out_4" ${DIFF} -in1 BSA_subset.tmp.tsv -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_1_subset_out.tsv )
set_tests_properties("TOPP_ProteomicsLFQ_1_subset_out_4" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_1_subset")


add_test("TOPP_ProteomicsLFQ_1" ${TOPP_BIN_PATH}/ProteomicsLFQ
         -in
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F1.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F2.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F1.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F2.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA3_F1.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA3_F2.mzML
         -ids
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F1.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F2.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F1.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F2.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA3_F1.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA3_F2.idXML
         -design
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA_design.tsv
         -Alignment:align_algorithm:max_rt_shift 0
         -fasta
         ${DATA_DIR_SHARE}/examples/TOPPAS/data/BSA_Identification/18Protein_SoCe_Tr_detergents_trace_target_decoy.fasta
         -targeted_only true
         -mass_recalibration false
         -out_cxml BSA.tmp.consensusXML
         -out_msstats BSA.tmp.csv
         -out BSA.tmp.mztab
         -out_triqler BSA.tmp.tsv
         -threads 1
         -proteinFDR 0.3
         -test
         )
add_test("TOPP_ProteomicsLFQ_1_out_1" ${DIFF} -whitelist "spectra_data" "map id=" "InferenceEngineVersion" -in1 BSA.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_1_out.consensusXML )
set_tests_properties("TOPP_ProteomicsLFQ_1_out_1" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_1")
add_test("TOPP_ProteomicsLFQ_1_out_2" ${DIFF}  -whitelist "software" "location" -in1 BSA.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_1_out.csv )
set_tests_properties("TOPP_ProteomicsLFQ_1_out_2" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_1")
add_test("TOPP_ProteomicsLFQ_1_out_3" ${DIFF} -whitelist "software" "location" "InferenceEngineVersion" "TOPPProteinInference q-value" -in1 BSA.tmp.mztab -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_1_out.mzTab )
set_tests_properties("TOPP_ProteomicsLFQ_1_out_3" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_1")
add_test("TOPP_ProteomicsLFQ_1_out_4" ${DIFF} -in1 BSA.tmp.tsv -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_1_out.tsv )
set_tests_properties("TOPP_ProteomicsLFQ_1_out_4" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_1")

# regression test (decoys in ID files)
add_test("TOPP_ProteomicsLFQ_2" ${TOPP_BIN_PATH}/ProteomicsLFQ
         -in
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F1.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F2.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F1.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F2.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA3_F1.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA3_F2.mzML
         -ids
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F1_decoys.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F2.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F1.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F2.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA3_F1.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA3_F2.idXML
         -design
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA_design.tsv
         -Alignment:align_algorithm:max_rt_shift 0
         -fasta
         ${DATA_DIR_SHARE}/examples/TOPPAS/data/BSA_Identification/18Protein_SoCe_Tr_detergents_trace_target_decoy.fasta
         -targeted_only true
         -mass_recalibration false
         -out_cxml BSA_plusTransfer.tmp.consensusXML
         -out_msstats BSA_plusTransfer.tmp.csv
         -out BSA_plusTransfer.tmp.mztab
         -out_triqler BSA_plusTransfer.tmp.tsv
         -threads 1
         -proteinFDR 0.3
         -test
         )
add_test("TOPP_ProteomicsLFQ_2_out_1" ${DIFF} -whitelist "spectra_data" "map id=" "InferenceEngineVersion" -in1 BSA_plusTransfer.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_2_out.consensusXML )
set_tests_properties("TOPP_ProteomicsLFQ_2_out_1" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_2")
add_test("TOPP_ProteomicsLFQ_2_out_2" ${DIFF}  -whitelist "software" "location" -in1 BSA_plusTransfer.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_2_out.csv )
set_tests_properties("TOPP_ProteomicsLFQ_2_out_2" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_2")
add_test("TOPP_ProteomicsLFQ_2_out_3" ${DIFF}  -whitelist "software" "location" "InferenceEngineVersion" "TOPPProteinInference q-value" -in1 BSA_plusTransfer.tmp.mztab -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_2_out.mzTab )
set_tests_properties("TOPP_ProteomicsLFQ_2_out_3" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_2")
add_test("TOPP_ProteomicsLFQ_2_out_4" ${DIFF}  -in1 BSA_plusTransfer.tmp.tsv -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_2_out.tsv )
set_tests_properties("TOPP_ProteomicsLFQ_2_out_4" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_2")

add_test("TOPP_ProteomicsLFQ_3" ${TOPP_BIN_PATH}/ProteomicsLFQ
         -in
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F1.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F2.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F1.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F2.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA3_F1.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA3_F2.mzML
         -ids
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F1.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F2.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F1.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F2.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA3_F1.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA3_F2.idXML
         -design
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA_design.tsv
         -Alignment:align_algorithm:max_rt_shift 0
         -Linking:min_nr_diffs_per_bin 10
         -fasta
         ${DATA_DIR_SHARE}/examples/TOPPAS/data/BSA_Identification/18Protein_SoCe_Tr_detergents_trace_target_decoy.fasta
         -targeted_only false
         -mass_recalibration false
         -out_cxml BSA_plusSeeds_plusTransfer.tmp.consensusXML
         -out_msstats BSA_plusSeeds_plusTransfer.tmp.csv
         -out BSA_plusSeeds_plusTransfer.tmp.mzTab
         -out_triqler BSA_plusSeeds_plusTransfer.tmp.tsv
         -threads 1
         -proteinFDR 0.3
         -test
         -debug 0
         )
add_test("TOPP_ProteomicsLFQ_3_out_1" ${DIFF} -whitelist "spectra_data" "map id=" "InferenceEngineVersion" -in1 BSA_plusSeeds_plusTransfer.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_3_out.consensusXML )
set_tests_properties("TOPP_ProteomicsLFQ_3_out_1" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_3")
add_test("TOPP_ProteomicsLFQ_3_out_2" ${DIFF}  -whitelist "software" "location" -in1 BSA_plusSeeds_plusTransfer.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_3_out.csv )
set_tests_properties("TOPP_ProteomicsLFQ_2_out_2" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_3")
add_test("TOPP_ProteomicsLFQ_3_out_3" ${DIFF}  -whitelist "software" "location" "InferenceEngineVersion" "TOPPProteinInference q-value" -in1 BSA_plusSeeds_plusTransfer.tmp.mzTab -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_3_out.mzTab )
set_tests_properties("TOPP_ProteomicsLFQ_2_out_3" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_3")
add_test("TOPP_ProteomicsLFQ_3_out_4" ${DIFF}  -in1 BSA_plusSeeds_plusTransfer.tmp.tsv -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_3_out.tsv )
set_tests_properties("TOPP_ProteomicsLFQ_3_out_4" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_3")

# spectral counting
add_test("TOPP_ProteomicsLFQ_4" ${TOPP_BIN_PATH}/ProteomicsLFQ
         -in
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F1.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F2.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F1.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F2.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA3_F1.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA3_F2.mzML
         -ids
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F1.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F2.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F1.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F2.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA3_F1.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA3_F2.idXML
         -design
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA_design.tsv
         -fasta
         ${DATA_DIR_SHARE}/examples/TOPPAS/data/BSA_Identification/18Protein_SoCe_Tr_detergents_trace_target_decoy.fasta
         -quantification_method spectral_counting
         -out_cxml BSA_sc.tmp.consensusXML
         -out BSA_sc.tmp.mztab
         -threads 1
         -proteinFDR 0.3
         -test
         )
add_test("TOPP_ProteomicsLFQ_4_out_1" ${DIFF} -whitelist "spectra_data" "map id=" "InferenceEngineVersion" -in1 BSA_sc.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_4_out.consensusXML )
set_tests_properties("TOPP_ProteomicsLFQ_4_out_1" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_4")
add_test("TOPP_ProteomicsLFQ_4_out_2" ${DIFF} -whitelist "software" "location" "InferenceEngineVersion" "TOPPProteinInference q-value" -in1 BSA_sc.tmp.mztab -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_4_out.mzTab )
set_tests_properties("TOPP_ProteomicsLFQ_4_out_2" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_4")

# test decoy report
add_test("TOPP_ProteomicsLFQ_5" ${TOPP_BIN_PATH}/ProteomicsLFQ
         -in
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F1.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F2.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F1.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F2.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA3_F1.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA3_F2.mzML
         -ids
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F1_decoys.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F2.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F1.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F2.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA3_F1.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA3_F2.idXML
         -design
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA_design.tsv
         -Alignment:align_algorithm:max_rt_shift 0
         -fasta
         ${DATA_DIR_SHARE}/examples/TOPPAS/data/BSA_Identification/18Protein_SoCe_Tr_detergents_trace_target_decoy.fasta
         -PeptideQuantification:quantify_decoys
         -targeted_only true
         -mass_recalibration false
         -out_cxml BSA_decoys.tmp.consensusXML
         -out_msstats BSA_decoys.tmp.csv
         -out BSA_decoys.tmp.mztab
         -out_triqler BSA_decoys.tmp.tsv
         -threads 1
         -proteinFDR 0.3
         -test
         )
add_test("TOPP_ProteomicsLFQ_5_out_1" ${DIFF} -whitelist "spectra_data" "map id=" "InferenceEngineVersion" -in1 BSA_decoys.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_5_out.consensusXML )
set_tests_properties("TOPP_ProteomicsLFQ_5_out_1" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_5")
add_test("TOPP_ProteomicsLFQ_5_out_2" ${DIFF}  -whitelist "software" "location" -in1 BSA_decoys.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_5_out.csv )
set_tests_properties("TOPP_ProteomicsLFQ_5_out_2" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_5")
add_test("TOPP_ProteomicsLFQ_5_out_3" ${DIFF} -whitelist "software" "location" "InferenceEngineVersion" "TOPPProteinInference q-value" -in1 BSA_decoys.tmp.mztab -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_5_out.mzTab )
set_tests_properties("TOPP_ProteomicsLFQ_5_out_3" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_5")
add_test("TOPP_ProteomicsLFQ_5_out_4" ${DIFF} -in1 BSA_decoys.tmp.tsv -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_5_out.tsv )
set_tests_properties("TOPP_ProteomicsLFQ_5_out_4" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_5")

# spectral counting (regression test with target + decoys)
add_test("TOPP_ProteomicsLFQ_6" ${TOPP_BIN_PATH}/ProteomicsLFQ
         -in
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F1.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F2.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F1.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F2.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA3_F1.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA3_F2.mzML
         -ids
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F1_msgf_idx_fdr_idpep_switched_filter.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F2_msgf_idx_fdr_idpep_switched_filter.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F1_msgf_idx_fdr_idpep_switched_filter.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F2_msgf_idx_fdr_idpep_switched_filter.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA3_F1_msgf_idx_fdr_idpep_switched_filter.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA3_F2_msgf_idx_fdr_idpep_switched_filter.idXML
         -design
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA_design.tsv
         -fasta
         ${DATA_DIR_SHARE}/examples/TOPPAS/data/BSA_Identification/18Protein_SoCe_Tr_detergents_trace_target_decoy.fasta
         -quantification_method spectral_counting
         -out_cxml BSA_sc_decoys.tmp.consensusXML
         -out BSA_sc_decoys.tmp.mztab
         -threads 1
         -proteinFDR 1.0 # disable FDR filtering to ensure we keep all decoys for testing purpose
         -test
         )
add_test("TOPP_ProteomicsLFQ_6_out_1" ${DIFF} -whitelist "spectra_data" "map id=" "InferenceEngineVersion" -in1 BSA_sc_decoys.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_6_out.consensusXML )
set_tests_properties("TOPP_ProteomicsLFQ_6_out_1" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_6")
add_test("TOPP_ProteomicsLFQ_6_out_2" ${DIFF} -whitelist "software" "location" "InferenceEngineVersion" "TOPPProteinInference q-value" -in1 BSA_sc_decoys.tmp.mztab -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_6_out.mzTab )
set_tests_properties("TOPP_ProteomicsLFQ_6_out_2" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_6")

add_test("TOPP_ProteomicsLFQ_7" ${TOPP_BIN_PATH}/ProteomicsLFQ
         -in
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F2.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F2.mzML
         -ids
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F2.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F2.idXML
         -design
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA_design.tsv
         -Alignment:align_algorithm:max_rt_shift 0
         -fasta
         ${DATA_DIR_SHARE}/examples/TOPPAS/data/BSA_Identification/18Protein_SoCe_Tr_detergents_trace_target_decoy.fasta
         -targeted_only true
         -mass_recalibration false
         -out_cxml BSA_sub.tmp.consensusXML
         -out_msstats BSA_sub.tmp.csv
         -out BSA_sub.tmp.mztab
         -out_triqler BSA_sub.tmp.tsv
         -threads 1
         -proteinFDR 0.8
         -test
         )
add_test("TOPP_ProteomicsLFQ_7_out_1" ${DIFF} -whitelist "spectra_data" "map id=" "InferenceEngineVersion" -in1 BSA_sub.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_7_out.consensusXML )
set_tests_properties("TOPP_ProteomicsLFQ_7_out_1" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_7")
add_test("TOPP_ProteomicsLFQ_7_out_2" ${DIFF}  -whitelist "software" "location" -in1 BSA_sub.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_7_out.csv )
set_tests_properties("TOPP_ProteomicsLFQ_7_out_2" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_7")
add_test("TOPP_ProteomicsLFQ_7_out_3" ${DIFF} -whitelist "software" "location" "InferenceEngineVersion" "TOPPProteinInference q-value" -in1 BSA_sub.tmp.mztab -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_7_out.mzTab )
set_tests_properties("TOPP_ProteomicsLFQ_7_out_3" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_7")
add_test("TOPP_ProteomicsLFQ_7_out_4" ${DIFF} -in1 BSA_sub.tmp.tsv -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_7_out.tsv )
set_tests_properties("TOPP_ProteomicsLFQ_7_out_4" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_7")

# test empty id
add_test("TOPP_ProteomicsLFQ_8" ${TOPP_BIN_PATH}/ProteomicsLFQ
         -in
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F2.mzML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F2.mzML
         -ids
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F2.idXML
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA2_F2_empty.idXML
         -design
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA_design.tsv
         -Alignment:align_algorithm:max_rt_shift 0
         -fasta
         ${DATA_DIR_SHARE}/examples/TOPPAS/data/BSA_Identification/18Protein_SoCe_Tr_detergents_trace_target_decoy.fasta
         -targeted_only true
         -mass_recalibration false
         -out_cxml BSA_sub_8.tmp.consensusXML
         -out_msstats BSA_sub_8.tmp.csv
         -out BSA_sub_8.tmp.mzTab
         -out_triqler BSA_sub_8.tmp.tsv
         -threads 1
         -proteinFDR 0.8
         -test
         )
add_test("TOPP_ProteomicsLFQ_8_out_1" ${DIFF} -whitelist "spectra_data" "map id=" "InferenceEngineVersion" -in1 BSA_sub_8.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_8_out.consensusXML )
set_tests_properties("TOPP_ProteomicsLFQ_8_out_1" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_8")
add_test("TOPP_ProteomicsLFQ_8_out_2" ${DIFF}  -whitelist "software" "location" -in1 BSA_sub_8.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_8_out.csv )
set_tests_properties("TOPP_ProteomicsLFQ_8_out_2" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_8")
add_test("TOPP_ProteomicsLFQ_8_out_3" ${DIFF} -whitelist "software" "location" "InferenceEngineVersion" "TOPPProteinInference q-value" -in1 BSA_sub_8.tmp.mzTab -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_8_out.mzTab )
set_tests_properties("TOPP_ProteomicsLFQ_8_out_3" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_8")
add_test("TOPP_ProteomicsLFQ_8_out_4" ${DIFF} -in1 BSA_sub_8.tmp.tsv -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_8_out.tsv )
set_tests_properties("TOPP_ProteomicsLFQ_8_out_4" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_8")

# test single file, no experimental design (regression)
add_test("TOPP_ProteomicsLFQ_9" ${TOPP_BIN_PATH}/ProteomicsLFQ
         -in
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F2.mzML
         -ids
         ${DATA_DIR_SHARE}/examples/FRACTIONS/BSA1_F2.idXML
         -Alignment:align_algorithm:max_rt_shift 0
         -fasta
         ${DATA_DIR_SHARE}/examples/TOPPAS/data/BSA_Identification/18Protein_SoCe_Tr_detergents_trace_target_decoy.fasta
         -targeted_only false
         -mass_recalibration false
         -out_cxml BSA_sub_9.tmp.consensusXML
         -out_msstats BSA_sub_9.tmp.csv
         -out BSA_sub_9.tmp.mzTab
         -out_triqler BSA_sub_9.tmp.tsv
         -threads 1
         -proteinFDR 0.8
         -test
         )
add_test("TOPP_ProteomicsLFQ_9_out_1" ${DIFF} -whitelist "spectra_data" "map id=" "InferenceEngineVersion" -in1 BSA_sub_9.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_9_out.consensusXML )
set_tests_properties("TOPP_ProteomicsLFQ_9_out_1" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_9")
add_test("TOPP_ProteomicsLFQ_9_out_2" ${DIFF}  -whitelist "software" "location" -in1 BSA_sub_9.tmp.csv -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_9_out.csv )
set_tests_properties("TOPP_ProteomicsLFQ_9_out_2" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_9")
add_test("TOPP_ProteomicsLFQ_9_out_3" ${DIFF} -whitelist "software" "location" "InferenceEngineVersion" "TOPPProteinInference q-value" -in1 BSA_sub_9.tmp.mzTab -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_9_out.mzTab )
set_tests_properties("TOPP_ProteomicsLFQ_9_out_3" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_9")
add_test("TOPP_ProteomicsLFQ_9_out_4" ${DIFF} -in1 BSA_sub_9.tmp.tsv -in2 ${DATA_DIR_TOPP}/ProteomicsLFQ_9_out.tsv )
set_tests_properties("TOPP_ProteomicsLFQ_9_out_4" PROPERTIES DEPENDS "TOPP_ProteomicsLFQ_9")


#------------------------------------------------------------------------------
# NucleicAcidSearchEngine:
add_test("TOPP_NucleicAcidSearchEngine_1" ${TOPP_BIN_PATH}/NucleicAcidSearchEngine -test -ini ${DATA_DIR_TOPP}/NucleicAcidSearchEngine_1.ini -in ${DATA_DIR_TOPP}/NucleicAcidSearchEngine_1.mzML -id_out NucleicAcidSearchEngine_11_out.tmp.idXML -out NucleicAcidSearchEngine_12_out.tmp.mzTab -db_out NucleicAcidSearchEngine_13_out.tmp.oms -digest_out NucleicAcidSearchEngine_1_digest.oms -database ${DATA_DIR_TOPP}/NucleicAcidSearchEngine_1.fasta)
add_test("TOPP_NucleicAcidSearchEngine_11_out" ${DIFF} -in1 NucleicAcidSearchEngine_11_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/NucleicAcidSearchEngine_11_out.idXML -whitelist "IdentificationRun date" "SearchParameters id=\"SP_0\" db=")
add_test("TOPP_NucleicAcidSearchEngine_12_out" ${DIFF} -in1 NucleicAcidSearchEngine_12_out.tmp.mzTab -in2 ${DATA_DIR_TOPP}/NucleicAcidSearchEngine_12_out.mzTab -whitelist "ms_run[1]-location")
# don't compare binary .oms (SQLite) output file (NucleicAcidSearchEngine_13_out.tmp.oms)
set_tests_properties("TOPP_NucleicAcidSearchEngine_11_out" PROPERTIES DEPENDS "TOPP_NucleicAcidSearchEngine_1")
set_tests_properties("TOPP_NucleicAcidSearchEngine_12_out" PROPERTIES DEPENDS "TOPP_NucleicAcidSearchEngine_1")
# use pre-digested sequence data:
add_test("TOPP_NucleicAcidSearchEngine_2" ${TOPP_BIN_PATH}/NucleicAcidSearchEngine -test -ini ${DATA_DIR_TOPP}/NucleicAcidSearchEngine_1.ini -in ${DATA_DIR_TOPP}/NucleicAcidSearchEngine_1.mzML -id_out NucleicAcidSearchEngine_21_out.tmp.idXML -out NucleicAcidSearchEngine_22_out.tmp.mzTab -db_out NucleicAcidSearchEngine_23_out.tmp.oms -digest NucleicAcidSearchEngine_1_digest.oms)
set_tests_properties("TOPP_NucleicAcidSearchEngine_2" PROPERTIES DEPENDS "TOPP_NucleicAcidSearchEngine_1")
add_test("TOPP_NucleicAcidSearchEngine_21_out" ${DIFF} -in1 NucleicAcidSearchEngine_21_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/NucleicAcidSearchEngine_11_out.idXML -whitelist "IdentificationRun date" "SearchParameters id=\"SP_0\" db=")
add_test("TOPP_NucleicAcidSearchEngine_22_out" ${DIFF} -in1 NucleicAcidSearchEngine_22_out.tmp.mzTab -in2 ${DATA_DIR_TOPP}/NucleicAcidSearchEngine_12_out.mzTab -whitelist "ms_run[1]-location")
# don't compare binary .oms (SQLite) output file (NucleicAcidSearchEngine_23_out.oms)
set_tests_properties("TOPP_NucleicAcidSearchEngine_21_out" PROPERTIES DEPENDS "TOPP_NucleicAcidSearchEngine_2")
set_tests_properties("TOPP_NucleicAcidSearchEngine_22_out" PROPERTIES DEPENDS "TOPP_NucleicAcidSearchEngine_2")

#------------------------------------------------------------------------------
# RNAMassCalculator:
add_test("TOPP_RNAMassCalculator_1" ${TOPP_BIN_PATH}/RNAMassCalculator -test -in_seq "AUCGGC" -charge -1 -2 -out RNAMassCalculator.tmp.txt)
add_test("TOPP_RNAMassCalculator_1_out" ${DIFF} -in1 RNAMassCalculator.tmp.txt -in2 ${DATA_DIR_TOPP}/RNAMassCalculator_1.txt)
set_tests_properties("TOPP_RNAMassCalculator_1_out" PROPERTIES DEPENDS "TOPP_RNAMassCalculator_1")

#------------------------------------------------------------------------------
# EPIFANY:
add_test("TOPP_Epifany_1" ${TOPP_BIN_PATH}/Epifany -test -in ${DATA_DIR_TOPP}/THIRDPARTY/FidoAdapter_1_input.idXML -algorithm:model_parameters:prot_prior 0.7 -algorithm:model_parameters:pep_spurious_emission 0.001 -algorithm:model_parameters:pep_emission 0.1 -out Epifany_1_out.tmp.idXML)
add_test("TOPP_Epifany_1_out" ${DIFF} -whitelist "InferenceEngineVersion" -in1 Epifany_1_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/Epifany_1_out.idXML)
set_tests_properties("TOPP_Epifany_1_out" PROPERTIES DEPENDS "TOPP_Epifany_1")
add_test("TOPP_Epifany_2" ${TOPP_BIN_PATH}/Epifany -test -in ${DATA_DIR_TOPP}/Epifany_2_input.consensusXML -algorithm:model_parameters:prot_prior 0.7 -algorithm:model_parameters:pep_spurious_emission 0.001 -algorithm:model_parameters:pep_emission 0.1 -out Epifany_2_out.tmp.consensusXML -out_type consensusXML)
add_test("TOPP_Epifany_2_out" ${DIFF} -whitelist "InferenceEngineVersion" -in1 Epifany_2_out.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/Epifany_2_out.consensusXML)
set_tests_properties("TOPP_Epifany_2_out" PROPERTIES DEPENDS "TOPP_Epifany_2")
add_test("TOPP_Epifany_3" ${TOPP_BIN_PATH}/Epifany -test -in ${DATA_DIR_TOPP}/Epifany_2_input.consensusXML -algorithm:keep_best_PSM_only false -algorithm:model_parameters:prot_prior 0.7 -algorithm:model_parameters:pep_spurious_emission 0.001 -algorithm:model_parameters:pep_emission 0.1 -out Epifany_3_out.tmp.consensusXML -out_type consensusXML)
add_test("TOPP_Epifany_3_out" ${DIFF} -whitelist "InferenceEngineVersion" -in1 Epifany_3_out.tmp.consensusXML -in2 ${DATA_DIR_TOPP}/Epifany_3_out.consensusXML)
set_tests_properties("TOPP_Epifany_3_out" PROPERTIES DEPENDS "TOPP_Epifany_3")
add_test("TOPP_Epifany_4" ${TOPP_BIN_PATH}/Epifany -test -in ${DATA_DIR_TOPP}/THIRDPARTY/FidoAdapter_1_input_singlerun.idXML -algorithm:model_parameters:prot_prior 0.7 -algorithm:model_parameters:pep_spurious_emission 0.001 -algorithm:model_parameters:pep_emission 0.1 -greedy_group_resolution remove_proteins_wo_evidence -out Epifany_4_out.tmp.idXML)
add_test("TOPP_Epifany_4_out" ${DIFF} -whitelist "InferenceEngineVersion" -in1 Epifany_4_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/Epifany_4_out_singlerun.idXML)
set_tests_properties("TOPP_Epifany_4_out" PROPERTIES DEPENDS "TOPP_Epifany_4")
add_test("TOPP_Epifany_5" ${TOPP_BIN_PATH}/Epifany -test -in ${DATA_DIR_TOPP}/THIRDPARTY/FidoAdapter_1_input_singlerun.idXML -algorithm:model_parameters:prot_prior 0.7 -algorithm:model_parameters:pep_spurious_emission 0.001 -algorithm:model_parameters:pep_emission 0.1 -greedy_group_resolution remove_proteins_wo_evidence -protein_fdr true -picked_fdr false -out Epifany_5_out.tmp.idXML)
add_test("TOPP_Epifany_5_out" ${DIFF} -whitelist "InferenceEngineVersion" -in1 Epifany_5_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/Epifany_4_out_singlerun_fdr.idXML)
set_tests_properties("TOPP_Epifany_5_out" PROPERTIES DEPENDS "TOPP_Epifany_5")

#------------------------------------------------------------------------------
# QCCalculator:
# test for qcML file export
add_test("TOPP_QCCalculator_1" ${TOPP_BIN_PATH}/QCCalculator -test -in ${DATA_DIR_TOPP}/QCCalculator_input.mzML -out QCCalculator_1.qcML)
add_test("TOPP_QCCalculator_1_out" ${DIFF} -in1 QCCalculator_1.qcML -in2 ${DATA_DIR_TOPP}/QCCalculator_1_output.qcML)
set_tests_properties("TOPP_QCCalculator_1_out" PROPERTIES DEPENDS "TOPP_QCCalculator_1")

# test for mzQC file export
add_test("TOPP_QCCalculator_2" ${TOPP_BIN_PATH}/QCCalculator -test -in ${DATA_DIR_TOPP}/QCCalculator_input.mzML -label label -name name -address address -description description -out QCCalculator_2.mzQC -feature ${DATA_DIR_TOPP}/FeatureFinderMetaboIdent_1_output.featureXML -id ${DATA_DIR_TOPP}/OpenPepXL_output.idXML)
add_test("TOPP_QCCalculator_2_out" ${DIFF} -whitelist "creationDate" "location" "version" -in1 QCCalculator_2.mzQC -in2 ${DATA_DIR_TOPP}/QCCalculator_2_output.mzQC)
set_tests_properties("TOPP_QCCalculator_2_out" PROPERTIES DEPENDS "TOPP_QCCalculator_2")

#------------------------------------------------------------------------------
# PSMFeatureExtractor
add_test("TOPP_PSMFeatureExtractor_1" ${TOPP_BIN_PATH}/PSMFeatureExtractor -test -in ${DATA_DIR_TOPP}/THIRDPARTY/MSGFPlusAdapter_1_out.idXML ${DATA_DIR_TOPP}/THIRDPARTY/XTandemAdapter_1_out.idXML -multiple_search_engines -skip_db_check -out_type idXML -out PSMFeatureExtractor_1_out.tmp.idXML)
add_test("TOPP_PSMFeatureExtractor_1_out" ${DIFF} -whitelist "InferenceEngineVersion" "MSGFPlusAdapter:" "date=" -in1 PSMFeatureExtractor_1_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/PSMFeatureExtractor_1_out.idXML)
set_tests_properties("TOPP_PSMFeatureExtractor_1_out" PROPERTIES DEPENDS "TOPP_PSMFeatureExtractor_1")

add_test("TOPP_PSMFeatureExtractor_2" ${TOPP_BIN_PATH}/PSMFeatureExtractor -test -in ${DATA_DIR_TOPP}/THIRDPARTY/MSGFPlusAdapter_1_out.idXML ${DATA_DIR_TOPP}/THIRDPARTY/XTandemAdapter_1_out.idXML -multiple_search_engines -skip_db_check -out_type mzid -out PSMFeatureExtractor_2_out.tmp.mzid)
add_test("TOPP_PSMFeatureExtractor_2_out" ${DIFF} -whitelist "InferenceEngineVersion" "creationDate=" "activityDate=" "version=" -in1 PSMFeatureExtractor_2_out.tmp.mzid -in2 ${DATA_DIR_TOPP}/PSMFeatureExtractor_2_out.mzid)
set_tests_properties("TOPP_PSMFeatureExtractor_2_out" PROPERTIES DEPENDS "TOPP_PSMFeatureExtractor_2")

add_test("TOPP_PSMFeatureExtractor_3" ${TOPP_BIN_PATH}/PSMFeatureExtractor -test -in ${DATA_DIR_TOPP}/THIRDPARTY/MSGFPlusAdapter_1_out.idXML -out_type idXML -out PSMFeatureExtractor_3_out.tmp.idXML)
add_test("TOPP_PSMFeatureExtractor_3_out" ${DIFF} -whitelist "InferenceEngineVersion" "date=" -in1 PSMFeatureExtractor_3_out.tmp.idXML -in2 ${DATA_DIR_TOPP}/PSMFeatureExtractor_3_out.idXML)
set_tests_properties("TOPP_PSMFeatureExtractor_3_out" PROPERTIES DEPENDS "TOPP_PSMFeatureExtractor_3")


#------------------------------------------------------------------------------
# TOPP tests reading nested json input
#------------------------------------------------------------------------------
# MassTraceExtractor test
add_test("TOPP_MassTraceExtractor_1_json" ${TOPP_BIN_PATH}/MassTraceExtractor -test -ini ${DATA_DIR_TOPP}/MassTraceExtractor.ini.json -in ${DATA_DIR_TOPP}/MassTraceExtractor_1_input.mzML -out MassTraceExtractor_1_json.tmp.featureXML)
add_test("TOPP_MassTraceExtractor_1_out_json" ${DIFF} -whitelist "id=" -in1 MassTraceExtractor_1_json.tmp.featureXML -in2 ${DATA_DIR_TOPP}/MassTraceExtractor_1_output.featureXML )
set_tests_properties("TOPP_MassTraceExtractor_1_out_json" PROPERTIES DEPENDS "TOPP_MassTraceExtractor_1_json")
add_test("TOPP_MassTraceExtractor_2_json" ${TOPP_BIN_PATH}/MassTraceExtractor -test -ini ${DATA_DIR_TOPP}/MassTraceExtractor_2.ini.json -in ${DATA_DIR_TOPP}/MassTraceExtractor_1_input.mzML -out MassTraceExtractor_2_json.tmp.featureXML)
add_test("TOPP_MassTraceExtractor_2_out_json" ${DIFF} -whitelist "id=" -in1 MassTraceExtractor_2_json.tmp.featureXML -in2 ${DATA_DIR_TOPP}/MassTraceExtractor_2_output.featureXML )
set_tests_properties("TOPP_MassTraceExtractor_2_out_json" PROPERTIES DEPENDS "TOPP_MassTraceExtractor_2_json")

#------------------------------------------------------------------------------
# TOPP tests reading flat json input
#------------------------------------------------------------------------------
# MassTraceExtractor test
add_test("TOPP_MassTraceExtractor_1_flat_json" ${TOPP_BIN_PATH}/MassTraceExtractor -test -ini ${DATA_DIR_TOPP}/MassTraceExtractor.ini.flat.json -in ${DATA_DIR_TOPP}/MassTraceExtractor_1_input.mzML -out MassTraceExtractor_1_json.tmp)
add_test("TOPP_MassTraceExtractor_1_out_flat_json" ${DIFF} -whitelist "id=" -in1 MassTraceExtractor_1_json.tmp -in2 ${DATA_DIR_TOPP}/MassTraceExtractor_1_output.featureXML )
set_tests_properties("TOPP_MassTraceExtractor_1_out_flat_json" PROPERTIES DEPENDS "TOPP_MassTraceExtractor_1_flat_json")
add_test("TOPP_MassTraceExtractor_2_flat_json" ${TOPP_BIN_PATH}/MassTraceExtractor -test -ini ${DATA_DIR_TOPP}/MassTraceExtractor_2.ini.flat.json -in ${DATA_DIR_TOPP}/MassTraceExtractor_1_input.mzML -out MassTraceExtractor_2_json.tmp)
add_test("TOPP_MassTraceExtractor_2_out_flat_json" ${DIFF} -whitelist "id=" -in1 MassTraceExtractor_2_json.tmp -in2 ${DATA_DIR_TOPP}/MassTraceExtractor_2_output.featureXML )
set_tests_properties("TOPP_MassTraceExtractor_2_out_flat_json" PROPERTIES DEPENDS "TOPP_MassTraceExtractor_2_flat_json")

#------------------------------------------------------------------------------
# TOPP tests reading json and writing ini
#------------------------------------------------------------------------------
# MassTraceExtractor test
add_test("TOPP_MassTraceExtractor_3_flat_json" ${TOPP_BIN_PATH}/MassTraceExtractor -test -ini ${DATA_DIR_TOPP}/MassTraceExtractor_3.ini.flat.json -write_ini MassTraceExtractor_3.ini)
#add_test("TOPP_MassTraceExtractor_3_json_out" ${CMAKE_COMMAND} -E compare_files --ignore-eol ${DATA_DIR_TOPP}/MassTraceExtractor_3_expected.ini)
add_test("TOPP_MassTraceExtractor_3_flat_json_out" ${DIFF} -whitelist "Version of the tool that generated this parameters file." -in1 MassTraceExtractor_3.ini -in2 ${DATA_DIR_TOPP}/MassTraceExtractor_3_expected.ini )
set_tests_properties("TOPP_MassTraceExtractor_3_flat_json_out" PROPERTIES DEPENDS "TOPP_MassTraceExtractor_3_flat_json")

add_test("TOPP_MassTraceExtractor_3_json" ${TOPP_BIN_PATH}/MassTraceExtractor -test -ini ${DATA_DIR_TOPP}/MassTraceExtractor_3.ini.json -write_ini MassTraceExtractor_3.ini)
#add_test("TOPP_MassTraceExtractor_3_json_out" ${CMAKE_COMMAND} -E compare_files --ignore-eol ${DATA_DIR_TOPP}/MassTraceExtractor_3_expected.ini)
add_test("TOPP_MassTraceExtractor_3_json_out" ${DIFF} -whitelist "Version of the tool that generated this parameters file." -in1 MassTraceExtractor_3.ini -in2 ${DATA_DIR_TOPP}/MassTraceExtractor_3_expected.ini )
set_tests_properties("TOPP_MassTraceExtractor_3_json_out" PROPERTIES DEPENDS "TOPP_MassTraceExtractor_3_json")
