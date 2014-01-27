################################## README ######################################
#
# This file contains the declarations of the pipeline tests executed by CTest
#
# In order to add a test, do the following:
#
# - Add a subfolder named 'WorkflowName' to source/TEST/TOPPAS
#
# - Put your workflow there, together with the input data
#   (please use the same directory structure as in the existing workflow tests)
#
# - Add the corresponding test declarations below.
#
# Please note:
#
# Test names (the first argument of 'add_test') must be distinct.
# If not, the first definition of the test is overwritten!
#
# Therefore, we use the following naming scheme for TOPPAS pipeline tests:
#
#
#                      TOPPAS_<workflow>[_out<num>]
#
#
# 'TOPPAS_'   This prefix is used to distinguish TOPPAS, TOPP, and OpenMS tests
#             in the nightly builds.
#
# <workflow>  The workflow name used to identify failed TOPP tool tests during
#             the nightly builds.
#
# '_out'      When a result file of a workflow is checked with FuzzyDiff, this
#             suffix is appended.
#
# <num>       If n output files of a workflow need to be checked, the checks are
#             numbered from 1 to n.
#
# When the result file of a test is checked with FuzzyDiff, we need to set the
# dependency to the test that generates the output file by
#
# set_tests_properties("TOPPAS_<workflow>_out<num>" PROPERTIES DEPENDS "TOPPAS_<workflow>")
#
################################################################################


### General settings
set(TOPP_BIN_PATH ${OPENMS_BINARY_DIR})
set(DATA_DIR_TOPP ${PROJECT_SOURCE_DIR}/TOPP)
set(DATA_DIR_TOPPAS ${PROJECT_SOURCE_DIR}/TOPPAS)
set(DIFF ${TOPP_BIN_PATH}/FuzzyDiff -ini ${DATA_DIR_TOPP}/FuzzyDiff.ini -whitelist "IdentificationRun " "SearchParameters " "xml-stylesheet")


################### Labelfree quantification with IDMapping ####################

# Test workflow
set(WF_NAME "LF_Quant")
add_test("TOPPAS_${WF_NAME}" ${TOPP_BIN_PATH}/ExecutePipeline -in ${DATA_DIR_TOPPAS}/${WF_NAME}/labelfree_quant_only.toppas -out_dir ${DATA_DIR_TOPPAS}/${WF_NAME})
set_tests_properties("TOPPAS_${WF_NAME}" PROPERTIES TIMEOUT 1800)

# Check output
set(OF_PATH "003-FeatureFinderCentroided/labelfree_1.featureXML")
add_test("TOPPAS_${WF_NAME}_out1" ${DIFF} -in1 ${DATA_DIR_TOPPAS}/${WF_NAME}/TOPPAS_out/${OF_PATH} -in2 ${DATA_DIR_TOPPAS}/${WF_NAME}/desired_output/${OF_PATH})
set_tests_properties("TOPPAS_${WF_NAME}_out1" PROPERTIES DEPENDS "TOPPAS_${WF_NAME}")

set(OF_PATH "003-FeatureFinderCentroided/labelfree_2.featureXML")
add_test("TOPPAS_${WF_NAME}_out2" ${DIFF} -in1 ${DATA_DIR_TOPPAS}/${WF_NAME}/TOPPAS_out/${OF_PATH} -in2 ${DATA_DIR_TOPPAS}/${WF_NAME}/desired_output/${OF_PATH})
set_tests_properties("TOPPAS_${WF_NAME}_out2" PROPERTIES DEPENDS "TOPPAS_${WF_NAME}")

set(OF_PATH "003-FeatureFinderCentroided/labelfree_3.featureXML")
add_test("TOPPAS_${WF_NAME}_out3" ${DIFF} -in1 ${DATA_DIR_TOPPAS}/${WF_NAME}/TOPPAS_out/${OF_PATH} -in2 ${DATA_DIR_TOPPAS}/${WF_NAME}/desired_output/${OF_PATH})
set_tests_properties("TOPPAS_${WF_NAME}_out3" PROPERTIES DEPENDS "TOPPAS_${WF_NAME}")

set(OF_PATH "010-ConsensusMapNormalizer/labelfree_1.consensusXML")
add_test("TOPPAS_${WF_NAME}_out4" ${DIFF} -in1 ${DATA_DIR_TOPPAS}/${WF_NAME}/TOPPAS_out/${OF_PATH} -in2 ${DATA_DIR_TOPPAS}/${WF_NAME}/desired_output/${OF_PATH})
set_tests_properties("TOPPAS_${WF_NAME}_out4" PROPERTIES DEPENDS "TOPPAS_${WF_NAME}")

set(OF_PATH "011-FeatureLinkerUnlabeledQT/labelfree_1.consensusXML")
add_test("TOPPAS_${WF_NAME}_out5" ${DIFF} -in1 ${DATA_DIR_TOPPAS}/${WF_NAME}/TOPPAS_out/${OF_PATH} -in2 ${DATA_DIR_TOPPAS}/${WF_NAME}/desired_output/${OF_PATH})
set_tests_properties("TOPPAS_${WF_NAME}_out5" PROPERTIES DEPENDS "TOPPAS_${WF_NAME}")


################# Labelfree quantification including OMSSA ID ##################

if (NOT (${OMSSA_BINARY} STREQUAL "OMSSA_BINARY-NOTFOUND"))

  # Test workflow
  set(WF_NAME "LF_Quant_and_OMSSA_ID")
  add_test("TOPPAS_${WF_NAME}" ${TOPP_BIN_PATH}/ExecutePipeline -in ${DATA_DIR_TOPPAS}/${WF_NAME}/labelfree_omssa.toppas -out_dir ${DATA_DIR_TOPPAS}/${WF_NAME})
  set_tests_properties("TOPPAS_${WF_NAME}" PROPERTIES TIMEOUT 1800)
  
  # Check output
  set(OF_PATH "004-FeatureFinderCentroided/labelfree_1.featureXML")
  add_test("TOPPAS_${WF_NAME}_out1" ${DIFF} -in1 ${DATA_DIR_TOPPAS}/${WF_NAME}/TOPPAS_out/${OF_PATH} -in2 ${DATA_DIR_TOPPAS}/${WF_NAME}/desired_output/${OF_PATH})
  set_tests_properties("TOPPAS_${WF_NAME}_out1" PROPERTIES DEPENDS "TOPPAS_${WF_NAME}")
  
  set(OF_PATH "004-FeatureFinderCentroided/labelfree_2.featureXML")
  add_test("TOPPAS_${WF_NAME}_out2" ${DIFF} -in1 ${DATA_DIR_TOPPAS}/${WF_NAME}/TOPPAS_out/${OF_PATH} -in2 ${DATA_DIR_TOPPAS}/${WF_NAME}/desired_output/${OF_PATH})
  set_tests_properties("TOPPAS_${WF_NAME}_out2" PROPERTIES DEPENDS "TOPPAS_${WF_NAME}")
  
  set(OF_PATH "004-FeatureFinderCentroided/labelfree_3.featureXML")
  add_test("TOPPAS_${WF_NAME}_out3" ${DIFF} -in1 ${DATA_DIR_TOPPAS}/${WF_NAME}/TOPPAS_out/${OF_PATH} -in2 ${DATA_DIR_TOPPAS}/${WF_NAME}/desired_output/${OF_PATH})
  set_tests_properties("TOPPAS_${WF_NAME}_out3" PROPERTIES DEPENDS "TOPPAS_${WF_NAME}")
  
  set(OF_PATH "011-IDFilter/labelfree_1.idXML")
  add_test("TOPPAS_${WF_NAME}_out4" ${DIFF} -in1 ${DATA_DIR_TOPPAS}/${WF_NAME}/TOPPAS_out/${OF_PATH} -in2 ${DATA_DIR_TOPPAS}/${WF_NAME}/desired_output/${OF_PATH})
  set_tests_properties("TOPPAS_${WF_NAME}_out4" PROPERTIES DEPENDS "TOPPAS_${WF_NAME}")
  
  set(OF_PATH "011-IDFilter/labelfree_2.idXML")
  add_test("TOPPAS_${WF_NAME}_out5" ${DIFF} -in1 ${DATA_DIR_TOPPAS}/${WF_NAME}/TOPPAS_out/${OF_PATH} -in2 ${DATA_DIR_TOPPAS}/${WF_NAME}/desired_output/${OF_PATH})
  set_tests_properties("TOPPAS_${WF_NAME}_out5" PROPERTIES DEPENDS "TOPPAS_${WF_NAME}")
  
  set(OF_PATH "011-IDFilter/labelfree_3.idXML")
  add_test("TOPPAS_${WF_NAME}_out6" ${DIFF} -in1 ${DATA_DIR_TOPPAS}/${WF_NAME}/TOPPAS_out/${OF_PATH} -in2 ${DATA_DIR_TOPPAS}/${WF_NAME}/desired_output/${OF_PATH})
  set_tests_properties("TOPPAS_${WF_NAME}_out6" PROPERTIES DEPENDS "TOPPAS_${WF_NAME}")
  
  set(OF_PATH "015-ConsensusMapNormalizer/labelfree_1.consensusXML")
  add_test("TOPPAS_${WF_NAME}_out7" ${DIFF} -in1 ${DATA_DIR_TOPPAS}/${WF_NAME}/TOPPAS_out/${OF_PATH} -in2 ${DATA_DIR_TOPPAS}/${WF_NAME}/desired_output/${OF_PATH})
  set_tests_properties("TOPPAS_${WF_NAME}_out7" PROPERTIES DEPENDS "TOPPAS_${WF_NAME}")
  
  set(OF_PATH "016-FeatureLinkerUnlabeledQT/labelfree_1.consensusXML")
  add_test("TOPPAS_${WF_NAME}_out8" ${DIFF} -in1 ${DATA_DIR_TOPPAS}/${WF_NAME}/TOPPAS_out/${OF_PATH} -in2 ${DATA_DIR_TOPPAS}/${WF_NAME}/desired_output/${OF_PATH})
  set_tests_properties("TOPPAS_${WF_NAME}_out8" PROPERTIES DEPENDS "TOPPAS_${WF_NAME}")

endif()

################################################################################

