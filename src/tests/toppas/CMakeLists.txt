project("OpenMS_TOPPAS_tests")
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

################################## README ######################################
#
# This file contains the declarations of the pipeline tests executed by CTest
#
# It is included when -D ENABLE_PIPELINE_TESTING is 'On'.
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

if(WITH_GUI) # ExecutePipeline requires Qt/GUI internally

  ### General settings
  set(TOPP_BIN_PATH ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
  set(DATA_DIR_TOPPAS ${PROJECT_SOURCE_DIR})
  set(DATA_DIR_TOPPAS_BIN ${PROJECT_BINARY_DIR})
  set(DIFF ${TOPP_BIN_PATH}/FuzzyDiff -test  -whitelist "IdentificationRun " "SearchParameters " "xml-stylesheet")


  if(HAS_XSERVER)
  #------------------------------------------------------------------------------
  # configure resource files (for TOPPAS/ExecutePipeline) - there is no other
  # way to specify the input
  configure_file(${DATA_DIR_TOPPAS}/ExecutePipeline_1.trf.in ${DATA_DIR_TOPPAS_BIN}/ExecutePipeline_1.trf)
  # ExecutePipeline tests (as substitute for TOPPAS) - the ResourceFiles are in
  # binary tree, as they have been configured from a .in file (see above)!
  add_test("TOPP_ExecutePipeline_1" ${TOPP_BIN_PATH}/ExecutePipeline -test -in ${DATA_DIR_TOPPAS}/ExecutePipeline_1.toppas -resource_file ${DATA_DIR_TOPPAS_BIN}/ExecutePipeline_1.trf -out_dir .)
  # do not test the output -- we just want the pipeline to run -- the tools
  # itself are tested separately
    
  ## Run the examples in <openms>/share/OpenMS/examples, to make sure they are up to date
  file(GLOB TOPPAS_EXAMPLE_PIPELINES ${CF_OPENMS_DATA_PATH}/examples/TOPPAS/*.toppas)
  foreach(I_PIPELINE_PATH ${TOPPAS_EXAMPLE_PIPELINES})
    get_filename_component(I_FILE "${I_PIPELINE_PATH}" NAME)
    if(I_FILE MATCHES "GenericWrapper_") # requires rscript as dependency...
      message(STATUS "Skipping  ${I_FILE}")
      continue()
    endif()
    message(STATUS "Adding TOPPAS example pipeline for testing: ${I_FILE} [WITH_GUI + ENABLE_PIPELINE_TESTING]")
    add_test("TOPP_ExecutePipeline_${I_FILE}" ${TOPP_BIN_PATH}/ExecutePipeline -test -in "${I_PIPELINE_PATH}" -out_dir .  -num_jobs 4)
    ## allow failure if search engines (comet.exe etc) are not found
    set_tests_properties("TOPP_ExecutePipeline_${I_FILE}" PROPERTIES SKIP_RETURN_CODE 14) ## EXTERNAL_PROGRAM_NOTFOUND
  endforeach()
  
  
    
  ################### Labelfree quantification with IDMapping ####################

  ## Example usage (no actual test or data currently available)
  
  ## Test workflow
  # set(WF_NAME "LF_Quant")
  # add_test("TOPPAS_${WF_NAME}" ${TOPP_BIN_PATH}/ExecutePipeline -in ${DATA_DIR_TOPPAS}/${WF_NAME}/labelfree_quant_only.toppas -out_dir ${DATA_DIR_TOPPAS}/${WF_NAME})
  # set_tests_properties("TOPPAS_${WF_NAME}" PROPERTIES TIMEOUT 1800)

  ## Check output
  # set(OF_PATH "003-FeatureFinderCentroided/labelfree_1.featureXML")
  # add_test("TOPPAS_${WF_NAME}_out1" ${DIFF} -in1 ${DATA_DIR_TOPPAS}/${WF_NAME}/TOPPAS_out/${OF_PATH} -in2 ${DATA_DIR_TOPPAS}/${WF_NAME}/desired_output/${OF_PATH})
  # set_tests_properties("TOPPAS_${WF_NAME}_out1" PROPERTIES DEPENDS "TOPPAS_${WF_NAME}")

  # set(OF_PATH "003-FeatureFinderCentroided/labelfree_2.featureXML")
  # add_test("TOPPAS_${WF_NAME}_out2" ${DIFF} -in1 ${DATA_DIR_TOPPAS}/${WF_NAME}/TOPPAS_out/${OF_PATH} -in2 ${DATA_DIR_TOPPAS}/${WF_NAME}/desired_output/${OF_PATH})
  # set_tests_properties("TOPPAS_${WF_NAME}_out2" PROPERTIES DEPENDS "TOPPAS_${WF_NAME}")

  # set(OF_PATH "003-FeatureFinderCentroided/labelfree_3.featureXML")
  # add_test("TOPPAS_${WF_NAME}_out3" ${DIFF} -in1 ${DATA_DIR_TOPPAS}/${WF_NAME}/TOPPAS_out/${OF_PATH} -in2 ${DATA_DIR_TOPPAS}/${WF_NAME}/desired_output/${OF_PATH})
  # set_tests_properties("TOPPAS_${WF_NAME}_out3" PROPERTIES DEPENDS "TOPPAS_${WF_NAME}")

  # set(OF_PATH "010-ConsensusMapNormalizer/labelfree_1.consensusXML")
  # add_test("TOPPAS_${WF_NAME}_out4" ${DIFF} -in1 ${DATA_DIR_TOPPAS}/${WF_NAME}/TOPPAS_out/${OF_PATH} -in2 ${DATA_DIR_TOPPAS}/${WF_NAME}/desired_output/${OF_PATH})
  # set_tests_properties("TOPPAS_${WF_NAME}_out4" PROPERTIES DEPENDS "TOPPAS_${WF_NAME}")

  # set(OF_PATH "011-FeatureLinkerUnlabeledQT/labelfree_1.consensusXML")
  # add_test("TOPPAS_${WF_NAME}_out5" ${DIFF} -in1 ${DATA_DIR_TOPPAS}/${WF_NAME}/TOPPAS_out/${OF_PATH} -in2 ${DATA_DIR_TOPPAS}/${WF_NAME}/desired_output/${OF_PATH})
  # set_tests_properties("TOPPAS_${WF_NAME}_out5" PROPERTIES DEPENDS "TOPPAS_${WF_NAME}")

  endif()


################################################################################


endif()
