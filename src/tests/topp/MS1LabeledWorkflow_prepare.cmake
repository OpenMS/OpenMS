# Prepares the input of the fractionated MS1LabeledWorkflow test: the experimental design names
# two fraction groups with two fractions each, so the one MS1 test run is copied under those four
# names. (The workflow matches design rows and input files by basename, so distinct names are needed.)
#
# Usage: cmake -DSOURCE=<mzML> -DTARGET_DIR=<dir> -P MS1LabeledWorkflow_prepare.cmake

if(NOT DEFINED SOURCE OR NOT DEFINED TARGET_DIR)
  message(FATAL_ERROR "SOURCE and TARGET_DIR must be given")
endif()

foreach(name A_F1 A_F2 B_F1 B_F2)
  configure_file("${SOURCE}" "${TARGET_DIR}/MS1LabeledWorkflow_${name}.mzML" COPYONLY)
endforeach()

file(MAKE_DIRECTORY "${TARGET_DIR}/MS1LabeledWorkflow_duplicate")
file(COPY "${SOURCE}" DESTINATION "${TARGET_DIR}/MS1LabeledWorkflow_duplicate")
get_filename_component(_source_basename "${SOURCE}" NAME)
# Only the second, identically named input is represented in the design. This used to select
# the first file silently during basename reconciliation.
file(WRITE "${TARGET_DIR}/MS1LabeledWorkflow_duplicate_design.tsv"
  "Fraction_Group\tFraction\tSpectra_Filepath\tLabel\tSample\n"
  "1\t1\t${TARGET_DIR}/MS1LabeledWorkflow_duplicate/${_source_basename}\t1\t1\n"
  "1\t1\t${TARGET_DIR}/MS1LabeledWorkflow_duplicate/${_source_basename}\t2\t2\n")
file(WRITE "${TARGET_DIR}/MS1LabeledWorkflow_ambiguous_design.tsv"
  "Fraction_Group\tFraction\tSpectra_Filepath\tLabel\tSample\n"
  "1\t1\t${SOURCE}\t1\t1\n"
  "1\t1\t${SOURCE}\t2\t2\n"
  "2\t1\t${TARGET_DIR}/MS1LabeledWorkflow_duplicate/${_source_basename}\t1\t3\n"
  "2\t1\t${TARGET_DIR}/MS1LabeledWorkflow_duplicate/${_source_basename}\t2\t4\n")
