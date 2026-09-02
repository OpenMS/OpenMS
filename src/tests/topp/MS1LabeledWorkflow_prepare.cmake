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
