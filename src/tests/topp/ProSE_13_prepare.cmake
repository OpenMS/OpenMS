# Prepares the inputs of the TOPP_ProSE_13_merged test: the one test run is copied under two names
# (`cmake -E copy` writes a single destination file), so that ProSE pools two distinct input files.
#
# Usage: cmake -DSOURCE=<mzML> -DTARGET_A=<mzML> -DTARGET_B=<mzML> -P ProSE_13_prepare.cmake

if(NOT DEFINED SOURCE OR NOT DEFINED TARGET_A OR NOT DEFINED TARGET_B)
  message(FATAL_ERROR "SOURCE, TARGET_A and TARGET_B must be given")
endif()

foreach(target ${TARGET_A} ${TARGET_B})
  file(COPY_FILE "${SOURCE}" "${target}")
endforeach()
