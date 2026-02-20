# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# algo_includes.cmake — aggregates sources for OpenMS_Algo library
# Algo modules: ML, ANALYSIS, PROCESSING, FEATUREFINDER, COMPARISON, QC, APPLICATIONS, OPENSWATH
#
# Note: ${CMAKE_CURRENT_LIST_DIR} is needed because this file is included from
# a subdirectory CMakeLists.txt whose CMAKE_CURRENT_SOURCE_DIR differs.

set(OpenMS_sources CACHE INTERNAL "")

## Source files
include(${CMAKE_CURRENT_LIST_DIR}/source/ML/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/ML/NNLS/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/ML/SVM/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/ML/CLUSTERING/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/ML/GRIDSEARCH/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/ML/CROSSVALIDATION/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/ML/INTERPOLATION/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/ML/ROCCURVE/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/ML/RANSAC/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/ML/REGRESSION/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/ANALYSIS/QUANTITATION/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/ANALYSIS/SEQUENCE/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/ANALYSIS/MAPMATCHING/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/ANALYSIS/DECHARGING/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/ANALYSIS/ID/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/ANALYSIS/MRM/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/ANALYSIS/NUXL/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/ANALYSIS/TARGETED/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/ANALYSIS/TOPDOWN/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/ANALYSIS/XLMS/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/PROCESSING/BASELINE/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/PROCESSING/FILTERING/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/PROCESSING/RESAMPLING/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/PROCESSING/MISC/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/PROCESSING/DEISOTOPING/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/PROCESSING/CALIBRATION/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/PROCESSING/FEATURE/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/PROCESSING/SMOOTHING/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/PROCESSING/NOISEESTIMATION/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/PROCESSING/ID/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/PROCESSING/CENTROIDING/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/PROCESSING/SPECTRAMERGING/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/PROCESSING/SCALING/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/FEATUREFINDER/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/COMPARISON/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/source/QC/sources.cmake)
if(NOT DISABLE_OPENSWATH)
  include(${CMAKE_CURRENT_LIST_DIR}/source/ANALYSIS/OPENSWATH/sources.cmake)
  include(${CMAKE_CURRENT_LIST_DIR}/source/ANALYSIS/OPENSWATH/DATAACCESS/sources.cmake)
endif()
include(${CMAKE_CURRENT_LIST_DIR}/source/APPLICATIONS/sources.cmake)

# IO (FORMAT) source files compiled here because they depend on Algo types.
# The files remain in their original FORMAT directories but are not listed in
# IO's sources.cmake — they appear only here.
list(APPEND OpenMS_sources
  source/FORMAT/FLASHDeconvFeatureFile.cpp    # uses PeakGroup, DeconvolvedSpectrum, FLASHHelperClasses
  source/FORMAT/FLASHDeconvSpectrumFile.cpp   # uses PeakGroup, DeconvolvedSpectrum, FLASHHelperClasses
)

set(OpenMS_Algo_sources ${OpenMS_sources} CACHE INTERNAL "OpenMS Algo source files")
set(OpenMS_sources CACHE INTERNAL "")

## Header files
set(OpenMS_sources_h CACHE INTERNAL "")

include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/ML/INTERPOLATION/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/ML/NNLS/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/ML/SVM/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/ML/RANSAC/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/ML/REGRESSION/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/ML/CLUSTERING/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/ML/GRIDSEARCH/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/ML/CROSSVALIDATION/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/ML/ROCCURVE/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/ANALYSIS/DECHARGING/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/ANALYSIS/ID/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/ANALYSIS/MAPMATCHING/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/ANALYSIS/QUANTITATION/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/ANALYSIS/SEQUENCE/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/ANALYSIS/MRM/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/ANALYSIS/NUXL/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/ANALYSIS/TARGETED/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/ANALYSIS/TOPDOWN/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/ANALYSIS/XLMS/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/COMPARISON/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/PROCESSING/BASELINE/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/PROCESSING/CALIBRATION/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/PROCESSING/MISC/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/PROCESSING/DEISOTOPING/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/PROCESSING/ID/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/PROCESSING/FEATURE/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/PROCESSING/NOISEESTIMATION/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/PROCESSING/SMOOTHING/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/PROCESSING/FILTERING/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/PROCESSING/RESAMPLING/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/PROCESSING/SCALING/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/FEATUREFINDER/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/PROCESSING/CENTROIDING/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/PROCESSING/SPECTRAMERGING/sources.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/QC/sources.cmake)
if(NOT DISABLE_OPENSWATH)
  include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/ANALYSIS/OPENSWATH/sources.cmake)
  include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/sources.cmake)
endif()
include(${CMAKE_CURRENT_LIST_DIR}/include/OpenMS/APPLICATIONS/sources.cmake)

# IO (FORMAT) headers for the source files above (parallel to source moves)
list(APPEND OpenMS_sources_h
  include/OpenMS/FORMAT/FLASHDeconvFeatureFile.h
  include/OpenMS/FORMAT/FLASHDeconvSpectrumFile.h
)

set(OpenMS_Algo_sources_h ${OpenMS_sources_h} CACHE INTERNAL "OpenMS Algo header files")
set(OpenMS_sources_h CACHE INTERNAL "")
