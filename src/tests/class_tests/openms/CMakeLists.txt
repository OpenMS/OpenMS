# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------


cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project("OpenMS_class_tests_openms")

#------------------------------------------------------------------------------
# Configure test file to get the TEST_DATA_PATH into the tests
set(CF_OPENMS_TEST_DATA_PATH "${PROJECT_SOURCE_DIR}/data/")
set (CONFIGURED_TEST_CONFIG_H "${PROJECT_BINARY_DIR}/include/OpenMS/test_config.h")
configure_file(${PROJECT_SOURCE_DIR}/include/OpenMS/test_config.h.in ${CONFIGURED_TEST_CONFIG_H})

#------------------------------------------------------------------------------
# get the test executables
include(executables.cmake)

#------------------------------------------------------------------------------
# Include directories for tests
set(OPENMS_CLASS_TESTS_INTERNAL_INCLUDE_DIRECTORIES "${PROJECT_BINARY_DIR}/include/")
# add OpenMS directories
set(OPENMS_CLASS_TESTS_EXTERNAL_INCLUDE_DIRECTORIES "${OpenMS_INCLUDE_DIRECTORIES};${Boost_INCLUDE_DIRS}")
include_directories(${OPENMS_CLASS_TESTS_INTERNAL_INCLUDE_DIRECTORIES})
include_directories(SYSTEM ${OPENMS_CLASS_TESTS_EXTERNAL_INCLUDE_DIRECTORIES})

#------------------------------------------------------------------------------
# disable optimization for tests for gcc like compilers
if (CMAKE_COMPILER_IS_INTELCXX OR CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG)
  set(_TMP_CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
  set(CMAKE_CXX_FLAGS_RELEASE "-O0")
endif()

#------------------------------------------------------------------------------
# Add the actual tests
foreach(_class_test ${TEST_executables})
  add_executable(${_class_test} source/${_class_test}.cpp)
  target_link_libraries(${_class_test} ${OpenMS_LIBRARIES})
  openms_add_executable_compiler_flags(${_class_test})
  add_test(${_class_test} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${_class_test})
  # only add OPENMP flags to gcc linker (except Mac OS X, due to compiler bug
  # see https://sourceforge.net/apps/trac/open-ms/ticket/280 for details)
  if (OPENMP_FOUND AND NOT MSVC AND NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(${_class_test} PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
  endif()
endforeach(_class_test)

#------------------------------------------------------------------------------
# some tests need special treatment

# StopWatch_test needs to be run alone to avoid timing errors
set_tests_properties(StopWatch_test PROPERTIES RUN_SERIAL 1)

# Percolator_subprocess_parity_test runs side-by-side comparisons against the
# external `percolator` binary. Tests that need it are gated via the
# PERCOLATOR_BINARY environment variable and skip silently when unset. Search
# PATH first, then fall back to the repo's bundled THIRDPARTY/<os>/<arch>/
# layout so local dev builds pick it up automatically.
set(_percolator_thirdparty_hints)
if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
  list(APPEND _percolator_thirdparty_hints
    "${OpenMS_SOURCE_DIR}/THIRDPARTY/Linux/${CMAKE_SYSTEM_PROCESSOR}/Percolator"
    "${CMAKE_SOURCE_DIR}/THIRDPARTY/Linux/${CMAKE_SYSTEM_PROCESSOR}/Percolator")
elseif (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  list(APPEND _percolator_thirdparty_hints
    "${OpenMS_SOURCE_DIR}/THIRDPARTY/MacOS/${CMAKE_SYSTEM_PROCESSOR}/Percolator"
    "${CMAKE_SOURCE_DIR}/THIRDPARTY/MacOS/${CMAKE_SYSTEM_PROCESSOR}/Percolator")
endif()
find_program(PERCOLATOR_BINARY_FOR_TEST percolator
  PATHS ENV PATH ${_percolator_thirdparty_hints})

if (PERCOLATOR_BINARY_FOR_TEST)
  set_tests_properties(Percolator_subprocess_parity_test PROPERTIES
    ENVIRONMENT "PERCOLATOR_BINARY=${PERCOLATOR_BINARY_FOR_TEST}")
  # Compute PercolatorAdapter path at configure time.
  # TOPP tools are placed in ${CMAKE_BINARY_DIR}/bin/ (the root build dir's bin/),
  # which is the same directory that the top-level CMakeLists sets as
  # CMAKE_RUNTIME_OUTPUT_DIRECTORY. Generator expressions like $<TARGET_FILE:...>
  # cannot be evaluated in ENVIRONMENT properties across sub-project boundaries.
  if (WIN32)
    set(_percolator_adapter_suffix ".exe")
  else()
    set(_percolator_adapter_suffix "")
  endif()
  set(_percolator_adapter_path "${CMAKE_BINARY_DIR}/bin/PercolatorAdapter${_percolator_adapter_suffix}")
  set_tests_properties(PercolatorAdapter_parity_test PROPERTIES
    ENVIRONMENT "PERCOLATOR_BINARY=${PERCOLATOR_BINARY_FOR_TEST};PERCOLATOR_ADAPTER_BINARY=${_percolator_adapter_path}")
  message(STATUS "Percolator_subprocess_parity_test: enabled (found ${PERCOLATOR_BINARY_FOR_TEST})")
  message(STATUS "PercolatorAdapter_parity_test: enabled (adapter at ${_percolator_adapter_path})")
else()
  message(STATUS "Percolator_subprocess_parity_test: subprocess sections skipped (percolator binary not on PATH nor in THIRDPARTY)")
  message(STATUS "PercolatorAdapter_parity_test: subprocess sections skipped (percolator binary not on PATH nor in THIRDPARTY)")
endif()

#------------------------------------------------------------------------------
# 2 - add link dependencies (TODO could be more finegrained for boost)
foreach(t ${Boost_dependent_tests})
  foreach(comp ${OpenMS_BOOST_COMPONENTS})
    target_link_libraries(${t} Boost::${comp})
  endforeach()
endforeach()

target_link_libraries(Base64_test ZLIB::ZLIB)
target_link_libraries(Libzip_test libzip::zip)
target_link_libraries(ZipRandomAccessFile_test libzip::zip)
target_link_libraries(LPWrapper_test ${LPTARGET})
target_link_libraries(SpectraSTSimilarityScore_test Eigen3::Eigen)
target_link_libraries(BinnedSpectrum_test Eigen3::Eigen)
target_link_libraries(BilinearInterpolation_test Eigen3::Eigen)
target_link_libraries(MatrixEigen_test Eigen3::Eigen)
target_link_libraries(MRMScoring_test Eigen3::Eigen)


if (WITH_HDF5)
  target_link_libraries(HDF5_test HDF5::HDF5)
endif()

target_link_libraries(Arrow_test
  ${OPENMS_ARROW_TARGET}
  ${OPENMS_PARQUET_TARGET}
)
target_link_libraries(MSExperimentArrowExport_test
  ${OPENMS_ARROW_TARGET}
  ${OPENMS_PARQUET_TARGET}
)
target_link_libraries(ConsensusMapArrowExport_test
  ${OPENMS_ARROW_TARGET}
  ${OPENMS_PARQUET_TARGET}
)
target_link_libraries(QPXFile_test
  ${OPENMS_ARROW_TARGET}
  ${OPENMS_PARQUET_TARGET}
)
target_link_libraries(ProteinIdentificationArrowIO_test
  ${OPENMS_ARROW_TARGET}
  ${OPENMS_PARQUET_TARGET}
)
target_link_libraries(FeatureMapArrowIO_test
  ${OPENMS_ARROW_TARGET}
  ${OPENMS_PARQUET_TARGET}
)
target_link_libraries(ConsensusMapArrowIO_test
  ${OPENMS_ARROW_TARGET}
  ${OPENMS_PARQUET_TARGET}
)
target_link_libraries(Libzip_test
  ${OPENMS_ARROW_TARGET}
  ${OPENMS_PARQUET_TARGET}
)
target_link_libraries(ZipRandomAccessFile_test
  ${OPENMS_ARROW_TARGET}
  ${OPENMS_PARQUET_TARGET}
)

if (TARGET TransitionParquetFile_test)
  target_link_libraries(TransitionParquetFile_test
    ${OPENMS_ARROW_TARGET}
    ${OPENMS_PARQUET_TARGET}
  )
endif()
if (TARGET OpenSwathOSWParquetReader_test)
  target_link_libraries(OpenSwathOSWParquetReader_test
    ${OPENMS_ARROW_TARGET}
    ${OPENMS_PARQUET_TARGET}
  )
endif()
if (TARGET OpenSwathOSWParquetWriter_test)
  target_link_libraries(OpenSwathOSWParquetWriter_test
    ${OPENMS_ARROW_TARGET}
    ${OPENMS_PARQUET_TARGET}
  )
endif()
if (TARGET OpenSwathOSWParquetRoundTrip_test)
  target_link_libraries(OpenSwathOSWParquetRoundTrip_test
    ${OPENMS_ARROW_TARGET}
    ${OPENMS_PARQUET_TARGET}
  )
endif()
if (TARGET OpenSwathExport_test)
  target_link_libraries(OpenSwathExport_test
    ${OPENMS_ARROW_TARGET}
    ${OPENMS_PARQUET_TARGET}
  )
endif()
if (TARGET ArrowSchemaRegistry_test)
  target_link_libraries(ArrowSchemaRegistry_test
    ${OPENMS_ARROW_TARGET}
    ${OPENMS_PARQUET_TARGET}
  )
endif()

if (TARGET BrukerTimsFile_test AND WITH_OPENTIMS)
  target_link_libraries(BrukerTimsFile_test opentims::opentims_cpp sqlite3)
endif()

if (TARGET BrukerTimsFile_test AND OPENTIMS_DDA_TEST_DATA)
  target_compile_definitions(BrukerTimsFile_test PRIVATE
    OPENTIMS_DDA_TEST_DATA="${OPENTIMS_DDA_TEST_DATA}"
    OPENTIMS_DIA_TEST_DATA="${OPENTIMS_DIA_TEST_DATA}"
    OPENTIMS_TEST_FASTA="${OPENTIMS_TEST_FASTA}")
  # This test loads real Bruker .d files (DDA + DIA) including frame aggregation,
  # a peptide search engine run, and three-way (greedy vs. hill) centroiding
  # comparisons (which add 4 extra full-fixture loads), so it needs a generous timeout.
  set_tests_properties(BrukerTimsFile_test PROPERTIES TIMEOUT 3600)
endif()

if (TARGET ThermoRawFile_test AND WITH_THERMO_RAW AND THERMO_RAW_TEST_DATA)
  target_compile_definitions(ThermoRawFile_test PRIVATE
    THERMO_RAW_TEST_DATA="${THERMO_RAW_TEST_DATA}")
  set_tests_properties(ThermoRawFile_test PROPERTIES TIMEOUT 1200)

  # Copy the managed bridge runtime (.dll + .runtimeconfig.json) next to the
  # test binary. The bridge looks first in the executable's own directory; on
  # Windows the test exe is in src/tests/class_tests/bin (not next to
  # libOpenMS), so the OpenMS-target POST_BUILD copy in src/openms/CMakeLists
  # is in a different directory and would not be found at test runtime.
  if (COMMAND openms_thermo_bridge_copy_runtime_files)
    openms_thermo_bridge_copy_runtime_files(TARGET ThermoRawFile_test)
  endif()
endif()

#------------------------------------------------------------------------------
# restore old compiler flags
if (CMAKE_COMPILER_IS_INTELCXX OR CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG)
  set(${CMAKE_CXX_FLAGS_RELEASE} ${_TMP_CMAKE_CXX_FLAGS_RELEASE})
endif()

#------------------------------------------------------------------------------
# add filenames to Visual Studio solution tree
set(sources_VS)
foreach(i ${TEST_executables})
  list(APPEND sources_VS "${i}.cpp")
endforeach(i)
source_group("" FILES ${sources_VS})
