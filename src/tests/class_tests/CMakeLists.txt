# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------


project("OpenMS_class_tests")

#------------------------------------------------------------------------------
# set new CMAKE_RUNTIME_OUTPUT_DIRECTORY for tests and remember old setting
set(_TMP_CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

#------------------------------------------------------------------------------
# The libOpenMS half of the class-test setup: a fixed UniqueIdGenerator seed and an
# OpenMS exception translator, registered at static init (the source says what breaks
# without it). Every class test below that links libOpenMS links this; the OpenSwathAlgo
# tests deliberately do not, which is why this cannot live in src/testframework -- that
# directory has to stay configurable without libOpenMS, and src/tests/header_sets builds
# it on its own to prove it.
#
# An OBJECT library, not STATIC, and that is the point: this translation unit exports no
# name that anyone references, so as a static archive member the linker would drop it and
# the seed would silently never be set -- the very nondeterminism it exists to prevent,
# and harder to notice than forgetting to list a source file. An object library puts its
# objects into the consumer unconditionally.
add_library(OpenMSTestSupport OBJECT source/OpenMSTestSupport.cpp)
target_link_libraries(OpenMSTestSupport PUBLIC OpenMSTestFramework OpenMS)
openms_add_library_compiler_flags(OpenMSTestSupport)
# The test projects below build test code at -O0, which also drops the NDEBUG that
# CMAKE_CXX_FLAGS_RELEASE carries. This translation unit ends up in the same
# executables and includes the same headers, so it has to agree with them: otherwise
# an NDEBUG-guarded inline would differ between the two and be an ODR mismatch with
# no diagnostic. It matters for more than optimization -- on MSVC, OPENMS_ASSERTIONS
# (config.h) is keyed on NDEBUG, and this file's OPENMS_ASSERTIONS branch decides
# whether precondition tests run at all. Set on the target, because the directory
# scope here also covers the OpenSwathAlgo tests, which keep the production flags.
# One quoted generator expression per flag: an unquoted one containing a space is split
# into two arguments by the CMake language, and a fragment that is no longer a complete
# expression is passed to the compiler verbatim rather than diagnosed.
if(CMAKE_COMPILER_IS_INTELCXX OR CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG)
  target_compile_options(OpenMSTestSupport PRIVATE
    "$<$<CONFIG:Release>:-O0>" "$<$<CONFIG:Release>:-UNDEBUG>")
endif()

#------------------------------------------------------------------------------
# add OpenSWATHAlgo specific tests
add_subdirectory(openswathalgo)

#------------------------------------------------------------------------------
# add OpenMS specific tests
add_subdirectory(openms)
add_subdirectory(openms_cli)

if(WITH_GUI)
  #------------------------------------------------------------------------------
  # add OpenMS_GUI specific tests
  add_subdirectory(openms_gui)
endif()

#------------------------------------------------------------------------------
# restore old CMAKE_RUNTIME_OUTPUT_DIRECTORY
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${_TMP_CMAKE_RUNTIME_OUTPUT_DIRECTORY})
