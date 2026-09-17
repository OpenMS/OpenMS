# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

project("OpenMSTestFramework")
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

#------------------------------------------------------------------------------
# The class-test framework (ClassTest.h macros, FuzzyStringComparator,
# MacrosTest.h) as its own STATIC library: standard-library only, no libOpenMS
# headers or symbols, so the tests of every OpenMS library (incl. ones that do
# not link libOpenMS, like OpenSwathAlgo) can use it. Library-specific behavior
# is registered by the test projects (see openms/source/OpenMSTestSupport.cpp).

set(OpenMSTestFramework_sources
  source/CONCEPT/ClassTest.cpp
  source/CONCEPT/FuzzyStringComparator.cpp
)

set(OpenMSTestFramework_sources_h
  include/OpenMS/CONCEPT/ClassTest.h
  include/OpenMS/CONCEPT/ClassTestUtils.h
  include/OpenMS/CONCEPT/FuzzyStringComparator.h
  include/OpenMS/CONCEPT/MacrosTest.h
)

source_group("Source Files\\CONCEPT" FILES ${OpenMSTestFramework_sources})
source_group("Header Files\\OpenMS\\CONCEPT" FILES ${OpenMSTestFramework_sources_h})

add_library(OpenMSTestFramework STATIC ${OpenMSTestFramework_sources} ${OpenMSTestFramework_sources_h})

target_include_directories(OpenMSTestFramework PUBLIC
                           "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
                           "$<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>")

target_compile_features(OpenMSTestFramework PUBLIC cxx_std_23)

openms_add_library_compiler_flags(OpenMSTestFramework)

#------------------------------------------------------------------------------
# In-repo dev tool, not shipped: EXCLUDE_FROM_ALL keeps the archive AND its
# headers out of default installs and the NSIS/zip installers (productbuild
# enumerates components by name, so it ignores the flag and drops both via an
# explicit list(REMOVE_ITEM) in cmake/package_mac_productbuild.cmake);
# `cmake --install . --component OpenMSTestFramework` (resp.
# `--component OpenMSTestFramework_headers`) still installs them. The headers are
# useless on their own anyway -- the target is not exported
# (OpenMSTargets.cmake) and not in the 'library' component, so an installation
# carries no library to link them against: the wheel builds install
# component-wise, so an exported target would dangle and break downstream
# find_package(OpenMS). Headers keep their historic include/OpenMS/CONCEPT/
# location so existing includes resolve. In-tree consumers link the target
# directly and read the headers from the source tree.
install(TARGETS OpenMSTestFramework
        ARCHIVE DESTINATION ${INSTALL_LIB_DIR} COMPONENT OpenMSTestFramework
        EXCLUDE_FROM_ALL)
install_headers("${OpenMSTestFramework_sources_h}" OpenMSTestFramework EXCLUDE_FROM_ALL)

openms_doc_path("${PROJECT_SOURCE_DIR}/include")
