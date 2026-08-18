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
# MacrosTest.h) as its own library, so that the tests of every OpenMS library
# (OpenSwathAlgo, OpenMS, OpenMS_GUI and future split-off libraries) can use it
# without the framework being part of libOpenMS itself.
#
# The framework depends on the C++ standard library ONLY -- no libOpenMS, no
# third-party headers. Library-specific behavior (exception naming, unique-ID
# seeding) is registered by the test projects, see
# src/tests/class_tests/openms/source/OpenMSTestSupport.cpp. Because the
# archive references no libOpenMS symbols, no link-order or --as-needed
# precautions are needed on any consumer's link line.
#
# The library is deliberately STATIC:
#  - the framework is per-executable global state, so there is no point in
#    sharing it between processes,
#  - its symbols are no longer part of the libOpenMS ABI, and consumers
#    (test executables, FuzzyDiff) simply embed it,
#  - no export/import decorations (OPENMS_DLLAPI) are needed on Windows.

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
# The framework is an in-repo development tool, not part of the shipped
# product: the archive is EXCLUDE_FROM_ALL, so default installs, the
# monolithic installers (NSIS/zip) and 'make install' skip it; it remains
# installable explicitly via `cmake --install . --component OpenMSTestFramework`.
# (package_mac_productbuild.cmake additionally drops the component from its
# per-component package list.)
#
# The headers keep their historic install location
# <prefix>/include/OpenMS/CONCEPT/ (they shipped inside the OpenMS header set
# before the framework was split out), so existing includes of
# <OpenMS/CONCEPT/ClassTest.h> still resolve against an installed OpenMS.
#
# The target is deliberately NOT part of the exported OpenMS package
# (OpenMSTargets.cmake) and NOT in the 'library' component: minimal install
# flows (e.g. the pyOpenMS wheel builds) build only 'ninja OpenMS OpenSwathAlgo'
# and install component-wise -- an exported test-framework target would make
# their installed OpenMSTargets.cmake reference an archive that was never
# built/installed and break every downstream find_package(OpenMS).
# All consumers (class tests, FuzzyDiff, tests of future split-off libraries)
# are in-tree and link the target directly.
install(TARGETS OpenMSTestFramework
        ARCHIVE DESTINATION ${INSTALL_LIB_DIR} COMPONENT OpenMSTestFramework
        EXCLUDE_FROM_ALL)
install_headers("${OpenMSTestFramework_sources_h}" OpenMSTestFramework)

openms_doc_path("${PROJECT_SOURCE_DIR}/include")
