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
# It still links libOpenMS for low-level utilities (StringUtils, Exception,
# DataValue, ...); once libOpenMS is split into smaller libraries, this
# dependency shrinks to the future core library.
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
  include/OpenMS/CONCEPT/FuzzyStringComparator.h
  include/OpenMS/CONCEPT/MacrosTest.h
)

source_group("Source Files\\CONCEPT" FILES ${OpenMSTestFramework_sources})
source_group("Header Files\\OpenMS\\CONCEPT" FILES ${OpenMSTestFramework_sources_h})

add_library(OpenMSTestFramework STATIC ${OpenMSTestFramework_sources} ${OpenMSTestFramework_sources_h})

target_include_directories(OpenMSTestFramework PUBLIC
                           "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
                           "$<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>")

# low-level utilities (StringUtils, Exception, DataValue, ...) come from libOpenMS
target_link_libraries(OpenMSTestFramework PUBLIC OpenMS)

target_compile_features(OpenMSTestFramework PUBLIC cxx_std_23)

openms_add_library_compiler_flags(OpenMSTestFramework)

#------------------------------------------------------------------------------
# install the archive and its headers; the headers keep their historic install
# location <prefix>/include/OpenMS/CONCEPT/, so existing includes of
# <OpenMS/CONCEPT/ClassTest.h> keep working against an installed OpenMS
install(TARGETS OpenMSTestFramework
        EXPORT ${OPENMS_EXPORT_SET}
        ARCHIVE DESTINATION ${INSTALL_LIB_DIR} COMPONENT library)
install_headers("${OpenMSTestFramework_sources_h}" OpenMSTestFramework)

#------------------------------------------------------------------------------
# register for export (OpenMSTargets), so find_package(OpenMS) consumers can
# link the test framework of an installed OpenMS
openms_register_export_target(OpenMSTestFramework)

openms_doc_path("${PROJECT_SOURCE_DIR}/include")
