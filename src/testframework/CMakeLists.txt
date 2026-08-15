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

# Low-level utilities (StringUtils, Exception, DataValue, ...) come from libOpenMS.
# Deliberately a *usage* (header) dependency instead of target_link_libraries(... OpenMS):
# a link edge would force CMake to order libOpenMS.so after libOpenMSTestFramework.a on
# every consumer's link line, which pushes it behind the static Arrow/Parquet archives
# that the Arrow-related class tests link. The linker would then satisfy the Arrow/Parquet
# symbols from those archives instead of from libOpenMS.so, pulling their objects (and the
# Thrift statics inside them) into each test executable -- 24x larger binaries and, before
# openms_hide_static_archive_symbols(), duplicated static objects that were destroyed twice
# at exit. Every consumer (class tests, FuzzyDiff) links OpenMS itself, and this target is
# not part of the installed package, so no link edge is needed here.
target_include_directories(OpenMSTestFramework PUBLIC
                           "$<TARGET_PROPERTY:OpenMS,INTERFACE_INCLUDE_DIRECTORIES>")
target_compile_definitions(OpenMSTestFramework PUBLIC
                           "$<TARGET_PROPERTY:OpenMS,INTERFACE_COMPILE_DEFINITIONS>")
add_dependencies(OpenMSTestFramework OpenMS)

# UniqueIdGenerator.h (seeding in ClassTest.cpp) includes boost/random headers.
# Boost::boost is header-only: include dirs only, no link-line entry, so this does
# not affect the link ordering discussed above.
target_link_libraries(OpenMSTestFramework PRIVATE Boost::boost)

target_compile_features(OpenMSTestFramework PUBLIC cxx_std_23)

openms_add_library_compiler_flags(OpenMSTestFramework)

#------------------------------------------------------------------------------
# install the archive and its headers under a dedicated component; the headers
# keep their historic install location <prefix>/include/OpenMS/CONCEPT/, so
# existing includes of <OpenMS/CONCEPT/ClassTest.h> keep working against a
# (fully) installed OpenMS.
#
# The target is deliberately NOT part of the exported OpenMS package
# (OpenMSTargets.cmake) and NOT in the 'library' component: minimal install
# flows (e.g. the pyOpenMS wheel builds) build only 'ninja OpenMS OpenSwathAlgo'
# and install component-wise -- an exported test-framework target would make
# their installed OpenMSTargets.cmake reference an archive that was never
# built/installed and break every downstream find_package(OpenMS).
# In-tree consumers (class tests, FuzzyDiff) link the target directly.
install(TARGETS OpenMSTestFramework
        ARCHIVE DESTINATION ${INSTALL_LIB_DIR} COMPONENT OpenMSTestFramework)
install_headers("${OpenMSTestFramework_sources_h}" OpenMSTestFramework)

openms_doc_path("${PROJECT_SOURCE_DIR}/include")
