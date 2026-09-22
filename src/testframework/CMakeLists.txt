# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

project("OpenMSTestFramework")
cmake_minimum_required(VERSION 3.24 FATAL_ERROR)

#------------------------------------------------------------------------------
# The class-test framework (ClassTest.h macros, FuzzyStringComparator,
# MacrosTest.h) as its own STATIC library: standard-library only, no libOpenMS
# headers or symbols, so the tests of every OpenMS library (incl. ones that do
# not link libOpenMS, like OpenSwathAlgo) can use it. Library-specific behavior
# is registered on top of it by OpenMSTestSupport below, which does link libOpenMS.

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

add_library(OpenMSTestFramework STATIC ${OpenMSTestFramework_sources})
target_sources(OpenMSTestFramework PUBLIC FILE_SET HEADERS
  BASE_DIRS "${PROJECT_SOURCE_DIR}/include"
  FILES ${OpenMSTestFramework_sources_h})
openms_validate_public_headers(${OpenMSTestFramework_sources_h})
if(OPENMS_VERIFY_INTERFACE_HEADER_SETS)
  set_target_properties(OpenMSTestFramework PROPERTIES VERIFY_INTERFACE_HEADER_SETS ON)
endif()

target_compile_features(OpenMSTestFramework PUBLIC cxx_std_23)

openms_add_library_compiler_flags(OpenMSTestFramework)

#------------------------------------------------------------------------------
# The libOpenMS half of the class-test setup: a fixed UniqueIdGenerator seed and an
# OpenMS exception translator, registered at static init (the source says what breaks
# without it). Every class test of a library that links libOpenMS needs it.
#
# An OBJECT library, not STATIC, and that is the point: this translation unit exports
# no name that anyone references, so as a static archive member the linker would drop
# it and the seed would silently never be set -- the very nondeterminism the file
# exists to prevent, and harder to notice than forgetting to list a source file.
# An object library puts its objects into the consumer unconditionally.
#
# Not installed and not exported; consumers link it in-tree. Exporting it later, to
# make a tool's test/ folder configurable out of tree, needs a STATIC library behind
# an INTERFACE target applying $<LINK_LIBRARY:WHOLE_ARCHIVE,...>, because an installed
# OBJECT library does not carry that guarantee to its consumers.
add_library(OpenMSTestSupport OBJECT source/OpenMSTestSupport.cpp)
target_link_libraries(OpenMSTestSupport PUBLIC OpenMSTestFramework OpenMS)
openms_add_library_compiler_flags(OpenMSTestSupport)

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
        EXCLUDE_FROM_ALL
        FILE_SET HEADERS DESTINATION ${INSTALL_INCLUDE_DIR} COMPONENT OpenMSTestFramework_headers
        EXCLUDE_FROM_ALL)

openms_doc_path("${PROJECT_SOURCE_DIR}/include")
