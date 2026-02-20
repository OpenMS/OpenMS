# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# OpenMS_Algo — algorithms, analysis, processing, ML, feature finding
# --------------------------------------------------------------------------

# Aggregate Algo source files using the reset-and-capture pattern
include(${CMAKE_CURRENT_LIST_DIR}/algo_includes.cmake)

# Merge header list
source_group("Header Files\\OpenMS" FILES ${OpenMS_Algo_sources_h})
list(APPEND OpenMS_Algo_sources ${OpenMS_Algo_sources_h})
list(REMOVE_DUPLICATES OpenMS_Algo_sources)

# --- Public dependencies ---
# OpenMS_IO is PUBLIC — transitively brings OpenMS_Core and all its deps
set(ALGO_DEP_LIBRARIES
  OpenMS_IO
)

if(NOT DISABLE_OPENSWATH)
  list(APPEND ALGO_DEP_LIBRARIES OpenSwathAlgo)
endif()

# --- Private dependencies ---
set(ALGO_DEP_PRIVATE_LIBRARIES
  ${LPTARGET}
  LibSVM::LibSVM
  Boost::boost
  Boost::regex
  Eigen3::Eigen
  Evergreen
  GTE
  Quadtree
  SQLiteCpp
)

if (WITH_CRAWDAD)
  list(APPEND ALGO_DEP_PRIVATE_LIBRARIES ${Crawdad_LIBRARY})
endif()

openms_add_library(TARGET_NAME  OpenMS_Algo
                   SOURCE_FILES  ${OpenMS_Algo_sources}
                   HEADER_FILES  ${OpenMS_Algo_sources_h}
                   INTERNAL_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/include
                                     ${CMAKE_CURRENT_BINARY_DIR}/include
                   LINK_LIBRARIES ${ALGO_DEP_LIBRARIES}
                   PRIVATE_LINK_LIBRARIES ${ALGO_DEP_PRIVATE_LIBRARIES}
                   DLL_EXPORT_PATH "")

# Phase 1: Use default visibility (see build_core.cmake for rationale)
set_target_properties(OpenMS_Algo PROPERTIES CXX_VISIBILITY_PRESET default)
set_target_properties(OpenMS_Algo PROPERTIES VISIBILITY_INLINES_HIDDEN 0)

# Phase 1: Allow deferred symbol resolution on macOS (see build_core.cmake)
if (APPLE)
  target_link_options(OpenMS_Algo PRIVATE "LINKER:-undefined,dynamic_lookup")
endif()

if (MSVC)
  target_compile_options(OpenMS_Algo PRIVATE "/we4100" "/we4189")
endif()
