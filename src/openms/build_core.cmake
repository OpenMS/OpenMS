# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# OpenMS_Core — fundamental data structures, kernel, chemistry, metadata
# --------------------------------------------------------------------------

# Aggregate Core source files using the reset-and-capture pattern
include(${CMAKE_CURRENT_LIST_DIR}/core_includes.cmake)

# Add configured headers and merge header list
source_group("Header Files\\OpenMS" FILES ${OpenMS_configured_headers})
list(APPEND OpenMS_Core_sources ${OpenMS_Core_sources_h} ${OpenMS_configured_headers})
list(REMOVE_DUPLICATES OpenMS_Core_sources)

# --- Public dependencies (propagated to consumers of OpenMS_Core) ---
set(CORE_DEP_LIBRARIES
  XercesC::XercesC
  Qt6::Core
  Qt6::Network
)

if (OPENMP_FOUND)
  list(APPEND CORE_DEP_LIBRARIES OpenMP::OpenMP_CXX)
endif()

if(APPLE)
  find_library(CoreFoundation_LIBRARY CoreFoundation)
  find_library(CoreServices_LIBRARY CoreServices)
  list(APPEND CORE_DEP_LIBRARIES ${CoreFoundation_LIBRARY} ${CoreServices_LIBRARY})
endif()

# --- Private dependencies (only needed to compile Core itself) ---
set(CORE_DEP_PRIVATE_LIBRARIES
  ${LPTARGET}
  BZip2::BZip2
  Boost::boost
  Boost::date_time
  Boost::regex
  Eigen3::Eigen
  eol-bspline
  Evergreen
  GTE
  IsoSpec
  Quadtree
  SIMDe
  ZLIB::ZLIB
)

openms_add_library(TARGET_NAME  OpenMS_Core
                   SOURCE_FILES  ${OpenMS_Core_sources}
                   HEADER_FILES  ${OpenMS_Core_sources_h}
                                 ${OpenMS_configured_headers}
                   INTERNAL_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/include
                                     ${CMAKE_CURRENT_BINARY_DIR}/include
                   LINK_LIBRARIES ${CORE_DEP_LIBRARIES}
                   PRIVATE_LINK_LIBRARIES ${CORE_DEP_PRIVATE_LIBRARIES}
                   DLL_EXPORT_PATH "")

# Phase 1: Use default visibility to avoid cross-library hidden symbol issues.
# The monolithic library had hidden visibility too, but all code was in one DSO
# so internal symbols were accessible. With the split, symbols that cross library
# boundaries (like template specializations) need to be visible.
# TODO Phase 4: Add OPENMS_DLLAPI to cross-boundary symbols and re-enable hidden.
set_target_properties(OpenMS_Core PROPERTIES CXX_VISIBILITY_PRESET default)
set_target_properties(OpenMS_Core PROPERTIES VISIBILITY_INLINES_HIDDEN 0)

if (MSVC)
  target_compile_options(OpenMS_Core PRIVATE "/we4100" "/we4189")
endif()
