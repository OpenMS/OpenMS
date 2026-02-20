# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# Monolithic OpenMS build — all sources compiled into a single shared library.
# Used on Windows (MSVC) where DLLs require all symbols resolved at link time
# with proper __declspec(dllexport/dllimport) annotations.
# --------------------------------------------------------------------------

# Source files are aggregated by core/io/algo_includes.cmake, which are
# included from the parent CMakeLists.txt (before this file).

# Combine all source and header files
set(OpenMS_all_sources
  ${OpenMS_Core_sources}
  ${OpenMS_IO_sources}
  ${OpenMS_Algo_sources}
)

set(OpenMS_all_sources_h
  ${OpenMS_Core_sources_h}
  ${OpenMS_IO_sources_h}
  ${OpenMS_Algo_sources_h}
  ${OpenMS_configured_headers}
)

source_group("Header Files\\OpenMS" FILES ${OpenMS_all_sources_h})
list(APPEND OpenMS_all_sources ${OpenMS_all_sources_h})
list(REMOVE_DUPLICATES OpenMS_all_sources)

# --- Public dependencies ---
set(OPENMS_DEP_LIBRARIES
  XercesC::XercesC
  Qt6::Core
  Qt6::Network
)

if(NOT DISABLE_OPENSWATH)
  list(APPEND OPENMS_DEP_LIBRARIES OpenSwathAlgo)
endif()

if (OPENMP_FOUND)
  list(APPEND OPENMS_DEP_LIBRARIES OpenMP::OpenMP_CXX)
endif()

if(APPLE)
  find_library(CoreFoundation_LIBRARY CoreFoundation)
  find_library(CoreServices_LIBRARY CoreServices)
  list(APPEND OPENMS_DEP_LIBRARIES ${CoreFoundation_LIBRARY} ${CoreServices_LIBRARY})
endif()

# --- Private dependencies ---
set(OPENMS_DEP_PRIVATE_LIBRARIES
  ${LPTARGET}
  LibSVM::LibSVM
  Boost::boost
  Boost::date_time
  Boost::iostreams
  Boost::regex
  BZip2::BZip2
  Eigen3::Eigen
  eol-bspline
  Evergreen
  GTE
  IsoSpec
  nlohmann_json::nlohmann_json
  Quadtree
  SIMDe
  SQLiteCpp
  ZLIB::ZLIB
)

if (WITH_HDF5)
  list(APPEND OPENMS_DEP_PRIVATE_LIBRARIES HDF5::HDF5)
endif()

if (WITH_PARQUET)
  list(APPEND OPENMS_DEP_PRIVATE_LIBRARIES
       ${OPENMS_ARROW_TARGET}
       ${OPENMS_ARROW_COMPUTE_TARGET}
       ${OPENMS_PARQUET_TARGET})
  if (OPENMS_ARROW_DATASET_TARGET)
    list(APPEND OPENMS_DEP_PRIVATE_LIBRARIES ${OPENMS_ARROW_DATASET_TARGET} LibXml2::LibXml2)
  endif()
endif()

if (ENABLE_TDL)
  list(APPEND OPENMS_DEP_PRIVATE_LIBRARIES tdl::tdl)
endif()

if (WITH_CRAWDAD)
  list(APPEND OPENMS_DEP_PRIVATE_LIBRARIES ${Crawdad_LIBRARY})
endif()

openms_add_library(TARGET_NAME  OpenMS
                   SOURCE_FILES  ${OpenMS_all_sources}
                   HEADER_FILES  ${OpenMS_all_sources_h}
                   INTERNAL_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/include
                                     ${CMAKE_CURRENT_BINARY_DIR}/include
                   LINK_LIBRARIES ${OPENMS_DEP_LIBRARIES}
                   PRIVATE_LINK_LIBRARIES ${OPENMS_DEP_PRIVATE_LIBRARIES}
                   DLL_EXPORT_PATH "OpenMS/")

if (ENABLE_TDL)
  target_compile_definitions(OpenMS PUBLIC ENABLE_TDL)
endif()

if (WITH_PARQUET)
  target_compile_definitions(OpenMS PUBLIC WITH_PARQUET=1)
endif()

if (OPENMS_ARROW_DATASET_TARGET)
  target_compile_definitions(OpenMS PRIVATE WITH_ARROW_DATASET)
endif()

if (MSVC)
  target_compile_options(OpenMS PRIVATE "/we4100" "/we4189")
endif()
