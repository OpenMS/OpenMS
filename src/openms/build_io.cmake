# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# OpenMS_IO — file format readers/writers (FORMAT module)
# --------------------------------------------------------------------------

# Aggregate IO source files using the reset-and-capture pattern
include(${CMAKE_CURRENT_LIST_DIR}/io_includes.cmake)

# Merge header list (no configured headers for IO — they belong to Core)
source_group("Header Files\\OpenMS\\FORMAT" FILES ${OpenMS_IO_sources_h})
list(APPEND OpenMS_IO_sources ${OpenMS_IO_sources_h})
list(REMOVE_DUPLICATES OpenMS_IO_sources)

# --- Public dependencies ---
# OpenMS_Core is PUBLIC so its include dirs and transitive deps propagate
set(IO_DEP_LIBRARIES
  OpenMS_Core
  XercesC::XercesC
)

if(NOT DISABLE_OPENSWATH)
  list(APPEND IO_DEP_LIBRARIES OpenSwathAlgo)
endif()

# --- Private dependencies ---
set(IO_DEP_PRIVATE_LIBRARIES
  Boost::boost
  Boost::date_time
  Boost::iostreams
  Boost::regex
  BZip2::BZip2
  Eigen3::Eigen
  SIMDe
  SQLiteCpp
  nlohmann_json::nlohmann_json
  ZLIB::ZLIB
)
# Note: BZip2 and ZLIB are here (not Core) because only FORMAT files use them.

if (WITH_HDF5)
  list(APPEND IO_DEP_PRIVATE_LIBRARIES HDF5::HDF5)
endif()

if (WITH_PARQUET)
  list(APPEND IO_DEP_PRIVATE_LIBRARIES
       ${OPENMS_ARROW_TARGET}
       ${OPENMS_ARROW_COMPUTE_TARGET}
       ${OPENMS_PARQUET_TARGET})
  if (OPENMS_ARROW_DATASET_TARGET)
    list(APPEND IO_DEP_PRIVATE_LIBRARIES ${OPENMS_ARROW_DATASET_TARGET} LibXml2::LibXml2)
  endif()
endif()

if (ENABLE_TDL)
  list(APPEND IO_DEP_PRIVATE_LIBRARIES tdl::tdl)
endif()

openms_add_library(TARGET_NAME  OpenMS_IO
                   SOURCE_FILES  ${OpenMS_IO_sources}
                   HEADER_FILES  ${OpenMS_IO_sources_h}
                   INTERNAL_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/include
                                     ${CMAKE_CURRENT_BINARY_DIR}/include
                   LINK_LIBRARIES ${IO_DEP_LIBRARIES}
                   PRIVATE_LINK_LIBRARIES ${IO_DEP_PRIVATE_LIBRARIES}
                   DLL_EXPORT_PATH "")

if (ENABLE_TDL)
  target_compile_definitions(OpenMS_IO PUBLIC ENABLE_TDL)
endif()

if (WITH_PARQUET)
  target_compile_definitions(OpenMS_IO PUBLIC WITH_PARQUET=1)
endif()

if (OPENMS_ARROW_DATASET_TARGET)
  target_compile_definitions(OpenMS_IO PRIVATE WITH_ARROW_DATASET)
endif()

# Phase 1: Use default visibility (see build_core.cmake for rationale)
set_target_properties(OpenMS_IO PROPERTIES CXX_VISIBILITY_PRESET default)
set_target_properties(OpenMS_IO PROPERTIES VISIBILITY_INLINES_HIDDEN 0)

# Phase 1: Allow deferred symbol resolution on macOS (see build_core.cmake)
if (APPLE)
  target_link_options(OpenMS_IO PRIVATE "LINKER:-undefined,dynamic_lookup")
endif()

if (MSVC)
  target_compile_options(OpenMS_IO PRIVATE "/we4100" "/we4189")
endif()
