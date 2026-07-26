# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche, Chris Bielow $
# $Authors: Andreas Bertsch, Chris Bielow, Stephan Aiche $
# --------------------------------------------------------------------------

project("OpenMS")
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

#------------------------------------------------------------------------------
# naming conventions:
#
# prefix a variable with 'CF_' if it is used to configure a file!
# e.g., CF_LibOpenMSExport

# Version should not have unescaped "_" because it is used in the Latex tutorials
set(CF_OPENMS_PACKAGE_VERSION "${OPENMS_PACKAGE_VERSION}" CACHE INTERNAL "OpenMS VERSION" FORCE)
set(CF_OPENMS_PACKAGE_VERSION_FULLSTRING "${OPENMS_PACKAGE_VERSION_FULLSTRING}" CACHE INTERNAL "OpenMS version with additional identifier" FORCE)
set(CF_OPENMS_GIT_BRANCH "${OPENMS_GIT_SHORT_REFSPEC}" CACHE INTERNAL "OpenMS Git branch for pyopenms setup.py" FORCE)

#------------------------------------------------------------------------------
# En/disable assertions
if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    set(CF_OPENMS_ASSERTIONS 1)
else()
    set(CF_OPENMS_ASSERTIONS 0)
endif()
set(CF_OPENMS_ASSERTIONS ${CF_OPENMS_ASSERTIONS} CACHE INTERNAL "Enables debug messages (precondition and postconditions are enabled, a bit slower) - this is NOT changing any compiler flags!" FORCE)


#------------------------------------------------------------------------------
# external third-party libs statically linked to OpenMS
#------------------------------------------------------------------------------
# The following library is in our public API, so create export targets
#install(EXPORT EvergreenTarget DESTINATION ${INSTALL_LIB_DIR}/cmake/OpenMS FILE OpenMSTargets.cmake COMPONENT cmake)
add_subdirectory(extern)

#------------------------------------------------------------------------------
# external libs (contrib or system)
#------------------------------------------------------------------------------
#include(${PROJECT_SOURCE_DIR}/cmake_findExternalLibs.cmake)

#------------------------------------------------------------------------------
# configure config.h
#------------------------------------------------------------------------------
include(${PROJECT_SOURCE_DIR}/configh.cmake)

#------------------------------------------------------------------------------
# big include file for headers and cpp files, that fills the OpenMS_sources variable
include (${PROJECT_SOURCE_DIR}/includes.cmake)

#------------------------------------------------------------------------------
# Vendored, namespace-isolated Percolator subtree.
# Built as a standalone static library (OpenMS_Percolator) with its OWN compile
# flags and linked PRIVATE into libOpenMS below. Keeping it a separate target
# isolates the vendored upstream code from the OpenMS target's treat-warning-as-
# error flags (/we4100, /we4189, -Werror) — folding the sources into libOpenMS
# made a clean build depend on the (generator-specific) flag ordering — and
# keeps the upstream sources out of the OpenMS target's file list in IDEs.
add_subdirectory(thirdparty/percolator)

#------------------------------------------------------------------------------
# all the dependency libraries are linked into libOpenMS.so
## Note: XercesC::XercesC is PRIVATE — Xerces is fully encapsulated behind the
## Xerces-free XMLHandler / XMLAttributes / StringManager public API and the
## internal SAX2HandlerAdapter, so it no longer appears in any installed header
## and downstream consumers neither include nor link it (see #8547).
## Note: Boost::boost (header-only) is PUBLIC because the installed public header
## config.h includes <boost/current_function.hpp> (via Macros.h), so consumers
## (e.g. pyOpenMS) need Boost's include dir. Only the header-only target is
## public — it propagates include dirs but no link libraries; the compiled Boost
## components (date_time/iostreams/regex) stay PRIVATE below.
set(OPENMS_DEP_LIBRARIES
  CURL::libcurl
  Boost::boost
)

## setup the arguments to 'target_link_libraries(OpenMS PRIVATE ${OPENMS_DEP_PRIVATE_LIBRARIES})'
## Note: Eigen3::Eigen is PRIVATE because Eigen types are not exposed in the public API
## Note: Evergreen is PRIVATE because it's only used internally by BayesianProteinInferenceAlgorithm
## Note: LibSVM::LibSVM is PRIVATE because LIBSVM types are not exposed in the public API;
##       they are fully encapsulated by the SimpleSVM wrapper (Pimpl) (see #8547).
set(OPENMS_DEP_PRIVATE_LIBRARIES
  XercesC::XercesC
  $<$<BOOL:${WITH_HDF5}>:HDF5::HDF5>
  ${LPTARGET}
  BZip2::BZip2
  Boost::date_time
  Boost::iostreams
  Boost::regex
  Eigen3::Eigen
  Evergreen
  GTE
  LibSVM::LibSVM
  Quadtree
  SIMDe
  SQLiteCpp
  ZLIB::ZLIB
  libzip::zip
  nlohmann_json::nlohmann_json
  IsoSpec
  eol-bspline
  OpenMS_Percolator
  )

option(WITH_ONNX "Build with ONNX Runtime support for PeptDeep ML bindings" OFF)

if(WITH_ONNX)
  if(OPENMS_USE_VCPKG)
    list(APPEND OPENMS_DEP_PRIVATE_LIBRARIES onnxruntime::onnxruntime)
    get_target_property(ORT_LIB_PATH onnxruntime::onnxruntime IMPORTED_LOCATION)
  else() 
    list(APPEND OPENMS_DEP_PRIVATE_LIBRARIES ONNXRuntime::ONNXRuntime)
    # Bypass IMPORTED_RUNTIME_ARTIFACTS by grabbing the file path directly.
    # This avoids the read-only TYPE property error completely.
    get_target_property(ORT_LIB_PATH ONNXRuntime::ONNXRuntime IMPORTED_LOCATION)
  endif()

  if(ORT_LIB_PATH)
      # Safely install ONNX shared libraries based on the operating system
      if(WIN32)
          install(PROGRAMS "${ORT_LIB_PATH}" DESTINATION "${INSTALL_BIN_DIR}")
      else()
          install(PROGRAMS "${ORT_LIB_PATH}" DESTINATION "${INSTALL_LIB_DIR}")
      endif()
  endif()

  # Download PeptDeep models once to the OpenMS model directory.
  # These ONNX files are embedded-weight models and do not require .onnx.data sidecars.
  set(OPENMS_PEPTDEEP_MODEL_URL "http://archive.openms.de/openms/models" CACHE STRING
      "Base URL for PeptDeep ONNX model downloads")
  set(OPENMS_PEPTDEEP_MODEL_DIR "${PROJECT_BINARY_DIR}/share/OpenMS/models" CACHE INTERNAL
      "Build-tree directory containing downloaded PeptDeep ONNX models")

  set(OPENMS_PEPTDEEP_MODELS
      peptdeep_ccs_dynamic.onnx
      peptdeep_ms2_dynamic.onnx
      peptdeep_rt_dynamic.onnx
      CACHE INTERNAL "PeptDeep ONNX models required by OpenMS"
  )

  set(OPENMS_PEPTDEEP_MODEL_SHA256S
      41816f325cdd838897b7e7dcaf291da5d4d3903a1532fc5dc9db7e0725015801
      4141ea233d663b6f5b98ad80b8923ee685dae5cfdbe12243e2e30cd9c44ecf00
      3f3b847ea37a9333ac163409f3d4b343d021471548091c61aed6cfa398ce62fd
      CACHE INTERNAL "SHA256 checksums for PeptDeep ONNX models"
  )

  file(MAKE_DIRECTORY "${OPENMS_PEPTDEEP_MODEL_DIR}")

  list(LENGTH OPENMS_PEPTDEEP_MODELS OPENMS_PEPTDEEP_MODEL_COUNT)
  list(LENGTH OPENMS_PEPTDEEP_MODEL_SHA256S OPENMS_PEPTDEEP_HASH_COUNT)

  if(NOT OPENMS_PEPTDEEP_MODEL_COUNT EQUAL OPENMS_PEPTDEEP_HASH_COUNT)
    message(FATAL_ERROR
      "OPENMS_PEPTDEEP_MODELS and OPENMS_PEPTDEEP_MODEL_SHA256S must have the same length")
  endif()

  math(EXPR OPENMS_PEPTDEEP_LAST_MODEL_INDEX "${OPENMS_PEPTDEEP_MODEL_COUNT} - 1")

  set(OPENMS_PEPTDEEP_MODEL_FILES)
  foreach(MODEL_INDEX RANGE 0 ${OPENMS_PEPTDEEP_LAST_MODEL_INDEX})
    list(GET OPENMS_PEPTDEEP_MODELS ${MODEL_INDEX} MODEL)
    list(GET OPENMS_PEPTDEEP_MODEL_SHA256S ${MODEL_INDEX} EXPECTED_SHA256)

    set(MODEL_PATH "${OPENMS_PEPTDEEP_MODEL_DIR}/${MODEL}")
    list(APPEND OPENMS_PEPTDEEP_MODEL_FILES "${MODEL_PATH}")

    set(MODEL_IS_VALID FALSE)

    if(EXISTS "${MODEL_PATH}")
      file(SHA256 "${MODEL_PATH}" ACTUAL_SHA256)

      if(ACTUAL_SHA256 STREQUAL EXPECTED_SHA256)
        set(MODEL_IS_VALID TRUE)
      else()
        message(WARNING
          "Existing PeptDeep model ${MODEL} has SHA256 checksum ${ACTUAL_SHA256}, "
          "expected ${EXPECTED_SHA256}; re-downloading.")
        file(REMOVE "${MODEL_PATH}")
      endif()
    endif()

    if(NOT MODEL_IS_VALID)
      message(STATUS "Fetching PeptDeep model: ${MODEL} from archive...")
      file(
        DOWNLOAD "${OPENMS_PEPTDEEP_MODEL_URL}/${MODEL}" "${MODEL_PATH}"
        STATUS DOWNLOAD_STATUS
        SHOW_PROGRESS
        EXPECTED_HASH "SHA256=${EXPECTED_SHA256}"
      )

      list(GET DOWNLOAD_STATUS 0 DOWNLOAD_CODE)
      list(GET DOWNLOAD_STATUS 1 DOWNLOAD_MESSAGE)

      if(NOT DOWNLOAD_CODE EQUAL 0)
        file(REMOVE "${MODEL_PATH}")
        message(FATAL_ERROR
          "Failed to download PeptDeep model ${MODEL}: ${DOWNLOAD_MESSAGE}")
      endif()
    endif()
  endforeach()

  # Ensure models are packaged during 'make install'.
  install(FILES ${OPENMS_PEPTDEEP_MODEL_FILES}
          DESTINATION "share/OpenMS/models")

endif()

# Parquet support requires Arrow core + compute (filter expressions) + Parquet.
list(APPEND OPENMS_DEP_PRIVATE_LIBRARIES
     ${OPENMS_ARROW_TARGET}
     ${OPENMS_ARROW_COMPUTE_TARGET}
     ${OPENMS_PARQUET_TARGET})
if (OPENMS_ARROW_DATASET_TARGET)
  list(APPEND OPENMS_DEP_PRIVATE_LIBRARIES ${OPENMS_ARROW_DATASET_TARGET})
  # Dataset pushdown needs libxml2 when Arrow Dataset is linked statically.
  list(APPEND OPENMS_DEP_PRIVATE_LIBRARIES LibXml2::LibXml2)
endif()

if (WITH_OPENTIMS)
  list(APPEND OPENMS_DEP_PRIVATE_LIBRARIES opentims::opentims_cpp)
endif()

if (WITH_THERMO_RAW)
  list(APPEND OPENMS_DEP_PRIVATE_LIBRARIES OpenMSThermoBridge::openms_thermo_bridge)
endif()

if (ENABLE_TDL)
set(OPENMS_DEP_PRIVATE_LIBRARIES
    ${OPENMS_DEP_PRIVATE_LIBRARIES}
    tdl::tdl
    )
endif()

# Xerces requires linking against CoreFoundation&CoreServices on macOS.
# Xerces is now a PRIVATE dependency, so these frameworks are private too.
# TODO check if this is still the case
if(APPLE)
  find_library(CoreFoundation_LIBRARY CoreFoundation )
  find_library(CoreServices_LIBRARY CoreServices )
  list(APPEND OPENMS_DEP_PRIVATE_LIBRARIES
                           ${CoreFoundation_LIBRARY}
                           ${CoreServices_LIBRARY})
endif()

if (OPENMP_FOUND)
  list(APPEND OPENMS_DEP_LIBRARIES OpenMP::OpenMP_CXX)
endif()

if (WITH_CRAWDAD) ## TODO check if still necessary
  list(APPEND OPENMS_DEP_LIBRARIES ${Crawdad_LIBRARY})
endif()

openms_add_library(TARGET_NAME  OpenMS
                   SOURCE_FILES  ${OpenMS_sources}
                   HEADER_FILES  ${OpenMS_sources_h}
                                 ${OpenMS_configured_headers}
                   INTERNAL_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/include
                                     ${CMAKE_CURRENT_BINARY_DIR}/include
                   PRIVATE_INCLUDES
                   EXTERNAL_INCLUDES
                   LINK_LIBRARIES OpenSwathAlgo
                                  ${OPENMS_DEP_LIBRARIES}
                   PRIVATE_LINK_LIBRARIES ${OPENMS_DEP_PRIVATE_LIBRARIES}
                   DLL_EXPORT_PATH "OpenMS/")

# The vendored percolator headers (#include "X.h", no subdir) are used only by
# libOpenMS's Percolator pimpl TU. Their directory is made visible via the
# PUBLIC include dir propagated from the linked OpenMS_Percolator target (see
# thirdparty/percolator/CMakeLists.txt), so no explicit include dir is needed.

if (OPENMS_ARROW_DATASET_TARGET)
  target_compile_definitions(OpenMS PRIVATE WITH_ARROW_DATASET)
endif()

# When Arrow is statically linked, hide its bundled third-party symbols from
# libOpenMS.so's dynamic symbol table to avoid conflicts with system libraries.
# NOTE: We intentionally do NOT hide libarrow.a or libparquet.a symbols because
# the OpenMS public API exposes arrow::Table etc., and downstream consumers
# (including test executables) need to resolve those symbols from libOpenMS.so.
if(ARROW_USE_STATIC AND NOT MSVC AND NOT APPLE)
  target_link_options(OpenMS PRIVATE "LINKER:--exclude-libs,libarrow_bundled_dependencies.a")
endif()

if (MSVC)
  ## treat warning of unused function parameter as error, similar to -Werror=unused-variable on GCC
  target_compile_options(OpenMS PRIVATE "/we4100")
  ## treat warning of unused local variable as error, similar to -Werror=unused-variable on GCC
  target_compile_options(OpenMS PRIVATE "/we4189")
endif()

if (ENABLE_TDL)
    target_compile_definitions(OpenMS PUBLIC ENABLE_TDL)
endif()

# Use std::filesystem instead of boost::filesystem in boost::process (avoids linking boost_filesystem)
# Requires Boost >= 1.81 (BOOST_PROCESS_USE_STD_FS was introduced in 1.78)
target_compile_definitions(OpenMS PRIVATE BOOST_PROCESS_USE_STD_FS)

if (WITH_OPENTIMS)
    target_compile_definitions(OpenMS PUBLIC WITH_OPENTIMS=1)
endif()

if (WITH_THERMO_RAW)
    target_compile_definitions(OpenMS PUBLIC WITH_THERMO_RAW=1)
    # Copy the thermo bridge managed runtime files (.dll + .runtimeconfig.json)
    # to the OpenMS library output directory so they can be found at runtime.
    if (COMMAND openms_thermo_bridge_copy_runtime_files)
      openms_thermo_bridge_copy_runtime_files(TARGET OpenMS)
    endif()
endif()

if(WITH_WNETALIGN)
  target_include_directories(OpenMS SYSTEM PRIVATE ${WNETALIGN_INCLUDE_DIRS})
  target_compile_definitions(OpenMS PRIVATE LEMON_ONLY_TEMPLATES OPENMS_HAS_WNETALIGN)
endif()

## Enable test-only hooks when class tests are being built.
if(ENABLE_CLASS_TESTING)
  target_compile_definitions(OpenMS PRIVATE OPENMS_ENABLE_TESTING_HOOKS=1)
endif()

#------------------------------------------------------------------------------
# since the share basically belongs to OpenMS core we control its installation
# here
# Note: that slash(/) is important here, otherwise the whole directory
#       (not its content) will be copied!
#install_directory(${OPENMS_HOST_DIRECTORY}/share/OpenMS/ ${INSTALL_SHARE_DIR} share)
install(DIRECTORY ${OPENMS_HOST_DIRECTORY}/share/OpenMS/
  DESTINATION ${INSTALL_SHARE_DIR}
  COMPONENT share
  FILE_PERMISSIONS      OWNER_WRITE OWNER_READ
                        GROUP_READ
                        WORLD_READ
  DIRECTORY_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                        GROUP_EXECUTE GROUP_READ
                        WORLD_EXECUTE WORLD_READ
  REGEX "^\\..*" EXCLUDE ## Exclude hidden files (svn, git, DSStore)
  REGEX ".*\\/\\..*" EXCLUDE ## Exclude hidden files in subdirectories
  REGEX "OpenMS/examples" EXCLUDE
)
install(FILES ${OPENMS_HOST_DIRECTORY}/License.txt
  DESTINATION ${INSTALL_SHARE_DIR}/LICENSES
  RENAME OpenMS-BSD-3-Clause.txt
  COMPONENT share
)
if(INSTALL_OPENMS_EXAMPLES)
  install_directory(${OPENMS_HOST_DIRECTORY}/share/OpenMS/examples ${INSTALL_SHARE_DIR} examples)
endif()

# Devs also need our custom find modules
install_directory(${OPENMS_HOST_DIRECTORY}/cmake/Modules/ ${INSTALL_CMAKE_DIR}/Modules cmake)
file(COPY ${OPENMS_HOST_DIRECTORY}/cmake/Modules DESTINATION ${OPENMS_HOST_BINARY_DIRECTORY})

openms_doc_path("${PROJECT_SOURCE_DIR}/include")

# Add the ONNX source files conditionally
if(WITH_ONNX)
  target_sources(OpenMS PRIVATE
    source/ML/ONNX/ONNXEnvironment.cpp
    source/ML/ONNX/ONNXPredictorBase.cpp
    source/ML/PEPTDEEP/PeptDeepCCSInference.cpp
    source/ML/PEPTDEEP/PeptDeepInput.cpp
    source/ML/PEPTDEEP/PeptDeepRTInference.cpp
    source/ML/PEPTDEEP/PeptDeepMS2Inference.cpp
  )
endif()
