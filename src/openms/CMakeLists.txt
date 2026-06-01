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
include(${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/percolator/sources.cmake)
list(APPEND OpenMS_sources ${OpenMS_percolator_sources})

# Suppress upstream-originated warnings on vendored files only.
# GCC/Clang use -Wno-*; MSVC needs /wd<code> with the corresponding warning
# numbers and chokes on the -Wno-* spellings (D8021).
if (MSVC)
  # /wd4100 + /wd4189 explicitly defeat the project-wide /we4100 + /we4189
  # ("unused parameter" and "unused-but-initialized local variable" elevated
  # to errors) on the vendored upstream code, which has plenty of both.
  set_source_files_properties(${OpenMS_percolator_sources}
    PROPERTIES
      COMPILE_FLAGS "/wd4018 /wd4100 /wd4189 /wd4244 /wd4267 /wd4996 /wd4458")
else()
  set_source_files_properties(${OpenMS_percolator_sources}
    PROPERTIES
      COMPILE_FLAGS "-Wno-sign-compare -Wno-unused-parameter -Wno-deprecated-declarations -Wno-unused-variable -Wno-overloaded-virtual -Wno-reorder -Wno-reorder-ctor -Wno-uninitialized -Wno-ignored-qualifiers -Wno-unused-but-set-variable")
  set_source_files_properties(${OpenMS_percolator_sources}
    PROPERTIES COMPILE_OPTIONS "-fvisibility=hidden")
endif()

#------------------------------------------------------------------------------
# all the dependency libraries are linked into libOpenMS.so
set(OPENMS_DEP_LIBRARIES
  LibSVM::LibSVM
  XercesC::XercesC
  CURL::libcurl
)

## setup the arguments to 'target_link_libraries(OpenMS PRIVATE ${OPENMS_DEP_PRIVATE_LIBRARIES})'
## Note: Eigen3::Eigen is PRIVATE because Eigen types are not exposed in the public API
## Note: Evergreen is PRIVATE because it's only used internally by BayesianProteinInferenceAlgorithm
set(OPENMS_DEP_PRIVATE_LIBRARIES
  $<$<BOOL:${WITH_HDF5}>:HDF5::HDF5>
  ${LPTARGET}
  BZip2::BZip2
  Boost::boost
  Boost::date_time
  Boost::iostreams
  Boost::regex
  Eigen3::Eigen
  Evergreen
  GTE
  Quadtree
  SIMDe
  SQLiteCpp
  ZLIB::ZLIB
  libzip::zip
  nlohmann_json::nlohmann_json
  IsoSpec
  eol-bspline
  )

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

# Xerces requires linking against CoreFoundation&CoreServices on macOS
# TODO check if this is still the case
if(APPLE)
  find_library(CoreFoundation_LIBRARY CoreFoundation )
  find_library(CoreServices_LIBRARY CoreServices )
  set(OPENMS_DEP_LIBRARIES ${OPENMS_DEP_LIBRARIES}
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

# Vendored percolator headers are resolved via #include "X.h" (no subdir) —
# add the dir PRIVATE so it's visible only to libOpenMS's own TUs.
target_include_directories(OpenMS PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/percolator)

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
install_directory(${OPENMS_HOST_DIRECTORY}/share/OpenMS/examples ${INSTALL_SHARE_DIR} examples)

# Devs also need our custom find modules
install_directory(${OPENMS_HOST_DIRECTORY}/cmake/Modules/ ${INSTALL_CMAKE_DIR}/Modules cmake)
file(COPY ${OPENMS_HOST_DIRECTORY}/cmake/Modules DESTINATION ${OPENMS_HOST_BINARY_DIRECTORY})

#------------------------------------------------------------------------------
# register relevant paths for the doxygen doc generation
openms_doc_path("${PROJECT_SOURCE_DIR}/include")
