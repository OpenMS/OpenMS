# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Julianus Pfeuffer $
# $Authors: Stephan Aiche, Julianus Pfeuffer $
# --------------------------------------------------------------------------



# --------------------------------------------------------------------------
# general definitions used for building OpenMS packages
set(CPACK_PACKAGE_NAME "OpenMS")
set(CPACK_PACKAGE_VENDOR "OpenMS.de")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "OpenMS - A framework for mass spectrometry")
set(CPACK_PACKAGE_VERSION "${OPENMS_PACKAGE_VERSION}")
set(CPACK_PACKAGE_VERSION_MAJOR "${OPENMS_PACKAGE_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${OPENMS_PACKAGE_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${OPENMS_PACKAGE_VERSION_PATCH}")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "${CPACK_PACKAGE_NAME}-${OPENMS_PACKAGE_VERSION}")
set(CPACK_PACKAGE_DESCRIPTION_FILE ${PROJECT_SOURCE_DIR}/cmake/OpenMSPackageDescriptionFile.txt)
set(CPACK_RESOURCE_FILE_LICENSE ${PROJECT_SOURCE_DIR}/License.txt)
set(CPACK_RESOURCE_FILE_WELCOME ${PROJECT_SOURCE_DIR}/cmake/OpenMSPackageResourceWelcomeFile.txt)
set(CPACK_RESOURCE_FILE_README ${PROJECT_SOURCE_DIR}/cmake/OpenMSPackageResourceReadme.txt)
set(CPACK_STRIP_FILES TRUE) # to save some space in the installers

set(OPENMS_LOGO_NAME openms_logo_large_transparent.png) ## The filename of the logo to be used for the OpenMS folder e.g. on the DMG
set(OPENMS_LOGO ${PROJECT_SOURCE_DIR}/cmake/MacOSX/${OPENMS_LOGO_NAME}) ## The logo to be used for the OpenMS folder e.g. on the DMG

set(OPENMS_LOGOSMALL_NAME openms_logo_corner_small.png) ## The filename of the logo to be used for the OpenMS folder e.g. on the PKG
set(OPENMS_LOGOSMALL ${PROJECT_SOURCE_DIR}/cmake/MacOSX/${OPENMS_LOGOSMALL_NAME}) ## The logo to be used for the OpenMS folder e.g. on the PKG

########################################################### Fixing dynamic dependencies
## Qt Plugins needed for the CL tools should have been installed before (such as QSqliteDriverPlugin)
## This currently works because our libs and TOPP tools include all dependencies. For macOS,
##  the app bundles need to have a different RUNTIME_DEPENDENCY_SET (TOPPView_DEPS, ...) due
##  to CMake assuming you want standalone bundles. But we want to share libs between them.

# Add CMAKE_PREFIX_PATH directories to search for runtime dependencies.
# We need this because TARGET_RUNTIME_DLLS (used during build stage copy) has limitations:
# - Incomplete CMake configs from third-party libraries that don't properly export targets
# - Private shared library dependencies not tracked in the dependency graph  
# - Manual find_library() results without imported targets
# By adding these paths, we ensure all dependencies are found during installation even if
# they were missed during the build-time copy step.
list(TRANSFORM CMAKE_PREFIX_PATH APPEND "/bin" OUTPUT_VARIABLE DEP_BIN_DIRS)
list(TRANSFORM CMAKE_PREFIX_PATH APPEND "/lib" OUTPUT_VARIABLE DEP_LIB_DIRS)

# Also include our own runtime directory where dependencies were copied during build.
# This serves as the primary source, with CMAKE_PREFIX_PATH directories as fallback.
list(APPEND DEP_BIN_DIRS $<TARGET_FILE_DIR:OpenMS>)

# Combine all search directories for comprehensive dependency resolution
set(RUNTIME_DEP_SEARCH_DIRS ${DEP_BIN_DIRS} ${DEP_LIB_DIRS})

# Ensure Arrow/Parquet shared libs are discoverable during install() dependency collection.
# This feeds RUNTIME_DEP_SEARCH_DIRS so install(RUNTIME_DEPENDENCY_SET ...) can locate Arrow libs
# that are not always in the main build tree. It is required for packaging (CPack installers and
# pyOpenMS wheels) to bundle Arrow runtime deps like libarrow_compute/libarrow_dataset, which are
# loaded by OpenMS/pyOpenMS at runtime. Without these search paths, the install step can miss those
# shared libs and the resulting binaries/wheels fail to link or import on user machines.
foreach(_arrow_dep IN ITEMS OPENMS_ARROW_TARGET OPENMS_ARROW_COMPUTE_TARGET OPENMS_PARQUET_TARGET OPENMS_ARROW_DATASET_TARGET)
  if(DEFINED ${_arrow_dep} AND NOT "${${_arrow_dep}}" STREQUAL "")
    if(TARGET ${${_arrow_dep}})
      list(APPEND RUNTIME_DEP_SEARCH_DIRS $<TARGET_FILE_DIR:${${_arrow_dep}}>)
    else()
      get_filename_component(_arrow_dep_dir "${${_arrow_dep}}" DIRECTORY)
      if(_arrow_dep_dir)
        list(APPEND RUNTIME_DEP_SEARCH_DIRS "${_arrow_dep_dir}")
      endif()
    endif()
  endif()
endforeach()


## Info on excluding dependencies:
# PRE_EXCLUDE_REGEXES: Excludes dependencies at the beginning of the dependency analysis. 
#                      This means that any dependency matching these patterns will be excluded before 
#                      CMake traverses its own dependencies. This can prevent CMake from spending 
#                      time analyzing entire dependency chains that you know you want to exclude.
# POST_EXCLUDE_REGEXES: Excludes dependencies after the complete dependency analysis is done. 
#                       This is applied to the final list of all found dependencies.

# On Windows we need to tell CMake where to look for.
# We also do not need API sets. So exclude them.


if(WIN32)
  # exclude dll's which are system dll's and should not be shipped (bloats installer and leads to incompatibilities)
  set(PRE_EXCLUDE
                  ## these two are direct systems deps of TOPPView etc. Exclude to save time
                  "api-ms" "ext-ms"
                  ## "HvsiFileTrust" "PdmUtilities" are detected as a dependency by CMake which cannot be resolved (and would lead to errors), so ignore it
                  "hvsi" "pdmutilities"  ## make all lower case, since this is what CMake extracts from the targets and the regex is case sensitive
                  ## MSVC runtime DLLs are handled separately by InstallRequiredSystemLibraries (in package_nsis.cmake).
                  ## Exclude them here to avoid "Multiple conflicting paths" errors when the same DLL
                  ## exists in multiple search directories (e.g. a Conda env and another prefix's bin/). CMake 4.x
                  ## treats such conflicts as fatal errors.
                  "vcruntime" "msvcp" "concrt" "vccorlib" "ucrtbase"
                  )
  ## exclude every Dll from c:\Windows\System32
  ## Note: CMake extracts Dll names and will have a list like
  ##-- Resolved runtime dependencies:
  ##--   C:/WINDOWS/system32/aclui.dll
  ##--   C:/WINDOWS/system32/activeds.dll
  ## ,i.e. folder names have weird cases (and the CMake regex engine is case sensitive)
  set(POST_EXCLUDE ".*[\\/][Ss][Yy][Ss][Tt][Ee][Mm]32[\\/].*")  # skip system32 DLLs completely (no matter how CMake names the path, could be 'System32', 'system32', 'SYSTEM32' etc)
elseif(APPLE)
  set(PRE_EXCLUDE "/usr/lib" "/System/")
  set(POST_EXCLUDE "")
else()
  set(PRE_EXCLUDE "")
  set(POST_EXCLUDE ".*/ld-linux-.*" ".*/linux-vdso.*" ".*/libm\\..*" ".*/libc\\..*" ".*/libpthread\\..*" ".*/libdl\\..*" ".*/libstdc\\+\\+\\..*" ".*/libgcc_s.*" ".*/libgomp\\..*" ".*/libQt6.*")
endif()

## The bundled libraries are writable by their owner only (0755), like any installed
## library. With GROUP_WRITE and WORLD_WRITE they were world-writable: the DEB's in
## /usr/lib, where any local user could replace code that the OpenMS tools load, also when
## root runs them, and the files of the Qt frameworks in the macOS package. The packaging
## steps that modify the staged libraries (code signing on macOS, the RUNPATH clean-up of
## the DEB) run as their owner.
install(RUNTIME_DEPENDENCY_SET OPENMS_DEPS
        DESTINATION ${INSTALL_LIB_DIR}
        PERMISSIONS
          OWNER_READ OWNER_WRITE OWNER_EXECUTE
          GROUP_READ GROUP_EXECUTE
          WORLD_READ WORLD_EXECUTE
        COMPONENT Dependencies
        PRE_EXCLUDE_REGEXES ${PRE_EXCLUDE}
        POST_EXCLUDE_REGEXES ${POST_EXCLUDE}
        DIRECTORIES ${RUNTIME_DEP_SEARCH_DIRS})

#install(RUNTIME_DEPENDENCY_SET TOPPView_DEPS) # I think without giving DESTINATION and COMPONENT it will be inferred
#install(RUNTIME_DEPENDENCY_SET TOPPAS_DEPS)
#...

########################################################### OpenMP runtime (macOS)
## Apple's Clang ships no OpenMP runtime, so the build takes libomp from Homebrew
## (tools/ci/deps-macos.sh; keg-only, found through OpenMP_ROOT). The runtime dependency
## set above bundles libomp.dylib. Its omp.h goes next to the OpenMS headers, some of which
## include it: OpenMSConfig.cmake points FindOpenMP of a consuming project at that header
## and the bundled library, so that a program built against the package or the SDK loads
## the bundled runtime only. Two OpenMP runtimes in one process abort ("OMP: Error #15").
set(_openms_omp_libraries "")
if(APPLE AND OPENMP_FOUND)
  set(_openms_omp_libraries ${OpenMP_CXX_LIBRARIES})
  foreach(_openms_omp_include_dir IN LISTS OpenMP_CXX_INCLUDE_DIRS)
    if(EXISTS "${_openms_omp_include_dir}/omp.h")
      install(FILES "${_openms_omp_include_dir}/omp.h"
              DESTINATION ${INSTALL_INCLUDE_DIR}
              COMPONENT thirdparty_headers)
      break()
    endif()
  endforeach()
endif()

########################################################### Third-party licenses
## The license texts of what only the packages bundle: Qt and, on macOS, the Homebrew
## formulae (cmake/third_party_licenses.cmake). Those of the vcpkg ports are installed for
## every installation (CMakeLists.txt). The DEB depends on the distribution's Qt (libQt6 is
## excluded above); the Windows and macOS packages bundle theirs.
set(_openms_bundled_qt_version "")
## On macOS Qt and libomp come from Homebrew (tools/ci/deps-macos.sh): the formulae of Qt are
## the kegs that hold its CMake packages, that of libomp the keg that holds the library.
set(_openms_homebrew_paths ${_openms_omp_libraries})
if(WITH_GUI AND Qt6Core_FOUND AND (WIN32 OR APPLE))
  set(_openms_bundled_qt_version "${Qt6Core_VERSION}")
  if(APPLE)
    foreach(_openms_qt_component IN ITEMS Core ${OpenMS_GUI_QT_COMPONENTS})
      list(APPEND _openms_homebrew_paths "${Qt6${_openms_qt_component}_DIR}")
    endforeach()
  endif()
endif()
set(_openms_homebrew_formulae "")
set(_openms_homebrew_prefix "")
if(APPLE)
  openms_homebrew_formulae_of(_openms_homebrew_formulae _openms_homebrew_prefix ${_openms_homebrew_paths})
endif()
openms_install_third_party_licenses(QT_VERSION "${_openms_bundled_qt_version}"
                                    HOMEBREW_PREFIX "${_openms_homebrew_prefix}"
                                    HOMEBREW_FORMULAE ${_openms_homebrew_formulae})
## On Linux, the runtime dependency set also bundles the libraries of the build machine that
## vcpkg's libraries need and the exclusions above keep (libgfortran, which vcpkg's LAPACK
## needs). Their licenses come from the Debian packages that installed them, so only the DEB
## does this.
if("${PACKAGE_TYPE}" STREQUAL "deb"
   AND OPENMS_USE_VCPKG AND VCPKG_INSTALLED_DIR AND VCPKG_TARGET_TRIPLET)
  file(GLOB _openms_vcpkg_shared_libraries LIST_DIRECTORIES false
       "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/lib/*.so"
       "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/lib/*.so.*")
  openms_install_system_library_licenses(ROOTS ${_openms_vcpkg_shared_libraries}
                                         EXCLUDE_REGEXES ${POST_EXCLUDE})
endif()

########################################################### SEARCHENGINES
set(THIRDPARTY_COMPONENT_GROUP)
## populates the THIRDPARTY_COMPONENT_GROUP list
if(EXISTS ${SEARCH_ENGINES_DIRECTORY})
  ## Automatically recurse over all subfolders in SEARCH_ENGINES_DIRECTORY
  file(GLOB THIRDPARTY_SUBDIRS RELATIVE ${SEARCH_ENGINES_DIRECTORY} ${SEARCH_ENGINES_DIRECTORY}/*)
  foreach(SUBDIR ${THIRDPARTY_SUBDIRS})
    ## Not shipped, although THIRDPARTY has them:
    ## - ProteoWizard (pwiz-bin, Windows only). OpenMS does not use it, and it carries the
    ##   libraries of several instrument vendors, each under its own license terms. Users get
    ##   msconvert from ProteoWizard itself.
    ## - X!Tandem (XTandem). No OpenMS tool runs it since XTandemAdapter was removed in 3.4.0,
    ##   and its Linux build embeds expat 2.0.1, which Critical CVEs affect.
    if(SUBDIR STREQUAL "pwiz-bin" OR SUBDIR STREQUAL "XTandem")
      continue()
    endif()
    if(IS_DIRECTORY ${SEARCH_ENGINES_DIRECTORY}/${SUBDIR})
      install_thirdparty_folder("${SUBDIR}")
      ## The license and notice files the tool ships with, for THIRD-PARTY-NOTICES.txt.
      openms_third_party_notice_files_of_tool(_openms_tool_notices
                                              "${SEARCH_ENGINES_DIRECTORY}/${SUBDIR}")
      foreach(_openms_tool_notice IN LISTS _openms_tool_notices)
        openms_add_third_party_notice("THIRDPARTY/${SUBDIR}/${_openms_tool_notice}"
                                      "${SEARCH_ENGINES_DIRECTORY}/${SUBDIR}/${_openms_tool_notice}")
      endforeach()
    endif()
  endforeach()
endif()
