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

set(OPENMS_LOGOSMALL_NAME openms_logo_corner_transparent.png) ## The filename of the logo to be used for the OpenMS folder e.g. on the PKG
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

# Normalize CMAKE_PREFIX_PATH entries to avoid mixed path separators (/ and \)
# This is especially important on Windows where paths might come from different sources
set(NORMALIZED_CMAKE_PREFIX_PATH "")
foreach(prefix_path ${CMAKE_PREFIX_PATH})
  if(prefix_path)
    file(TO_CMAKE_PATH "${prefix_path}" normalized_path)
    list(APPEND NORMALIZED_CMAKE_PREFIX_PATH "${normalized_path}")
  endif()
endforeach()

# Use normalized paths and append /bin and /lib subdirectories
list(TRANSFORM NORMALIZED_CMAKE_PREFIX_PATH APPEND "/bin" OUTPUT_VARIABLE DEP_BIN_DIRS)
list(TRANSFORM NORMALIZED_CMAKE_PREFIX_PATH APPEND "/lib" OUTPUT_VARIABLE DEP_LIB_DIRS)

# Print the paths for debugging
message(STATUS "Dependency search paths (bin): ${DEP_BIN_DIRS}")
message(STATUS "Dependency search paths (lib): ${DEP_LIB_DIRS}")

# Also include our own runtime directory where dependencies were copied during build.
# This serves as the primary source, with CMAKE_PREFIX_PATH directories as fallback.
list(APPEND DEP_BIN_DIRS $<TARGET_FILE_DIR:OpenMS>)

# Combine all search directories for comprehensive dependency resolution
set(RUNTIME_DEP_SEARCH_DIRS ${DEP_BIN_DIRS} ${DEP_LIB_DIRS})


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
                  ) ## TODO I have found that sometimes vcruntime140_1.dll cannot be resolved (but vcruntime140.dll is found). Maybe add it here too?
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

# TODO check if we can reduce the permissions
install(RUNTIME_DEPENDENCY_SET OPENMS_DEPS
        DESTINATION ${INSTALL_LIB_DIR}
        PERMISSIONS
          OWNER_READ OWNER_WRITE OWNER_EXECUTE
          GROUP_READ GROUP_WRITE GROUP_EXECUTE
          WORLD_READ WORLD_WRITE WORLD_EXECUTE
        COMPONENT Dependencies
        PRE_EXCLUDE_REGEXES ${PRE_EXCLUDE}
        POST_EXCLUDE_REGEXES ${POST_EXCLUDE}
        DIRECTORIES ${RUNTIME_DEP_SEARCH_DIRS})

#install(RUNTIME_DEPENDENCY_SET TOPPView_DEPS) # I think without giving DESTINATION and COMPONENT it will be inferred
#install(RUNTIME_DEPENDENCY_SET TOPPAS_DEPS)
#...

########################################################### SEARCHENGINES
set(THIRDPARTY_COMPONENT_GROUP)
## populates the THIRDPARTY_COMPONENT_GROUP list
if(EXISTS ${SEARCH_ENGINES_DIRECTORY})
  ## Automatically recurse over all subfolders in SEARCH_ENGINES_DIRECTORY
  file(GLOB THIRDPARTY_SUBDIRS RELATIVE ${SEARCH_ENGINES_DIRECTORY} ${SEARCH_ENGINES_DIRECTORY}/*)
  foreach(SUBDIR ${THIRDPARTY_SUBDIRS})
    if(IS_DIRECTORY ${SEARCH_ENGINES_DIRECTORY}/${SUBDIR})
      install_thirdparty_folder("${SUBDIR}")
    endif()
  endforeach()
endif()
