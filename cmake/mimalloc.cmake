# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# --------------------------------------------------------------------------
# $Maintainer: Oliver Kohlbacher $
# $Authors: Oliver Kohlbacher $
# --------------------------------------------------------------------------
# WITH_MIMALLOC: the TOPP tools allocate through mimalloc instead of the NT heap (Windows, MSVC).
#
# mimalloc's documented dynamic override for Windows (its readme, "Dynamic Override on Windows"): the
# executable links mimalloc's import library ahead of every other library (/INCLUDE:mi_version keeps
# the import), and mimalloc.dll and mimalloc-redirect.dll sit next to it. At load time the redirect
# patches the heap entry points of the UCRT, so the allocations of OpenMS.dll and of the static
# dependencies go to mimalloc as well, not only the executable's own. Measured with ProSE (HeLa,
# 3 files, Ryzen 9 7900X, MSVC 14.44): 12 threads 57 -> 42 s (-26 %, the NT heap's page faults and
# system time under contention), 1 thread -5 %, outputs identical.
#
# mimalloc comes from its upstream release through FetchContent, pinned by hash: the vcpkg port cannot
# provide the override under OpenMS's x64-windows-static-md triplet (its 'override' feature supports
# "!windows | !static": a static mimalloc has no DLL to redirect through), and a per-port linkage
# override in the triplet would change the ABI hash of every cached package. pyOpenMS keeps the NT
# heap: the override needs the executable's import table.

option(WITH_MIMALLOC "Windows/MSVC x64 only: the TOPP tools allocate through mimalloc (mimalloc.dll and mimalloc-redirect.dll next to the executables)." OFF)

## openms_link_mimalloc(<executable>): a no-op unless WITH_MIMALLOC is in effect. Tested through a global
## property, not `TARGET mimalloc`: other packages (e.g. vcpkg's mimalloc wrapper, pulled in by arrow) may
## define a target of that name.
function(openms_link_mimalloc target)
  get_property(_active GLOBAL PROPERTY OPENMS_MIMALLOC_ACTIVE)
  if(NOT _active)
    return()
  endif()
  target_link_libraries(${target} mimalloc) # build order
  ## Link options precede the libraries on the link line, so mimalloc.dll becomes the first import of
  ## the executable, as the override requires (check with mimalloc's minject -l <exe>).
  target_link_options(${target} PRIVATE /INCLUDE:mi_version $<TARGET_LINKER_FILE:mimalloc>)
endfunction()

if(NOT WITH_MIMALLOC)
  return()
endif()
## x64 only: the redirect DLL name and the undecorated /INCLUDE:mi_version symbol are x64-specific
## (mimalloc-redirect-arm64.dll, mimalloc-redirect32.dll, _mi_version on x86).
if(NOT (WIN32 AND MSVC AND CMAKE_CXX_COMPILER_ARCHITECTURE_ID STREQUAL "x64"))
  message(STATUS "WITH_MIMALLOC is Windows/MSVC x64-only; ignored on this platform.")
  return()
endif()
if(NOT BUILD_TOPP_TOOLS)
  message(STATUS "WITH_MIMALLOC: no TOPP tools in this build, nothing to link.")
  return()
endif()

include(FetchContent)
FetchContent_Declare(mimalloc
  URL https://github.com/microsoft/mimalloc/archive/refs/tags/v3.6.1.tar.gz
  URL_HASH SHA256=2672e6a99d4579fb4bb32d9de0d63aa424ce600d31fbfaf06eb74d940d2e07a6
  ## A subdirectory without a CMakeLists.txt: FetchContent_MakeAvailable() then only downloads, and the
  ## add_subdirectory() below adds mimalloc with EXCLUDE_FROM_ALL, which drops its install rules
  ## (headers, import library, CMake package); OpenMS installs the two DLLs and the license itself.
  SOURCE_SUBDIR cmake
)
if(TARGET mimalloc)
  message(FATAL_ERROR "WITH_MIMALLOC: a target 'mimalloc' already exists (e.g. from vcpkg's mimalloc wrapper); "
                      "the fetched mimalloc cannot be added next to it.")
endif()
## Settings for the mimalloc subproject only: normal variables in the top-level scope, restored afterwards
## (same pattern as HiGHS in cmake_findExternalLibs.cmake).
set(_openms_mimalloc_vars MI_BUILD_SHARED MI_BUILD_STATIC MI_BUILD_OBJECT MI_BUILD_TESTS MI_OVERRIDE MI_WIN_REDIRECT)
foreach(_var IN LISTS _openms_mimalloc_vars)
  if(DEFINED ${_var})
    set(_openms_mimalloc_saved_${_var} "${${_var}}")
  endif()
endforeach()
set(MI_BUILD_SHARED ON)
set(MI_BUILD_STATIC OFF)
set(MI_BUILD_OBJECT OFF)
set(MI_BUILD_TESTS OFF)
set(MI_OVERRIDE ON)
set(MI_WIN_REDIRECT ON)
FetchContent_MakeAvailable(mimalloc)
add_subdirectory(${mimalloc_SOURCE_DIR} ${mimalloc_BINARY_DIR} EXCLUDE_FROM_ALL)
foreach(_var IN LISTS _openms_mimalloc_vars)
  if(DEFINED _openms_mimalloc_saved_${_var})
    set(${_var} "${_openms_mimalloc_saved_${_var}}")
    unset(_openms_mimalloc_saved_${_var})
  else()
    unset(${_var})
  endif()
endforeach()
set_property(GLOBAL PROPERTY OPENMS_MIMALLOC_ACTIVE TRUE)

## mimalloc.dll is built into CMAKE_RUNTIME_OUTPUT_DIRECTORY next to the TOPP tools, and mimalloc's own
## post-build step copies mimalloc-redirect.dll (x64) there. Installed next to the executables; the license
## goes to LICENSES/vendored/mimalloc (openms_install_vendored_licenses()).
install(FILES $<TARGET_FILE:mimalloc> $<TARGET_FILE_DIR:mimalloc>/mimalloc-redirect.dll
        DESTINATION ${INSTALL_BIN_DIR} COMPONENT library)
message(STATUS "WITH_MIMALLOC: the TOPP tools allocate through mimalloc 3.6.1 (${mimalloc_SOURCE_DIR})")
