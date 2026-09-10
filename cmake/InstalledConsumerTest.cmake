# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
#
# Installs the development components of a finished OpenMS build into a
# scratch prefix below the build tree, configures src/tests/external against
# that installation, builds it and runs its CTest.
#
# Runs in CMake script mode. It is registered as a CTest test by
# src/tests/CMakeLists.txt, which passes every setting explicitly:
#
#   cmake -DBUILD_DIR=<openms build dir> -DSOURCE_DIR=<openms source dir>
#         -DCONFIG=<Release|Debug|...> -DGENERATOR=<CMAKE_GENERATOR>
#         -DCTEST_COMMAND=<ctest> -DINSTALL_CMAKE_DIR=... -DINSTALL_LIB_DIR=...
#         -DINSTALL_DATA_PATH=<compiled-in install data dir>
#         [-DCMAKE_TOOLCHAIN_FILE=... -DVCPKG_INSTALLED_DIR=... ...]
#         -P cmake/InstalledConsumerTest.cmake
#
# Nothing is read back from CMakeCache.txt: the values are the ordinary
# configure-time variables of the OpenMS build, forwarded by the caller.

cmake_minimum_required(VERSION 3.21)

foreach(_required IN ITEMS BUILD_DIR SOURCE_DIR CONFIG GENERATOR CTEST_COMMAND INSTALL_CMAKE_DIR INSTALL_LIB_DIR)
  if(NOT DEFINED ${_required} OR "${${_required}}" STREQUAL "")
    message(FATAL_ERROR "InstalledConsumerTest: ${_required} must be passed with -D")
  endif()
endforeach()

set(scratch  "${BUILD_DIR}/installed-consumer")
set(prefix   "${scratch}/prefix")
set(consumer "${scratch}/build")

# On Linux and macOS the library probes the data directory compiled in from
# CMAKE_INSTALL_PREFIX before anything else. An OpenMS installation already
# present there would therefore shadow the scratch installation this test
# makes, and the consumer's installed-data check could not tell the two apart.
string(FIND "${INSTALL_DATA_PATH}" "${prefix}/" _data_below_prefix)
if(NOT CMAKE_HOST_WIN32 AND INSTALL_DATA_PATH
   AND EXISTS "${INSTALL_DATA_PATH}/CHEMISTRY/unimod.xml"
   AND NOT _data_below_prefix EQUAL 0)
  message(FATAL_ERROR "InstalledConsumerTest: an OpenMS installation exists at the "
    "compiled-in install prefix (${INSTALL_DATA_PATH}). It would shadow the "
    "scratch installation under ${prefix}. Remove or move that installation, "
    "or configure with a CMAKE_INSTALL_PREFIX that has no OpenMS installed.")
endif()

# Print a command, run it, and fail the test on a non-zero exit code.
function(run)
  list(JOIN ARGN " " _shown)
  message(STATUS "+ ${_shown}")
  execute_process(COMMAND ${ARGN} RESULT_VARIABLE _rc)
  if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "command failed with exit code ${_rc}: ${_shown}")
  endif()
endfunction()

#------------------------------------------------------------------------------
# 0. Step 4 moves the source data directory aside while the consumer runs. If an
#    earlier run was interrupted between its two renames, the data is still in
#    the backup location and the "share" component below could not be installed.
#    Put it back first.
set(_source_data "${SOURCE_DIR}/share/OpenMS")
set(_saved_data  "${SOURCE_DIR}/share/OpenMS.installed-consumer-backup")
if(EXISTS "${_saved_data}")
  if(EXISTS "${_source_data}")
    message(FATAL_ERROR "Both ${_source_data} and ${_saved_data} exist; remove the stale backup first.")
  endif()
  message(WARNING "Restoring ${_source_data} from the backup left by an interrupted run")
  file(RENAME "${_saved_data}" "${_source_data}")
endif()

#------------------------------------------------------------------------------
# 1. Install the development package into a scratch prefix.
#    Application, documentation and redistribution components are left out:
#    this test is a development-package consumer.
file(REMOVE_RECURSE "${scratch}")
foreach(_component IN ITEMS library OpenMS_headers OpenSwathAlgo_headers cmake share)
  run("${CMAKE_COMMAND}" --install "${BUILD_DIR}" --config "${CONFIG}"
      --prefix "${prefix}" --component ${_component})
endforeach()

set(package_dir "${prefix}/${INSTALL_CMAKE_DIR}")
if(NOT EXISTS "${package_dir}/OpenMSConfig.cmake")
  message(FATAL_ERROR "Installed OpenMSConfig.cmake is missing from ${package_dir}")
endif()

#------------------------------------------------------------------------------
# 2. Runtime environment. The installed shared library keeps its private
#    dependencies (Xerces, CURL, ...) private, so both the native linker and the
#    consumer executable have to find them through the loader path, exactly as
#    an external user's environment would provide them.
set(_runtime_dirs "${prefix}/${INSTALL_LIB_DIR}" "${prefix}/bin")
set(_dependency_prefixes ${CMAKE_PREFIX_PATH})
if(VCPKG_INSTALLED_DIR AND VCPKG_TARGET_TRIPLET)
  list(APPEND _dependency_prefixes "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}")
endif()
string(TOLOWER "${CONFIG}" _config_lower)
foreach(_dep IN LISTS _dependency_prefixes)
  if(_config_lower STREQUAL "debug")
    list(APPEND _runtime_dirs "${_dep}/debug/bin" "${_dep}/debug/lib")
  endif()
  list(APPEND _runtime_dirs "${_dep}/bin" "${_dep}/lib")
endforeach()

set(_native_dirs)
foreach(_dir IN LISTS _runtime_dirs)
  file(TO_NATIVE_PATH "${_dir}" _native)
  list(APPEND _native_dirs "${_native}")
endforeach()
if(CMAKE_HOST_WIN32)
  list(JOIN _native_dirs ";" _path_prefix)
else()
  list(JOIN _native_dirs ":" _path_prefix)
endif()
foreach(_var IN ITEMS PATH LD_LIBRARY_PATH DYLD_LIBRARY_PATH)
  if(DEFINED ENV{${_var}} AND NOT "$ENV{${_var}}" STREQUAL "")
    if(CMAKE_HOST_WIN32)
      set(ENV{${_var}} "${_path_prefix};$ENV{${_var}}")
    else()
      set(ENV{${_var}} "${_path_prefix}:$ENV{${_var}}")
    endif()
  else()
    set(ENV{${_var}} "${_path_prefix}")
  endif()
endforeach()
# The consumer must resolve data through the installation, never through an
# environment override left over from the OpenMS build.
unset(ENV{OPENMS_DATA_PATH})

#------------------------------------------------------------------------------
# 3. Configure the consumer against the installation only.
string(TOUPPER "${CONFIG}" _config_upper)
set(_configure
  "${CMAKE_COMMAND}" -S "${SOURCE_DIR}/src/tests/external" -B "${consumer}"
  -G "${GENERATOR}"
  "-DOpenMS_DIR=${package_dir}"
  "-DOPENMS_EXPECTED_PREFIX=${prefix}"
  "-DCMAKE_BUILD_TYPE=${CONFIG}"
  # The executable lives in <prefix>/bin so OpenMS finds its data exe-relative
  # (../share/OpenMS), which is how an installed consumer resolves it.
  "-DCMAKE_RUNTIME_OUTPUT_DIRECTORY=${prefix}/bin"
  "-DCMAKE_RUNTIME_OUTPUT_DIRECTORY_${_config_upper}=${prefix}/bin"
  -DCMAKE_FIND_USE_PACKAGE_REGISTRY=OFF
  -DCMAKE_FIND_USE_SYSTEM_PACKAGE_REGISTRY=OFF
  -DCMAKE_DISABLE_FIND_PACKAGE_CURL=ON
  -DCMAKE_DISABLE_FIND_PACKAGE_XercesC=ON
  -DVCPKG_MANIFEST_MODE=OFF)

# Reuse the compiler and dependency installation of the OpenMS build. Each of
# these is forwarded by the registering CMakeLists only when it is set there.
# They go through an initial-cache file (-C) rather than -D arguments: the
# command line is a CMake list, and a list-valued setting such as
# CMAKE_PREFIX_PATH or CMAKE_OSX_ARCHITECTURES would be split at its
# semicolons into separate arguments. Bracket arguments keep every character
# of the value, including semicolons, quotes and backslashes.
set(_initial_cache "${scratch}/consumer-initial-cache.cmake")
set(_cache_content "# generated by InstalledConsumerTest.cmake\n")
foreach(_var IN ITEMS
    CMAKE_TOOLCHAIN_FILE VCPKG_INSTALLED_DIR VCPKG_TARGET_TRIPLET VCPKG_HOST_TRIPLET
    CMAKE_PREFIX_PATH CMAKE_C_COMPILER CMAKE_CXX_COMPILER CMAKE_MAKE_PROGRAM
    CMAKE_MSVC_RUNTIME_LIBRARY CMAKE_OSX_ARCHITECTURES CMAKE_OSX_DEPLOYMENT_TARGET
    CMAKE_OSX_SYSROOT)
  if(DEFINED ${_var} AND NOT "${${_var}}" STREQUAL "")
    string(APPEND _cache_content "set(${_var} [==[${${_var}}]==] CACHE STRING \"forwarded from the OpenMS build\")\n")
  endif()
endforeach()
file(WRITE "${_initial_cache}" "${_cache_content}")
list(APPEND _configure -C "${_initial_cache}")
if(GENERATOR_PLATFORM)
  list(APPEND _configure -A "${GENERATOR_PLATFORM}")
endif()
if(GENERATOR_TOOLSET)
  list(APPEND _configure -T "${GENERATOR_TOOLSET}")
endif()

run(${_configure})
run("${CMAKE_COMMAND}" --build "${consumer}" --config "${CONFIG}" --parallel 2)

#------------------------------------------------------------------------------
# 4. Run the consumer with the source data directory moved aside. OpenMS probes
#    the compiled-in build-tree path (OPENMS_DATA_PATH) before the exe-relative
#    ../share/OpenMS, and skips the compiled-in install prefix on Windows, so
#    while share/OpenMS exists the consumer would silently read the source tree
#    and a missing install rule could not be detected. The rename is undone
#    whether or not the tests pass, which is why this test is registered
#    RUN_SERIAL: other tests read that directory.
file(RENAME "${_source_data}" "${_saved_data}")
# The consumer tests take seconds. The bound stays well below the outer CTest
# TIMEOUT (1800 s) so that a hanging consumer test is killed here, where the
# rename is undone, and not by the outer ctest, which would leave share/OpenMS
# moved aside.
execute_process(
  COMMAND "${CTEST_COMMAND}" --test-dir "${consumer}" -C "${CONFIG}" --output-on-failure --no-tests=error
  TIMEOUT 600
  RESULT_VARIABLE _rc)
file(RENAME "${_saved_data}" "${_source_data}")
if(NOT _rc EQUAL 0)
  message(FATAL_ERROR "Installed consumer tests failed: ${_rc}")
endif()
