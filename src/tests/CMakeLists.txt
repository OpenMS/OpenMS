# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche, Chris Bielow $
# --------------------------------------------------------------------------


project("OpenMS_tests")

set(_openms_topp_testing_default ON)
if(DEFINED BUILD_TOPP_TOOLS)
  set(_openms_topp_testing_default ${BUILD_TOPP_TOOLS})
endif()
option(ENABLE_TOPP_TESTING "Enables tests for TOPP. Should be disabled only on time constraints (e.g. chunking during continuous integration)." ${_openms_topp_testing_default})
option(ENABLE_CLASS_TESTING "Enables tests for library classes. Should be disabled only on time constraints (e.g. chunking during continuous integration)." ON)
option(ENABLE_PIPELINE_TESTING "Enables the additional testing of various TOPPAS pipelines when 'make test' is called." ${_openms_topp_testing_default})
unset(_openms_topp_testing_default)

if(DEFINED BUILD_TOPP_TOOLS
   AND NOT BUILD_TOPP_TOOLS
   AND (ENABLE_TOPP_TESTING OR ENABLE_PIPELINE_TESTING))
  message(FATAL_ERROR
    "ENABLE_TOPP_TESTING and ENABLE_PIPELINE_TESTING require BUILD_TOPP_TOOLS to be ON.")
endif()

#------------------------------------------------------------------------------
# we only test if we have no package target
if("${PACKAGE_TYPE}" STREQUAL "none")
  ## configure the regular class ..
  if(ENABLE_CLASS_TESTING)
    add_subdirectory(class_tests)
  endif()
  # TOPP Test project ..
  if(ENABLE_TOPP_TESTING)
    add_subdirectory(topp)
  endif()
  # and if requested also TOPPAS pipeline tests
  if(ENABLE_PIPELINE_TESTING)
    add_subdirectory(toppas)
  endif()
elseif( ENABLE_CLASS_TESTING OR ENABLE_TOPP_TESTING OR ENABLE_PIPELINE_TESTING)
  message(WARNING "PACKAGE_TYPE specified along with enabled testing. Tests are only enabled when PACKAGE_TYPE is none. Ignoring tests.")
endif("${PACKAGE_TYPE}" STREQUAL "none")

#------------------------------------------------------------------------------
# Installed-consumer test: install the development package into a scratch
# prefix below the build tree and build/run src/tests/external against it
# (see cmake/InstalledConsumerTest.cmake). Off by default: it performs an
# install and moves share/OpenMS aside while it runs. The CI presets enable it.
option(OPENMS_TEST_INSTALLED_CONSUMER
  "Install into a scratch prefix and build/run the external CMake consumer against it" OFF)

if(OPENMS_TEST_INSTALLED_CONSUMER AND "${PACKAGE_TYPE}" STREQUAL "none")
  set(_consumer_args
    "-DBUILD_DIR=${CMAKE_BINARY_DIR}"
    "-DSOURCE_DIR=${CMAKE_SOURCE_DIR}"
    "-DCONFIG=$<CONFIG>"
    "-DGENERATOR=${CMAKE_GENERATOR}"
    "-DGENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}"
    "-DGENERATOR_TOOLSET=${CMAKE_GENERATOR_TOOLSET}"
    "-DCTEST_COMMAND=${CMAKE_CTEST_COMMAND}"
    "-DINSTALL_CMAKE_DIR=${INSTALL_CMAKE_DIR}"
    "-DINSTALL_LIB_DIR=${INSTALL_LIB_DIR}"
    # the data directory compiled into the library (CMAKE_INSTALL_PREFIX based)
    "-DINSTALL_DATA_PATH=${CF_OPENMS_INSTALL_DATA_PATH}")
  # Forward the toolchain and dependency settings of this build so the consumer
  # compiles with the same compiler against the same vcpkg/contrib installation.
  foreach(_var IN ITEMS
      CMAKE_TOOLCHAIN_FILE VCPKG_INSTALLED_DIR VCPKG_TARGET_TRIPLET VCPKG_HOST_TRIPLET
      CMAKE_PREFIX_PATH CMAKE_C_COMPILER CMAKE_CXX_COMPILER CMAKE_MAKE_PROGRAM
      CMAKE_MSVC_RUNTIME_LIBRARY CMAKE_OSX_ARCHITECTURES CMAKE_OSX_DEPLOYMENT_TARGET
      CMAKE_OSX_SYSROOT)
    if(DEFINED ${_var} AND NOT "${${_var}}" STREQUAL "")
      # keep list-valued settings such as CMAKE_PREFIX_PATH as one argument
      string(REPLACE ";" "$<SEMICOLON>" _value "${${_var}}")
      list(APPEND _consumer_args "-D${_var}=${_value}")
    endif()
  endforeach()
  add_test(NAME TestExternalCode_installed
    COMMAND ${CMAKE_COMMAND} ${_consumer_args}
            -P "${CMAKE_SOURCE_DIR}/cmake/InstalledConsumerTest.cmake")
  # Moves share/OpenMS aside while it runs, which every other test reads.
  set_tests_properties(TestExternalCode_installed PROPERTIES
    RUN_SERIAL TRUE
    TIMEOUT 1800
    LABELS "installed_consumer")
  unset(_consumer_args)
endif()
