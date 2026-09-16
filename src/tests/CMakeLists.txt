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
# Installed-consumer tests: install the development package into scratch prefixes
# below the build tree, then configure, build and run consumer projects against
# them the way an external project does: find_package(OpenMS) with OpenMS_DIR
# pointing at the installation. Each step is one plain cmake/ctest invocation,
# chained through CTest fixtures. Off by default because it performs installs;
# the CI presets enable it.
#
# The package is layered (cmake/install_macros.cmake), so two installations are
# tested:
#  * the full development installation (core and CLI layer; the GUI layer too in
#    a WITH_GUI build) with src/tests/external, the example consumer project;
#  * a core-only installation, the components the pyOpenMS wheels are built
#    against, with src/tests/external/core_only: the package has to work without
#    the CLI and GUI layers, report them absent and refuse a required CLI
#    component.
option(OPENMS_TEST_INSTALLED_CONSUMER
  "Install into scratch prefixes and build/run the external CMake consumers against them" OFF)

if(OPENMS_TEST_INSTALLED_CONSUMER AND "${PACKAGE_TYPE}" STREQUAL "none")
  # Only the toolchain of this build is forwarded to the consumers, as a consumer
  # building with the same compiler and dependency provider would pass it;
  # everything about OpenMS comes from the installed OpenMSConfig.cmake.
  set(_configure_args)
  foreach(_var IN ITEMS
      CMAKE_TOOLCHAIN_FILE VCPKG_INSTALLED_DIR VCPKG_TARGET_TRIPLET VCPKG_HOST_TRIPLET
      CMAKE_PREFIX_PATH CMAKE_C_COMPILER CMAKE_CXX_COMPILER CMAKE_MAKE_PROGRAM
      CMAKE_MSVC_RUNTIME_LIBRARY CMAKE_OSX_ARCHITECTURES CMAKE_OSX_DEPLOYMENT_TARGET
      CMAKE_OSX_SYSROOT)
    if(DEFINED ${_var} AND NOT "${${_var}}" STREQUAL "")
      # keep list-valued settings such as CMAKE_PREFIX_PATH as one argument
      string(REPLACE ";" "$<SEMICOLON>" _value "${${_var}}")
      list(APPEND _configure_args "-D${_var}=${_value}")
    endif()
  endforeach()
  if(CMAKE_GENERATOR_PLATFORM)
    list(APPEND _configure_args -A "${CMAKE_GENERATOR_PLATFORM}")
  endif()
  if(CMAKE_GENERATOR_TOOLSET)
    list(APPEND _configure_args -T "${CMAKE_GENERATOR_TOOLSET}")
  endif()

  # The tests of one installation: <name>_install_<component> for each of the
  # install components given after the consumer build directory (fixture
  # <name>_installed), then <name>_configure, <name>_build and <name>_run of
  # the consumer project in source_dir.
  function(_openms_add_installed_consumer_tests name source_dir prefix consumer)
    set(_install_tests)
    foreach(_component IN LISTS ARGN)
      add_test(NAME ${name}_install_${_component}
        COMMAND ${CMAKE_COMMAND} --install "${CMAKE_BINARY_DIR}" --config $<CONFIG>
                --prefix "${prefix}" --component ${_component})
      list(APPEND _install_tests ${name}_install_${_component})
    endforeach()
    set_tests_properties(${_install_tests} PROPERTIES FIXTURES_SETUP ${name}_installed)

    add_test(NAME ${name}_configure
      COMMAND ${CMAKE_COMMAND} -S "${source_dir}" -B "${consumer}"
              -G "${CMAKE_GENERATOR}" "-DCMAKE_BUILD_TYPE=$<CONFIG>"
              "-DOpenMS_DIR=${prefix}/${INSTALL_CMAKE_DIR}" ${_configure_args})
    set_tests_properties(${name}_configure PROPERTIES
      FIXTURES_REQUIRED ${name}_installed
      FIXTURES_SETUP ${name}_configured)

    add_test(NAME ${name}_build
      COMMAND ${CMAKE_COMMAND} --build "${consumer}" --config $<CONFIG>)
    set_tests_properties(${name}_build PROPERTIES
      FIXTURES_REQUIRED ${name}_configured
      FIXTURES_SETUP ${name}_built)
    add_test(NAME ${name}_run
      COMMAND ${CMAKE_CTEST_COMMAND} --test-dir "${consumer}" -C $<CONFIG>
              --output-on-failure --no-tests=error)
    set_tests_properties(${name}_run PROPERTIES FIXTURES_REQUIRED ${name}_built)

    set_tests_properties(${_install_tests} ${name}_configure ${name}_build ${name}_run
      PROPERTIES LABELS installed_consumer)
  endfunction()

  # 1. The full development installation (no applications, documentation or
  #    packaging extras): what a library consumer that also builds TOPP-style
  #    tools has. The GUI layer comes without its headers, as TestExternalCodeGUI
  #    only links the library.
  set(_full_components
    library library_cli OpenMS_headers OpenMS_CLI_headers OpenSwathAlgo_headers cmake cmake_cli share)
  if(WITH_GUI)
    list(APPEND _full_components library_gui cmake_gui)
  endif()
  _openms_add_installed_consumer_tests(TestExternalCode
    "${CMAKE_SOURCE_DIR}/src/tests/external"
    "${CMAKE_BINARY_DIR}/installed-consumer/prefix"
    "${CMAKE_BINARY_DIR}/installed-consumer/build"
    ${_full_components})

  # 2. The core-only installation: exactly the components the pyOpenMS wheel
  #    builds install (.github/workflows/pyopenms-wheels-cibuildwheel.yml) ...
  set(_core_prefix "${CMAKE_BINARY_DIR}/installed-consumer/prefix-core")
  _openms_add_installed_consumer_tests(TestExternalCodeCore
    "${CMAKE_SOURCE_DIR}/src/tests/external/core_only"
    "${_core_prefix}"
    "${CMAKE_BINARY_DIR}/installed-consumer/build-core"
    library OpenMS_headers OpenSwathAlgo_headers thirdparty_headers cmake share)
  #    ... which has to refuse a consumer that requires the CLI component. That
  #    configure step fails, and the test passes on the reason the package gives.
  add_test(NAME TestExternalCodeCore_reject_cli
    COMMAND ${CMAKE_COMMAND} -S "${CMAKE_SOURCE_DIR}/src/tests/external/core_only"
            -B "${CMAKE_BINARY_DIR}/installed-consumer/build-core-reject"
            -G "${CMAKE_GENERATOR}" "-DCMAKE_BUILD_TYPE=$<CONFIG>"
            -DOPENMS_CORE_ONLY_REQUIRE_CLI=ON
            "-DOpenMS_DIR=${_core_prefix}/${INSTALL_CMAKE_DIR}" ${_configure_args})
  set_tests_properties(TestExternalCodeCore_reject_cli PROPERTIES
    FIXTURES_REQUIRED TestExternalCodeCore_installed
    PASS_REGULAR_EXPRESSION "does not include the TOPP tool framework library"
    LABELS installed_consumer)

  unset(_full_components)
  unset(_core_prefix)
  unset(_configure_args)
endif()
