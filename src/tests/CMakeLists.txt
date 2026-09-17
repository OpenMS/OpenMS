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
# Installed-consumer test: install the development package into a scratch prefix
# below the build tree, then configure, build and run src/tests/external against
# it the way an external project does: find_package(OpenMS) with OpenMS_DIR
# pointing at the installation. Each step is one plain cmake/ctest invocation,
# chained through CTest fixtures. Off by default because it performs an
# install; the CI presets enable it.
option(OPENMS_TEST_INSTALLED_CONSUMER
  "Install into a scratch prefix and build/run the external CMake consumer against it" OFF)

if(OPENMS_TEST_INSTALLED_CONSUMER AND "${PACKAGE_TYPE}" STREQUAL "none")
  set(_prefix   "${CMAKE_BINARY_DIR}/installed-consumer/prefix")
  set(_consumer "${CMAKE_BINARY_DIR}/installed-consumer/build")

  # 1. Install the development components (no applications, documentation or
  #    packaging extras): this is what a library consumer has.
  set(_install_tests)
  foreach(_component IN ITEMS library OpenMS_headers OpenMS_CLI_headers OpenSwathAlgo_headers cmake share)
    add_test(NAME TestExternalCode_install_${_component}
      COMMAND ${CMAKE_COMMAND} --install "${CMAKE_BINARY_DIR}" --config $<CONFIG>
              --prefix "${_prefix}" --component ${_component})
    list(APPEND _install_tests TestExternalCode_install_${_component})
  endforeach()
  set_tests_properties(${_install_tests} PROPERTIES FIXTURES_SETUP openms_installed)

  # 2. Configure the consumer against the installation. Only the toolchain of
  #    this build is forwarded, as a consumer building with the same compiler
  #    and dependency provider would pass it; everything about OpenMS comes
  #    from the installed OpenMSConfig.cmake.
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
  add_test(NAME TestExternalCode_configure
    COMMAND ${CMAKE_COMMAND} -S "${CMAKE_SOURCE_DIR}/src/tests/external" -B "${_consumer}"
            -G "${CMAKE_GENERATOR}" "-DCMAKE_BUILD_TYPE=$<CONFIG>"
            "-DOpenMS_DIR=${_prefix}/${INSTALL_CMAKE_DIR}" ${_configure_args})
  set_tests_properties(TestExternalCode_configure PROPERTIES
    FIXTURES_REQUIRED openms_installed
    FIXTURES_SETUP consumer_configured)

  # 3. Build it and run its own tests.
  add_test(NAME TestExternalCode_build
    COMMAND ${CMAKE_COMMAND} --build "${_consumer}" --config $<CONFIG>)
  set_tests_properties(TestExternalCode_build PROPERTIES
    FIXTURES_REQUIRED consumer_configured
    FIXTURES_SETUP consumer_built)
  add_test(NAME TestExternalCode_run
    COMMAND ${CMAKE_CTEST_COMMAND} --test-dir "${_consumer}" -C $<CONFIG>
            --output-on-failure --no-tests=error)
  set_tests_properties(TestExternalCode_run PROPERTIES FIXTURES_REQUIRED consumer_built)

  set(_consumer_tests TestExternalCode_configure TestExternalCode_build TestExternalCode_run)

  # 4. Configure and build the same consumer once more against the same installation,
  #    with Eigen made unfindable. No installed header includes Eigen and a shared
  #    OpenMS carries no Eigen3::Eigen in its exported link interface, so finding and
  #    linking the package has to succeed without it. Either regression fails here: a
  #    find_dependency(Eigen3) reintroduced into OpenMSConfig.cmake reports the package
  #    as not found, and an Eigen3::Eigen left dangling in OpenMSTargets.cmake is
  #    diagnosed when the consumer's build system is generated.
  #
  #    Shared builds only, which is what the OpenMS build produces (BUILD_SHARED_LIBS is
  #    set unconditionally in the top-level CMakeLists.txt). A static OpenMS keeps its
  #    private dependencies in the export, so its OpenMSConfig.cmake does look Eigen up
  #    -- disabling Eigen would then make find_package(OpenMS REQUIRED) fail by design,
  #    which is correct behaviour rather than a regression. The matching assertion in
  #    src/tests/external/CMakeLists.txt makes the same distinction consumer-side, by
  #    asking the imported target for its TYPE.
  if(BUILD_SHARED_LIBS)
    set(_consumer_noeigen "${CMAKE_BINARY_DIR}/installed-consumer/build-no-eigen")
    add_test(NAME TestExternalCode_noeigen_configure
      COMMAND ${CMAKE_COMMAND} -S "${CMAKE_SOURCE_DIR}/src/tests/external" -B "${_consumer_noeigen}"
              -G "${CMAKE_GENERATOR}" "-DCMAKE_BUILD_TYPE=$<CONFIG>"
              "-DOpenMS_DIR=${_prefix}/${INSTALL_CMAKE_DIR}"
              -DCMAKE_DISABLE_FIND_PACKAGE_Eigen3=TRUE ${_configure_args})
    set_tests_properties(TestExternalCode_noeigen_configure PROPERTIES
      FIXTURES_REQUIRED openms_installed
      FIXTURES_SETUP consumer_noeigen_configured)
    add_test(NAME TestExternalCode_noeigen_build
      COMMAND ${CMAKE_COMMAND} --build "${_consumer_noeigen}" --config $<CONFIG>)
    set_tests_properties(TestExternalCode_noeigen_build PROPERTIES
      FIXTURES_REQUIRED consumer_noeigen_configured)
    list(APPEND _consumer_tests TestExternalCode_noeigen_configure TestExternalCode_noeigen_build)
    unset(_consumer_noeigen)
  endif()

  set_tests_properties(${_install_tests} ${_consumer_tests}
    PROPERTIES LABELS installed_consumer)
  unset(_consumer_tests)
  unset(_install_tests)
  unset(_configure_args)
  unset(_consumer)
  unset(_prefix)
endif()
