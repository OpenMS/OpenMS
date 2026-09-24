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
  add_test(NAME HeaderSets
    COMMAND ${CMAKE_COMMAND}
      "-DBINARY_DIR=${CMAKE_BINARY_DIR}/header-set-tests"
      "-DGENERATOR=${CMAKE_GENERATOR}"
      "-DGENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}"
      "-DGENERATOR_TOOLSET=${CMAKE_GENERATOR_TOOLSET}"
      ${_configure_args} -P "${CMAKE_SOURCE_DIR}/src/tests/header_sets/run.cmake")
  set_tests_properties(HeaderSets PROPERTIES LABELS "installed_consumer;header_sets")

  if(CMAKE_GENERATOR_PLATFORM)
    list(APPEND _configure_args -A "${CMAKE_GENERATOR_PLATFORM}")
  endif()
  if(CMAKE_GENERATOR_TOOLSET)
    list(APPEND _configure_args -T "${CMAKE_GENERATOR_TOOLSET}")
  endif()

  # The tests of one installation: <name>_clean (fixture <name>_cleaned), then
  # <name>_install_<component> for each install component in COMPONENTS (fixture
  # <name>_installed), then <name>_configure, <name>_build and <name>_run of the
  # consumer project in source_dir.
  #
  # The prefix is emptied before every run because `cmake --install --component`
  # only ever adds files: a component dropped from COMPONENTS would otherwise
  # leave the files of an earlier run of the same build tree behind, and since
  # OpenMSConfig.cmake reports a layer from the presence of its target file, the
  # consumer would keep finding a layer the installation no longer installs.
  #
  # RELOCATE_TO <dir> renames the finished installation to <dir> and points the
  # consumer there instead of at prefix, which is how the package is checked for
  # being relocatable. It adds <dir> to the directories <name>_clean empties, so
  # the rename has a free destination and the moved prefix holds this
  # installation and nothing else, and <name>_move (fixture <name>_moved), which
  # performs the rename. CONFIGURE_ARGS are passed to the consumer's configure
  # step on top of the toolchain arguments above.
  function(_openms_add_installed_consumer_tests name source_dir prefix consumer)
    cmake_parse_arguments(PARSE_ARGV 4 _ic "" "RELOCATE_TO" "COMPONENTS;CONFIGURE_ARGS")
    if(_ic_UNPARSED_ARGUMENTS)
      message(FATAL_ERROR "_openms_add_installed_consumer_tests(${name}): unexpected arguments '${_ic_UNPARSED_ARGUMENTS}'")
    endif()
    # Without this an empty COMPONENTS silently produces a chain that installs
    # nothing: the foreach below adds no tests, set_tests_properties() on an
    # empty list is a no-op, and <name>_configure then runs against whatever a
    # previous run happened to leave in the prefix.
    if(NOT _ic_COMPONENTS)
      message(FATAL_ERROR "_openms_add_installed_consumer_tests(${name}): COMPONENTS must name at least one install component")
    endif()

    # the prefix the consumer is pointed at, which is the one it ends up in
    set(_found_prefix "${prefix}")
    set(_clean_dirs "${prefix}")
    if(_ic_RELOCATE_TO)
      set(_found_prefix "${_ic_RELOCATE_TO}")
      list(APPEND _clean_dirs "${_ic_RELOCATE_TO}")
    endif()
    add_test(NAME ${name}_clean COMMAND ${CMAKE_COMMAND} -E rm -rf ${_clean_dirs})
    set_tests_properties(${name}_clean PROPERTIES FIXTURES_SETUP ${name}_cleaned)
    set(_step_tests ${name}_clean)

    set(_install_tests)
    foreach(_component IN LISTS _ic_COMPONENTS)
      add_test(NAME ${name}_install_${_component}
        COMMAND ${CMAKE_COMMAND} --install "${CMAKE_BINARY_DIR}" --config $<CONFIG>
                --prefix "${prefix}" --component ${_component})
      list(APPEND _install_tests ${name}_install_${_component})
    endforeach()
    # Every chain installs from this one build tree, and each `cmake --install
    # --component <c>` writes ${CMAKE_BINARY_DIR}/install_manifest_<c>.txt, so
    # concurrent installs of the same component interleave in that file and leave
    # it describing whichever scratch prefix finished last. They also contend for
    # install throughput rather than gaining from it. One lock across all chains.
    set_tests_properties(${_install_tests} PROPERTIES
      FIXTURES_REQUIRED ${name}_cleaned
      FIXTURES_SETUP ${name}_installed
      RESOURCE_LOCK openms_installed_consumer_install)
    list(APPEND _step_tests ${_install_tests})

    set(_configure_fixture ${name}_installed)
    if(_ic_RELOCATE_TO)
      add_test(NAME ${name}_move
        COMMAND ${CMAKE_COMMAND} -E rename "${prefix}" "${_ic_RELOCATE_TO}")
      set_tests_properties(${name}_move PROPERTIES
        FIXTURES_REQUIRED ${name}_installed
        FIXTURES_SETUP ${name}_moved)
      list(APPEND _step_tests ${name}_move)
      set(_configure_fixture ${name}_moved)
    endif()

    add_test(NAME ${name}_configure
      COMMAND ${CMAKE_COMMAND} -S "${source_dir}" -B "${consumer}"
              -G "${CMAKE_GENERATOR}" "-DCMAKE_BUILD_TYPE=$<CONFIG>"
              "-DOpenMS_DIR=${_found_prefix}/${INSTALL_CMAKE_DIR}"
              ${_ic_CONFIGURE_ARGS} ${_configure_args})
    set_tests_properties(${name}_configure PROPERTIES
      FIXTURES_REQUIRED ${_configure_fixture}
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

    set_tests_properties(${_step_tests} ${name}_configure ${name}_build ${name}_run
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
  set(_full_prefix "${CMAKE_BINARY_DIR}/installed-consumer/prefix")
  _openms_add_installed_consumer_tests(TestExternalCode
    "${CMAKE_SOURCE_DIR}/src/tests/external"
    "${_full_prefix}"
    "${CMAKE_BINARY_DIR}/installed-consumer/build"
    COMPONENTS ${_full_components})

  # 1b. Configure the same consumer against the same installation once more, with Eigen
  #     made unfindable: no installed header includes Eigen, and a shared OpenMS carries
  #     no Eigen3::Eigen in its exported link interface, so find_package(OpenMS REQUIRED)
  #     has to succeed without it. Both regressions this guards against are reported
  #     before anything is compiled: a find_dependency(Eigen3) reintroduced into
  #     OpenMSConfig.cmake is a hard error, because CMAKE_DISABLE_FIND_PACKAGE_Eigen3
  #     refuses a REQUIRED lookup, and an Eigen3::Eigen left dangling in an exported
  #     target file fails while the consumer's build system is generated ("but the target
  #     was not found"), which is part of this very invocation.
  #
  #     Configure only, deliberately. For a shared installation the build system generated
  #     here is identical to the one TestExternalCode_configure produces -- same compile
  #     commands, same link commands -- because disabling the Eigen *package* changes
  #     CMake-level discovery, not compiler include paths, and TestExternalCode_build
  #     already compiles and links exactly that. An installed header that includes
  #     <Eigen/Core> again fails in that build too, but only so long as nothing hands
  #     consumers a path Eigen resolves from. That proviso is not free: a soft lookup
  #     reinstated in OpenMSConfig.cmake could append Eigen's include directories to the
  #     imported target, leaving the link interface clean and the normal build green while
  #     every consumer without Eigen breaks. src/tests/external/CMakeLists.txt asserts
  #     against precisely that, by testing whether <Eigen/Core> resolves from any include
  #     directory the package exposes, so that combination is caught here at configure
  #     time instead of needing a second consumer build to surface it.
  #
  #     Shared builds only, which is what the OpenMS build produces (BUILD_SHARED_LIBS is
  #     set unconditionally in the top-level CMakeLists.txt). A static OpenMS keeps its
  #     private dependencies in the export, so its OpenMSConfig.cmake does look Eigen up
  #     -- disabling Eigen would then make find_package(OpenMS REQUIRED) fail by design,
  #     which is correct behaviour rather than a regression. The matching assertion in
  #     src/tests/external/CMakeLists.txt makes the same distinction consumer-side, by
  #     asking the imported target for its TYPE.
  if(BUILD_SHARED_LIBS)
    add_test(NAME TestExternalCode_noeigen_configure
      COMMAND ${CMAKE_COMMAND} -S "${CMAKE_SOURCE_DIR}/src/tests/external"
              -B "${CMAKE_BINARY_DIR}/installed-consumer/build-no-eigen"
              -G "${CMAKE_GENERATOR}" "-DCMAKE_BUILD_TYPE=$<CONFIG>"
              "-DOpenMS_DIR=${_full_prefix}/${INSTALL_CMAKE_DIR}"
              -DCMAKE_DISABLE_FIND_PACKAGE_Eigen3=TRUE ${_configure_args})
    set_tests_properties(TestExternalCode_noeigen_configure PROPERTIES
      FIXTURES_REQUIRED TestExternalCode_installed
      LABELS installed_consumer)
  endif()

  # 2. The core-only installation: exactly the components the pyOpenMS wheel
  #    builds install (.github/workflows/pyopenms-wheels-cibuildwheel.yml) ...
  set(_core_prefix "${CMAKE_BINARY_DIR}/installed-consumer/prefix-core")
  _openms_add_installed_consumer_tests(TestExternalCodeCore
    "${CMAKE_SOURCE_DIR}/src/tests/external/core_only"
    "${_core_prefix}"
    "${CMAKE_BINARY_DIR}/installed-consumer/build-core"
    COMPONENTS library OpenMS_headers OpenSwathAlgo_headers thirdparty_headers cmake share)
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

  # 3. A relocated installation: the prefix is renamed once the install is
  #    complete and the consumer is pointed at the new location. That is what
  #    unpacking a binary package produces -- the prefix a user installs to is
  #    never the one the release was built with -- and it is what
  #    OpenMSConfig.cmake deriving everything from its own location buys (#4603).
  #    Installed into a staging prefix of its own and then *renamed*, rather than
  #    copied, so that the prefix it was installed to is really gone: a path baked
  #    in at install time would still resolve against a copy's original.
  #    The destination has a space in its name on purpose, because <prefix> is
  #    "C:\Program Files\OpenMS" for every user of the Windows installer.
  #
  #    The core-only component set without share: the shared data relocates
  #    through the same PACKAGE_PREFIX_DIR arithmetic as the directory variables
  #    the consumer does check, so copying ~200 MB a third time would cost more
  #    than it proves. The consumer accordingly stays off the shared data.
  set(_relocate_prefix "${CMAKE_BINARY_DIR}/installed-consumer/prefix-staged")
  set(_relocate_moved "${CMAKE_BINARY_DIR}/installed-consumer/relocated prefix")
  _openms_add_installed_consumer_tests(TestExternalCodeRelocated
    "${CMAKE_SOURCE_DIR}/src/tests/external/relocated"
    "${_relocate_prefix}"
    "${CMAKE_BINARY_DIR}/installed-consumer/build-relocated"
    COMPONENTS library OpenMS_headers OpenSwathAlgo_headers thirdparty_headers cmake
    RELOCATE_TO "${_relocate_moved}"
    CONFIGURE_ARGS "-DOPENMS_RELOCATED_FROM=${_relocate_prefix}"
                   "-DOPENMS_RELOCATED_TO=${_relocate_moved}")

  # 4. The package of a build tree across reconfigurations: the GUI target file
  #    of a WITH_GUI=ON configuration must not survive a WITH_GUI=OFF
  #    reconfiguration of the same build directory, since OpenMSConfig.cmake
  #    reports a layer from the presence of its file. Reconfiguring OpenMS
  #    itself twice is far too slow for a test, so src/tests/package_layers, a
  #    stub project using the real macros and package template, stands in.
  set(_layers_source "${CMAKE_SOURCE_DIR}/src/tests/package_layers")
  set(_layers_build  "${CMAKE_BINARY_DIR}/package-layers/build")
  add_test(NAME PackageLayers_configure_with_gui
    COMMAND ${CMAKE_COMMAND} -S "${_layers_source}" -B "${_layers_build}"
            -G "${CMAKE_GENERATOR}" -DWITH_GUI=ON ${_configure_args})
  set_tests_properties(PackageLayers_configure_with_gui PROPERTIES
    FIXTURES_SETUP package_layers_with_gui)
  add_test(NAME PackageLayers_probe_with_gui
    COMMAND ${CMAKE_COMMAND} -S "${_layers_source}/probe" -B "${_layers_build}-probe-with-gui"
            -G "${CMAKE_GENERATOR}" "-DOpenMS_DIR=${_layers_build}" -DEXPECT_GUI=ON ${_configure_args})
  set_tests_properties(PackageLayers_probe_with_gui PROPERTIES
    FIXTURES_REQUIRED package_layers_with_gui
    FIXTURES_SETUP package_layers_probed_with_gui)
  add_test(NAME PackageLayers_reconfigure_without_gui
    COMMAND ${CMAKE_COMMAND} -S "${_layers_source}" -B "${_layers_build}" -DWITH_GUI=OFF)
  set_tests_properties(PackageLayers_reconfigure_without_gui PROPERTIES
    FIXTURES_REQUIRED package_layers_probed_with_gui
    FIXTURES_SETUP package_layers_without_gui)
  add_test(NAME PackageLayers_probe_without_gui
    COMMAND ${CMAKE_COMMAND} -S "${_layers_source}/probe" -B "${_layers_build}-probe-without-gui"
            -G "${CMAKE_GENERATOR}" "-DOpenMS_DIR=${_layers_build}" -DEXPECT_GUI=OFF ${_configure_args})
  set_tests_properties(PackageLayers_probe_without_gui PROPERTIES
    FIXTURES_REQUIRED package_layers_without_gui)
  add_test(NAME PackageLayers_reject_gui_without_gui
    COMMAND ${CMAKE_COMMAND} -S "${_layers_source}/probe" -B "${_layers_build}-probe-reject-gui"
            -G "${CMAKE_GENERATOR}" "-DOpenMS_DIR=${_layers_build}" -DREQUIRE_GUI=ON ${_configure_args})
  set_tests_properties(PackageLayers_reject_gui_without_gui PROPERTIES
    FIXTURES_REQUIRED package_layers_without_gui
    PASS_REGULAR_EXPRESSION "built without the GUI library")
  set_tests_properties(PackageLayers_configure_with_gui PackageLayers_probe_with_gui
    PackageLayers_reconfigure_without_gui PackageLayers_probe_without_gui
    PackageLayers_reject_gui_without_gui PROPERTIES LABELS "installed_consumer;package_layers")

  unset(_layers_source)
  unset(_layers_build)
  unset(_full_components)
  unset(_full_prefix)
  unset(_core_prefix)
  unset(_relocate_prefix)
  unset(_relocate_moved)
  unset(_configure_args)
endif()
