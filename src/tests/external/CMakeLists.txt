cmake_minimum_required(VERSION 3.22 FATAL_ERROR)

### example CMakeLists.txt to develop programs using OpenMS
project("OpenMS_ExternalCodeTest")

include(CTest)

message("Received CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")
if (DEFINED OPENMS_CONTRIB_LIBS)
  message("Received OPENMS_CONTRIB_LIBS: ${OPENMS_CONTRIB_LIBS}")
endif()

## list all your executables here (a corresponding .cpp file should exist, e.g. TestExternalCode.cpp)
set(my_executables
	TestExternalCode
)

## list all classes here, which are required by your executables
## (all these classes will be linked into a library)
set(my_sources
	ExampleLibraryFile.cpp
)

## Find the OpenMS package. It provides the imported targets OpenMS::OpenMS and
## OpenMS::OpenSwathAlgo, OpenMS::OpenMS_CLI when the installation includes the TOPP
## tool framework (OpenMS_CLI_FOUND; a TOPP-style tool requests COMPONENTS CLI) and,
## for an installation built with WITH_GUI that includes the GUI library,
## OpenMS::OpenMS_GUI, whose Qt6 modules are found when the GUI component is requested.
## Select an installation with -DOpenMS_DIR=<prefix>/lib/cmake/OpenMS
## (<prefix>/CMake on Windows) or by adding <prefix> to CMAKE_PREFIX_PATH.
find_package(OpenMS CONFIG REQUIRED COMPONENTS CLI OPTIONAL_COMPONENTS GUI)
message(STATUS "Found OpenMS ${OpenMS_VERSION} at ${OpenMS_DIR}")

#------------------------------------------------------------------------------
## Package contract. This project doubles as the installed-consumer test of the
## OpenMS build (OPENMS_TEST_INSTALLED_CONSUMER) for the full development
## installation (core_only/ covers a core-only one), so it checks what the package
## promises before using it. A plain consumer project does not need this block.

## the CLI layer of the full installation was found ...
if(NOT OpenMS_CLI_FOUND)
  message(FATAL_ERROR "OpenMSConfig.cmake found the required CLI component without setting OpenMS_CLI_FOUND")
endif()
## ... the libraries are namespaced imported targets ...
foreach(_target IN ITEMS OpenMS::OpenMS OpenMS::OpenSwathAlgo OpenMS::OpenMS_CLI)
  if(NOT TARGET ${_target})
    message(FATAL_ERROR "OpenMSConfig.cmake did not provide the target ${_target}")
  endif()
  get_target_property(_imported ${_target} IMPORTED)
  if(NOT _imported)
    message(FATAL_ERROR "${_target} is not an imported target")
  endif()
endforeach()
## ... the un-namespaced names of earlier releases remain as aliases ...
foreach(_alias IN ITEMS OpenMS OpenSwathAlgo OpenMS_CLI)
  if(NOT TARGET ${_alias})
    message(FATAL_ERROR "OpenMSConfig.cmake did not provide the compatibility alias ${_alias}")
  endif()
endforeach()
## ... and the bundled third-party libraries exist under the OpenMS:: namespace only
foreach(_leaked IN ITEMS IsoSpec Evergreen eol-bspline opentims_cpp openms_thermo_bridge)
  if(TARGET ${_leaked})
    message(FATAL_ERROR "OpenMSConfig.cmake leaked the un-namespaced third-party target ${_leaked}")
  endif()
endforeach()

## Qt6 is not part of the core library's link interface (see TestExternalCodeQt below)
get_target_property(_core_interface OpenMS::OpenMS INTERFACE_LINK_LIBRARIES)
if(_core_interface MATCHES "Qt6::")
  message(FATAL_ERROR "OpenMS::OpenMS leaks Qt6 into its link interface: ${_core_interface}")
endif()

## the directory variables describe this installation, not the prefix of a dependency
if(NOT EXISTS "${OPENMS_DATA_DIR}/CHEMISTRY/unimod.xml")
  message(FATAL_ERROR "OPENMS_DATA_DIR='${OPENMS_DATA_DIR}' does not hold the OpenMS shared data")
endif()
if(NOT IS_DIRECTORY "${OPENMS_LIB_DIR}")
  message(FATAL_ERROR "OPENMS_LIB_DIR='${OPENMS_LIB_DIR}' is not a directory")
endif()
file(REAL_PATH "${OPENMS_LIB_DIR}" _lib_dir)
get_target_property(_configurations OpenMS::OpenMS IMPORTED_CONFIGURATIONS)
if(NOT _configurations)
  message(FATAL_ERROR "OpenMS::OpenMS has no imported configurations")
endif()
foreach(_configuration IN LISTS _configurations)
  get_target_property(_location OpenMS::OpenMS IMPORTED_LOCATION_${_configuration})
  get_filename_component(_location_dir "${_location}" DIRECTORY)
  file(REAL_PATH "${_location_dir}" _location_dir)
  if(NOT _location_dir STREQUAL _lib_dir)
    message(FATAL_ERROR "OPENMS_LIB_DIR='${OPENMS_LIB_DIR}' but the ${_configuration} OpenMS library is ${_location}")
  endif()
endforeach()

## GUI component: requested as optional above, so a runner whose Qt6 is not findable
## from this project (unlike from the OpenMS build) skips the check instead of failing.
if(OpenMS_WITH_GUI)
  if(OpenMS_GUI_FOUND)
    foreach(_target IN ITEMS OpenMS::OpenMS_GUI OpenMS_GUI Qt6::Core Qt6::Widgets)
      if(NOT TARGET ${_target})
        message(FATAL_ERROR "the OpenMS GUI component did not provide the target ${_target}")
      endif()
    endforeach()
    ## linking the GUI library resolves its Qt6 link interface, which is what the component supplies
    add_executable(TestExternalCodeGUI TestExternalCodeGUI.cpp)
    target_link_libraries(TestExternalCodeGUI OpenMS::OpenMS_GUI)
    add_test(TestExternalCodeGUI TestExternalCodeGUI)
  else()
    message(STATUS "Qt6 not found from this project - skipping the GUI component check.")
  endif()
elseif(OpenMS_GUI_FOUND)
  message(FATAL_ERROR "OpenMSConfig.cmake reports the GUI component for an installation without the GUI library")
endif()

#------------------------------------------------------------------------------
## library with additional classes from above
add_library(my_custom_lib STATIC ${my_sources})
## Projects written for earlier releases link the un-namespaced name. This library
## keeps that spelling on purpose, so the compatibility alias stays covered.
target_link_libraries(my_custom_lib OpenMS)
## disable auto-linking of boost (mainly needed for Windows. CMake takes care of choosing the right library)
target_compile_definitions(my_custom_lib PUBLIC BOOST_ALL_NO_LIB)

## add targets for the executables
foreach(i ${my_executables})
  # create the executable
  add_executable(${i} ${i}.cpp)
  ## link executables against OpenMS
  target_link_libraries(${i} OpenMS::OpenMS my_custom_lib)
endforeach(i)

## Downstream link contract: since the Qt removal, OpenMSConfig.cmake no longer
## declares Qt6 as a (transitive) dependency of the core library. An external consumer that wants Qt must
## find and link it itself. This target uses a Qt type and links Qt6::Core *explicitly* after finding Qt6
## here -- it builds only because this project supplies Qt, not because OpenMS::OpenMS leaks it
## transitively. Guarded so the external-code test still configures cleanly on runners without Qt6.
find_package(Qt6 COMPONENTS Core QUIET)
if (Qt6_FOUND)
  add_executable(TestExternalCodeQt TestExternalCodeQt.cpp)
  target_link_libraries(TestExternalCodeQt OpenMS::OpenMS my_custom_lib Qt6::Core)
  add_test(TestExternalCodeQt TestExternalCodeQt)
else()
  message(STATUS "Qt6 not found - skipping TestExternalCodeQt (transitive-Qt link-contract check).")
endif()

## A TOPP-style tool of an external project: derives from TOPPBase, so it links the
## tool framework OpenMS::OpenMS_CLI (which carries OpenMS::OpenMS transitively).
add_executable(TestExternalCodeCLI TestExternalCodeCLI.cpp)
target_link_libraries(TestExternalCodeCLI OpenMS::OpenMS_CLI)
add_test(TestExternalCodeCLI TestExternalCodeCLI -greeting "hello from OpenMS_CLI")
set_tests_properties(TestExternalCodeCLI PROPERTIES PASS_REGULAR_EXPRESSION "hello from OpenMS_CLI")

## Enable testing - for Nightly Build log
add_test(TestExternalCode TestExternalCode)

## Windows has no RPATH: executables find OpenMS.dll (and the DLLs of the
## dependencies OpenMS was built against) through PATH. OpenMSConfig.cmake
## tells us where the OpenMS libraries were installed.
if(WIN32)
  foreach(_test IN ITEMS TestExternalCode TestExternalCodeCLI TestExternalCodeQt TestExternalCodeGUI)
    if(TEST ${_test})
      set_property(TEST ${_test} APPEND PROPERTY
        ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${OPENMS_LIB_DIR}")
    endif()
  endforeach()
endif()
