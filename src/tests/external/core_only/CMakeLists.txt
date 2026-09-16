cmake_minimum_required(VERSION 3.22 FATAL_ERROR)

### Installed-consumer test of a core-only OpenMS installation: the install
### components library, OpenMS_headers, OpenSwathAlgo_headers, thirdparty_headers,
### share and cmake, without the CLI and GUI layers (see cmake/install_macros.cmake).
### This is what the pyOpenMS wheels are built against. The package of such an
### installation has to work on its own and to report exactly what is installed.
### Driven by the TestExternalCodeCore_* tests of src/tests/CMakeLists.txt
### (OPENMS_TEST_INSTALLED_CONSUMER). For an example of an external project see
### ../CMakeLists.txt.
project("OpenMS_CoreOnlyConsumerTest")

include(CTest)

option(OPENMS_CORE_ONLY_REQUIRE_CLI
  "Require the CLI component, which a core-only installation has to refuse" OFF)

if(OPENMS_CORE_ONLY_REQUIRE_CLI)
  ## has to fail: the driving test passes on the reason OpenMSConfig.cmake gives
  find_package(OpenMS CONFIG REQUIRED COMPONENTS CLI)
  message(FATAL_ERROR "find_package(OpenMS COMPONENTS CLI) succeeded on a core-only installation")
endif()

## The optional components are absent here; the package itself is found.
find_package(OpenMS CONFIG REQUIRED OPTIONAL_COMPONENTS CLI GUI)
message(STATUS "Found OpenMS ${OpenMS_VERSION} at ${OpenMS_DIR}")

#------------------------------------------------------------------------------
## Package contract of the core-only installation

## the core layer is there, with the compatibility aliases ...
foreach(_target IN ITEMS OpenMS::OpenMS OpenMS::OpenSwathAlgo OpenMS OpenSwathAlgo)
  if(NOT TARGET ${_target})
    message(FATAL_ERROR "the core-only installation did not provide the target ${_target}")
  endif()
endforeach()
## ... and nothing of the other layers: neither their targets ...
foreach(_target IN ITEMS OpenMS::OpenMS_CLI OpenMS_CLI OpenMS::OpenMS_GUI OpenMS_GUI)
  if(TARGET ${_target})
    message(FATAL_ERROR "the core-only installation provided the target ${_target}")
  endif()
endforeach()
## ... nor a claim to have them
foreach(_var IN ITEMS OpenMS_CLI_FOUND OpenMS_GUI_FOUND OpenMS_WITH_GUI)
  if(${_var})
    message(FATAL_ERROR "the core-only installation reports ${_var}=${${_var}}")
  endif()
endforeach()

#------------------------------------------------------------------------------
## A program of the core layer builds and runs against it: the sources of the
## full consumer test, linked against the core library only.
add_library(my_custom_lib STATIC ../ExampleLibraryFile.cpp)
target_link_libraries(my_custom_lib OpenMS::OpenMS)
## disable auto-linking of boost (mainly needed for Windows. CMake takes care of choosing the right library)
target_compile_definitions(my_custom_lib PUBLIC BOOST_ALL_NO_LIB)

add_executable(TestExternalCode ../TestExternalCode.cpp)
target_link_libraries(TestExternalCode OpenMS::OpenMS my_custom_lib)
add_test(TestExternalCode TestExternalCode)

## Windows has no RPATH: the executable finds OpenMS.dll through PATH.
if(WIN32)
  set_property(TEST TestExternalCode APPEND PROPERTY
    ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${OPENMS_LIB_DIR}")
endif()
