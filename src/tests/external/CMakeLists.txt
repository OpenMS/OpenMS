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

## find OpenMS package and register target "OpenMS" (our library)
## Select an installation with -DOpenMS_DIR=<prefix>/lib/cmake/OpenMS
## (<prefix>/CMake on Windows) or by adding <prefix> to CMAKE_PREFIX_PATH.
find_package(OpenMS CONFIG REQUIRED)
message(STATUS "Found OpenMS at ${OpenMS_DIR}")

# check whether the OpenMS package was found
if (OpenMS_FOUND)

  ## library with additional classes from above
  add_library(my_custom_lib STATIC ${my_sources})
  target_link_libraries(my_custom_lib OpenMS)
  ## disable auto-linking of boost (mainly needed for Windows. CMake takes care of choosing the right library)
  target_compile_definitions(my_custom_lib PUBLIC BOOST_ALL_NO_LIB)

  ## add targets for the executables
  foreach(i ${my_executables})
    # create the executable
    add_executable(${i} ${i}.cpp)
    ## link executables against OpenMS
    target_link_libraries(${i} OpenMS my_custom_lib)
  endforeach(i)

  ## Downstream link contract: since the Qt removal, OpenMSConfig.cmake no longer
  ## declares Qt6 as a (transitive) dependency. An external consumer that wants Qt must find and link
  ## it itself. This target uses a Qt type and links Qt6::Core *explicitly* after finding Qt6 here --
  ## it builds only because this project supplies Qt, not because OpenMS leaks it transitively. Guarded
  ## so the external-code test still configures cleanly on runners without Qt6.
  find_package(Qt6 COMPONENTS Core QUIET)
  if (Qt6_FOUND)
    add_executable(TestExternalCodeQt TestExternalCodeQt.cpp)
    target_link_libraries(TestExternalCodeQt OpenMS my_custom_lib Qt6::Core)
    add_test(TestExternalCodeQt TestExternalCodeQt)
  else()
    message(STATUS "Qt6 not found - skipping TestExternalCodeQt (transitive-Qt link-contract check).")
  endif()

else(OpenMS_FOUND)
  message(FATAL_ERROR "OpenMSConfig.cmake file not found!")
endif(OpenMS_FOUND)

## Enable testing - for Nightly Build log
add_test(TestExternalCode TestExternalCode)

## Windows has no RPATH: executables find OpenMS.dll (and the DLLs of the
## dependencies OpenMS was built against) through PATH. OpenMSConfig.cmake
## tells us where the OpenMS libraries were installed.
if(WIN32)
  set_property(TEST TestExternalCode APPEND PROPERTY
    ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${OPENMS_LIB_DIR}")
  if(TARGET TestExternalCodeQt)
    set_property(TEST TestExternalCodeQt APPEND PROPERTY
      ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${OPENMS_LIB_DIR}")
  endif()
endif()
