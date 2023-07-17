## Use CMake 3.17 to correctly find OpenMP on macOS
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

### example CMakeLists.txt to develop programs using OpenMS
project("OpenMS_ExternalCodeTest")

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
## Note: This is customized to fit the nightly test scenario. In a
##       regular build find_package(OpenMS) should be sufficient.
## This will also try to find all dependencies of OpenMS, so make sure
##  they are findable by e.g. providing CMAKE_PREFIX_PATH
find_package(OpenMS PATHS "$ENV{OPENMS_BUILD_TREE}" NO_CMAKE_PACKAGE_REGISTRY)

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

else(OpenMS_FOUND)
  message(FATAL_ERROR "OpenMSConfig.cmake file not found!")
endif(OpenMS_FOUND)

## Enable testing - for Nightly Build log
include(Dart)
add_test(TestExternalCode TestExternalCode)
