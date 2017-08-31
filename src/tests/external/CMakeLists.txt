### example CMakeLists.txt to develop programs using OpenMS
project("OpenMS_ExternalCodeTest")
cmake_minimum_required(VERSION 2.8.3)

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
find_package(OpenMS PATHS "$ENV{OPENMS_BUILD_TREE}" NO_CMAKE_PACKAGE_REGISTRY)

# check whether the OpenMS package was found
if (OpenMS_FOUND)

  ## include directories for OpenMS headers (and contrib)
  include_directories(${OpenMS_INCLUDE_DIRECTORIES})

  ## append precompiler macros specific to OpenMS
  ## Warning: this could be harmful to your project. Check this
  ## if problems occur
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPENMS_ADDCXX_FLAGS}")
  add_definitions(${OPENMS_DEFINITIONS})

  ## library with additional classes from above
  add_library(my_custom_lib STATIC ${my_sources})
  target_link_libraries(my_custom_lib OpenMS)

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
