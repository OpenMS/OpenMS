cmake_minimum_required(VERSION 3.19 FATAL_ERROR)

### example CMakeLists.txt to develop C++ programs using OpenMS
project("Example_Project_using_OpenMS")

## list all your executables here (a corresponding .cpp file should exist, e.g. Main.cpp)
set(my_executables
  Main
)

## list all classes here, which are required by your executables
## (all these classes will be linked into a library)
set(my_sources
  ExampleLibraryFile.cpp
)

## find OpenMS: provides the imported targets OpenMS::OpenMS (the library) and OpenMS::OpenSwathAlgo.
## If this fails, point CMake at an installation or build tree with -DOpenMS_DIR=<prefix>/lib/cmake/OpenMS
## (<prefix>/CMake on Windows, or the OpenMS build directory), or add <prefix> to CMAKE_PREFIX_PATH.
find_package(OpenMS CONFIG REQUIRED)
message(STATUS "Found OpenMS ${OpenMS_VERSION} at ${OpenMS_DIR}")

## library with additional classes from above
add_library(my_custom_lib STATIC ${my_sources})
target_link_libraries(my_custom_lib PUBLIC OpenMS::OpenMS)

## add targets for the executables
foreach(i ${my_executables})
  add_executable(${i} ${i}.cpp)
  ## link executables against OpenMS
  target_link_libraries(${i} PRIVATE OpenMS::OpenMS my_custom_lib)
endforeach(i)
