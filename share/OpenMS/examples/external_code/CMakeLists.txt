### example CMakeLists.txt to develop C++ programs using OpenMS
project("Example_Project_using_OpenMS")
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
 
## list all your executables here (a corresponding .cpp file should exist, e.g. Main.cpp)
set(my_executables
  Main
)
 
## list all classes here, which are required by your executables
## (all these classes will be linked into a library)
set(my_sources
  ExampleLibraryFile.cpp
)
 
## find OpenMS configuration and register target "OpenMS" (our library)
find_package(OpenMS)
## if the above fails you can try calling cmake with -D OpenMS_DIR=/path/to/OpenMS/
## or modify the find_package() call accordingly
## find_package(OpenMS PATHS "</path/to/OpenMS//")
 
# check whether the OpenMS package was found
if (OpenMS_FOUND)
  message(STATUS "\nFound OpenMS at ${OpenMS_DIR}\n")
    
  ## library with additional classes from above
  add_library(my_custom_lib STATIC ${my_sources})
 
  ## add targets for the executables
  foreach(i ${my_executables})
    add_executable(${i} ${i}.cpp)
    ## link executables against OpenMS
    target_link_libraries(${i} OpenMS my_custom_lib)
  endforeach(i)
 
 
else(OpenMS_FOUND)
  message(FATAL_ERROR "OpenMSConfig.cmake file not found!")
endif(OpenMS_FOUND)
