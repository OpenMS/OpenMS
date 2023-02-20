cmake_minimum_required(VERSION 3.10)
project(blosc_hdf5)
include(ExternalProject)

set(BLOSC_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/blosc")
set(BLOSC_INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}/blosc")
set(BLOSC_CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${BLOSC_INSTALL_DIR})

message("BLOSC_PREFIX='${BLOSC_PREFIX}'")
message("BLOSC_INSTALL_DIR='${BLOSC_INSTALL_DIR}'")
message("BLOSC_CMAKE_ARGS='${BLOSC_CMAKE_ARGS}'")

ExternalProject_Add(project_blosc
  PREFIX ${BLOSC_PREFIX}
  GIT_TAG v1.21.3
  GIT_REPOSITORY https://github.com/Blosc/c-blosc.git    
  INSTALL_DIR ${BLOSC_INSTALL_DIR}
  CMAKE_ARGS ${BLOSC_CMAKE_ARGS}
)

# sources
set(SOURCES src/blosc_filter.c)
set(PLUGIN_SOURCES src/blosc_filter.c src/blosc_plugin.c )

# dependencies
if(MSVC)
    # FindHDF5.cmake does not find Windows installations. Try to
    # use an environment variable instead until the official "find"
    # file can be updated for Windows.
    #
    # Note that you have to set this environment variable by hand.
    file(TO_CMAKE_PATH "$ENV{HDF5_DIR}" HDF5_HINT)
    set(HDF5_DIR ${HDF5_HINT} CACHE STRING "Path to HDF5 CMake config directory.")
    find_package(HDF5 REQUIRED HINTS ${HDF5_DIR})
else(MSVC)
    find_package(HDF5 REQUIRED)
endif(MSVC)
include_directories(${HDF5_INCLUDE_DIRS})


# add blosc libraries
add_library(blosc_static STATIC IMPORTED)
set_property(TARGET blosc_static PROPERTY IMPORTED_LOCATION ${BLOSC_INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}blosc${CMAKE_STATIC_LIBRARY_SUFFIX})
add_dependencies(blosc_static project_blosc)
include_directories(${BLOSC_INSTALL_DIR}/include)

add_library(blosc_filter_static STATIC ${SOURCES})
set_target_properties(
  blosc_filter_static PROPERTIES OUTPUT_NAME blosc_filter)
target_link_libraries(blosc_filter_static blosc_static ${HDF5_LIBRARIES})
