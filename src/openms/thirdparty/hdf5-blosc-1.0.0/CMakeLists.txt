cmake_minimum_required(VERSION 3.20)
project(blosc_hdf5)
include(ExternalProject)

set(BLOSC_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/blosc")
set(BLOSC_INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}/blosc")
set(BLOSC_CMAKE_ARGS 
  -DCMAKE_INSTALL_PREFIX=${BLOSC_INSTALL_DIR} 
  -DBUILD_SHARED=OFF 
  -DBUILD_TESTS=OFF
  -DBUILD_FUZZERS=OFF
  -DBUILD_BENCHMARKS=OFF
  )

message("BLOSC_PREFIX='${BLOSC_PREFIX}'")
message("BLOSC_INSTALL_DIR='${BLOSC_INSTALL_DIR}'")
message("BLOSC_CMAKE_ARGS='${BLOSC_CMAKE_ARGS}'")

# download at build time
ExternalProject_Add(project_blosc
  PREFIX ${BLOSC_PREFIX}
  GIT_TAG v1.21.3
  GIT_REPOSITORY https://github.com/Blosc/c-blosc.git 
  INSTALL_DIR ${BLOSC_INSTALL_DIR}
  CMAKE_ARGS ${BLOSC_CMAKE_ARGS}
)


# add blosc-c library
add_library(blosc_static STATIC IMPORTED)
file(MAKE_DIRECTORY ${BLOSC_INSTALL_DIR}/include) # cmake hack: create include folder otherwise cmake complains about it at configure time
set_target_properties(blosc_static PROPERTIES 
  IMPORTED_LOCATION ${BLOSC_INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}blosc${CMAKE_STATIC_LIBRARY_SUFFIX}
  INTERFACE_INCLUDE_DIRECTORIES ${BLOSC_INSTALL_DIR}/include
  )

# add blosc filter library
# sources
set(SOURCES src/blosc_filter.c src/blosc_filter.h)
add_library(blosc_filter_static STATIC ${SOURCES})
add_dependency(blosc_filter_static project_blosc)
set_target_properties(blosc_filter_static PROPERTIES OUTPUT_NAME blosc_filter)
target_link_libraries(blosc_filter_static blosc_static HDF5::HDF5)
target_include_directories(blosc_filter_static PUBLIC src/)

