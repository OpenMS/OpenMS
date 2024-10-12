cmake_minimum_required(VERSION 3.20)
project(blosc_hdf5)

set(BLOSC_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/blosc")
set(BLOSC_INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}/blosc")
set(BLOSC_CMAKE_ARGS 
  -DCMAKE_INSTALL_PREFIX=${BLOSC_INSTALL_DIR} 
  -DBUILD_SHARED=OFF 
  -DBUILD_TESTS=OFF
  -DBUILD_FUZZERS=OFF
  -DBUILD_BENCHMARKS=OFF
  )

# download at build time
include(FetchContent)
FetchContent_Declare(project_blosc
  PREFIX ${BLOSC_PREFIX}
  GIT_TAG v1.21.3
  GIT_REPOSITORY https://github.com/Blosc/c-blosc.git 
  INSTALL_DIR ${BLOSC_INSTALL_DIR}
  CMAKE_ARGS ${BLOSC_CMAKE_ARGS}
)
FetchContent_MakeAvailable(project_blosc)

# add blosc filter library
# sources
set(SOURCES src/blosc_filter.c src/blosc_filter.h)
add_library(blosc_filter_static STATIC ${SOURCES})
set_target_properties(blosc_filter_static PROPERTIES OUTPUT_NAME blosc_filter)
target_link_libraries(blosc_filter_static blosc_static HDF5::HDF5)
target_include_directories(blosc_filter_static PUBLIC src/)

