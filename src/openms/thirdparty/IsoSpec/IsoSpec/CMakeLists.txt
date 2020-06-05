project("IsoSpec")
cmake_minimum_required(VERSION 2.6)

set(my_sources
cwrapper.cpp
allocator.cpp  
dirtyAllocator.cpp 
isoSpec++.cpp 
isoMath.cpp 
marginalTrek++.cpp 
operators.cpp
element_tables.cpp
misc.cpp
)
SET_SOURCE_FILES_PROPERTIES( ${my_sources} PROPERTIES LANGUAGE CXX )

# export include directory
set(ISOSPEC_INCLUDE_DIRECTORY
  ${PROJECT_SOURCE_DIR}
  CACHE
  INTERNAL "isospec include directories" FORCE
)

## platform dependent compiler flags:
include(CheckCXXCompilerFlag)
if (NOT WIN32) # we only want fPIC on non-windows systems (fPIC is implicitly true there)
  CHECK_CXX_COMPILER_FLAG("-fPIC" WITH_FPIC)
  if (WITH_FPIC)
    add_definitions(-fPIC)
  endif()

  # Only add those definitions on non-windows systems
  add_definitions(-std=c++11 -Wall -pedantic -Wextra)
endif()

## add_library(IsoSpec SHARED ${my_sources})
## add_executable(mytest test.cpp)
## target_link_libraries(mytest IsoSpec)

configure_file(IsoSpecConfig.cmake.in "${PROJECT_BINARY_DIR}/IsoSpecConfig.cmake" @ONLY)

## export(TARGETS IsoSpec FILE IsoSpecLibrary.cmake)

