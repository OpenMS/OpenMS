# - Finds OpenMP support (for GNU, INTEL, MSVC compiler and cmake 2.6.0 to cmake 2.6.2)
#
# inspired by:
# Copyright 2008, 2009 <André Rigland Brodtkorb> Andre.Brodtkorb@ifi.uio.no
#
# Redistribution AND use is allowed according to the terms of the New 
# BSD license. 
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#
#
# The following variables are set:
#   OpenMP_CXX_FLAGS - flags to add to the CXX compiler for OpenMP support
#   OPENMP_FOUND - true if openmp is detected

include(CheckCSourceCompiles)
include(CheckCXXSourceCompiles)
include(FindPackageHandleStandardArgs)

# sample openmp source code to test
set(OpenMP_C_TEST_SOURCE 
"
#include <omp.h>
int main() { 
#ifdef _OPENMP
  return 0; 
#else
  breaks_on_purpose
#endif
}
")
# use the same source for CXX as C for now
set(OpenMP_CXX_TEST_SOURCE ${OpenMP_C_TEST_SOURCE})


if (MSVC)
	set(OpenMP_CXX_FLAGS_INTERNAL "/openmp")
elseif(CMAKE_COMPILER_IS_GNUCXX)
	set(FLAG "-fopenmp")
	set(SAFE_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
  set(CMAKE_REQUIRED_FLAGS "${FLAG}")
 	check_cxx_source_compiles("${OpenMP_C_TEST_SOURCE}" OpenMP_FLAG_DETECTED_GXX)
  set(CMAKE_REQUIRED_FLAGS "${SAFE_CMAKE_REQUIRED_FLAGS}")
  if(OpenMP_FLAG_DETECTED_GXX)
    set(OpenMP_CXX_FLAGS_INTERNAL "${FLAG}")
  endif()
else() ## just test Intel
	set(FLAG "-openmp")
	set(SAFE_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
  set(CMAKE_REQUIRED_FLAGS "${FLAG}")
 	check_cxx_source_compiles("${OpenMP_C_TEST_SOURCE}" OpenMP_FLAG_DETECTED_INTEL)
  set(CMAKE_REQUIRED_FLAGS "${SAFE_CMAKE_REQUIRED_FLAGS}")
  if(OpenMP_FLAG_DETECTED_INTEL)
    set(OpenMP_CXX_FLAGS_INTERNAL "${FLAG}")
  endif()
endif()

set(OpenMP_CXX_FLAGS "${OpenMP_CXX_FLAGS_INTERNAL}"
  CACHE STRING "C++ compiler flags for OpenMP parallization")
	
# handle the standard arguments for find_package
find_package_handle_standard_args(OpenMP DEFAULT_MSG 
  OpenMP_CXX_FLAGS )

mark_as_advanced(
  OpenMP_CXX_FLAGS
)
