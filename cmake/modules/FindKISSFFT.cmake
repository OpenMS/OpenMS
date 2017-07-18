# - Try to find KISSFFT
#
#  BZIP2_INCLUDE_DIR - the BZip2 include directory
#  BZIP2_LIBRARY - Link these to use BZip2


include(${CMAKE_CURRENT_LIST_DIR}/SelectLibraryConfigurations.cmake)

set(KISSFFT_ROOT_DIR "" CACHE PATH "KISSFFT root directory")

find_path(KISSFFT_INCLUDE_DIR kiss_fft.h
  HINTS ${KISSFFT_ROOT_DIR}/include
)
find_library(KISSFFT_LIBRARY_RELEASE NAMES libkfft.a
    HINTS ${KISSFFT_ROOT_DIR}/lib
)
select_library_configurations(KISSFFT)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(KISSFFT DEFAULT_MSG
  KISSFFT_INCLUDE_DIR
  KISSFFT_LIBRARY
)

mark_as_advanced(
  KISSFFT_INCLUDE_DIR
  KISSFFT_LIBRARY
  )
