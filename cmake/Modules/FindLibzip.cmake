# FindLibzip.cmake
# Finds libzip and creates the libzip::zip imported target.
#
# We intentionally avoid find_package(libzip CONFIG) because libzip <= 1.10
# ships a single CMake targets file that defines both the library and CLI tool
# targets (zipcmp, zipmerge, ziptool). If the tools package is not installed,
# the config file raises a FATAL_ERROR (which QUIET cannot suppress).
#
# Instead we always do a manual header+library search, which works on all
# platforms (apt, brew, conda, vcpkg, contrib).

find_path(LIBZIP_INCLUDE_DIR
  NAMES zip.h
)

find_library(LIBZIP_LIBRARY
  NAMES zip libzip
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Libzip
  REQUIRED_VARS LIBZIP_LIBRARY LIBZIP_INCLUDE_DIR
)

mark_as_advanced(LIBZIP_INCLUDE_DIR LIBZIP_LIBRARY)

if(Libzip_FOUND AND NOT TARGET libzip::zip)
  add_library(libzip::zip UNKNOWN IMPORTED)
  set_target_properties(libzip::zip PROPERTIES
    IMPORTED_LOCATION "${LIBZIP_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${LIBZIP_INCLUDE_DIR}"
  )
endif()
