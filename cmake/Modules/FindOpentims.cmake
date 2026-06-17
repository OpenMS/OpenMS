# - Try to find opentims (Bruker TimsTOF .d file reading library)
# Once done this will define:
#
#  Opentims_FOUND        - system has opentims
#  Opentims_INCLUDE_DIR  - the opentims include directory (parent of opentims++/)
#  Opentims_LIBRARIES    - link these to use opentims
#
# And the imported target:
#  opentims::opentims_cpp
#
# Standard CMake mechanisms (all honoured automatically):
#  opentims_DIR          - directory containing opentimsConfig.cmake
#  opentims_ROOT         - root of the opentims installation (CMake 3.12+ / CMP0074)
#  CMAKE_PREFIX_PATH     - list of installation prefixes to search
#
# Legacy OpenMS hint (kept for backwards compatibility):
#  Opentims_ROOT_DIR     - same role as opentims_ROOT

# Prefer a CMake config package produced by opentims's own install rules.
# opentims_ROOT / opentims_DIR / CMAKE_PREFIX_PATH are checked automatically
# by CMake (CMP0074); Opentims_ROOT_DIR is the legacy OpenMS-specific hint.
find_package(opentims CONFIG QUIET
  HINTS
    "${Opentims_ROOT_DIR}"
    "${Opentims_ROOT_DIR}/lib/cmake/opentims"
)
if(opentims_FOUND)
  if(NOT TARGET opentims::opentims_cpp)
    message(FATAL_ERROR "opentims config package found but target opentims::opentims_cpp is missing")
  endif()
  set(Opentims_FOUND TRUE)
  get_target_property(Opentims_INCLUDE_DIR opentims::opentims_cpp INTERFACE_INCLUDE_DIRECTORIES)
  # IMPORTED_LOCATION is not set on interface-only or alias targets; try
  # configuration-specific locations before giving up.
  get_target_property(Opentims_LIBRARIES opentims::opentims_cpp IMPORTED_LOCATION)
  if(Opentims_LIBRARIES MATCHES "NOTFOUND")
    foreach(_cfg RELEASE DEBUG RELWITHDEBINFO MINSIZEREL)
      get_target_property(_loc opentims::opentims_cpp "IMPORTED_LOCATION_${_cfg}")
      if(_loc AND NOT _loc MATCHES "NOTFOUND")
        set(Opentims_LIBRARIES "${_loc}")
        break()
      endif()
    endforeach()
    if(Opentims_LIBRARIES MATCHES "NOTFOUND")
      set(Opentims_LIBRARIES "")
    endif()
  endif()
  return()
endif()

# Fall back to manual search when no config package is present.
find_path(Opentims_INCLUDE_DIR
  NAMES opentims++/opentims.h
  HINTS
    "${opentims_ROOT}/include"
    "${opentims_ROOT}"
    "${Opentims_ROOT_DIR}/include"
    "${Opentims_ROOT_DIR}"
  DOC "opentims include directory (parent of opentims++/)"
)

find_library(Opentims_LIBRARY
  NAMES opentims_cpp
  HINTS
    "${opentims_ROOT}/lib"
    "${opentims_ROOT}"
    "${Opentims_ROOT_DIR}/lib"
    "${Opentims_ROOT_DIR}"
  DOC "opentims library"
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Opentims
  REQUIRED_VARS Opentims_LIBRARY Opentims_INCLUDE_DIR
)

if(Opentims_FOUND)
  set(Opentims_LIBRARIES "${Opentims_LIBRARY}")

  if(NOT TARGET opentims::opentims_cpp)
    # NOTE: This manual-import path does not know whether the library was built
    # with OPENTIMS_LINK_SQLITE_STATICALLY (which requires the caller to supply
    # sqlite3 headers via Opentims_INCLUDE_DIR and link sqlite3).
    # INTERFACE_COMPILE_DEFINITIONS is intentionally set to an empty list so
    # that downstream consumers (e.g. cmake_findExternalLibs.cmake) can
    # distinguish "property not set at all" (CONFIG target) from "property set
    # but empty" (manual target) when deciding whether to inject sqlite3.
    # If OPENTIMS_LINK_SQLITE_STATICALLY is required, override this property or
    # pass -DOpentims_ROOT_DIR pointing to an installation with a config package.
    add_library(opentims::opentims_cpp UNKNOWN IMPORTED)
    set_target_properties(opentims::opentims_cpp PROPERTIES
      IMPORTED_LOCATION             "${Opentims_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${Opentims_INCLUDE_DIR}"
      INTERFACE_COMPILE_DEFINITIONS ""
    )
  endif()
endif()

mark_as_advanced(Opentims_INCLUDE_DIR Opentims_LIBRARY)
