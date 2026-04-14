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
# Hints:
#  Set Opentims_ROOT_DIR to the root of a custom installation.

# Prefer a CMake config package produced by opentims's own install rules.
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
  get_target_property(Opentims_LIBRARIES   opentims::opentims_cpp IMPORTED_LOCATION)
  return()
endif()

# Fall back to manual search when no config package is present.
find_path(Opentims_INCLUDE_DIR
  NAMES opentims++/opentims.h
  HINTS
    "${Opentims_ROOT_DIR}/include"
    "${Opentims_ROOT_DIR}"
  DOC "opentims include directory (parent of opentims++/)"
)

find_library(Opentims_LIBRARY
  NAMES opentims_cpp
  HINTS
    "${Opentims_ROOT_DIR}/lib"
    "${Opentims_ROOT_DIR}"
  DOC "opentims static library"
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Opentims
  REQUIRED_VARS Opentims_LIBRARY Opentims_INCLUDE_DIR
)

if(Opentims_FOUND)
  set(Opentims_LIBRARIES "${Opentims_LIBRARY}")

  if(NOT TARGET opentims::opentims_cpp)
    add_library(opentims::opentims_cpp STATIC IMPORTED)
    set_target_properties(opentims::opentims_cpp PROPERTIES
      IMPORTED_LOCATION             "${Opentims_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${Opentims_INCLUDE_DIR}"
    )
  endif()
endif()

mark_as_advanced(Opentims_INCLUDE_DIR Opentims_LIBRARY)
