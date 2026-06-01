find_path(SIMDe_INCLUDE_DIR simde/simde-align.h)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(
  SIMDe
  REQUIRED_VARS SIMDe_INCLUDE_DIR
)

if(NOT TARGET SIMDe)
  add_library(SIMDe INTERFACE)
  target_include_directories(SIMDe SYSTEM INTERFACE ${SIMDe_INCLUDE_DIR})
endif()

mark_as_advanced(SIMDe_INCLUDE_DIR)
