# Remove stale bundled libraries, retaining current files and unversioned aliases.
if(NOT DIR OR NOT CURRENT_OPENMS OR NOT CURRENT_SWATH)
  message(FATAL_ERROR "DIR and current library filenames are required")
endif()
file(GLOB _old "${DIR}/libOpenMS.so" "${DIR}/libOpenMS.so.*"
  "${DIR}/libOpenMS.dylib" "${DIR}/libOpenMS.*.dylib"
  "${DIR}/libOpenSwathAlgo.so" "${DIR}/libOpenSwathAlgo.so.*"
  "${DIR}/libOpenSwathAlgo.dylib" "${DIR}/libOpenSwathAlgo.*.dylib")
foreach(_file IN LISTS _old)
  get_filename_component(_name "${_file}" NAME)
  if(_name STREQUAL CURRENT_OPENMS OR _name STREQUAL CURRENT_SWATH)
    continue()
  endif()
  if(IS_SYMLINK "${_file}" AND _name MATCHES "^(libOpenMS|libOpenSwathAlgo)([.]so|[.]dylib)$")
    continue()
  endif()
  file(REMOVE "${_file}")
endforeach()
