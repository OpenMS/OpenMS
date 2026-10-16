# Remove old versioned OpenMS shared libraries from the pyopenms build dir
# so that a version bump in a reused build tree doesn't leave stale copies.
if(NOT DIR)
  return()
endif()
file(GLOB _old_openms "${DIR}/libOpenMS.so.*" "${DIR}/libOpenMS.*.dylib")
file(GLOB _old_swath  "${DIR}/libOpenSwathAlgo.so.*" "${DIR}/libOpenSwathAlgo.*.dylib")
set(_old ${_old_openms} ${_old_swath})
if(_old)
  file(REMOVE ${_old})
endif()
