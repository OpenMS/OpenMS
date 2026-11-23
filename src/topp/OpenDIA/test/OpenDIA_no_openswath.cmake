# OpenDIA uses OpenMS readers and low-level classes only, no OpenSWATH code (design rule). Fails if the OpenDIA core
# library sources (include/, source/) reference OpenSWATH headers or the OpenSwath:: namespace. The command-line tool
# may use OpenSWATH only at its explicitly marked edge (OpenDIA -decoys predict).
file(GLOB core "${SOURCE_DIR}/include/*.h" "${SOURCE_DIR}/source/*.cpp")
set(bad "")
foreach(f IN LISTS core)
  file(STRINGS "${f}" hits REGEX "OPENSWATH|OpenSwath::|PeakMapExtractor|TransitionPQPFile|DIAHelper")
  if(hits)
    list(APPEND bad "${f}: ${hits}")
  endif()
endforeach()
if(bad)
  string(REPLACE ";" "\n" bad "${bad}")
  message(FATAL_ERROR "OpenSWATH code in the OpenDIA core:\n${bad}")
endif()
message(STATUS "OpenDIA core: no OpenSWATH references in ${SOURCE_DIR}/include and /source")
