# Usage: cmake -P check_psm_floor.cmake <idXML-file> <floor>
#
# Counts <PeptideHit ...> occurrences in an idXML file and fails if the count
# is below <floor>. Intended as a regression tripwire on PSM yields from the
# Bruker DDA opentims integration tests (1% FDR outputs). Floor-only: legitimate
# improvements that raise the count should NOT fail; only drops do.

if(NOT DEFINED CMAKE_ARGV3 OR NOT DEFINED CMAKE_ARGV4)
  message(FATAL_ERROR "Usage: cmake -P check_psm_floor.cmake <file> <floor>")
endif()

set(file "${CMAKE_ARGV3}")
set(floor "${CMAKE_ARGV4}")

if(NOT EXISTS "${file}")
  message(FATAL_ERROR "PSM-count input not found: ${file}")
endif()

# Use file(READ) + string(REGEX MATCHALL) rather than file(STRINGS ... REGEX).
# file(STRINGS) silently drops long lines/entries and produces wrong counts on
# idXML (observed: 186 vs the true 4537) — MATCHALL over the full file content
# is the reliable path.
file(READ "${file}" content)
string(REGEX MATCHALL "<PeptideHit " hits "${content}")
list(LENGTH hits count)

get_filename_component(basename "${file}" NAME)
message(STATUS "1%FDR PSMs in ${basename}: ${count} (floor=${floor})")

if(count LESS floor)
  message(FATAL_ERROR
    "PSM-count regression: ${basename} has ${count} PeptideHits, below floor ${floor}. "
    "Either a real regression or (if algorithm changes legitimately shifted yield) "
    "update the floor.")
endif()
