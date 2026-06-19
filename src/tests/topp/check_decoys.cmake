# Usage: cmake -P check_decoys.cmake <idXML-file> <present|absent>
#
# Counts decoy accessions (accession="DECOY_...) in an idXML file — these occur
# in both <ProteinHit> and <PeptideEvidence> entries — and asserts they are
# present or absent. Used to pin ProSE's decoy-reporting contract without giant
# golden reference files:
#   - protein FDR finalizes a result  -> decoys absent
#   - protein FDR off (intermediate)  -> decoys present (target+decoy retained)

if(NOT DEFINED CMAKE_ARGV3 OR NOT DEFINED CMAKE_ARGV4)
  message(FATAL_ERROR "Usage: cmake -P check_decoys.cmake <file> <present|absent>")
endif()

set(file "${CMAKE_ARGV3}")
set(mode "${CMAKE_ARGV4}")

if(NOT EXISTS "${file}")
  message(FATAL_ERROR "decoy-check input not found: ${file}")
endif()

# Use file(READ) + string(REGEX MATCHALL) rather than file(STRINGS ... REGEX):
# file(STRINGS) silently drops long lines and miscounts on idXML.
file(READ "${file}" content)
string(REGEX MATCHALL "accession=\"DECOY_" decoys "${content}")
list(LENGTH decoys count)

get_filename_component(basename "${file}" NAME)
message(STATUS "DECOY_ accessions in ${basename}: ${count} (expect ${mode})")

if(mode STREQUAL "present")
  if(count LESS 1)
    message(FATAL_ERROR
      "Decoy-reporting regression: expected decoys RETAINED in ${basename} (intermediate / "
      "protein FDR off) but found none. PSM-level FDR must not strip decoys.")
  endif()
elseif(mode STREQUAL "absent")
  if(count GREATER 0)
    message(FATAL_ERROR
      "Decoy-reporting regression: expected NO decoys in ${basename} (protein FDR finalized "
      "the result) but found ${count}.")
  endif()
else()
  message(FATAL_ERROR "Unknown mode '${mode}' (use present|absent)")
endif()
