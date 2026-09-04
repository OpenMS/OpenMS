# Usage: cmake -DINPUT_FILE=<idXML> -DREFERENCE_FILE=<idXML> -P compare_sequence_attributes.cmake
#
# Requires both files to carry the same sequence="..." attributes in the same order. This is
# the round-trip check for a peptidoform spelling: a reader that degrades a modification
# (drops it, renames it, turns a named modification into a mass bracket) changes the
# attribute even when every mass and score still compares clean. FuzzyDiff cannot do this
# job: it compares lines positionally, and UserParam order is process-dependent.
foreach(required INPUT_FILE REFERENCE_FILE)
  if (NOT DEFINED ${required})
    message(FATAL_ERROR "${required} is required")
  endif()
endforeach()

file(READ "${INPUT_FILE}" input_contents)
file(READ "${REFERENCE_FILE}" reference_contents)
string(REGEX MATCHALL "sequence=\"[^\"]*\"" input_sequences "${input_contents}")
string(REGEX MATCHALL "sequence=\"[^\"]*\"" reference_sequences "${reference_contents}")
list(LENGTH input_sequences input_count)
list(LENGTH reference_sequences reference_count)

if (reference_count EQUAL 0)
  message(FATAL_ERROR "${REFERENCE_FILE} carries no sequence attributes")
endif()

if (NOT input_count EQUAL reference_count)
  message(FATAL_ERROR
    "${INPUT_FILE} carries ${input_count} sequence attributes, ${REFERENCE_FILE} carries ${reference_count}")
endif()

math(EXPR last "${reference_count} - 1")
foreach(i RANGE ${last})
  list(GET input_sequences ${i} input_sequence)
  list(GET reference_sequences ${i} reference_sequence)
  if (NOT input_sequence STREQUAL reference_sequence)
    message(FATAL_ERROR
      "sequence attribute ${i} differs:\n  input:     ${input_sequence}\n  reference: ${reference_sequence}")
  endif()
endforeach()

message(STATUS "${input_count} sequence attributes identical to the reference")
