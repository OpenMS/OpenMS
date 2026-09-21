# Usage: cmake -DINPUT_FILE=<file> [-DEXPECTED=<literal;...>] [-DFORBIDDEN=<literal;...>] -P check_file_contains.cmake
#
# Requires every EXPECTED literal to occur in the file and every FORBIDDEN literal to be
# absent. Literals rather than regexes, so bracket-heavy strings such as ProForma
# peptidoforms need no escaping. A literal cannot contain ';' (the list separator).
if (NOT DEFINED INPUT_FILE)
  message(FATAL_ERROR "INPUT_FILE is required")
endif()
if (NOT DEFINED EXPECTED AND NOT DEFINED FORBIDDEN)
  message(FATAL_ERROR "EXPECTED or FORBIDDEN is required")
endif()

file(READ "${INPUT_FILE}" contents)

foreach(literal IN LISTS EXPECTED)
  string(FIND "${contents}" "${literal}" position)
  if (position EQUAL -1)
    message(FATAL_ERROR "${INPUT_FILE} does not contain: ${literal}")
  endif()
endforeach()

foreach(literal IN LISTS FORBIDDEN)
  string(FIND "${contents}" "${literal}" position)
  if (NOT position EQUAL -1)
    message(FATAL_ERROR "${INPUT_FILE} must not contain: ${literal}")
  endif()
endforeach()

message(STATUS "${INPUT_FILE} passed the content checks")
