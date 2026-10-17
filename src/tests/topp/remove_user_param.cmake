# Usage: cmake -DINPUT_FILE=<xml> -DOUTPUT_FILE=<xml> -DNAME=<UserParam name> -P remove_user_param.cmake
#
# Writes a copy of an OpenMS XML file with every single-line <UserParam ... name="NAME" .../>
# element deleted. It derives a "same file minus one element" input at test time, so a
# compatibility property can be asserted without a second copy of a large fixture.
foreach(required INPUT_FILE OUTPUT_FILE NAME)
  if (NOT DEFINED ${required})
    message(FATAL_ERROR "${required} is required")
  endif()
endforeach()

file(READ "${INPUT_FILE}" contents)
string(REGEX MATCHALL "[^\n]*<UserParam [^\n]*name=\"${NAME}\"[^\n]*\n" removed "${contents}")
list(LENGTH removed removed_count)
if (removed_count EQUAL 0)
  message(FATAL_ERROR "${INPUT_FILE} carries no UserParam named ${NAME}")
endif()
string(REGEX REPLACE "[^\n]*<UserParam [^\n]*name=\"${NAME}\"[^\n]*\n" "" contents "${contents}")
file(WRITE "${OUTPUT_FILE}" "${contents}")

message(STATUS "Removed ${removed_count} UserParam(s) named ${NAME}")
