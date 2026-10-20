# Usage: cmake "-DCOMMAND=<;-separated argv>" -DLOG=<file> -P run_logged.cmake
#
# Runs the command, writes its standard output and error to LOG (for checks of what it logged) and fails if the command
# does not return 0.
if (NOT DEFINED COMMAND OR NOT DEFINED LOG)
  message(FATAL_ERROR "COMMAND and LOG are required")
endif()
execute_process(COMMAND ${COMMAND} RESULT_VARIABLE code OUTPUT_VARIABLE out ERROR_VARIABLE err)
file(WRITE "${LOG}" "${out}${err}")
if (NOT code EQUAL 0)
  message(FATAL_ERROR "Exit code '${code}' (output in ${LOG})")
endif()
