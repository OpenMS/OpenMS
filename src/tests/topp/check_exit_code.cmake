# Run a command and require an EXACT exit code.
#
# ctest's WILL_FAIL only distinguishes "zero" from "non-zero", and PASS_REGULAR_EXPRESSION
# ignores the exit code entirely -- neither can tell a clean error return from a crash, which
# is exactly the distinction some regression tests need (e.g. a tool that must return
# ILLEGAL_PARAMETERS rather than abort). Use this when the specific code is the assertion.
#
# Required: -DCOMMAND=<;-separated argv> -DEXPECTED_CODE=<int>
# Optional: -DEXPECTED_OUTPUT_REGEX=<regex> -- stdout+stderr must also match, so a test can
# insist on the specific failure it is about rather than on any error return.
if(NOT DEFINED COMMAND)
  message(FATAL_ERROR "check_exit_code.cmake: -DCOMMAND=... is required")
endif()
if(NOT DEFINED EXPECTED_CODE)
  message(FATAL_ERROR "check_exit_code.cmake: -DEXPECTED_CODE=... is required")
endif()

execute_process(COMMAND ${COMMAND} RESULT_VARIABLE actual_code OUTPUT_VARIABLE out ERROR_VARIABLE err)

if(NOT actual_code STREQUAL EXPECTED_CODE)
  # A signal name (e.g. "SIGABRT") rather than a number means the process died, which is what
  # this check most often exists to catch.
  message(FATAL_ERROR
    "Expected exit code ${EXPECTED_CODE}, got '${actual_code}'.\n"
    "--- stdout ---\n${out}\n--- stderr ---\n${err}")
endif()
if(DEFINED EXPECTED_OUTPUT_REGEX AND NOT "${out}${err}" MATCHES "${EXPECTED_OUTPUT_REGEX}")
  message(FATAL_ERROR
    "Output does not match '${EXPECTED_OUTPUT_REGEX}'.\n"
    "--- stdout ---\n${out}\n--- stderr ---\n${err}")
endif()
message(STATUS "Exit code ${actual_code} as expected")
