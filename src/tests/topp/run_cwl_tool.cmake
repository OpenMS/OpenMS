# Run a generated CWL description through a CWL runner and require it to succeed.
#
# A CWL runner does not pass the parameters on the command line: the generated descriptions write
# them to a JSON file and invoke the tool as '<Tool> -ini cwl_inputs.json' (see ParamCWLFile.cpp).
# Whether a tool can read back what a runner wrote is therefore only covered by actually running
# one -- the 'cwltool --validate' in cmake/cwl_generation.cmake checks that the description itself
# is well-formed, which a description the tool cannot consume passes just as happily.
#
# Required: -DCWLTOOL=<cwltool> -DTOOL_BIN_PATH=<directory holding the TOPP tools>
#           -DCWL=<description> -DJOB=<job file> -DOUT_DIR=<output directory>
#           -DEXPECTED_OUTPUT=<file that has to exist afterwards>
foreach(required CWLTOOL TOOL_BIN_PATH CWL JOB OUT_DIR EXPECTED_OUTPUT)
  if(NOT DEFINED ${required})
    message(FATAL_ERROR "run_cwl_tool.cmake: -D${required}=... is required")
  endif()
endforeach()

# The generated descriptions use the bare tool name as their baseCommand, so the runner can only
# find the tool via PATH.
if(WIN32)
  set(path_separator ";")
else()
  set(path_separator ":")
endif()
set(ENV{PATH} "${TOOL_BIN_PATH}${path_separator}$ENV{PATH}")

file(REMOVE "${EXPECTED_OUTPUT}")
file(MAKE_DIRECTORY "${OUT_DIR}")

execute_process(COMMAND ${CWLTOOL} --no-container --outdir ${OUT_DIR} ${CWL} ${JOB}
                RESULT_VARIABLE exit_code OUTPUT_VARIABLE out ERROR_VARIABLE err)

if(NOT exit_code EQUAL 0)
  message(FATAL_ERROR
    "'${CWLTOOL} ${CWL} ${JOB}' failed with '${exit_code}'.\n"
    "--- stdout ---\n${out}\n--- stderr ---\n${err}")
endif()
if(NOT EXISTS "${EXPECTED_OUTPUT}")
  message(FATAL_ERROR
    "The runner reported success, but '${EXPECTED_OUTPUT}' was not produced.\n"
    "--- stdout ---\n${out}\n--- stderr ---\n${err}")
endif()
message(STATUS "CWL run produced ${EXPECTED_OUTPUT}")
