file(READ "${INPUT}" actual_xml)
file(READ "${REFERENCE}" reference_xml)
string(REGEX MATCHALL "<feature id=" actual_features "${actual_xml}")
string(REGEX MATCHALL "<feature id=" reference_features "${reference_xml}")
list(LENGTH actual_features actual_count)
list(LENGTH reference_features reference_count)

if(actual_count LESS_EQUAL reference_count)
  message(FATAL_ERROR
    "OpenSwathWorkflow default reported ${actual_count} candidates, but the capped reference has ${reference_count}")
endif()
