file(READ "${INPUT}" default_ini)
string(REGEX MATCH "name=\"stop_report_after_feature\" value=\"-1\"" default_report_limit "${default_ini}")

if(NOT default_report_limit)
  message(FATAL_ERROR "Expected OpenSwathWorkflow to report all features by default (stop_report_after_feature = -1)")
endif()
