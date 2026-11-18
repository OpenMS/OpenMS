if (NOT DEFINED INPUT_FILE)
  message(FATAL_ERROR "INPUT_FILE is required")
endif()

file(READ "${INPUT_FILE}" contents)

# The search of TOPP_OpenNuXL_1 with -NuXL:localization consistency (see OpenNuXL_1_output.idXML for 'fragments'):
# - FDQKIIDK: the shifted b ions start at K3 and the shifted y ions end there. 'fragments' rejects every ladder site
#   and puts U-H2O on the C-terminal K7, where the shifted lysine immonium ions tie.
# - LAGPMQAWGQPTFEGTRPTGR: 'fragments' follows a weak shifted y3 to T18, the ladder points to G8.
foreach(localized IN ITEMS "FDQK(NuXL:U-H2O1)IIDK" "LAGPMQAWG(NuXL:U)QPTFEGTRPTGR")
  string(FIND "${contents}" "sequence=\"${localized}\"" position)
  if (position EQUAL -1)
    message(FATAL_ERROR "Cross-link not localized as ${localized}")
  endif()
endforeach()

# Every cross-link of the fixture has fragment evidence: localized by consistency, with its range and delta
string(REGEX MATCHALL "name=\"NuXL:isXL\" value=\"1\"" xls "${contents}")
string(REGEX MATCHALL "name=\"NuXL:localization_method\" value=\"consistency\"" methods "${contents}")
string(REGEX MATCHALL "name=\"NuXL:localization_range\" value=\"[0-9]+-[0-9]+\"" ranges "${contents}")
string(REGEX MATCHALL "name=\"NuXL:localization_delta\"" deltas "${contents}")
list(LENGTH xls xl_count)
list(LENGTH methods method_count)
list(LENGTH ranges range_count)
list(LENGTH deltas delta_count)
if (xl_count EQUAL 0 OR NOT method_count EQUAL xl_count OR NOT range_count EQUAL xl_count OR NOT delta_count EQUAL xl_count)
  message(FATAL_ERROR
    "Expected each of the ${xl_count} cross-links localized by consistency with range and delta, found "
    "${method_count} methods, ${range_count} ranges, ${delta_count} deltas")
endif()

message(STATUS "Validated consistency localization of ${xl_count} cross-links")
