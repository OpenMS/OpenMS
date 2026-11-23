# TEMPORARY (scratch branch only): print the Measure_error FragmentArray lines of the
# generated OpenPepXL mzIdentML so the golden file can be regenerated from a real run.
if(NOT EXISTS "${MZID}")
  message(FATAL_ERROR "missing ${MZID}")
endif()
file(STRINGS "${MZID}" lines REGEX "measure_ref=\"Measure_error\"")
list(LENGTH lines n)
message("MZID_DUMP_BEGIN ${n}")
foreach(l IN LISTS lines)
  message("MZID_DUMP|${l}")
endforeach()
message("MZID_DUMP_END")
message(FATAL_ERROR "intentional failure so that ctest prints the output above")
