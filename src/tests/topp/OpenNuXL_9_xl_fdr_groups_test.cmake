if (NOT DEFINED INPUT_FILE)
  message(FATAL_ERROR "INPUT_FILE is required")
endif()

file(READ "${INPUT_FILE}" contents)

# -report:xl_fdr_groups on writes the cross-links of -out_xls. The test data have short and other cross-links but no
# decoys among them, and FalseDiscoveryRate would give every target of a group without decoys q-value 0. So the
# cross-links keep one FDR and are not annotated with NuXL:fdr_group. (The grouped q-values are covered by
# NuXLFDR_test.)
string(REGEX MATCHALL "<PeptideHit " hits "${contents}")
list(LENGTH hits hit_count)
if (hit_count EQUAL 0)
  message(FATAL_ERROR "No cross-links in -out_xls")
endif()

string(FIND "${contents}" "NuXL:fdr_group" group)
if (NOT group EQUAL -1)
  message(FATAL_ERROR "Cross-links annotated with NuXL:fdr_group although the test data have no cross-link decoys")
endif()
