if (NOT DEFINED INPUT_FILE)
  message(FATAL_ERROR "INPUT_FILE is required")
endif()

file(READ "${INPUT_FILE}" contents)

# NuXL:rescoring native annotates the delta score and the ion likelihood features on every hit and passes them to
# Percolator. The test spectra are too few to train the ion likelihood, which then reports 0 instead of omitting
# the features.
set(features NuXL:delta_score NuXL:delta_score_class NuXL:ion_llr NuXL:ion_llr_shifted NuXL:ion_explained NuXL:ion_topk_observed)

string(REGEX MATCH "name=\"extra_features\" value=\"[^\"]*\"" extra_features "${contents}")
if (extra_features STREQUAL "")
  message(FATAL_ERROR "No extra_features in the search parameters")
endif()

string(REGEX MATCHALL "<PeptideHit " hits "${contents}")
list(LENGTH hits hit_count)
if (hit_count EQUAL 0)
  message(FATAL_ERROR "No PeptideHits")
endif()

foreach(feature IN LISTS features)
  string(FIND "${extra_features}" "${feature}," in_list)
  string(FIND "${extra_features}" "${feature}\"" at_end)
  if (in_list EQUAL -1 AND at_end EQUAL -1)
    message(FATAL_ERROR "${feature} is not passed to Percolator (extra_features)")
  endif()
  string(REGEX MATCHALL "name=\"${feature}\" value=" annotated "${contents}")
  list(LENGTH annotated annotated_count)
  if (NOT annotated_count EQUAL hit_count)
    message(FATAL_ERROR "${feature} is annotated on ${annotated_count} of ${hit_count} PeptideHits")
  endif()
endforeach()

# the matched ions are kept only until the ion likelihood is trained
string(FIND "${contents}" "NuXL:ion_evidence" evidence)
if (NOT evidence EQUAL -1)
  message(FATAL_ERROR "NuXL:ion_evidence is left in the output")
endif()
