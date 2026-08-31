if (NOT DEFINED INPUT_FILE)
  message(FATAL_ERROR "INPUT_FILE is required")
endif()

file(READ "${INPUT_FILE}" contents)

# The first localized fixture hit carries U-H2O at zero-based residue 8.
string(REGEX MATCH
  "sequence=\"AEADNLDDK\\[\\+306\\.0253048409[0-9]*\\]K\""
  localized_sequence "${contents}")
if (localized_sequence STREQUAL "")
  message(FATAL_ERROR
    "Localized U-H2O adduct was not written as a +306.0253048409 Da bracket at residue 8")
endif()

# These two fixture hits have best_localization_position == -1 and must stay bare.
foreach(unlocalized_sequence IN ITEMS VEDVHPISERPQYFSGDGK SGYYAQNIGESSLRTIHLAQLR)
  string(FIND "${contents}" "sequence=\"${unlocalized_sequence}\"" sequence_position)
  if (sequence_position EQUAL -1)
    message(FATAL_ERROR
      "Unlocalized sequence ${unlocalized_sequence} was changed or given a fabricated localization")
  endif()
endforeach()

string(REGEX MATCHALL
  "NuXL:best_localization_position\" value=\"-1\""
  unlocalized_hits "${contents}")
list(LENGTH unlocalized_hits unlocalized_count)
if (NOT unlocalized_count EQUAL 2)
  message(FATAL_ERROR
    "Expected two unlocalized fixture PSMs, found ${unlocalized_count}")
endif()

message(STATUS
  "Validated localized +306.0253048409 Da bracket and ${unlocalized_count} bare unlocalized PSMs")
