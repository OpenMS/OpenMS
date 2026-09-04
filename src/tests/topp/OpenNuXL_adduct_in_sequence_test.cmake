if (NOT DEFINED INPUT_FILE)
  message(FATAL_ERROR "INPUT_FILE is required")
endif()

file(READ "${INPUT_FILE}" contents)

# The first localized fixture hit carries U-H2O at zero-based residue 8, written as a named modification.
string(FIND "${contents}" "sequence=\"AEADNLDDK(NuXL:U-H2O1)K\"" localized_sequence_position)
if (localized_sequence_position EQUAL -1)
  message(FATAL_ERROR
    "Localized U-H2O adduct was not written as K(NuXL:U-H2O1) at residue 8")
endif()

# Its definition (name, origin K, formula) travels in the search parameters.
string(FIND "${contents}" "name=\"modification_definitions\"" definitions_position)
if (definitions_position EQUAL -1)
  message(FATAL_ERROR "SearchParameters carry no modification_definitions UserParam")
endif()
string(FIND "${contents}" "1|NuXL:U-H2O1|NuXL:U-H2O1 (K)||K|none|C9H11N2O8P1|" definition_position)
if (definition_position EQUAL -1)
  message(FATAL_ERROR
    "modification_definitions does not define NuXL:U-H2O1 on K with formula C9H11N2O8P1")
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
  "Validated localized K(NuXL:U-H2O1) with its definition and ${unlocalized_count} bare unlocalized PSMs")
