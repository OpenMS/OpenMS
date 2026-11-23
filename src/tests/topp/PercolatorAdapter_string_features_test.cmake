# Stores the search-engine-specific features of the PercolatorAdapter_1 fixture as string meta
# values, as adapters do that read search engine scores from text (e.g. SageAdapter).
if (NOT DEFINED INPUT_FILE OR NOT DEFINED OUTPUT_FILE)
  message(FATAL_ERROR "INPUT_FILE and OUTPUT_FILE are required")
endif()

file(READ "${INPUT_FILE}" contents)
# the fixture's 'extra_features'
foreach(feature IN ITEMS "COMET:deltCn" "COMET:deltLCn" "COMET:lnExpect" "MS:1002252" "MS:1002255"
                         "COMET:lnNumSP" "COMET:lnRankSP" "COMET:IonFrac")
  string(REPLACE "<UserParam type=\"float\" name=\"${feature}\""
                 "<UserParam type=\"string\" name=\"${feature}\"" retyped "${contents}")
  if (retyped STREQUAL contents)
    message(FATAL_ERROR "Feature '${feature}' is not a float meta value in ${INPUT_FILE}")
  endif()
  set(contents "${retyped}")
endforeach()
file(WRITE "${OUTPUT_FILE}" "${contents}")
