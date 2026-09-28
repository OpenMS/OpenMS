if (NOT DEFINED OPERATION)
  message(FATAL_ERROR "OPERATION is required")
endif()

if (OPERATION STREQUAL "prepare")
  if (NOT DEFINED INPUT_FILE OR NOT DEFINED OUTPUT_FILE OR NOT DEFINED FEATURE_DIR OR NOT DEFINED CAPTURE_FILE)
    message(FATAL_ERROR "INPUT_FILE, OUTPUT_FILE, FEATURE_DIR and CAPTURE_FILE are required for prepare")
  endif()

  file(READ "${INPUT_FILE}" contents)
  set(hit_opening [=[<PeptideHit score="0" sequence="SHC(Carbamidomethyl)IAEVEK" charge="3" aa_before="K" aa_after="D" protein_refs="PH_0" >]=])
  string(FIND "${contents}" "${hit_opening}" hit_position)
  if (hit_position EQUAL -1)
    message(FATAL_ERROR "Could not find the fixture's first PeptideHit")
  endif()

  # A small but unambiguous offset keeps the target near its observed precursor
  # while forcing FeatureFinderIdentificationAlgorithm's CalcMass branch.
  set(stamped_hit "${hit_opening}\n\t\t\t\t<UserParam type=\"float\" name=\"CalcMass\" value=\"358.184682617188\"/>")
  string(REPLACE "${hit_opening}" "${stamped_hit}" contents "${contents}")
  file(WRITE "${OUTPUT_FILE}" "${contents}")
  file(REMOVE_RECURSE "${FEATURE_DIR}")
  file(REMOVE "${CAPTURE_FILE}")
elseif (OPERATION STREQUAL "run")
  foreach(required_var TOOL MZML_FILE ID_FILE DESIGN_FILE FASTA_FILE FEATURE_DIR CAPTURE_FILE)
    if (NOT DEFINED ${required_var})
      message(FATAL_ERROR "${required_var} is required for run")
    endif()
  endforeach()

  execute_process(
    COMMAND "${TOOL}"
            -in "${MZML_FILE}"
            -ids "${ID_FILE}"
            -design "${DESIGN_FILE}"
            -fasta "${FASTA_FILE}"
            -targeted_only true
            -mass_recalibration false
            -threads 1
            -proteinFDR 0.3
            -test
            -debug 1
            -feat_dir "${FEATURE_DIR}"
            -detect_only
    RESULT_VARIABLE tool_result
    OUTPUT_VARIABLE tool_stdout
    ERROR_VARIABLE tool_stderr)
  file(WRITE "${CAPTURE_FILE}" "${tool_stdout}${tool_stderr}")
  if (NOT tool_result EQUAL 0)
    message(FATAL_ERROR "ProteomicsLFQ failed with exit status ${tool_result}:\n${tool_stdout}${tool_stderr}")
  endif()
elseif (OPERATION STREQUAL "check")
  if (NOT DEFINED CAPTURE_FILE OR NOT EXISTS "${CAPTURE_FILE}")
    message(FATAL_ERROR "The ProteomicsLFQ output capture was not produced")
  endif()
  file(READ "${CAPTURE_FILE}" contents)
  string(FIND "${contents}" "Adding unknown modification" branch_message)
  if (branch_message EQUAL -1)
    message(FATAL_ERROR "FeatureFinderIdentification did not consume the input CalcMass")
  endif()
else()
  message(FATAL_ERROR "Unknown OPERATION: ${OPERATION}")
endif()
