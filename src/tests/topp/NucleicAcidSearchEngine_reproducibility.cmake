# Run a real tied-localization spectrum with different OpenMP schedules.
file(MAKE_DIRECTORY "${OUTPUT_DIR}")
set(previous "")
foreach(run RANGE 1 3)
  if(run EQUAL 1)
    set(threads 1)
  else()
    set(threads 4)
  endif()
  set(prefix "${OUTPUT_DIR}/NASE_ties_${run}")
  execute_process(COMMAND "${NASE}"
    -test -ini "${DATA_DIR}/NucleicAcidSearchEngine_ties.ini"
    -in "${DATA_DIR}/NucleicAcidSearchEngine_ties.mzML"
    -database "${DATA_DIR}/NucleicAcidSearchEngine_ties.fasta"
    -out "${prefix}.mzTab" -id_out "${prefix}.idXML"
    -bedrmod_out "${prefix}.bed" -threads ${threads}
    -bedrmod_chebi_mapping "${DATA_DIR}/NucleicAcidSearchEngine_ties.csv"
    RESULT_VARIABLE result OUTPUT_VARIABLE log ERROR_VARIABLE errors)
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "NASE failed: ${log} ${errors}")
  endif()
  file(READ "${prefix}.idXML" ids)
  # Isolate the known tied spectrum so other spectra cannot satisfy this check.
  string(FIND "${ids}" "spectrum_reference=\"controllerType=0 controllerNumber=1 scan=43600\"" scan_start)
  if(scan_start EQUAL -1)
    message(FATAL_ERROR "Missing regression spectrum")
  endif()
  string(SUBSTRING "${ids}" ${scan_start} -1 scan_tail)
  string(FIND "${scan_tail}" "</PeptideIdentification>" scan_end)
  if(scan_end EQUAL -1)
    message(FATAL_ERROR "Malformed idXML output")
  endif()
  string(SUBSTRING "${scan_tail}" 0 ${scan_end} tied_scan)
  # Do not merely verify that an arbitrary single winner repeats: both equally
  # scoring modification placements must survive candidate collection.
  foreach(sequence "CCCCCCG[mU?]U[ac4C]CUCCCGp" "CCCCCCG[mU?]UC[ac4C]UCCCGp")
    string(FIND "${tied_scan}" "value=\"${sequence}\"" found)
    if(found EQUAL -1)
      message(FATAL_ERROR "Lost tied candidate ${sequence} with ${threads} threads")
    endif()
  endforeach()
  file(STRINGS "${prefix}.bed" data_rows REGEX "^[^#]")
  list(LENGTH data_rows row_count)
  if(row_count EQUAL 0)
    message(FATAL_ERROR "No BED rows were exported")
  endif()
  if(NOT previous STREQUAL "")
    execute_process(COMMAND "${CMAKE_COMMAND}" -E compare_files "${previous}" "${prefix}.bed"
                    RESULT_VARIABLE different)
    if(NOT different EQUAL 0)
      message(FATAL_ERROR "BED output changed between thread counts or repeated runs")
    endif()
  endif()
  set(previous "${prefix}.bed")
endforeach()
