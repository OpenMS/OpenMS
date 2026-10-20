# Checks the spectrum columns that ProSE's -qpx_spectra adds to a QPX psm.parquet, from its text dump
# (ParquetDiff -out_tsv). Usage:
#
#   cmake -DTSV=<dump.tsv> -DMZML=<the searched mzML> -DNATIVE_ID_PREFIX=<e.g. spectrum=>
#         [-DPREDICTIONS=ON|OFF] [-DYAML=<-summary_out of the run> -DPEPTIDE=<peptidoform> -DCHARGE=<z>
#          -DREFERENCE=<tsv: nce, ion, mz, intensity> -DOUT_PREDICTED=<file> -DOUT_REFERENCE=<file>]
#         -P check_prose_qpx_spectra.cmake
#
# - The five columns exist.
# - Every row has its spectrum's peaks: peak_mz_array and peak_intensity_array have as many values as the spectrum
#   with the row's native ID (NATIVE_ID_PREFIX + scan number) has in the mzML (its defaultArrayLength).
# - PREDICTIONS ON (default): every row of precursor charge 1 to 6 has a non-empty prediction, with as many ions as
#   m/z and intensities. OFF (OpenMS without ONNX Runtime): every prediction column is null.
# - With PEPTIDE: the predicted ions of the first row of that peptidoform and CHARGE are written to OUT_PREDICTED, and
#   the rows of REFERENCE for the NCE that the summary (YAML) reports for the file to OUT_REFERENCE, as
#   "ion<TAB>m/z<TAB>intensity" lines, for FuzzyDiff to compare.
foreach(_required TSV MZML NATIVE_ID_PREFIX)
  if (NOT DEFINED ${_required})
    message(FATAL_ERROR "${_required} is required")
  endif()
endforeach()
if (NOT DEFINED PREDICTIONS)
  set(PREDICTIONS ON)
endif()

# Peak counts of the spectra, by native ID.
file(STRINGS "${MZML}" _spectra REGEX "<spectrum [^>]*defaultArrayLength=")
foreach(_s IN LISTS _spectra)
  if (NOT _s MATCHES "id=\"([^\"]*)\"")
    continue()
  endif()
  set(_id "${CMAKE_MATCH_1}")
  if (NOT _s MATCHES "defaultArrayLength=\"([0-9]+)\"")
    continue()
  endif()
  string(MAKE_C_IDENTIFIER "${_id}" _key)
  set(_peaks_${_key} "${CMAKE_MATCH_1}")
endforeach()

# The dump, one list element per line. Brackets and semicolons would change how CMake splits lists, so the list
# brackets of the dump become '<' '>' first.
file(READ "${TSV}" _dump)
string(REPLACE ";" "," _dump "${_dump}")
string(REPLACE "[" "<" _dump "${_dump}")
string(REPLACE "]" ">" _dump "${_dump}")
string(REPLACE "\n" ";" _lines "${_dump}")
list(POP_FRONT _lines _header)
string(REPLACE "\t" ";" _header "${_header}")

foreach(_column peptidoform charge scan peak_mz_array peak_intensity_array predicted_ion_type_array predicted_mz_array predicted_intensity_array)
  list(FIND _header "${_column}" _col_${_column})
  if (_col_${_column} EQUAL -1)
    message(FATAL_ERROR "${TSV}: no column '${_column}'")
  endif()
endforeach()

# The number of values of a dumped list "<a, b, c>"; -1 for null.
function(list_length value out)
  if (value STREQUAL "null")
    set(${out} -1 PARENT_SCOPE)
    return()
  endif()
  string(REGEX REPLACE "^<(.*)>$" "\\1" _inner "${value}")
  if (_inner STREQUAL "")
    set(${out} 0 PARENT_SCOPE)
    return()
  endif()
  string(REPLACE ", " ";" _items "${_inner}")
  list(LENGTH _items _n)
  set(${out} ${_n} PARENT_SCOPE)
endfunction()

set(_rows 0)
set(_predicted_rows 0)
set(_found_peptide OFF)
foreach(_line IN LISTS _lines)
  if (_line STREQUAL "")
    continue()
  endif()
  math(EXPR _rows "${_rows} + 1")
  string(REPLACE "\t" ";" _cells "${_line}")
  list(GET _cells ${_col_peptidoform} _peptidoform)
  list(GET _cells ${_col_charge} _charge)
  list(GET _cells ${_col_scan} _scan)
  list(GET _cells ${_col_peak_mz_array} _peak_mz)
  list(GET _cells ${_col_peak_intensity_array} _peak_intensity)
  list(GET _cells ${_col_predicted_ion_type_array} _ions)
  list(GET _cells ${_col_predicted_mz_array} _mz)
  list(GET _cells ${_col_predicted_intensity_array} _intensity)

  # peaks
  string(REGEX REPLACE "^<(.*)>$" "\\1" _scan "${_scan}")
  string(MAKE_C_IDENTIFIER "${NATIVE_ID_PREFIX}${_scan}" _key)
  if (NOT DEFINED _peaks_${_key})
    message(FATAL_ERROR "${TSV}: row ${_rows} (${_peptidoform}): no spectrum '${NATIVE_ID_PREFIX}${_scan}' in ${MZML}")
  endif()
  list_length("${_peak_mz}" _n_mz)
  list_length("${_peak_intensity}" _n_intensity)
  if (NOT _n_mz EQUAL _peaks_${_key} OR NOT _n_intensity EQUAL _peaks_${_key})
    message(FATAL_ERROR "${TSV}: row ${_rows} (${_peptidoform}, scan ${_scan}): ${_n_mz} m/z and ${_n_intensity} "
                        "intensities, but the spectrum has ${_peaks_${_key}} peaks")
  endif()

  # predictions
  list_length("${_ions}" _n_ions)
  list_length("${_mz}" _n_pred_mz)
  list_length("${_intensity}" _n_pred_intensity)
  if (PREDICTIONS)
    if (_charge GREATER_EQUAL 1 AND _charge LESS_EQUAL 6)
      if (_n_ions LESS 1 OR NOT _n_pred_mz EQUAL _n_ions OR NOT _n_pred_intensity EQUAL _n_ions)
        message(FATAL_ERROR "${TSV}: row ${_rows} (${_peptidoform}, charge ${_charge}): ${_n_ions} predicted ions, "
                            "${_n_pred_mz} m/z, ${_n_pred_intensity} intensities")
      endif()
      math(EXPR _predicted_rows "${_predicted_rows} + 1")
    endif()
  elseif (NOT _n_ions EQUAL -1 OR NOT _n_pred_mz EQUAL -1 OR NOT _n_pred_intensity EQUAL -1)
    message(FATAL_ERROR "${TSV}: row ${_rows} (${_peptidoform}): a prediction without ONNX Runtime")
  endif()

  if (DEFINED PEPTIDE AND NOT _found_peptide AND _peptidoform STREQUAL PEPTIDE AND _charge EQUAL CHARGE)
    set(_found_peptide ON)
    string(REGEX REPLACE "^<(.*)>$" "\\1" _ions "${_ions}")
    string(REGEX REPLACE "^<(.*)>$" "\\1" _mz "${_mz}")
    string(REGEX REPLACE "^<(.*)>$" "\\1" _intensity "${_intensity}")
    string(REPLACE ", " ";" _ions "${_ions}")
    string(REPLACE ", " ";" _mz "${_mz}")
    string(REPLACE ", " ";" _intensity "${_intensity}")
    set(_predicted "")
    foreach(_ion _m _i IN ZIP_LISTS _ions _mz _intensity)
      string(APPEND _predicted "${_ion}\t${_m}\t${_i}\n")
    endforeach()
    file(WRITE "${OUT_PREDICTED}" "${_predicted}")
  endif()
endforeach()

if (_rows EQUAL 0)
  message(FATAL_ERROR "${TSV}: no rows")
endif()
if (PREDICTIONS AND _predicted_rows EQUAL 0)
  message(FATAL_ERROR "${TSV}: no row of charge 1 to 6")
endif()

if (DEFINED PEPTIDE)
  if (NOT _found_peptide)
    message(FATAL_ERROR "${TSV}: no row of ${PEPTIDE} with charge ${CHARGE}")
  endif()
  # The NCE of the predictions, as the summary reports it for the (one) file.
  file(STRINGS "${YAML}" _nce REGEX "^ +peptdeep_nce: ")
  list(LENGTH _nce _n_nce)
  if (NOT _n_nce EQUAL 1)
    message(FATAL_ERROR "${YAML}: expected one 'peptdeep_nce', found ${_n_nce}")
  endif()
  string(REGEX REPLACE "^ +peptdeep_nce: " "" _nce "${_nce}")
  set(_reference "")
  file(STRINGS "${REFERENCE}" _reference_lines)
  foreach(_r IN LISTS _reference_lines)
    string(REPLACE "\t" ";" _r "${_r}")
    list(POP_FRONT _r _r_nce)
    if (_r_nce EQUAL _nce)
      string(REPLACE ";" "\t" _r "${_r}")
      string(APPEND _reference "${_r}\n")
    endif()
  endforeach()
  if (_reference STREQUAL "")
    message(FATAL_ERROR "${REFERENCE}: no prediction for NCE ${_nce}")
  endif()
  file(WRITE "${OUT_REFERENCE}" "${_reference}")
endif()

message(STATUS "${TSV}: ${_rows} rows with their spectra's peaks, ${_predicted_rows} with a prediction")
