# Inputs of TOPP_ProSE_18 (databases ProSE refuses, or searches with a warning, before the search) and of
# TOPP_ProSE_11_rescore_mixed, in OUT_DIR, from the target database FASTA:
# - ProSE_18_twice.fasta: FASTA twice (every entry repeated: searched once, with a warning);
# - ProSE_18_conflict.fasta: FASTA and its first accession again, with another sequence (refused);
# - ProSE_18_few_decoys.fasta: FASTA (no decoys) and one DECOY_ entry (decoys=auto refuses it).
file(READ "${FASTA}" _content)
file(WRITE "${OUT_DIR}/ProSE_18_twice.fasta" "${_content}${_content}")
string(REGEX MATCH "^>([^ \r\n]+)" _first "${_content}")
file(WRITE "${OUT_DIR}/ProSE_18_conflict.fasta" "${_content}\n>${CMAKE_MATCH_1}\nPEPTIDEKPEPTIDEK\n")
file(WRITE "${OUT_DIR}/ProSE_18_few_decoys.fasta" "${_content}\n>DECOY_ProSE_18\nKEDITPEPKEDITPEP\n")
# - ProSE_18_small.mzML (with -DMZML): the first 60 spectra of MZML, as a plain mzML (without the index), an input with
#   fewer than 100 PSMs, too few for a Percolator model of its own.
if(DEFINED MZML)
  file(READ "${MZML}" _mzml)
  string(FIND "${_mzml}" "<spectrum " _first)
  string(SUBSTRING "${_mzml}" 0 ${_first} _header)
  string(REGEX REPLACE "<indexedmzML[^>]*>" "" _header "${_header}")
  string(REGEX REPLACE "<spectrumList count=\"[0-9]+\"" "<spectrumList count=\"60\"" _header "${_header}")
  set(_end ${_first})
  foreach(_i RANGE 1 60)
    string(SUBSTRING "${_mzml}" ${_end} -1 _rest)
    string(FIND "${_rest}" "</spectrum>" _close)
    math(EXPR _end "${_end} + ${_close} + 11")
  endforeach()
  math(EXPR _length "${_end} - ${_first}")
  string(SUBSTRING "${_mzml}" ${_first} ${_length} _spectra)
  file(WRITE "${OUT_DIR}/ProSE_18_small.mzML" "${_header}${_spectra}\n</spectrumList>\n</run>\n</mzML>\n")
endif()
