# Inputs of TOPP_ProSE_18 (databases ProSE refuses, or searches with a warning, before the search), in OUT_DIR, from the
# target database FASTA:
# - ProSE_18_twice.fasta: FASTA twice (every entry repeated: searched once, with a warning);
# - ProSE_18_conflict.fasta: FASTA and its first accession again, with another sequence (refused);
# - ProSE_18_few_decoys.fasta: FASTA (no decoys) and one DECOY_ entry (decoys=auto refuses it).
file(READ "${FASTA}" _content)
file(WRITE "${OUT_DIR}/ProSE_18_twice.fasta" "${_content}${_content}")
string(REGEX MATCH "^>([^ \r\n]+)" _first "${_content}")
file(WRITE "${OUT_DIR}/ProSE_18_conflict.fasta" "${_content}\n>${CMAKE_MATCH_1}\nPEPTIDEKPEPTIDEK\n")
file(WRITE "${OUT_DIR}/ProSE_18_few_decoys.fasta" "${_content}\n>DECOY_ProSE_18\nKEDITPEPKEDITPEP\n")
