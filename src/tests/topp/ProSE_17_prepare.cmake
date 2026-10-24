# Inputs of TOPP_ProSE_17 (databases ProSE refuses, or searches with a warning, before the search), in OUT_DIR, from the
# target database FASTA:
# - ProSE_17_twice.fasta: FASTA twice (every entry repeated: searched once, with a warning);
# - ProSE_17_conflict.fasta: FASTA and its first accession again, with another sequence (refused);
# - ProSE_17_few_decoys.fasta: FASTA (no decoys) and one DECOY_ entry (decoys=auto refuses it).
file(READ "${FASTA}" _content)
file(WRITE "${OUT_DIR}/ProSE_17_twice.fasta" "${_content}${_content}")
string(REGEX MATCH "^>([^ \r\n]+)" _first "${_content}")
file(WRITE "${OUT_DIR}/ProSE_17_conflict.fasta" "${_content}\n>${CMAKE_MATCH_1}\nPEPTIDEKPEPTIDEK\n")
file(WRITE "${OUT_DIR}/ProSE_17_few_decoys.fasta" "${_content}\n>DECOY_ProSE_17\nKEDITPEPKEDITPEP\n")
