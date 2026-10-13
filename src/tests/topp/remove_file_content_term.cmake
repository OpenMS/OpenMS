# Copies INPUT_FILE to OUTPUT_FILE without the lines in the mzML <fileContent> block that contain TERM.
# The rest of the file is copied unchanged.
if(NOT DEFINED INPUT_FILE OR NOT DEFINED OUTPUT_FILE OR NOT DEFINED TERM)
  message(FATAL_ERROR "INPUT_FILE, OUTPUT_FILE and TERM must be set")
endif()
file(READ "${INPUT_FILE}" content)
string(FIND "${content}" "</fileContent>" end)
if(end EQUAL -1)
  message(FATAL_ERROR "No fileContent block in ${INPUT_FILE}")
endif()
string(SUBSTRING "${content}" 0 ${end} head)
string(SUBSTRING "${content}" ${end} -1 tail)
string(REGEX REPLACE "[^\n]*${TERM}[^\n]*\n" "" head "${head}")
file(WRITE "${OUTPUT_FILE}" "${head}${tail}")
