## MS2 Search-Engines go here...
## MACRO OPENMS_FINDBINARY:
## fills ${varname} with the path to the binary given in ${binaryname}
## @param varname      Name of the variable which will hold the result string (e.g. OMSSA_BINARY)
## @param binaryname   List of binary names which are searched
## @param name         Human readable version of binaryname for messages
macro (OPENMS_FINDBINARY varname binaryname name)
  find_program(${varname} ${binaryname} PATHS ENV PATH)
  if (${${varname}} STREQUAL "${varname}-NOTFOUND")
    message(STATUS "  - ${name} not found")
  else()
    get_filename_component(found_executable_name ${${varname}} NAME)
    message(STATUS "  + ${name} binary found at ${found_executable_name} -> Enabling corresponding tests.")
  endif()
endmacro (OPENMS_FINDBINARY)

macro (openms_check_tandem_version binary valid)
  if(NOT (${XTANDEM_BINARY} STREQUAL "XTANDEM_BINARY-NOTFOUND"))
    set(${valid} FALSE)
    execute_process(COMMAND ${XTANDEM_BINARY}
      RESULT_VARIABLE _tandem_result
      OUTPUT_VARIABLE _tandem_output
      ERROR_VARIABLE _tandem_error
      INPUT_FILE ${DATA_DIR_TOPP}/SEARCHENGINES/tandem_break.txt
    )

    # we are looking for something like (2013.09.01.1)
    string(REGEX MATCH "\([0-9]+[.][0-9]+[.][0-9]+([.][0-9]+)\)"
          _tandem_version ${_tandem_output})

    if("${_tandem_version}" VERSION_LESS "2013.09.01")
      message(STATUS "  - X!Tandem too old. Please provide a X!Tandem version >= 2013.09.01 to enable the tests.")
    else()
      set(${valid} TRUE)
    endif()
  endif()
endmacro (openms_check_tandem_version)

message(STATUS "Searching for MS2 search engines ...")

#------------------------------------------------------------------------------
# OMSSA
OPENMS_FINDBINARY(OMSSA_BINARY "omssacl" "OMSSA")

#------------------------------------------------------------------------------
# X!Tandem
OPENMS_FINDBINARY(XTANDEM_BINARY "tandem;tandem.exe" "X!Tandem")
openms_check_tandem_version(${XTANDEM_BINARY} xtandem_valid)

#------------------------------------------------------------------------------
# MyriMatch
OPENMS_FINDBINARY(MYRIMATCH_BINARY "myrimatch" "Myrimatch")

#------------------------------------------------------------------------------
## optional tests
if (NOT (${OMSSA_BINARY} STREQUAL "OMSSA_BINARY-NOTFOUND"))
  add_test("TOPP_OMSSAAdapter_1" ${TOPP_BIN_PATH}/OMSSAAdapter -test -ini ${DATA_DIR_TOPP}/SEARCHENGINES/OMSSAAdapter_1.ini -database ${DATA_DIR_TOPP}/SEARCHENGINES/OMSSAAdapter_1.fasta -in ${DATA_DIR_TOPP}/SEARCHENGINES/OMSSAAdapter_1.mzML -out OMSSAAdapter_1_out.tmp -omssa_executable "${OMSSA_BINARY}")
  add_test("TOPP_OMSSAAdapter_1_out" ${DIFF} -in1 OMSSAAdapter_1_out.tmp -in2 ${DATA_DIR_TOPP}/SEARCHENGINES/OMSSAAdapter_1_out.idXML -whitelist "IdentificationRun date" "SearchParameters id=\"SP_0\" db=")
  set_tests_properties("TOPP_OMSSAAdapter_1_out" PROPERTIES DEPENDS "TOPP_OMSSAAdapter_1")

  # test charge check
  add_test("TOPP_OMSSAAdapter_2" ${TOPP_BIN_PATH}/OMSSAAdapter -test -ini ${DATA_DIR_TOPP}/SEARCHENGINES/OMSSAAdapter_1.ini -database ${DATA_DIR_TOPP}/SEARCHENGINES/OMSSAAdapter_1.fasta -in ${DATA_DIR_TOPP}/SEARCHENGINES/OMSSAAdapter_1.mzML -out OMSSAAdapter_1_out.tmp -omssa_executable "${OMSSA_BINARY}" -min_precursor_charge 4 -max_precursor_charge 3)
  set_tests_properties("TOPP_OMSSAAdapter_2" PROPERTIES WILL_FAIL 1) ## has invalid charge range
endif()

#------------------------------------------------------------------------------
if (NOT (${XTANDEM_BINARY} STREQUAL "XTANDEM_BINARY-NOTFOUND") AND xtandem_valid)
  add_test("TOPP_XTandemAdapter_1" ${TOPP_BIN_PATH}/XTandemAdapter -test -ini ${DATA_DIR_TOPP}/SEARCHENGINES/XTandemAdapter_1.ini -database ${DATA_DIR_TOPP}/SEARCHENGINES/OMSSAAdapter_1.fasta -in ${DATA_DIR_TOPP}/SEARCHENGINES/OMSSAAdapter_1.mzML -out XTandemAdapter_1_out.tmp -xtandem_executable "${XTANDEM_BINARY}")
  add_test("TOPP_XTandemAdapter_1_out" ${DIFF} -in1 XTandemAdapter_1_out.tmp -in2 ${DATA_DIR_TOPP}/SEARCHENGINES/XTandemAdapter_2_out.idXML -whitelist "IdentificationRun date" "SearchParameters id=\"SP_0\" db=")
  set_tests_properties("TOPP_XTandemAdapter_1_out" PROPERTIES DEPENDS "TOPP_XTandemAdapter_1")

  # test charge check
  add_test("TOPP_XTandemAdapter_2" ${TOPP_BIN_PATH}/XTandemAdapter -test -ini ${DATA_DIR_TOPP}/SEARCHENGINES/XTandemAdapter_1.ini -database ${DATA_DIR_TOPP}/SEARCHENGINES/OMSSAAdapter_1.fasta -in ${DATA_DIR_TOPP}/SEARCHENGINES/OMSSAAdapter_1.mzML -out XTandemAdapter_1_out.tmp -xtandem_executable "${XTANDEM_BINARY}" -min_precursor_charge 4 -max_precursor_charge 3)
  set_tests_properties("TOPP_XTandemAdapter_2" PROPERTIES WILL_FAIL 1) ## has invalid charge range
endif()

#------------------------------------------------------------------------------
if (NOT (${MYRIMATCH_BINARY} STREQUAL "MYRIMATCH_BINARY-NOTFOUND"))
  add_test("TOPP_MyriMatchAdapter_1" ${TOPP_BIN_PATH}/MyriMatchAdapter -test -ini ${DATA_DIR_TOPP}/SEARCHENGINES/MyriMatchAdapter_1.ini -database ${DATA_DIR_TOPP}/SEARCHENGINES/OMSSAAdapter_1.fasta -in ${DATA_DIR_TOPP}/SEARCHENGINES/OMSSAAdapter_1.mzML -out MyriMatchAdapter_1_out.tmp -myrimatch_executable "${MYRIMATCH_BINARY}")
  add_test("TOPP_MyriMatchAdapter_1_out" ${DIFF} -in1 MyriMatchAdapter_1_out.tmp -in2 ${DATA_DIR_TOPP}/SEARCHENGINES/MyriMatchAdapter_1_out.idXML -whitelist "IdentificationRun date" "SearchParameters id=\"SP_0\" db=")
  set_tests_properties("TOPP_MyriMatchAdapter_1_out" PROPERTIES DEPENDS "TOPP_MyriMatchAdapter_1")

  # test charge check
  add_test("TOPP_MyriMatchAdapter_2" ${TOPP_BIN_PATH}/MyriMatchAdapter -test -ini ${DATA_DIR_TOPP}/SEARCHENGINES/MyriMatchAdapter_1.ini -database ${DATA_DIR_TOPP}/SEARCHENGINES/OMSSAAdapter_1.fasta -in ${DATA_DIR_TOPP}/SEARCHENGINES/OMSSAAdapter_1.mzML -out MyriMatchAdapter_1_out.tmp -myrimatch_executable "${MYRIMATCH_BINARY}" -min_precursor_charge 4 -max_precursor_charge 3)
  set_tests_properties("TOPP_MyriMatchAdapter_2" PROPERTIES WILL_FAIL 1) ## has invalid charge range
endif()


# TODO the following tests are waiting for better implementations of InspectAdapter and associated classes
#add_test("TOPP_InspectAdapter_3" ${TOPP_BIN_PATH}/InspectAdapter -ini ${DATA_DIR_TOPP}/InspectAdapter_1_parameters.ini -trie_dbs ${DATA_DIR_TOPP}/Inspect_FASTAFile_test2.trie -in ${DATA_DIR_TOPP}/InspectAdapter.out -dbs ${DATA_DIR_TOPP}/Inspect_FASTAFile_test.fasta -out InspectAdapter_4_output.tmp -inspect_out)
#add_test("TOPP_InspectAdapter_3_out1" ${DIFF} -whitelist "?xml-stylesheet" "IdentificationRun date" -in1 InspectAdapter_4_output.tmp -in2 ${DATA_DIR_TOPP}/InspectAdapter_4_output.idXML )
#set_tests_properties("TOPP_InspectAdapter_3_out1" PROPERTIES DEPENDS "TOPP_InspectAdapter_3")

### SequestAdapter tests
# TODO disabled until tool is reactivated
#add_test("TOPP_SequestAdapter_1" ${TOPP_BIN_PATH}/SequestAdapter -ini ${DATA_DIR_TOPP}/SequestAdapter_1_parameters.ini -in ${DATA_DIR_TOPP}/Sequest.mzXML -mz_files ${DATA_DIR_TOPP}/Sequest.mzXML -modifications_xml_file ${DATA_DIR_TOPP}/Sequest_PTMs.xml -out SequestAdapter_2_output.tmp -sequest_in -temp_data_directory ${DATA_DIR_TOPP} -db ${DATA_DIR_TOPP}/Sequest_FASTAFile_test.fasta)
#add_test("TOPP_SequestAdapter_1_out1" ${DIFF} -in1 SequestAdapter_2_output.tmp -in2 ${DATA_DIR_TOPP}/SequestAdapter_2_output.sequest_in)
#add_test("TOPP_SequestAdapter_2" ${TOPP_BIN_PATH}/SequestAdapter -ini ${DATA_DIR_TOPP}/SequestAdapter_1_parameters.ini -in ${DATA_DIR_TOPP}/Sequest.mzData -mz_files ${DATA_DIR_TOPP}/Sequest.mzXML -modifications_xml_file ${DATA_DIR_TOPP}/Sequest_PTMs.xml -out SequestAdapter_3_output.tmp -sequest_in -temp_data_directory ${DATA_DIR_TOPP} -db ${DATA_DIR_TOPP}/Sequest_FASTAFile_test.fasta)
#add_test("TOPP_SequestAdapter_2_out1" ${DIFF} -in1 SequestAdapter_3_output.tmp -in2 ${DATA_DIR_TOPP}/SequestAdapter_2_output.sequest_in)

# TODO the following tests are waiting for better implementations of InspectAdapter and
# associated classes
#add_test("TOPP_SequestAdapter_3" ${TOPP_BIN_PATH}/SequestAdapter -ini ${DATA_DIR_TOPP}/SequestAdapter_2_parameters.ini -mz_files ${DATA_DIR_TOPP}/Sequest.mzXML -modifications_xml_file ${DATA_DIR_TOPP}/Sequest_PTMs.xml -in ${DATA_DIR_TOPP}/Sequest.mzXML -out SequestAdapter_4_output.tmp -temp_data_directory ${DATA_DIR_TOPP} -db ${DATA_DIR_TOPP}/Sequest_FASTAFile_test.fasta)
#add_test("TOPP_SequestAdapter_3_out1" ${DIFF} -in1 SequestAdapter_4_output.tmp -in2 ${DATA_DIR_TOPP}/SequestAdapter_4_output.idXML)

### SpecLibSearcher tests
#add_test("TOPP_SpecLibSearcher_1" ${TOPP_BIN_PATH}/SpecLibSearcher -test -ini ${DATA_DIR_TOPP}/SpecLibSearcher_1_parameters.ini -in ${DATA_DIR_TOPP}/SpecLibSearcher_1.MzData -lib $(DATA_DIR_TOPP)/SpecLibSearcher_1.MSP -out SpecLibSearcher_1.tmp)
#add_test("TOPP_SpecLibSearcher_1_out1" ${DIFF} -in1 SpecLibSearcher_1.tmp  -in2 $(DATA_DIR_TOPP)/SpecLibSearcher_1.idXML)
### PepNovoAdapter tests
#The PepNovoAdapter now only works as a frontend and cannot be run without an installation of PepNovo.Therefore no test possible
#add_test("TOPP_PepNovoAdapter_1" ${TOPP_BIN_PATH}/PepNovoAdapter -ini ${DATA_DIR_TOPP}/PepNovoAdapter_1_parameters.ini -in ${DATA_DIR_TOPP}/PepNovo.mzXML -pepnovo_in -out PepNovoAdapter_3_output.tmp -dta_list ${DATA_DIR_TOPP}/tmp/dta_list.txt -model_directory ${DATA_DIR_TOPP}/tmp/ -temp_data_directory ${DATA_DIR_TOPP}/tmp/)
#add_test("TOPP_PepNovoAdapter_1_out1" ${DIFF} -in1 ${DATA_DIR_TOPP}/tmp/PepNovo_PTMs_.txt -in2 ${DATA_DIR_TOPP}/tmp/PepNovo_PTMs.txt)
#TODO ANDREAS - We have to clean up the /tmp/ directory to run this test multiple times
#add_test("TOPP_PepNovoAdapter_2" ${TOPP_BIN_PATH}/PepNovoAdapter -ini ${DATA_DIR_TOPP}/PepNovoAdapter_1_parameters.ini -in ${DATA_DIR_TOPP}/PepNovo.mzData -pepnovo_in -out PepNovoAdapter_4_output.tmp -temp_data_directory ${DATA_DIR_TOPP})
#add_test("TOPP_PepNovoAdapter_2_out1" ${DIFF} -in1 ${DATA_DIR_TOPP}/PepNovo_PTMs_.txt -in2 ${DATA_DIR_TOPP}/PepNovo_PTMs.txt)
#add_test("TOPP_PepNovoAdapter_3" ${TOPP_BIN_PATH}/PepNovoAdapter -ini ${DATA_DIR_TOPP}/PepNovoAdapter_5_parameters.ini -in ${DATA_DIR_TOPP}/PepNovoAdapter_5_output.pepnovo_out -out PepNovoAdapter_5_output.tmp -pepnovo_out -dta_list ${DATA_DIR_TOPP}/tmp/dta_list.txt -model_directory ${DATA_DIR_TOPP}/tmp/ -temp_data_directory ${DATA_DIR_TOPP}/tmp/ -modifications_xml_file ${DATA_DIR_TOPP}/PepNovo_PTMs.xml -mz_files ${DATA_DIR_TOPP}/PepNovo.mzXML)
#add_test("TOPP_PepNovoAdapter_3_out1" ${DIFF} -whitelist "?xml-stylesheet" "date_group_1" -in1 PepNovoAdapter_5_output.tmp -in2 ${DATA_DIR_TOPP}/PepNovoAdapter_5_output.idXML)
