## Third-party tools (e.g. MS2 search engines) go here...
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
    execute_process(COMMAND "${XTANDEM_BINARY}"
      RESULT_VARIABLE _tandem_result
      OUTPUT_VARIABLE _tandem_output
      ERROR_VARIABLE _tandem_output  ## write to the same variable, in case Tandem decides to use std::cerr one day
      INPUT_FILE ${DATA_DIR_TOPP}/THIRDPARTY/tandem_break.txt  ## provide some input, otherwise tandem.exe will block and not finish
    )

    # we are looking for something like (2013.09.01.1)
    string(REGEX MATCH "\([0-9]+[.][0-9]+[.][0-9]+([.][0-9]+)\)"
          _tandem_version "${_tandem_output}")

    if("${_tandem_version}" VERSION_LESS "2013.09.01")
      message(STATUS "  - X! Tandem too old (${_tandem_version}). Please provide an X! Tandem version >= 2013.09.01 to enable the tests.")
    else()
      message(STATUS "  + X! Tandem version: ${_tandem_version}.")
      set(${valid} TRUE)
    endif()
  endif()
endmacro (openms_check_tandem_version)

message(STATUS "Searching for third party tools...")

#------------------------------------------------------------------------------
# Comet
OPENMS_FINDBINARY(COMET_BINARY "comet.exe" "Comet")

#------------------------------------------------------------------------------
# OMSSA
OPENMS_FINDBINARY(OMSSA_BINARY "omssacl" "OMSSA")

#------------------------------------------------------------------------------
# X!Tandem
OPENMS_FINDBINARY(XTANDEM_BINARY "tandem;tandem.exe" "X! Tandem")
openms_check_tandem_version(${XTANDEM_BINARY} xtandem_valid)

#------------------------------------------------------------------------------
# MyriMatch
OPENMS_FINDBINARY(MYRIMATCH_BINARY "myrimatch" "Myrimatch")

#------------------------------------------------------------------------------
# MS-GF+
OPENMS_FINDBINARY(MSGFPLUS_BINARY "MSGFPlus.jar" "MS-GF+")

#------------------------------------------------------------------------------
# percolator
OPENMS_FINDBINARY(PERCOLATOR_BINARY "percolator" "Percolator")

#------------------------------------------------------------------------------
# Fido
OPENMS_FINDBINARY(FIDO_BINARY "Fido" "Fido")
OPENMS_FINDBINARY(FIDOCHOOSEPARAMS_BINARY "FidoChooseParameters" "FidoChooseParameters")

#------------------------------------------------------------------------------
# Sirius
OPENMS_FINDBINARY(SIRIUS_BINARY "sirius;sirius-console-64.exe" "Sirius")

#------------------------------------------------------------------------------
# Spectrast
OPENMS_FINDBINARY(SPECTRAST_BINARY "spectrast" "SpectraST")

#------------------------------------------------------------------------------
## optional tests
if (NOT (${OMSSA_BINARY} STREQUAL "OMSSA_BINARY-NOTFOUND"))
  add_test("TOPP_OMSSAAdapter_1" ${TOPP_BIN_PATH}/OMSSAAdapter -test -ini ${DATA_DIR_TOPP}/THIRDPARTY/OMSSAAdapter_1.ini -database ${DATA_DIR_TOPP}/THIRDPARTY/proteins.fasta -in ${DATA_DIR_TOPP}/THIRDPARTY/spectra.mzML -out OMSSAAdapter_1_out.tmp -omssa_executable "${OMSSA_BINARY}")
  add_test("TOPP_OMSSAAdapter_1_out" ${DIFF} -in1 OMSSAAdapter_1_out.tmp -in2 ${DATA_DIR_TOPP}/THIRDPARTY/OMSSAAdapter_1_out.idXML -whitelist "IdentificationRun date" "SearchParameters id=\"SP_0\" db=" "UserParam type=\"stringList\" name=\"spectra_data\" value=")
  set_tests_properties("TOPP_OMSSAAdapter_1_out" PROPERTIES DEPENDS "TOPP_OMSSAAdapter_1")

  # test charge check
  add_test("TOPP_OMSSAAdapter_2" ${TOPP_BIN_PATH}/OMSSAAdapter -test -ini ${DATA_DIR_TOPP}/THIRDPARTY/OMSSAAdapter_1.ini -database ${DATA_DIR_TOPP}/THIRDPARTY/proteins.fasta -in ${DATA_DIR_TOPP}/THIRDPARTY/spectra.mzML -out OMSSAAdapter_1_out.tmp -omssa_executable "${OMSSA_BINARY}" -min_precursor_charge 4 -max_precursor_charge 3)
  set_tests_properties("TOPP_OMSSAAdapter_2" PROPERTIES WILL_FAIL 1) # has invalid charge range
endif()

#------------------------------------------------------------------------------
if (NOT (${XTANDEM_BINARY} STREQUAL "XTANDEM_BINARY-NOTFOUND") AND xtandem_valid)
  add_test("TOPP_XTandemAdapter_1" ${TOPP_BIN_PATH}/XTandemAdapter -test -ini ${DATA_DIR_TOPP}/THIRDPARTY/XTandemAdapter_1.ini -database ${DATA_DIR_TOPP}/THIRDPARTY/proteins.fasta -in ${DATA_DIR_TOPP}/THIRDPARTY/spectra.mzML -out XTandemAdapter_1_out.tmp -xtandem_executable "${XTANDEM_BINARY}")
  add_test("TOPP_XTandemAdapter_1_out" ${DIFF} -in1 XTandemAdapter_1_out.tmp -in2 ${DATA_DIR_TOPP}/THIRDPARTY/XTandemAdapter_1_out.idXML -whitelist "IdentificationRun date" "SearchParameters id=\"SP_0\" db=" "UserParam type=\"stringList\" name=\"spectra_data\" value=")
  set_tests_properties("TOPP_XTandemAdapter_1_out" PROPERTIES DEPENDS "TOPP_XTandemAdapter_1")

  # test output result option (set it to 'valid')
  add_test("TOPP_XTandemAdapter_2" ${TOPP_BIN_PATH}/XTandemAdapter -test -ini ${DATA_DIR_TOPP}/THIRDPARTY/XTandemAdapter_1.ini -database ${DATA_DIR_TOPP}/THIRDPARTY/proteins.fasta -in ${DATA_DIR_TOPP}/THIRDPARTY/spectra.mzML -out XTandemAdapter_2_out.tmp -output_results valid -xtandem_executable "${XTANDEM_BINARY}" -max_valid_expect 1e-14)
  add_test("TOPP_XTandemAdapter_2_out" ${DIFF} -in1 XTandemAdapter_2_out.tmp -in2 ${DATA_DIR_TOPP}/THIRDPARTY/XTandemAdapter_2_out.idXML -whitelist "IdentificationRun date" "SearchParameters id=\"SP_0\" db=" "UserParam type=\"stringList\" name=\"spectra_data\" value=")
  set_tests_properties("TOPP_XTandemAdapter_2_out" PROPERTIES DEPENDS "TOPP_XTandemAdapter_2")
endif()

#------------------------------------------------------------------------------
if (NOT (${MYRIMATCH_BINARY} STREQUAL "MYRIMATCH_BINARY-NOTFOUND"))
  add_test("TOPP_MyriMatchAdapter_1" ${TOPP_BIN_PATH}/MyriMatchAdapter -test -ini ${DATA_DIR_TOPP}/THIRDPARTY/MyriMatchAdapter_1.ini -database ${DATA_DIR_TOPP}/THIRDPARTY/proteins.fasta -in ${DATA_DIR_TOPP}/THIRDPARTY/spectra.mzML -out MyriMatchAdapter_1_out.tmp -myrimatch_executable "${MYRIMATCH_BINARY}")
  add_test("TOPP_MyriMatchAdapter_1_out" ${DIFF} -in1 MyriMatchAdapter_1_out.tmp -in2 ${DATA_DIR_TOPP}/THIRDPARTY/MyriMatchAdapter_1_out.idXML -whitelist "IdentificationRun date" "SearchParameters id=\"SP_0\" db=" "UserParam type=\"stringList\" name=\"spectra_data\" value=")
  set_tests_properties("TOPP_MyriMatchAdapter_1_out" PROPERTIES DEPENDS "TOPP_MyriMatchAdapter_1")

  # test charge check
  add_test("TOPP_MyriMatchAdapter_2" ${TOPP_BIN_PATH}/MyriMatchAdapter -test -ini ${DATA_DIR_TOPP}/THIRDPARTY/MyriMatchAdapter_1.ini -database ${DATA_DIR_TOPP}/THIRDPARTY/proteins.fasta -in ${DATA_DIR_TOPP}/THIRDPARTY/spectra.mzML -out MyriMatchAdapter_1_out.tmp -myrimatch_executable "${MYRIMATCH_BINARY}" -min_precursor_charge 4 -max_precursor_charge 3)
  set_tests_properties("TOPP_MyriMatchAdapter_2" PROPERTIES WILL_FAIL 1) # has invalid charge range
endif()

#------------------------------------------------------------------------------
if (NOT (${MSGFPLUS_BINARY} STREQUAL "MSGFPLUS_BINARY-NOTFOUND"))
  add_test("TOPP_MSGFPlusAdapter_1" ${TOPP_BIN_PATH}/MSGFPlusAdapter -test -ini ${DATA_DIR_TOPP}/THIRDPARTY/MSGFPlusAdapter_1.ini -database ${DATA_DIR_TOPP}/THIRDPARTY/proteins.fasta -in ${DATA_DIR_TOPP}/THIRDPARTY/spectra.mzML -out MSGFPlusAdapter_1_out1.tmp -mzid_out MSGFPlusAdapter_1_out2.tmp.mzid -executable "${MSGFPLUS_BINARY}")
  add_test("TOPP_MSGFPlusAdapter_1_out1" ${DIFF} -in1 MSGFPlusAdapter_1_out1.tmp -in2 ${DATA_DIR_TOPP}/THIRDPARTY/MSGFPlusAdapter_1_out.idXML -whitelist "IdentificationRun date" "SearchParameters id=\"SP_0\" db=" "UserParam type=\"stringList\" name=\"spectra_data\" value=")
  set_tests_properties("TOPP_MSGFPlusAdapter_1_out1" PROPERTIES DEPENDS "TOPP_MSGFPlusAdapter_1")
  add_test("TOPP_MSGFPlusAdapter_1_out2" ${DIFF} -in1 MSGFPlusAdapter_1_out2.tmp.mzid -in2 ${DATA_DIR_TOPP}/THIRDPARTY/MSGFPlusAdapter_1_out.mzid -whitelist "creationDate=" "SearchDatabase numDatabaseSequences=\"10\" location=" "SpectraData location=" "AnalysisSoftware")
  set_tests_properties("TOPP_MSGFPlusAdapter_1_out2" PROPERTIES DEPENDS "TOPP_MSGFPlusAdapter_1")
endif()

#------------------------------------------------------------------------------
if (NOT (${COMET_BINARY} STREQUAL "COMET_BINARY-NOTFOUND"))
  ### NOT needs to be added after the binarys have been included
  add_test("TOPP_CometAdapter_1" ${TOPP_BIN_PATH}/CometAdapter -test -ini ${DATA_DIR_TOPP}/THIRDPARTY/CometAdapter_1.ini -database ${DATA_DIR_TOPP}/THIRDPARTY/proteins.fasta -in ${DATA_DIR_TOPP}/THIRDPARTY/spectra_comet.mzML -out CometAdapter_1_out1.tmp -pin_out CometAdapter_1_out2.tmp.csv -comet_executable "${COMET_BINARY}")
  add_test("TOPP_CometAdapter_1_out1" ${DIFF} -in1 CometAdapter_1_out1.tmp -in2 ${DATA_DIR_TOPP}/THIRDPARTY/CometAdapter_1_out.idXML -whitelist "IdentificationRun date" "SearchParameters id=\"SP_0\" db=" "UserParam type=\"stringList\" name=\"spectra_data\" value=")
  set_tests_properties("TOPP_CometAdapter_1_out1" PROPERTIES DEPENDS "TOPP_CometAdapter_1")
  ### Second test for optional pin file needs to be added, not sure how to do FuzzyDiff on the csv style pin file, whitelisting the first id column
endif()

#------------------------------------------------------------------------------
if (NOT (${PERCOLATOR_BINARY} STREQUAL "PERCOLATOR_BINARY-NOTFOUND"))
  ### NOT needs to be added after the binarys have been included
  add_test("TOPP_PercolatorAdapter_1" ${TOPP_BIN_PATH}/PercolatorAdapter -test -ini ${DATA_DIR_TOPP}/THIRDPARTY/PercolatorAdapter_1.ini -in ${DATA_DIR_TOPP}/THIRDPARTY/PercolatorAdapter_1.idXML -out PercolatorAdapter_1_out1.tmp -percolator_executable "${PERCOLATOR_BINARY}")
  add_test("TOPP_PercolatorAdapter_1_out1" ${DIFF} -in1 PercolatorAdapter_1_out1.tmp -in2 ${DATA_DIR_TOPP}/THIRDPARTY/PercolatorAdapter_1_out.idXML -whitelist "IdentificationRun date" "SearchParameters id=\"SP_0\" db=" "UserParam type=\"stringList\" name=\"spectra_data\" value=")
  set_tests_properties("TOPP_PercolatorAdapter_1_out1" PROPERTIES DEPENDS "TOPP_PercolatorAdapter_1")
  add_test("TOPP_PercolatorAdapter_2" ${TOPP_BIN_PATH}/PercolatorAdapter -test -osw_level ms1 -in_osw ${DATA_DIR_TOPP}/THIRDPARTY/PercolatorAdapter_2.osw -osw_out PercolatorAdapter_2_out1.osw -percolator_executable "${PERCOLATOR_BINARY}")
  add_test("TOPP_PercolatorAdapter_3" ${TOPP_BIN_PATH}/PercolatorAdapter -test -osw_level ms2 -in_osw PercolatorAdapter_2_out1.osw -osw_out PercolatorAdapter_3_out1.osw -percolator_executable "${PERCOLATOR_BINARY}")
  set_tests_properties("TOPP_PercolatorAdapter_3" PROPERTIES DEPENDS "TOPP_PercolatorAdapter_2")
  add_test("TOPP_PercolatorAdapter_4" ${TOPP_BIN_PATH}/PercolatorAdapter -test -osw_level transition -in_osw PercolatorAdapter_3_out1.osw -osw_out PercolatorAdapter_4_out1.osw -percolator_executable "${PERCOLATOR_BINARY}")
  set_tests_properties("TOPP_PercolatorAdapter_4" PROPERTIES DEPENDS "TOPP_PercolatorAdapter_3")
  ### TOPP_PercolatorAdapter_2-4 do not validate output, but checks whether OSW files can be read and written to. 
endif()

#------------------------------------------------------------------------------
if (NOT (${FIDOCHOOSEPARAMS_BINARY} STREQUAL "FIDOCHOOSEPARAMS_BINARY-NOTFOUND"))
  add_test("TOPP_FidoAdapter_1" ${TOPP_BIN_PATH}/FidoAdapter -test -in ${DATA_DIR_TOPP}/THIRDPARTY/FidoAdapter_1_input.idXML -out FidoAdapter_1_output.tmp -fidocp_executable "${FIDOCHOOSEPARAMS_BINARY}")
  add_test("TOPP_FidoAdapter_1_out" ${DIFF} -in1 FidoAdapter_1_output.tmp -in2 ${DATA_DIR_TOPP}/THIRDPARTY/FidoAdapter_1_output.idXML -whitelist "IdentificationRun date")
  set_tests_properties("TOPP_FidoAdapter_1_out" PROPERTIES DEPENDS "TOPP_FidoAdapter_1")

  add_test("TOPP_FidoAdapter_2" ${TOPP_BIN_PATH}/FidoAdapter -test -in ${DATA_DIR_TOPP}/THIRDPARTY/FidoAdapter_1_input.idXML -out FidoAdapter_2_output.tmp -fidocp_executable "${FIDOCHOOSEPARAMS_BINARY}" -separate_runs)
  add_test("TOPP_FidoAdapter_2_out" ${DIFF} -in1 FidoAdapter_2_output.tmp -in2 ${DATA_DIR_TOPP}/THIRDPARTY/FidoAdapter_2_output.idXML -whitelist "IdentificationRun date")
  set_tests_properties("TOPP_FidoAdapter_2_out" PROPERTIES DEPENDS "TOPP_FidoAdapter_2")

  add_test("TOPP_FidoAdapter_3" ${TOPP_BIN_PATH}/FidoAdapter -test -in ${DATA_DIR_TOPP}/THIRDPARTY/FidoAdapter_1_input.idXML -out FidoAdapter_3_output.tmp -fidocp_executable "${FIDOCHOOSEPARAMS_BINARY}" -group_level -all_PSMs)
  add_test("TOPP_FidoAdapter_3_out" ${DIFF} -in1 FidoAdapter_3_output.tmp -in2 ${DATA_DIR_TOPP}/THIRDPARTY/FidoAdapter_3_output.idXML -whitelist "IdentificationRun date")
  set_tests_properties("TOPP_FidoAdapter_3_out" PROPERTIES DEPENDS "TOPP_FidoAdapter_3")

  add_test("TOPP_FidoAdapter_4" ${TOPP_BIN_PATH}/FidoAdapter -test -in ${DATA_DIR_TOPP}/THIRDPARTY/FidoAdapter_4_input.idXML -out FidoAdapter_4_output.tmp -fidocp_executable "${FIDOCHOOSEPARAMS_BINARY}")
  add_test("TOPP_FidoAdapter_4_out" ${DIFF} -in1 FidoAdapter_4_output.tmp -in2 ${DATA_DIR_TOPP}/THIRDPARTY/FidoAdapter_4_output.idXML -whitelist "IdentificationRun date")
  set_tests_properties("TOPP_FidoAdapter_4_out" PROPERTIES DEPENDS "TOPP_FidoAdapter_4")

  add_test("TOPP_FidoAdapter_5" ${TOPP_BIN_PATH}/FidoAdapter -test -greedy_group_resolution -in ${DATA_DIR_TOPP}/THIRDPARTY/FidoAdapter_5_input.idXML -out FidoAdapter_5_output.tmp -fidocp_executable "${FIDOCHOOSEPARAMS_BINARY}")
  add_test("TOPP_FidoAdapter_5_out" ${DIFF} -in1 FidoAdapter_5_output.tmp -in2 ${DATA_DIR_TOPP}/THIRDPARTY/FidoAdapter_5_output.idXML -whitelist "IdentificationRun date")
  set_tests_properties("TOPP_FidoAdapter_5_out" PROPERTIES DEPENDS "TOPP_FidoAdapter_5")
endif()

if (NOT (${FIDO_BINARY} STREQUAL "FIDO_BINARY-NOTFOUND"))
  add_test("TOPP_FidoAdapter_6" ${TOPP_BIN_PATH}/FidoAdapter -test -in ${DATA_DIR_TOPP}/THIRDPARTY/FidoAdapter_1_input.idXML -out FidoAdapter_6_output.tmp -fido_executable "${FIDO_BINARY}" -prob:protein 0.9 -prob:peptide 0.01 -prob:spurious 0.0)
  add_test("TOPP_FidoAdapter_6_out" ${DIFF} -in1 FidoAdapter_6_output.tmp -in2 ${DATA_DIR_TOPP}/THIRDPARTY/FidoAdapter_1_output.idXML -whitelist "IdentificationRun date")
  set_tests_properties("TOPP_FidoAdapter_6_out" PROPERTIES DEPENDS "TOPP_FidoAdapter_6")
endif()


# TODO the following tests are waiting for better implementations of InspectAdapter and associated classes
#add_test("TOPP_InspectAdapter_3" ${TOPP_BIN_PATH}/InspectAdapter -ini ${DATA_DIR_TOPP}/InspectAdapter_1_parameters.ini -trie_dbs ${DATA_DIR_TOPP}/Inspect_FASTAFile_test2.trie -in ${DATA_DIR_TOPP}/InspectAdapter.out -dbs ${DATA_DIR_TOPP}/Inspect_FASTAFile_test.fasta -out InspectAdapter_4_output.tmp -inspect_out)
#add_test("TOPP_InspectAdapter_3_out1" ${DIFF} -whitelist "?xml-stylesheet" "IdentificationRun date" -in1 InspectAdapter_4_output.tmp -in2 ${DATA_DIR_TOPP}/InspectAdapter_4_output.idXML )
#set_tests_properties("TOPP_InspectAdapter_3_out1" PROPERTIES DEPENDS "TOPP_InspectAdapter_3")

### PepNovoAdapter tests
#The PepNovoAdapter now only works as a frontend and cannot be run without an installation of PepNovo.Therefore no test possible
#add_test("TOPP_PepNovoAdapter_1" ${TOPP_BIN_PATH}/PepNovoAdapter -ini ${DATA_DIR_TOPP}/PepNovoAdapter_1_parameters.ini -in ${DATA_DIR_TOPP}/PepNovo.mzXML -pepnovo_in -out PepNovoAdapter_3_output.tmp -dta_list ${DATA_DIR_TOPP}/tmp/dta_list.txt -model_directory ${DATA_DIR_TOPP}/tmp/ -temp_data_directory ${DATA_DIR_TOPP}/tmp/)
#add_test("TOPP_PepNovoAdapter_1_out1" ${DIFF} -in1 ${DATA_DIR_TOPP}/tmp/PepNovo_PTMs_.txt -in2 ${DATA_DIR_TOPP}/tmp/PepNovo_PTMs.txt)
#TODO ANDREAS - We have to clean up the /tmp/ directory to run this test multiple times
#add_test("TOPP_PepNovoAdapter_2" ${TOPP_BIN_PATH}/PepNovoAdapter -ini ${DATA_DIR_TOPP}/PepNovoAdapter_1_parameters.ini -in ${DATA_DIR_TOPP}/PepNovo.mzData -pepnovo_in -out PepNovoAdapter_4_output.tmp -temp_data_directory ${DATA_DIR_TOPP})
#add_test("TOPP_PepNovoAdapter_2_out1" ${DIFF} -in1 ${DATA_DIR_TOPP}/PepNovo_PTMs_.txt -in2 ${DATA_DIR_TOPP}/PepNovo_PTMs.txt)
#add_test("TOPP_PepNovoAdapter_3" ${TOPP_BIN_PATH}/PepNovoAdapter -ini ${DATA_DIR_TOPP}/PepNovoAdapter_5_parameters.ini -in ${DATA_DIR_TOPP}/PepNovoAdapter_5_output.pepnovo_out -out PepNovoAdapter_5_output.tmp -pepnovo_out -dta_list ${DATA_DIR_TOPP}/tmp/dta_list.txt -model_directory ${DATA_DIR_TOPP}/tmp/ -temp_data_directory ${DATA_DIR_TOPP}/tmp/ -modifications_xml_file ${DATA_DIR_TOPP}/PepNovo_PTMs.xml -mz_files ${DATA_DIR_TOPP}/PepNovo.mzXML)
#add_test("TOPP_PepNovoAdapter_3_out1" ${DIFF} -whitelist "?xml-stylesheet" "date_group_1" -in1 PepNovoAdapter_5_output.tmp -in2 ${DATA_DIR_TOPP}/PepNovoAdapter_5_output.idXML)

#------------------------------------------------------------------------------
if (NOT (${SIRIUS_BINARY} STREQUAL "SIRIUS_BINARY-NOTFOUND"))
  add_test("TOPP_SiriusAdapter_1" ${TOPP_BIN_PATH}/SiriusAdapter -test -executable "${SIRIUS_BINARY}" -in ${DATA_DIR_TOPP}/THIRDPARTY/SiriusAdapter_1_input.mzML -out_sirius SiriusAdapter_1_output.tmp -auto_charge -profile qtof -database all)
  add_test("TOPP_SiriusAdapter_1_out" ${DIFF} -in1 SiriusAdapter_1_output.tmp -in2 ${DATA_DIR_TOPP}/THIRDPARTY/SiriusAdapter_1_output.mzTab -whitelist "MTD")
  set_tests_properties("TOPP_SiriusAdapter_1_out" PROPERTIES DEPENDS "TOPP_SiriusAdapter_1")
  # Note that with FingerID, only one spectrum (number 4) should produce a result
  if (ENABLE_FINGERID_TEST)
  add_test("TOPP_SiriusAdapter_2" ${TOPP_BIN_PATH}/SiriusAdapter -test -executable "${SIRIUS_BINARY}" -in ${DATA_DIR_TOPP}/THIRDPARTY/SiriusAdapter_1_input.mzML -out_sirius SiriusAdapter_2_output.tmp -out_fingerid SiriusAdapter_2_foutput.tmp -auto_charge -profile qtof -database all)
  add_test("TOPP_SiriusAdapter_2_out" ${DIFF} -in1 SiriusAdapter_2_foutput.tmp -in2 ${DATA_DIR_TOPP}/THIRDPARTY/SiriusAdapter_2_foutput.mzTab -whitelist "MTD")
  set_tests_properties("TOPP_SiriusAdapter_2_out" PROPERTIES DEPENDS "TOPP_SiriusAdapter_2")
  endif()
endif()

# made library with spectrast -cNtestLib -cP0.0 CometAdapter_1_out.pep.xml 
#------------------------------------------------------------------------------
if (NOT (${SPECTRAST_BINARY} STREQUAL "SPECTRAST_BINARY-NOTFOUND") AND FALSE)
  add_test("TOPP_SpectrastSearchAdapter_0_prepare" ${TOPP_BIN_PATH}/FileConverter -test -force_TPP_compatibility -in ${DATA_DIR_TOPP}/THIRDPARTY/spectra_spectrast.mzXML -out SpectrastAdapter_1_hack.mzML)
  add_test("TOPP_SpectrastSearchAdapter_1" ${TOPP_BIN_PATH}/SpectraSTSearchAdapter -test -library_file ${DATA_DIR_TOPP}/THIRDPARTY/testLib.splib -spectra_files SpectrastAdapter_1_hack.mzML -output_files SpectrastAdapter_1_out1.tmp.pep.xml -executable "${SPECTRAST_BINARY}")
  add_test("TOPP_SpectrastSearchAdapter_1_out" ${DIFF} -in1 SpectrastAdapter_1_out1.tmp.pep.xml -in2 ${DATA_DIR_TOPP}/THIRDPARTY/SpectrastAdapter_1_output.pep.xml -whitelist "msms_pipeline_analysis date" "?xml-stylesheet" "summary base_name")
  set_tests_properties("TOPP_SpectrastSearchAdapter_1" PROPERTIES DEPENDS "TOPP_SpectrastSearchAdapter_0_prepare")
  set_tests_properties("TOPP_SpectrastSearchAdapter_1_out" PROPERTIES DEPENDS "TOPP_SpectrastSearchAdapter_1")
endif()

