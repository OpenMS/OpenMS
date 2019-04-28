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
# MaRaCluster
OPENMS_FINDBINARY(MARACLUSTER_BINARY "maracluster" "MaRaCluster")

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
# MSFragger
OPENMS_FINDBINARY(MSFRAGGER_BINARY "MSFragger.jar" "MSFragger")

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
# Novor
OPENMS_FINDBINARY(NOVOR_BINARY "novor.jar" "Novor")

#------------------------------------------------------------------------------
# Crux
OPENMS_FINDBINARY(CRUX_BINARY "crux;crux.exe" "Crux")

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
if (NOT (${CRUX_BINARY} STREQUAL "CRUX_BINARY-NOTFOUND"))
  ### NOT needs to be added after the binarys have been included
  add_test("TOPP_CruxAdapter_1" ${TOPP_BIN_PATH}/CruxAdapter -test -ini ${DATA_DIR_TOPP}/THIRDPARTY/CruxAdapter_1.ini -database ${DATA_DIR_TOPP}/THIRDPARTY/proteins.fasta -in ${DATA_DIR_TOPP}/THIRDPARTY/spectra_comet.mzML -out CruxAdapter_1_out1.tmp -crux_executable "${CRUX_BINARY}" -run_percolator false -decoy-format peptide-reverse)
  add_test("TOPP_CruxAdapter_1_out1" ${DIFF} -in1 CruxAdapter_1_out1.tmp -in2 ${DATA_DIR_TOPP}/THIRDPARTY/CruxAdapter_1_out.idXML -whitelist "IdentificationRun date" "SearchParameters id=\"SP_0\" db=" "UserParam type=\"stringList\" name=\"spectra_data\" value=")
  set_tests_properties("TOPP_CruxAdapter_1_out1" PROPERTIES DEPENDS "TOPP_CruxAdapter_1")
endif()

#------------------------------------------------------------------------------
if (NOT (${COMET_BINARY} STREQUAL "COMET_BINARY-NOTFOUND"))
  ### NOT needs to be added after the binarys have been included
  add_test("TOPP_CometAdapter_1" ${TOPP_BIN_PATH}/CometAdapter -test -ini ${DATA_DIR_TOPP}/THIRDPARTY/CometAdapter_1.ini -database ${DATA_DIR_TOPP}/THIRDPARTY/proteins.fasta -in ${DATA_DIR_TOPP}/THIRDPARTY/spectra_comet.mzML -out CometAdapter_1_out1.tmp -pin_out CometAdapter_1_out2.tmp.csv -comet_executable "${COMET_BINARY}")
  add_test("TOPP_CometAdapter_1_out1" ${DIFF} -in1 CometAdapter_1_out1.tmp -in2 ${DATA_DIR_TOPP}/THIRDPARTY/CometAdapter_1_out.idXML -whitelist "IdentificationRun date" "SearchParameters id=\"SP_0\" db=" "UserParam type=\"stringList\" name=\"spectra_data\" value=")
  set_tests_properties("TOPP_CometAdapter_1_out1" PROPERTIES DEPENDS "TOPP_CometAdapter_1")
  ### Second test for optional pin file needs to be added, not sure how to do FuzzyDiff on the csv style pin file, whitelisting the first id column
  add_test("TOPP_CometAdapter_2_prepare" ${TOPP_BIN_PATH}/FileConverter -test -in ${DATA_DIR_TOPP}/THIRDPARTY/CometAdapter_2_in.mzML -out CometAdapter_2_prepared.mzML -force_TPP_compatibility)
  add_test("TOPP_CometAdapter_2" ${TOPP_BIN_PATH}/CometAdapter -test -database ${DATA_DIR_TOPP}/THIRDPARTY/CometAdapter_2_in.fasta -in CometAdapter_2_prepared.mzML -out CometAdapter_2_out1.tmp -pin_out CometAdapter_2_out2.tmp.csv -comet_executable "${COMET_BINARY}" -precursor_mass_tolerance 3 -precursor_error_units Da -ini ${DATA_DIR_TOPP}/THIRDPARTY/CometAdapter_1.ini)
  add_test("TOPP_CometAdapter_2_out1" ${DIFF} -in1 CometAdapter_2_out1.tmp -in2 ${DATA_DIR_TOPP}/THIRDPARTY/CometAdapter_2_out.idXML -whitelist "IdentificationRun date" "SearchParameters id=\"SP_0\" db=" "UserParam type=\"stringList\" name=\"spectra_data\" value=")
  set_tests_properties("TOPP_CometAdapter_2" PROPERTIES DEPENDS "TOPP_CometAdapter_2_prepare")
  set_tests_properties("TOPP_CometAdapter_2_out1" PROPERTIES DEPENDS "TOPP_CometAdapter_2")
  ### Second test for optional pin file needs to be added, not sure how to do FuzzyDiff on the csv style pin file, whitelisting the first id column
  add_test("TOPP_CometAdapter_3" ${TOPP_BIN_PATH}/CometAdapter -test -ini ${DATA_DIR_TOPP}/THIRDPARTY/CometAdapter_3.ini -database ${DATA_DIR_TOPP}/THIRDPARTY/CometAdapter_3.fasta -in ${DATA_DIR_TOPP}/THIRDPARTY/CometAdapter_3.mzML -out CometAdapter_3_out1.tmp -pin_out CometAdapter_3_out2.tmp.csv -comet_executable "${COMET_BINARY}")
  add_test("TOPP_CometAdapter_3_out1" ${DIFF} -in1 CometAdapter_3_out1.tmp -in2 ${DATA_DIR_TOPP}/THIRDPARTY/CometAdapter_3_out.idXML -whitelist "IdentificationRun date" "SearchParameters id=\"SP_0\" db=" "UserParam type=\"stringList\" name=\"spectra_data\" value=")
  set_tests_properties("TOPP_CometAdapter_3_out1" PROPERTIES DEPENDS "TOPP_CometAdapter_3")
endif()

#------------------------------------------------------------------------------
if (NOT (${MARACLUSTER_BINARY} STREQUAL "MARACLUSTER_BINARY-NOTFOUND"))
  ### NOT needs to be added after the binarys have been included
  add_test("TOPP_MaRaClusterAdapter_1" ${TOPP_BIN_PATH}/MaRaClusterAdapter -test -ini ${DATA_DIR_TOPP}/THIRDPARTY/MaRaClusterAdapter_1.ini -in ${DATA_DIR_TOPP}/THIRDPARTY/MaRaClusterAdapter_1_in_1.mzML ${DATA_DIR_TOPP}/THIRDPARTY/MaRaClusterAdapter_1_in_2.mzML -consensus_out MaRaClusterAdapter_1_out_1.tmp.mzML -maracluster_executable "${MARACLUSTER_BINARY}")
  add_test("TOPP_MaRaClusterAdapter_1_out_1" ${DIFF} -in1 MaRaClusterAdapter_1_out_1.tmp.part1.mzML -in2 ${DATA_DIR_TOPP}/THIRDPARTY/MaRaClusterAdapter_1_out_1.part1.mzML -whitelist "IdentificationRun date" "SearchParameters id=\"SP_0\" db=" "UserParam type=\"stringList\" name=\"spectra_data\" value=" "sourceFile id=" "fileChecksum" "cvParam cvRef=\"MS\" accession=\"MS:1000569\" name=\"SHA-1\"" "software id=\"MaRaCluster\" version=" "cvParam cvRef=\"MS\" accession=\"MS:1000747\" name=\"completion time\"" "cv id=\"MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Ontology\" version=" "software id=\"pwiz_3.0" "processingMethod order=\"0\" softwareRef=")
  set_tests_properties("TOPP_MaRaClusterAdapter_1_out_1" PROPERTIES DEPENDS "TOPP_MaRaClusterAdapter_1")
  add_test("TOPP_MaRaClusterAdapter_2" ${TOPP_BIN_PATH}/MaRaClusterAdapter -test -ini ${DATA_DIR_TOPP}/THIRDPARTY/MaRaClusterAdapter_2.ini -in ${DATA_DIR_TOPP}/THIRDPARTY/MaRaClusterAdapter_1_in_1.mzML ${DATA_DIR_TOPP}/THIRDPARTY/MaRaClusterAdapter_1_in_2.mzML -id_in ${DATA_DIR_TOPP}/THIRDPARTY/MaRaClusterAdapter_1_in_3.idXML -out MaRaClusterAdapter_2_out_1.tmp.idXML -maracluster_executable "${MARACLUSTER_BINARY}")
  add_test("TOPP_MaRaClusterAdapter_2_out_1" ${DIFF} -in1 MaRaClusterAdapter_2_out_1.tmp.idXML -in2 ${DATA_DIR_TOPP}/THIRDPARTY/MaRaClusterAdapter_2_out_1.idXML -whitelist "IdentificationRun date" "SearchParameters id=\"SP_0\" db=" "UserParam type=\"stringList\" name=\"spectra_data\" value=" "UserParam type=\"string\" name=\"file_origin\" value=")
  set_tests_properties("TOPP_MaRaClusterAdapter_2_out_1" PROPERTIES DEPENDS "TOPP_MaRaClusterAdapter_2")
endif()

#------------------------------------------------------------------------------
if (NOT (${PERCOLATOR_BINARY} STREQUAL "PERCOLATOR_BINARY-NOTFOUND"))
  ### NOT needs to be added after the binarys have been included
  add_test("TOPP_PercolatorAdapter_1" ${TOPP_BIN_PATH}/PercolatorAdapter -test -ini ${DATA_DIR_TOPP}/THIRDPARTY/PercolatorAdapter_1.ini -in ${DATA_DIR_TOPP}/THIRDPARTY/PercolatorAdapter_1.idXML -out PercolatorAdapter_1_out1.tmp -out_type idXML -percolator_executable "${PERCOLATOR_BINARY}")
  add_test("TOPP_PercolatorAdapter_1_out1" ${DIFF} -in1 PercolatorAdapter_1_out1.tmp -in2 ${DATA_DIR_TOPP}/THIRDPARTY/PercolatorAdapter_1_out.idXML -whitelist "IdentificationRun date" "SearchParameters id=\"SP_0\" db=" "UserParam type=\"stringList\" name=\"spectra_data\" value=")
  set_tests_properties("TOPP_PercolatorAdapter_1_out1" PROPERTIES DEPENDS "TOPP_PercolatorAdapter_1")
  add_test("TOPP_PercolatorAdapter_2" ${TOPP_BIN_PATH}/PercolatorAdapter -test -osw_level ms1 -in_osw ${DATA_DIR_TOPP}/THIRDPARTY/PercolatorAdapter_2.osw -out PercolatorAdapter_2_out1.osw -out_type osw -percolator_executable "${PERCOLATOR_BINARY}")
  add_test("TOPP_PercolatorAdapter_3" ${TOPP_BIN_PATH}/PercolatorAdapter -test -osw_level ms2 -in_osw PercolatorAdapter_2_out1.osw -out PercolatorAdapter_3_out1.osw -out_type osw -percolator_executable "${PERCOLATOR_BINARY}")
  set_tests_properties("TOPP_PercolatorAdapter_3" PROPERTIES DEPENDS "TOPP_PercolatorAdapter_2")
  add_test("TOPP_PercolatorAdapter_4" ${TOPP_BIN_PATH}/PercolatorAdapter -test -osw_level transition -in_osw PercolatorAdapter_3_out1.osw -out PercolatorAdapter_4_out1.osw -out_type osw -percolator_executable "${PERCOLATOR_BINARY}")
  set_tests_properties("TOPP_PercolatorAdapter_4" PROPERTIES DEPENDS "TOPP_PercolatorAdapter_3")
  ### TOPP_PercolatorAdapter_2-4 do not validate output, but checks whether OSW files can be read and written to.
endif()

#------------------------------------------------------------------------------
option(WITH_MASCOT_TEST "Runs the Mascot Online test (do not turn this on unless you know what you are doing)" OFF)
if (WITH_MASCOT_TEST)
  add_test("TOPP_MascotAdapterOnline_1" ${TOPP_BIN_PATH}/MascotAdapterOnline -test -ini ${DATA_DIR_TOPP}/THIRDPARTY/MascotAdapterOnline_1.ini -Mascot_parameters:database SwissProt -in ${DATA_DIR_TOPP}/THIRDPARTY/spectra_comet.mzML -out MascotAdapterOnline_1_out1.tmp)
  add_test("TOPP_MascotAdapterOnline_1_out1" ${DIFF} -in1 MascotAdapterOnline_1_out1.tmp -in2 ${DATA_DIR_TOPP}/THIRDPARTY/MascotAdapterOnline_1_out.idXML -whitelist "IdentificationRun date" "UserParam type=\"string\" name=\"SearchNumber\" value=" "db=\"SwissProt\" db_version=")
  set_tests_properties("TOPP_MascotAdapterOnline_1_out1" PROPERTIES DEPENDS "TOPP_MascotAdapterOnline_1")
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

#------------------------------------------------------------------------------
# MSFragger
option(WITH_MSFRAGGER_TEST "Runs the MSFragger test (which needs at least 2.6GB free RAM during testing)" ON)
if (WITH_MSFRAGGER_TEST)
  if (NOT (${MSFRAGGER_BINARY} STREQUAL "MSFRAGGER_BINARY-NOTFOUND"))
    add_test("TOPP_MSFraggerAdapter_7" ${TOPP_BIN_PATH}/MSFraggerAdapter -test -java_heapmemory 2600  -in ${DATA_DIR_TOPP}/THIRDPARTY/spectra.mzML -executable "${MSFRAGGER_BINARY}" -database ${DATA_DIR_TOPP}/THIRDPARTY/proteins.fasta -out MSFraggerAdapter_7_out_tmp.idXML -opt_out MSFraggerAdapter_7_opt_out_tmp.pepXML -varmod:enable_common -digest:num_enzyme_termini semi)
    add_test("TOPP_MSFraggerAdapter_7_out" ${DIFF} -in1 MSFraggerAdapter_7_out_tmp.idXML -in2 ${DATA_DIR_TOPP}/THIRDPARTY/MSFraggerAdapter_7_out.idXML -whitelist "date" "search_database" "db") # Because MSFragger links the search database in a temporary directory
    add_test("TOPP_MSFraggerAdapter_7_opt_out" ${DIFF} -in1 MSFraggerAdapter_7_opt_out_tmp.pepXML -in2 ${DATA_DIR_TOPP}/THIRDPARTY/MSFraggerAdapter_7_opt_out.pepXML -whitelist "date" "search_database" "db")
    set_tests_properties("TOPP_MSFraggerAdapter_7_out" PROPERTIES DEPENDS "TOPP_MSFraggerAdapter_7")
    set_tests_properties("TOPP_MSFraggerAdapter_7_opt_out" PROPERTIES DEPENDS "TOPP_MSFraggerAdapter_7")

    add_test("TOPP_MSFraggerAdapter_8" ${TOPP_BIN_PATH}/MSFraggerAdapter -test -java_heapmemory 2600 -in ${DATA_DIR_TOPP}/THIRDPARTY/spectra_comet.mzML -executable "${MSFRAGGER_BINARY}" -database ${DATA_DIR_TOPP}/THIRDPARTY/proteins.fasta -out MSFraggerAdapter_8_out_tmp.idXML -varmod:enable_common -digest:num_enzyme_termini semi)
    add_test("TOPP_MSFraggerAdapter_8_out" ${DIFF} -in1 MSFraggerAdapter_8_out_tmp.idXML -in2 ${DATA_DIR_TOPP}/THIRDPARTY/MSFraggerAdapter_8_out.idXML -whitelist "date" "search_database" "db") # Because MSFragger links the search database in a temporary directory
    set_tests_properties("TOPP_MSFraggerAdapter_8_out" PROPERTIES DEPENDS "TOPP_MSFraggerAdapter_8")
  endif()
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
  # Note: Following test are performed without adduct/id information, since these are obtained by the MetaboliteAdductDecharger/AccurateMassSearch
  
  if (ENABLE_SIRIUS_TEST)
  # add dependencies for one test at a time to reduce memory and cpu consumption
  # test mzMl as input
  add_test("TOPP_SiriusAdapter_1" ${TOPP_BIN_PATH}/SiriusAdapter -test -executable "${SIRIUS_BINARY}" -in ${DATA_DIR_TOPP}/THIRDPARTY/SiriusAdapter_1_input.mzML -out_sirius SiriusAdapter_1_output.tmp -sirius:auto_charge -sirius:profile qtof -sirius:database all -sirius:cores 1)
  add_test("TOPP_SiriusAdapter_1_out" ${DIFF} -in1 SiriusAdapter_1_output.tmp -in2 ${DATA_DIR_TOPP}/THIRDPARTY/SiriusAdapter_1_output.mzTab -whitelist "MTD")
  set_tests_properties("TOPP_SiriusAdapter_1_out" PROPERTIES DEPENDS "TOPP_SiriusAdapter_1")
  # test mzML and featureXML with feature_only
  add_test("TOPP_SiriusAdapter_2" ${TOPP_BIN_PATH}/SiriusAdapter -test -executable "${SIRIUS_BINARY}" -in ${DATA_DIR_TOPP}/THIRDPARTY/SiriusAdapter_2_input.mzML -in_featureinfo ${DATA_DIR_TOPP}/THIRDPARTY/SiriusAdapter_2_input.featureXML -out_sirius SiriusAdapter_2_output.tmp -preprocessing:feature_only -preprocessing:filter_by_num_masstraces 3 -sirius:auto_charge -sirius:profile qtof -sirius:database all -sirius:cores 1)
  add_test("TOPP_SiriusAdapter_2_out" ${DIFF} -in1 SiriusAdapter_2_output.tmp -in2 ${DATA_DIR_TOPP}/THIRDPARTY/SiriusAdapter_2_output.mzTab -whitelist "MTD")
  set_tests_properties("TOPP_SiriusAdapter_2" PROPERTIES DEPENDS "TOPP_SiriusAdapter_1")
  set_tests_properties("TOPP_SiriusAdapter_2_out" PROPERTIES DEPENDS "TOPP_SiriusAdapter_2")
  # test mzML and featureXML without feature_only (filter_by_num_masstraces should automatically set to 1)
  add_test("TOPP_SiriusAdapter_3" ${TOPP_BIN_PATH}/SiriusAdapter -test -executable "${SIRIUS_BINARY}" -in ${DATA_DIR_TOPP}/THIRDPARTY/SiriusAdapter_3_input.mzML -in_featureinfo ${DATA_DIR_TOPP}/THIRDPARTY/SiriusAdapter_3_input.featureXML -out_sirius SiriusAdapter_3_output.tmp -preprocessing:filter_by_num_masstraces 3 -sirius:auto_charge -sirius:profile qtof -sirius:database all -sirius:cores 1)
  add_test("TOPP_SiriusAdapter_3_out" ${DIFF} -in1 SiriusAdapter_3_output.tmp -in2 ${DATA_DIR_TOPP}/THIRDPARTY/SiriusAdapter_3_output.mzTab -whitelist "MTD")
  set_tests_properties("TOPP_SiriusAdapter_3" PROPERTIES DEPENDS "TOPP_SiriusAdapter_2")
  set_tests_properties("TOPP_SiriusAdapter_3_out" PROPERTIES DEPENDS "TOPP_SiriusAdapter_3")
  # test internal .ms output (converter mode) 
  add_test("TOPP_SiriusAdapter_5" ${TOPP_BIN_PATH}/SiriusAdapter -test -executable "${SIRIUS_BINARY}" -in ${DATA_DIR_TOPP}/THIRDPARTY/SiriusAdapter_3_input.mzML -in_featureinfo ${DATA_DIR_TOPP}/THIRDPARTY/SiriusAdapter_3_input.featureXML -out_ms SiriusAdapter_5_output.tmp -converter_mode) 
  add_test("TOPP_SiriusAdapter_5_out" ${DIFF} -in1 SiriusAdapter_5_output.tmp -in2 ${DATA_DIR_TOPP}/THIRDPARTY/SiriusAdapter_5_output.ms)
  set_tests_properties("TOPP_SiriusAdapter_5" PROPERTIES DEPENDS "TOPP_SiriusAdapter_3")
  set_tests_properties("TOPP_SiriusAdapter_5_out" PROPERTIES DEPENDS "TOPP_SiriusAdapter_5")
  # test internal .ms output negative 
  add_test("TOPP_SiriusAdapter_6" ${TOPP_BIN_PATH}/SiriusAdapter -test -executable "${SIRIUS_BINARY}" -in ${DATA_DIR_TOPP}/THIRDPARTY/SiriusAdapter_4_input.mzML -in_featureinfo ${DATA_DIR_TOPP}/THIRDPARTY/SiriusAdapter_4_input.featureXML -out_ms SiriusAdapter_6_output.tmp -converter_mode) 
  add_test("TOPP_SiriusAdapter_6_out" ${DIFF} -in1 SiriusAdapter_6_output.tmp -in2 ${DATA_DIR_TOPP}/THIRDPARTY/SiriusAdapter_6_output.ms)
  set_tests_properties("TOPP_SiriusAdapter_6" PROPERTIES DEPENDS "TOPP_SiriusAdapter_5")
  set_tests_properties("TOPP_SiriusAdapter_6_out" PROPERTIES DEPENDS "TOPP_SiriusAdapter_6")
  # test mzML and featureXML negative
  add_test("TOPP_SiriusAdapter_7" ${TOPP_BIN_PATH}/SiriusAdapter -test -executable "${SIRIUS_BINARY}" -in ${DATA_DIR_TOPP}/THIRDPARTY/SiriusAdapter_4_input.mzML -in_featureinfo ${DATA_DIR_TOPP}/THIRDPARTY/SiriusAdapter_4_input.featureXML -out_sirius SiriusAdapter_7_output.tmp -preprocessing:feature_only -sirius:profile qtof -sirius:database all -sirius:cores 1)
  add_test("TOPP_SiriusAdapter_7_out" ${DIFF} -in1 SiriusAdapter_7_output.tmp -in2 ${DATA_DIR_TOPP}/THIRDPARTY/SiriusAdapter_7_output.mzTab -whitelist "MTD")
  set_tests_properties("TOPP_SiriusAdapter_7" PROPERTIES DEPENDS "TOPP_SiriusAdapter_6")
  set_tests_properties("TOPP_SiriusAdapter_7_out" PROPERTIES DEPENDS "TOPP_SiriusAdapter_7")
 
  # use AccurateMassSearch data
  add_test("UTILS_AssayGeneratorMetabo_7" ${TOPP_BIN_PATH}/AssayGeneratorMetabo -test -executable "${SIRIUS_BINARY}" -in ${DATA_DIR_TOPP}/AssayGeneratorMetabo_input.mzML -in_id ${DATA_DIR_TOPP}/AssayGeneratorMetabo_ams_input.featureXML -out AssayGeneratorMetabo_ams_sirius_output.tmp.tsv -fragment_annotation sirius -use_exact_mass -transition_threshold 3.0 -min_transitions 2 -max_transitions 3 -preprocessing:filter_by_num_masstraces 1 -preprocessing:precursor_mz_tolerance 10 -preprocessing:precursor_mz_tolerance_unit ppm -preprocessing:feature_only -sirius:profile qtof -sirius:compound_timeout 100)
  add_test("UTILS_AssayGeneratorMetabo_7_out1" ${DIFF} -in1 AssayGeneratorMetabo_ams_sirius_output.tmp.tsv -in2 ${DATA_DIR_TOPP}/AssayGeneratorMetabo_ams_sirius_output.tsv)
  set_tests_properties("UTILS_AssayGeneratorMetabo_7_out1" PROPERTIES DEPENDS "UTILS_AssayGeneratorMetabo_7")

  # use AccurateMassSearch data
  add_test("UTILS_AssayGeneratorMetabo_8" ${TOPP_BIN_PATH}/AssayGeneratorMetabo -test -executable "${SIRIUS_BINARY}" -in ${DATA_DIR_TOPP}/AssayGeneratorMetabo_input.mzML -in_id ${DATA_DIR_TOPP}/AssayGeneratorMetabo_ams_input.featureXML -out AssayGeneratorMetabo_ams_sirius_ukn_output.tmp.tsv -fragment_annotation sirius -use_exact_mass -transition_threshold 3.0 -min_transitions 2 -max_transitions 3 -preprocessing:filter_by_num_masstraces 1 -preprocessing:precursor_mz_tolerance 10 -preprocessing:precursor_mz_tolerance_unit ppm -preprocessing:feature_only -sirius:profile qtof -sirius:compound_timeout 100 -use_known_unknowns)
  add_test("UTILS_AssayGeneratorMetabo_8_out1" ${DIFF} -in1 AssayGeneratorMetabo_ams_sirius_ukn_output.tmp.tsv -in2 ${DATA_DIR_TOPP}/AssayGeneratorMetabo_ams_sirius_ukn_output.tsv)
  set_tests_properties("UTILS_AssayGeneratorMetabo_8" PROPERTIES DEPENDS "UTILS_AssayGeneratorMetabo_7")
  set_tests_properties("UTILS_AssayGeneratorMetabo_8_out1" PROPERTIES DEPENDS "UTILS_AssayGeneratorMetabo_8")

  # use AccurateMassSearch data
  add_test("UTILS_AssayGeneratorMetabo_9" ${TOPP_BIN_PATH}/AssayGeneratorMetabo -test -executable "${SIRIUS_BINARY}" -in ${DATA_DIR_TOPP}/AssayGeneratorMetabo_intsort_input.mzML -in_id ${DATA_DIR_TOPP}/AssayGeneratorMetabo_intsort_input.featureXML -out AssayGeneratorMetabo_ams_sirius_intsort_output.tmp.tsv -fragment_annotation sirius -use_exact_mass -transition_threshold 3.0 -min_transitions 2 -max_transitions 3 -preprocessing:filter_by_num_masstraces 1 -preprocessing:precursor_mz_tolerance 10 -preprocessing:precursor_mz_tolerance_unit ppm -preprocessing:feature_only -sirius:profile qtof -sirius:compound_timeout 100)
  add_test("UTILS_AssayGeneratorMetabo_9_out1" ${DIFF} -in1 AssayGeneratorMetabo_ams_sirius_intsort_output.tmp.tsv -in2 ${DATA_DIR_TOPP}/AssayGeneratorMetabo_ams_sirius_intsort_output.tsv)
  set_tests_properties("UTILS_AssayGeneratorMetabo_9" PROPERTIES DEPENDS "UTILS_AssayGeneratorMetabo_8")
  set_tests_properties("UTILS_AssayGeneratorMetabo_9_out1" PROPERTIES DEPENDS "UTILS_AssayGeneratorMetabo_9")
  endif()

  # Note that with FingerID, output for compound 79 without feature only
  if (ENABLE_FINGERID_TEST)
  add_test("TOPP_SiriusAdapter_4" ${TOPP_BIN_PATH}/SiriusAdapter -test -executable "${SIRIUS_BINARY}" -in ${DATA_DIR_TOPP}/THIRDPARTY/SiriusAdapter_2_input.mzML -in_featureinfo ${DATA_DIR_TOPP}/THIRDPARTY/SiriusAdapter_2_input.featureXML  -out_sirius SiriusAdapter_4_output.tmp -out_fingerid SiriusAdapter_4_foutput.tmp -sirius:auto_charge -sirius:profile qtof -sirius:database all)
  add_test("TOPP_SiriusAdapter_4_out" ${DIFF} -in1 SiriusAdapter_4_foutput.tmp -in2 ${DATA_DIR_TOPP}/THIRDPARTY/SiriusAdapter_4_foutput.mzTab -whitelist "MTD")
  set_tests_properties("TOPP_SiriusAdapter_4_out" PROPERTIES DEPENDS "TOPP_SiriusAdapter_4")
  endif()

endif()


#------------------------------------------------------------------------------
if (NOT (${NOVOR_BINARY} STREQUAL "NOVOR_BINARY-NOTFOUND"))
  add_test("TOPP_NovorAdapter_1" ${TOPP_BIN_PATH}/NovorAdapter -test -java_memory 512 -executable "${NOVOR_BINARY}" -in ${DATA_DIR_TOPP}/THIRDPARTY/NovorAdapter_in.mzML -out NovorAdapter_1_out.tmp -variable_modifications "Acetyl (K)" -fixed_modifications "Carbamidomethyl (C)" -forbiddenResidues "I")
  add_test("TOPP_NovorAdapter_1_out" ${DIFF} -in1 NovorAdapter_1_out.tmp -in2 ${DATA_DIR_TOPP}/THIRDPARTY/NovorAdapter_1_out.idXML -whitelist "IdentificationRun date")
  set_tests_properties("TOPP_NovorAdapter_1_out" PROPERTIES DEPENDS "TOPP_NovorAdapter_1")
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

