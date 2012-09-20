## DB test settings
set(DB_TEST OFF CACHE BOOL "If true, the DB tests are enabled.")
set(DB_TEST_HOST "localhost" CACHE STRING "Test database server name (only used if DB_TEST is true).")
set(DB_TEST_PORT "3307" CACHE STRING "Test database server port (only used if DB_TEST is true).")
set(DB_TEST_DB "OPENMS_TEST_DB" CACHE STRING "Test database name (only used for DB_TEST).")
set(DB_TEST_USER "openms_test_user" CACHE STRING "Test database user name (only used for DB_TEST).")
set(DB_TEST_PW "openms_test_password" CACHE STRING "Test database user password (only used for DB_TEST).")
if (DB_TEST)
	# output
	message(STATUS "DB testing enabled")
	message(STATUS "DB testing - creating credentials files ...")
	# create assorted credentials files
	configure_file(${PROJECT_SOURCE_DIR}/source/TEST/DB_credentials.txt.in ${PROJECT_BINARY_DIR}/source/TEST/DB_credentials.txt)
	configure_file(${PROJECT_SOURCE_DIR}/source/TEST/TOPP/DBImporter.ini.in ${PROJECT_BINARY_DIR}/source/TEST/TOPP/DBImporter.ini)
	configure_file(${PROJECT_SOURCE_DIR}/source/TEST/TOPP/DBExporter.ini.in ${PROJECT_BINARY_DIR}/source/TEST/TOPP/DBExporter.ini)
endif()


## MS2 Search-Engines go here...
## MACRO OPENMS_FINDBINARY:
## fills ${varname} with the path to the binary given in ${binaryname}
## @param varname      Name of the variable which will hold the result string (e.g. OMSSA_BINARY)
## @param binaryname   List of binary names which are searched
## @param name         Human readable version of binaryname for messages
MACRO (OPENMS_FINDBINARY varname binaryname name)
  find_program(${varname} ${binaryname} PATHS ENV PATH)
  if (${${varname}} STREQUAL "${varname}-NOTFOUND")
    message(STATUS "  - ${name} not found.")
  else()
    message(STATUS "  + ${name} binary found at ${binaryname}. Enabling corresponding tests.")
  endif()  
ENDMACRO (OPENMS_FINDBINARY)

message(STATUS "Searching for MS2 search engines ...")

## OMSSA
OPENMS_FINDBINARY(OMSSA_BINARY "omssacl" "OMSSA")

## X!Tandem
OPENMS_FINDBINARY(XTANDEM_BINARY "tandem" "X!Tandem")

## MyriMatch
OPENMS_FINDBINARY(MYRIMATCH_BINARY "myrimatch" "Myrimatch")
