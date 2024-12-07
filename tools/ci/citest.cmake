##
## this script REQUIRES the previous execution of the cibuild.cmake script!!
##

# define build name&co for easier identification on CDash
set(CTEST_BUILD_NAME "$ENV{BUILD_NAME}")
set(CTEST_SITE "$ENV{CI_PROVIDER}")
set(CTEST_SOURCE_DIRECTORY "$ENV{SOURCE_DIRECTORY}")
set(CTEST_BINARY_DIRECTORY "${CTEST_SOURCE_DIRECTORY}/bld")

# cdash server (fu-berlin) SSL certificate sometimes is revoked. Keeps CI running.
set (CTEST_CURL_OPTIONS       CURLOPT_SSL_VERIFYHOST_OFF CURLOPT_SSL_VERIFYPEER_OFF )

message(STATUS "CTEST_SOURCE_DIRECTORY: ${CTEST_SOURCE_DIRECTORY}")
message(STATUS "CTEST_BINARY_DIRECTORY: ${CTEST_BINARY_DIRECTORY}")

# ignore failing GzipIfstream_test which seems to be related to the used
# zlib version
set(CTEST_CUSTOM_TESTS_IGNORE
	GzipIfstream_test
)

# customize reporting of errors in CDash
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_ERRORS 1000)
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS 1000)

ctest_start(APPEND)

## run tests
## for pyopenms build, only run pyopenms tests
if("$ENV{PYOPENMS}" STREQUAL "ON")
  ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}" INCLUDE "pyopenms" PARALLEL_LEVEL 3 RETURN_VALUE _test_errors)
else()
  ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}" PARALLEL_LEVEL 3 RETURN_VALUE _test_errors)
endif()
## send test results to CDash
ctest_submit(PARTS Test Done)

string(REPLACE "+" "%2B" BUILD_NAME_SAFE ${CTEST_BUILD_NAME})
string(REPLACE "." "%2E" BUILD_NAME_SAFE ${BUILD_NAME_SAFE})
string(REPLACE "/" "%2F" BUILD_NAME_SAFE ${BUILD_NAME_SAFE})

if (_test_errors)
  message(FATAL_ERROR "There were errors: Please check the test results at: https://cdash.seqan.de/index.php?project=OpenMS&begin=2023-01-01&end=2030-01-01&filtercount=1&field1=buildname&compare1=63&value1=${BUILD_NAME_SAFE}")
else()
  message("Testing successful: Please check the test results at: https://cdash.seqan.de/index.php?project=OpenMS&begin=2023-01-01&end=2030-01-01&filtercount=1&field1=buildname&compare1=63&value1=${BUILD_NAME_SAFE}")
endif()
