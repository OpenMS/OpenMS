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

# Define patterns that should NOT be classified as warnings
set(CTEST_CUSTOM_WARNING_EXCEPTION
    # Ignore the generic "non-zero return value" warnings - let the actual errors be classified properly
    ".*WARNING.*non-zero return value.*cmake.*"
    )

# Define patterns that should be classified as errors (these override warning classification)
set(CTEST_CUSTOM_ERROR_MATCH
    "subprocess-exited-with-error"
    "FAILED:.*pyopenms"
    "× Getting requirements to build wheel did not run successfully"
    "error: subprocess-exited-with-error"
    ".*FAILED:.*pyOpenMS.*"
    ".*pip wheel.*"
    "Exception:.*"
    )

ctest_start(APPEND)

## run tests
ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}" PARALLEL_LEVEL 3 RETURN_VALUE _test_errors)

## send test results to CDash
ctest_submit(PARTS Test Done CAPTURE_CMAKE_ERROR _submit_result)

if(NOT _submit_result EQUAL 0)
  execute_process(COMMAND ${CMAKE_COMMAND} -E echo "::warning file=citest.cmake,line=37::CTest submission failed, CDASH server is not available. Continuing execution.")
  message(WARNING "CTest submission failed, no detailed logs will be available.")
  if (_test_errors)
    message(FATAL_ERROR "There were errors: aborting")
  endif()
else()
  string(REPLACE "+" "%2B" BUILD_NAME_SAFE ${CTEST_BUILD_NAME})
  string(REPLACE "." "%2E" BUILD_NAME_SAFE ${BUILD_NAME_SAFE})
  string(REPLACE "/" "%2F" BUILD_NAME_SAFE ${BUILD_NAME_SAFE})
  if (_test_errors)
  message(FATAL_ERROR "There were errors: Please check the test results at: https://cdash.seqan.de/index.php?project=OpenMS&begin=2023-01-01&end=2030-01-01&filtercount=1&field1=buildname&compare1=63&value1=${BUILD_NAME_SAFE}")
  else()
  message("Testing successful: Please check the test results at: https://cdash.seqan.de/index.php?project=OpenMS&begin=2023-01-01&end=2030-01-01&filtercount=1&field1=buildname&compare1=63&value1=${BUILD_NAME_SAFE}")
  endif()
endif()

