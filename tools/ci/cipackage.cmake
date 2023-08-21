##
## this script REQUIRES the previous execution of the cibuild.cmake script!!
##

# define build name&co for easier identification on CDash
set(CTEST_BUILD_NAME "$ENV{BUILD_NAME}_Package")
set(CTEST_SITE "$ENV{CI_PROVIDER}")
set(CTEST_SOURCE_DIRECTORY "$ENV{SOURCE_DIRECTORY}")
set(CTEST_BINARY_DIRECTORY "${CTEST_SOURCE_DIRECTORY}/bld")

# cdash server (fu-berlin) SSL certificate sometimes is revoked. Keeps CI running.
set (CTEST_CURL_OPTIONS       CURLOPT_SSL_VERIFYHOST_OFF CURLOPT_SSL_VERIFYPEER_OFF )

message(STATUS "CTEST_SOURCE_DIRECTORY: ${CTEST_SOURCE_DIRECTORY}")
message(STATUS "CTEST_BINARY_DIRECTORY: ${CTEST_BINARY_DIRECTORY}")

ctest_start(Nightly GROUP Package)

# build the dist target for packages
if("$ENV{ENABLE_STYLE_TESTING}" STREQUAL "OFF")
  if(NOT "$ENV{PYOPENMS}" STREQUAL "ON")
    ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" TARGET "dist" NUMBER_ERRORS _build_errors)
    ctest_submit(PARTS Build)
  endif()
else()
  set(_build_errors 0)
endif()

if (_build_errors)
  message(FATAL_ERROR "Please check the build results at: https://cdash.openms.de/index.php?project=OpenMS&begin=2023-01-01&end=2030-01-01&filtercount=1&field1=buildname&compare1=63&value1=${CTEST_BUILD_NAME}")
else()
  message("Please check the build results at: https://cdash.openms.de/index.php?project=OpenMS&begin=2023-01-01&end=2030-01-01&filtercount=1&field1=buildname&compare1=63&value1=${CTEST_BUILD_NAME}")
endif()
