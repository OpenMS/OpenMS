##
## this script REQUIRES the previous execution of the cibuild.cmake script!!
##

# TODO (in all cmake scripts) explicitly require certain variables to be in the ENV (e.g. by using a macro). Otherwise, you might have random hard-to-debug failures if not set.

# define build name&co for easier identification on CDash
set(CTEST_BUILD_NAME "$ENV{BUILD_NAME}_Package")
set(CTEST_SITE "$ENV{CI_PROVIDER}")
set(CTEST_SOURCE_DIRECTORY "$ENV{SOURCE_DIRECTORY}")
set(CTEST_BINARY_DIRECTORY "${CTEST_SOURCE_DIRECTORY}/bld")
set(CTEST_CMAKE_GENERATOR "Ninja")

set(CONFIGURE_OPTIONS)
set(VARS_TO_LOAD PACKAGE_TYPE SEARCH_ENGINES_DIRECTORY SIGNING_IDENTITY SIGNING_EMAIL)
foreach(VAR ${VARS_TO_LOAD})
  if(DEFINED ENV{${VAR}})
    list(APPEND CONFIGURE_OPTIONS "-D${VAR}=$ENV{${VAR}}")
  endif()
endforeach()

if(APPLE)
  list(APPEND CONFIGURE_OPTIONS "-DCPACK_COMMAND_HDIUTIL='/usr/bin/sudo /usr/bin/hdiutil'")
endif()

# cdash server (fu-berlin) SSL certificate sometimes is revoked. Keeps CI running.
set (CTEST_CURL_OPTIONS       CURLOPT_SSL_VERIFYHOST_OFF CURLOPT_SSL_VERIFYPEER_OFF )

message(STATUS "CTEST_SOURCE_DIRECTORY: ${CTEST_SOURCE_DIRECTORY}")
message(STATUS "CTEST_BINARY_DIRECTORY: ${CTEST_BINARY_DIRECTORY}")

#TODO set the Model to Continuous or Nightly depending on branch
ctest_start(Nightly GROUP Package)

# build the dist target for packages
# we assume the configuration is correct for this. So please no KNIME, PYOPENMS, STYLE or COVERAGE.

ctest_configure(OPTIONS "${CONFIGURE_OPTIONS}" RETURN_VALUE _reconfig_package_ret_val)
ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" TARGET "dist" NUMBER_ERRORS _build_errors)
ctest_submit(PARTS Build)

string(REPLACE "+" "%2B" BUILD_NAME_SAFE ${CTEST_BUILD_NAME})
string(REPLACE "." "%2E" BUILD_NAME_SAFE ${BUILD_NAME_SAFE})
string(REPLACE "/" "%2F" BUILD_NAME_SAFE ${BUILD_NAME_SAFE})

if (_build_errors)
  message(FATAL_ERROR "There were errors: Please check the build results at: https://cdash.seqan.de/index.php?project=OpenMS&begin=2023-01-01&end=2030-01-01&filtercount=1&field1=buildname&compare1=63&value1=${BUILD_NAME_SAFE}")
else()
  message("Packaging successful: Please check the build results at: https://cdash.seqan.de/index.php?project=OpenMS&begin=2023-01-01&end=2030-01-01&filtercount=1&field1=buildname&compare1=63&value1=${BUILD_NAME_SAFE}")
endif()
