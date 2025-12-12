##
## this script is invoked by lnx-cibuild.sh during the main "script:" section in .travis.yml
##

# define build name&co for easier identification on CDash
set(CTEST_BUILD_NAME "$ENV{BUILD_NAME}")

set(CTEST_SITE "$ENV{CI_PROVIDER}")
set(CTEST_SOURCE_DIRECTORY "$ENV{SOURCE_DIRECTORY}")
set(CTEST_BINARY_DIRECTORY "${CTEST_SOURCE_DIRECTORY}/bld")
set(CTEST_CONFIGURATION_TYPE "$ENV{BUILD_TYPE}")
set(CTEST_BUILD_CONFIGURATION "$ENV{BUILD_TYPE}")
set(CTEST_CMAKE_GENERATOR "$ENV{CMAKE_GENERATOR}")
# custom build flags
set(CTEST_BUILD_FLAGS "$ENV{BUILD_FLAGS}")

# cdash server (fu-berlin) SSL certificate sometimes is revoked. Keeps CI running.
set (CTEST_CURL_OPTIONS       CURLOPT_SSL_VERIFYHOST_OFF CURLOPT_SSL_VERIFYPEER_OFF )
# customize reporting of errors in CDash
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_ERRORS 1000)
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS 1000)

set (CTEST_CUSTOM_WARNING_EXCEPTION
    # Suppress warnings imported from qt
    ".*qsharedpointer_impl.h:595:43.*"
    ".*src/openms/extern/.*"
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
    )

message(STATUS "CTEST_SOURCE_DIRECTORY: ${CTEST_SOURCE_DIRECTORY}")
message(STATUS "CTEST_BINARY_DIRECTORY: ${CTEST_BINARY_DIRECTORY}")

set(OWN_OPTIONS "")
if($ENV{CMAKE_GENERATOR} MATCHES ".*Visual Studio.*")
  set(OWN_OPTIONS "-DCMAKE_CXX_RELEASE_FLAGS='/MD /Od /Ob0 /DNDEBUG /EHsc'")
endif()

# run the classical CTest suite
ctest_start(Continuous) # TODO think about adding GROUP GitHub-Actions to separate visually

# Gather update information.
find_package(Git)
set(CTEST_UPDATE_VERSION_ONLY ON)
set(CTEST_UPDATE_COMMAND "${GIT_EXECUTABLE}")
ctest_update()

ctest_configure (BUILD "${CTEST_BINARY_DIRECTORY}" OPTIONS "${OWN_OPTIONS}" RETURN_VALUE _configure_ret)
ctest_submit(PARTS Update Configure CAPTURE_CMAKE_ERROR _submit_result)

if(NOT _submit_result EQUAL 0)
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E echo "::warning file=cibuild.cmake,line=168::CTest submission failed, CDASH server is not available. Continuing execution."
    )
    message(WARNING "CTest submission failed, no detailed logs will be available.")
endif()


# we only build when we do non-style testing and we may have special targets like pyopenms
if("$ENV{ENABLE_STYLE_TESTING}" STREQUAL "OFF")
  # Initialize _submit_result to track any submission failures
  set(_submit_result 0)
  
  ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" NUMBER_ERRORS _build_errors)
  ctest_submit(PARTS Build CAPTURE_CMAKE_ERROR _this_submit_result)
  if(NOT _this_submit_result EQUAL 0)
    set(_submit_result 1)
  endif()

  if("$ENV{PYOPENMS}" STREQUAL "ON")
    ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" TARGET "pyopenms" APPEND NUMBER_ERRORS _py_build_errors)
    ctest_submit(PARTS Build CAPTURE_CMAKE_ERROR _this_submit_result)
    if(NOT _this_submit_result EQUAL 0)
      set(_submit_result 1)
    endif()
    math(EXPR _build_errors "${_build_errors} + ${_py_build_errors}")
  endif()

  # Only build compile_pxds if PYOPENMS is not ON (since it's already a subtarget of pyopenms)
  if("$ENV{COMPILE_PXDS}" STREQUAL "ON" AND "$ENV{PYOPENMS}" STREQUAL "OFF")
    ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" TARGET "compile_pxds" APPEND NUMBER_ERRORS _pdxs_build_errors)
    ctest_submit(PARTS Build CAPTURE_CMAKE_ERROR _this_submit_result)
    if(NOT _this_submit_result EQUAL 0)
      set(_submit_result 1)
    endif()
    math(EXPR _build_errors "${_build_errors} + ${_pdxs_build_errors}")
  endif()

  # Generate and validate the CWL files if "ENABLE_CWL_GENERATION" is set
  if("$ENV{ENABLE_CWL_GENERATION}" STREQUAL "ON")
    ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" TARGET "generate_cwl_files" APPEND NUMBER_ERRORS _cwl_build_errors)
    ctest_submit(PARTS Build CAPTURE_CMAKE_ERROR _this_submit_result)
    if(NOT _this_submit_result EQUAL 0)
      set(_submit_result 1)
    endif()
    math(EXPR _build_errors "${_build_errors} + ${_cwl_build_errors}")
  endif()
else()
  set(_build_errors 0)
  set(_submit_result 0)
endif()

if(NOT _submit_result EQUAL 0)
  execute_process(COMMAND ${CMAKE_COMMAND} -E echo "::warning file=cibuild.cmake,line=193::CTest submission failed, CDASH server is not available. Continuing execution.")
  message(WARNING "CTest submission failed, no detailed logs will be available.")
  if (_test_errors)
    message(FATAL_ERROR "There were errors: aborting")
  endif()
else()
  string(REPLACE "+" "%2B" BUILD_NAME_SAFE ${CTEST_BUILD_NAME})
  string(REPLACE "." "%2E" BUILD_NAME_SAFE ${BUILD_NAME_SAFE})
  string(REPLACE "/" "%2F" BUILD_NAME_SAFE ${BUILD_NAME_SAFE})
  
  if (_build_errors)
    message(FATAL_ERROR "There were errors: Please check the build results at: https://cdash.seqan.de/index.php?project=OpenMS&begin=2023-01-01&end=2030-01-01&filtercount=1&field1=buildname&compare1=63&value1=${BUILD_NAME_SAFE}")
  else()
    message("Build successful: Please check the build results at: https://cdash.seqan.de/index.php?project=OpenMS&begin=2023-01-01&end=2030-01-01&filtercount=1&field1=buildname&compare1=63&value1=${BUILD_NAME_SAFE}")
  endif()
endif()
