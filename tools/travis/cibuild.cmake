##
## this script is invoked by lnx-cibuild.sh during the main "script:" section in .travis.yml
##

# define build name&co for easier identification on CDash
set(CTEST_BUILD_NAME "$ENV{BUILD_NAME}")

set(CTEST_SITE "travis-ci-build-server")
set(CTEST_SOURCE_DIRECTORY "$ENV{SOURCE_DIRECTORY}")
set(CTEST_BINARY_DIRECTORY "${CTEST_SOURCE_DIRECTORY}/_build")

message(STATUS "CTEST_SOURCE_DIRECTORY: ${CTEST_SOURCE_DIRECTORY}")
message(STATUS "CTEST_BINARY_DIRECTORY: ${CTEST_BINARY_DIRECTORY}")

set(INITIAL_CACHE
"CMAKE_PREFIX_PATH=$ENV{OS_PREFIX_PATH}
OPENMS_CONTRIB_LIBS=$ENV{SOURCE_DIRECTORY}/contrib
BOOST_USE_STATIC=Off
CMAKE_BUILD_TYPE=$ENV{BUILD_TYPE}
ENABLE_TUTORIALS=Off
ENABLE_GCC_WERROR=On
PYOPENMS=$ENV{PYOPENMS}
MT_ENABLE_OPENMP=$ENV{OPENMP}
PYTHON_EXECUTABLE:FILEPATH=$ENV{PYTHON_EXE}
PY_NUM_THREADS=4
PY_NUM_MODULES=4
PY_NO_OPTIMIZATION=ON
PY_SINGLE_THREADED=ON
ENABLE_STYLE_TESTING=$ENV{ENABLE_STYLE_TESTING}
ENABLE_TOPP_TESTING=$ENV{ENABLE_TOPP_TESTING}
ENABLE_CLASS_TESTING=$ENV{ENABLE_CLASS_TESTING}
WITH_GUI=$ENV{WITH_GUI}
ADDRESS_SANITIZER=$ENV{ADDRESS_SANITIZER}"
)

# create cache
file(WRITE "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt" ${INITIAL_CACHE})

# ignore failing GzipIfstream_test which seems to be related to the used
# zlib version
set(CTEST_CUSTOM_TESTS_IGNORE
	GzipIfstream_test
)

# customize reporting of errors in CDash
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_ERRORS 1000)
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS 1000)

set (CTEST_CUSTOM_WARNING_EXCEPTION
    # Suppress warnings imported from qt
    ".*qsharedpointer_impl.h:595:43.*"
    )

# try to speed up the builds so we don't get killed
set(CTEST_BUILD_FLAGS -j5)

## speed up compile time on GCC
if (CMAKE_COMPILER_IS_GNUCXX)
	add_definitions(-O0)
endif()

# we want makefiles
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")

# run the classical CTest suite without update
# travis-ci handles this for us
ctest_start     (Continuous)
ctest_configure (BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE _configure_ret)

# we only build when we do non-style testing and we may have special targets like pyopenms
if("$ENV{ENABLE_STYLE_TESTING}" STREQUAL "OFF")
  if("$ENV{PYOPENMS}" STREQUAL "ON")
    ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" TARGET "pyopenms" NUMBER_ERRORS _build_errors)
  else()
    ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" NUMBER_ERRORS _build_errors)
  endif()
else()
	set(_build_errors 0)
endif()

## build lib&executables, run tests
## for pyopenms build, only run pyopenms tests
if("$ENV{PYOPENMS}" STREQUAL "ON")
  ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}" INCLUDE "pyopenms" PARALLEL_LEVEL 3)
else()
  ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}" PARALLEL_LEVEL 3)
endif()
## send to CDash
ctest_submit()

# indicate errors
if(${_build_errors} GREATER 0 OR NOT ${_configure_ret} EQUAL 0)
  file(WRITE "$ENV{SOURCE_DIRECTORY}/failed" "build_failed")
endif()
