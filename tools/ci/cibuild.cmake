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

# cdash server (fu-berlin) SSL certificate sometimes is revoked. Keeps CI running.
# set (CTEST_CURL_OPTIONS       CURLOPT_SSL_VERIFYHOST_OFF CURLOPT_SSL_VERIFYPEER_OFF )

message(STATUS "CTEST_SOURCE_DIRECTORY: ${CTEST_SOURCE_DIRECTORY}")
message(STATUS "CTEST_BINARY_DIRECTORY: ${CTEST_BINARY_DIRECTORY}")

set(INITIAL_CACHE
"
GIT_TRACKING=OFF
OPENMS_GIT_SHORT_SHA1=ci
OPENMS_GIT_SHORT_REFSPEC=$ENV{RUN_NAME}
OPENMS_GIT_LC_DATE=1970-01-01
Boost_DEBUG=OFF
CMAKE_PREFIX_PATH=$ENV{OS_PREFIX_PATH}
OPENMS_CONTRIB_LIBS=$ENV{CONTRIB_BUILD_DIRECTORY}
BOOST_USE_STATIC=$ENV{USE_STATIC_BOOST}
CMAKE_BUILD_TYPE=$ENV{BUILD_TYPE}
ENABLE_TUTORIALS=Off
ENABLE_GCC_WERROR=Off
PYOPENMS=$ENV{PYOPENMS}
MT_ENABLE_OPENMP=$ENV{OPENMP}
PYTHON_EXECUTABLE:FILEPATH=$ENV{PYTHON_EXE}
PY_NUM_THREADS=4
PY_NUM_MODULES=4
PY_NO_OPTIMIZATION=ON
PY_MEMLEAK_DISABLE=ON
PY_SINGLE_THREADED=ON
PY_NO_OUTPUT=On
ENABLE_STYLE_TESTING=$ENV{ENABLE_STYLE_TESTING}
ENABLE_TOPP_TESTING=$ENV{ENABLE_TOPP_TESTING}
ENABLE_CLASS_TESTING=$ENV{ENABLE_CLASS_TESTING}
WITH_GUI=$ENV{WITH_GUI}
ADDRESS_SANITIZER=$ENV{ADDRESS_SANITIZER}
WITH_THERMORAWFILEPARSER_TEST=Off"
)

set(OWN_OPTIONS "")
if($ENV{CMAKE_GENERATOR} MATCHES ".*Visual Studio.*")
  set(INITIAL_CACHE 
"${INITIAL_CACHE}
CMAKE_GENERATOR_PLATFORM=x64"
  )
  set(OWN_OPTIONS "-DCMAKE_CXX_RELEASE_FLAGS='/MD /Od /Ob0 /DNDEBUG /EHsc'")
endif()

if(DEFINED ENV{CMAKE_CCACHE_EXE})
  set(INITIAL_CACHE 
"${INITIAL_CACHE}
CMAKE_CXX_COMPILER_LAUNCHER=$ENV{CMAKE_CCACHE_EXE}
CMAKE_C_COMPILER_LAUNCHER=$ENV{CMAKE_CCACHE_EXE}"
  )
endif()

# create cache
file(WRITE "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt" ${INITIAL_CACHE})

# customize reporting of errors in CDash
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_ERRORS 1000)
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS 1000)

set (CTEST_CUSTOM_WARNING_EXCEPTION
    # Suppress warnings imported from qt
    ".*qsharedpointer_impl.h:595:43.*"
    )

# try to speed up the builds so we don't get killed
set(CTEST_BUILD_FLAGS "$ENV{BUILD_FLAGS}")

## speed up compile time on GCC
if (CMAKE_COMPILER_IS_GNUCXX)
	add_compile_options(-O0)
endif()

set(CTEST_CMAKE_GENERATOR "$ENV{CMAKE_GENERATOR}")

# run the classical CTest suite
ctest_start(Continuous) # TODO think about adding GROUP GitHub-Actions to separate visually

# Gather update information.
find_package(Git)
set(CTEST_UPDATE_VERSION_ONLY ON)
set(CTEST_UPDATE_COMMAND "${GIT_EXECUTABLE}")
ctest_update()

ctest_configure (BUILD "${CTEST_BINARY_DIRECTORY}" OPTIONS "${OWN_OPTIONS}" RETURN_VALUE _configure_ret)
ctest_submit(PARTS Update Configure)

# we only build when we do non-style testing and we may have special targets like pyopenms
if("$ENV{ENABLE_STYLE_TESTING}" STREQUAL "OFF")
  if("$ENV{PYOPENMS}" STREQUAL "ON")
    ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" TARGET "pyopenms" NUMBER_ERRORS _build_errors)
  else()
    if(WIN32)
       ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" TARGET "GUI" NUMBER_ERRORS _build_errors)
       ctest_submit(PARTS Build)
       ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" APPEND TARGET "TOPP" NUMBER_ERRORS _build_errors)
       ctest_submit(PARTS Build)
       ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" APPEND TARGET "UTILS" NUMBER_ERRORS _build_errors)
       ctest_submit(PARTS Build)
       set(CTEST_BUILD_FLAGS "-j1") # Super duper hack, since no one wants to debug excessive memory usage on win
       ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" APPEND NUMBER_ERRORS _build_errors)
       ctest_submit(PARTS Build)
    else()
      ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" NUMBER_ERRORS _build_errors)
      ctest_submit(PARTS Build)
    endif()
  endif()
else()
  set(_build_errors 0)
endif()

message("Please check the build results at: https://cdash.openms.de/index.php?project=OpenMS&begin=2023-01-01&end=2030-01-01&filtercount=1&field1=buildname&compare1=63&value1=${CTEST_BUILD_NAME}")
