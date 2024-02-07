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
    )

message(STATUS "CTEST_SOURCE_DIRECTORY: ${CTEST_SOURCE_DIRECTORY}")
message(STATUS "CTEST_BINARY_DIRECTORY: ${CTEST_BINARY_DIRECTORY}")

# function to set a cache variable from an env variable if it exists only
macro(add_env_var_to_cache_if_exists VAR_NAME)
  if (DEFINED ENV{${VAR_NAME}})
    message("tools/ci/cibuild.cmake: Found ${VAR_NAME} with value $ENV{${VAR_NAME}}")
    set(INITIAL_CACHE
"${INITIAL_CACHE}
${VAR_NAME}=$ENV{${VAR_NAME}}"
    )
  endif()
endmacro()

# same but for multiple variables
macro(add_env_vars_to_cache_if_exists VAR_NAMES)
  foreach(VAR_NAME ${VAR_NAMES})
    add_env_var_to_cache_if_exists(${VAR_NAME})
  endforeach()
endmacro()

macro(add_env_var_to_cache_with_default VAR_NAME DEFAULT)
  if (DEFINED ENV{${VAR_NAME}})
    set(INITIAL_CACHE
"${INITIAL_CACHE}
${VAR_NAME}=$ENV{${VAR_NAME}}"
    )
  else()
      set(INITIAL_CACHE
"${INITIAL_CACHE}
${VAR_NAME}=$ENV{${DEFAULT}}"
    )
  endif()
endmacro()

set(INITIAL_CACHE "")

if(DEFINED ENV{CMAKE_CCACHE_EXE})
  set(INITIAL_CACHE 
"${INITIAL_CACHE}
CMAKE_CXX_COMPILER_LAUNCHER=$ENV{CMAKE_CCACHE_EXE}
CMAKE_C_COMPILER_LAUNCHER=$ENV{CMAKE_CCACHE_EXE}"
  )
endif()

set(VARS_TO_LOAD
  "ADDRESS_SANITIZER"
  "CMAKE_PREFIX_PATH"
  "CMAKE_BUILD_TYPE"
  "CMAKE_GENERATOR_PLATFORM"
  "Boost_DEBUG"
  "BOOST_USE_STATIC"
  "OPENMS_CONTRIB_LIBS"
  "ENABLE_CLASS_TESTING"
  "ENABLE_GCC_WERROR"
  "ENABLE_STYLE_TESTING"
  "ENABLE_TOPP_TESTING"
  "ENABLE_TUTORIALS"
  "ENABLE_UPDATE_CHECK"
  "MT_ENABLE_OPENMP"
  "SEARCH_ENGINES_DIRECTORY"
  "PACKAGE_TYPE"
  "PYOPENMS"
  "PY_MEMLEAK_DISABLE"
  "PY_NO_OPTIMIZATION"
  "PY_NO_OUTPUT"
  "PY_NUM_MODULES"
  "PY_NUM_THREADS"
  "WITH_GUI"
  "WITH_THERMORAWFILEPARSER_TEST"
 )

message("tools/ci/cibuild.cmake: Loading the following vars from ENV if available: ${VARS_TO_LOAD}")
add_env_vars_to_cache_if_exists("${VARS_TO_LOAD}")

# Unused now! If you want to set a variable to a non-default, you have to set it in your environment.
# To make build scripts/workflows clearer from the outside.
set(OLD_VALUES
"
GIT_TRACKING=OFF
OPENMS_GIT_SHORT_SHA1=ci
OPENMS_GIT_SHORT_REFSPEC=$ENV{RUN_NAME}
OPENMS_GIT_LC_DATE=1970-01-01
ENABLE_UPDATE_CHECK:BOOL=OFF
Boost_DEBUG=OFF
CMAKE_PREFIX_PATH=$ENV{OS_PREFIX_PATH}
OPENMS_CONTRIB_LIBS=$ENV{CONTRIB_BUILD_DIRECTORY}
BOOST_USE_STATIC=$ENV{USE_STATIC_BOOST}
CMAKE_BUILD_TYPE=$ENV{BUILD_TYPE}
PACKAGE_TYPE=$ENV{PACKAGE_TYPE}
SEARCH_ENGINES_DIRECTORY=$ENV{SEARCH_ENGINES_DIRECTORY}
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

# create cache
file(WRITE "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt" ${INITIAL_CACHE})

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
ctest_submit(PARTS Update Configure)

# we only build when we do non-style testing and we may have special targets like pyopenms
if("$ENV{ENABLE_STYLE_TESTING}" STREQUAL "OFF")
  if("$ENV{PYOPENMS}" STREQUAL "ON")
    ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" TARGET "pyopenms" NUMBER_ERRORS _build_errors)
  else()
    ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" NUMBER_ERRORS _build_errors)
    ctest_submit(PARTS Build)
  endif()
else()
  set(_build_errors 0)
endif()

string(REPLACE "+" "%2B" BUILD_NAME_SAFE ${CTEST_BUILD_NAME})
string(REPLACE "." "%2E" BUILD_NAME_SAFE ${BUILD_NAME_SAFE})
string(REPLACE "/" "%2F" BUILD_NAME_SAFE ${BUILD_NAME_SAFE})

if (_build_errors)
  message(FATAL_ERROR "There were errors: Please check the build results at: https://cdash.openms.de/index.php?project=OpenMS&begin=2023-01-01&end=2030-01-01&filtercount=1&field1=buildname&compare1=63&value1=${BUILD_NAME_SAFE}")
else()
  message("Build successful: Please check the build results at: https://cdash.openms.de/index.php?project=OpenMS&begin=2023-01-01&end=2030-01-01&filtercount=1&field1=buildname&compare1=63&value1=${BUILD_NAME_SAFE}")
endif()
