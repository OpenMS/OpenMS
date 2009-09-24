## general purpose nightly testing script for OpenMS
## call this using
## e.g. ctest -S ../Nightly.cmake -VV > ../CTEST_debug.log
##
## Customize to fit your system
##
## TODO
##   * add initial checkout
##   * maybe we should read arguments from cmd line to be more generic

CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

# set variables describing the build
SET (CTEST_SITE       "scratchy.mi.fu-berlin.de")
SET (CTEST_BUILD_NAME "Vista-x64-VS9-x64 Debug")

# set variables describing the build environments
SET (CTEST_SOURCE_DIRECTORY "C:/dev/NIGHTLY/OpenMS")
SET (CTEST_BINARY_DIRECTORY "C:/dev/NIGHTLY/OpenMS_build_debug")
SET (CTEST_BINARY_TEST_DIRECTORY "C:/dev/NIGHTLY/OpenMS_build_debug/source/TEST/")

# general ctest/cmake/other programs
SET (CTEST_CMAKE_COMMAND   "cmake")
SET (CTEST_UPDATE_COMMAND  "svn")

# define generator
SET (CTEST_CMAKE_GENERATOR "Visual Studio 9 2008 Win64" )

SET (CTEST_COMMAND "ctest.exe -D Nightly -C Debug")

SET(CTEST_BUILD_CONFIGURATION Debug) 
  
# clear the binary directory to avoid problems
CTEST_EMPTY_BINARY_DIRECTORY (${CTEST_BINARY_DIRECTORY})

# this is the initial cache to use for the binary tree, be careful to escape
# any quotes inside of this string if you use it
FILE(WRITE "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt" "
CONTRIB_CUSTOM_DIR:PATH=C:/dev/contrib_build
CMAKE_GENERATOR:INTERNAL=Visual Studio 9 2008 Win64
OPENMS_BUILD_TYPE:STRING=Debug
CMAKE_CONFIGURATION_TYPES:STRING=Debug
")

# do the dashboard/testings steps
CTEST_START     (Nightly)
CTEST_UPDATE    (SOURCE "${CTEST_SOURCE_DIRECTORY}")
CTEST_CONFIGURE (BUILD "${CTEST_BINARY_DIRECTORY}")
CTEST_BUILD     (BUILD "${CTEST_BINARY_DIRECTORY}")
SET (CTEST_PROJECT_NAME "OpenMS_tests")
CTEST_BUILD     (BUILD "${CTEST_BINARY_TEST_DIRECTORY}")
SET (CTEST_PROJECT_NAME "OpenMS")
CTEST_TEST      (BUILD "${CTEST_BINARY_DIRECTORY}")
CTEST_SUBMIT    ()