#!/bin/bash

# helper function to convert build name into a proper, CDash-compatible format
function cdashify()
{
  cdash_version=$(echo $1 | sed 's/\//_/g')
  echo ${cdash_version}
}

# get some infos from git to embed it in the build name
export SOURCE_DIRECTORY=`pwd`
mkdir _build

# additional variables
export CMAKE_GENERATOR="Unix Makefiles"
export CONTRIB_BUILD_DIRECTORY="$SOURCE_DIRECTORY/contrib"
export OPENMS_CONTRIB_LIBS="$SOURCE_DIRECTORY/contrib"
export USE_STATIC_BOOST="Off"

# assemble a proper build name
_build_name=""

# extend with specific information
if [ "$TRAVIS_PULL_REQUEST" != "false" ]; then
  _build_name="pr"$(cdashify ${TRAVIS_PULL_REQUEST})
else
  _build_name=$(cdashify ${TRAVIS_BRANCH})
fi

# append OS
_build_name=${_build_name}"-macOS"

# append compiler info
_build_name=${_build_name}"-"${CXX}

# add style to build name if requested
if [ "${ENABLE_STYLE_TESTING}" = "ON" ]; then
  _build_name=${_build_name}"-codingStyle"
fi

# add class-testing to build name if requested
if [ "${ENABLE_CLASS_TESTING}" = "ON" ]; then
  _build_name=${_build_name}"-testClass"
fi
# add TOPP-testing to build name if requested
if [ "${ENABLE_TOPP_TESTING}" = "ON" ]; then
  _build_name=${_build_name}"-testTOPP"
fi

# add GUI to build name if requested
if [ "${WITH_GUI}" = "OFF" ]; then
  _build_name=${_build_name}"-noGUI"
fi

# we will use this in the CMake script
export BUILD_NAME=${_build_name}

# we need an X-server for building the documentation and some tests
# so we start xvfb
#export DISPLAY=:99.0

#sh -e /usr/X11/bin/xvfb start

# add third-party binaries (e.g. search engines) to PATH
#export PATH=${SOURCE_DIRECTORY}/_thirdparty/XTandem:$PATH
#export PATH=${SOURCE_DIRECTORY}/_thirdparty/MSGFPlus:$PATH

# if we perform style tests, add cppcheck to path
#if [ $ENABLE_STYLE_TESTING == "ON" ]; then
#  export PATH=${SOURCE_DIRECTORY}/cppcheck:$PATH
#fi

# set os dependent folder for preinstalled libraries (in this case to homebrew libraries)
export OS_PREFIX_PATH=/usr/local
ctest --output-on-failure -V -S tools/travis/cibuild.cmake

# tell the user where he can find the results
echo "Please check the build results at: http://cdash.seqan.de/index.php?project=OpenMS&date="$(date +"%y-%m-%d")"#Continuous"
echo "This build has the name: ${BUILD_NAME}"

# we indicate build failures if CTest experienced any errors
if [ -f ${SOURCE_DIRECTORY}/failed ]; then
  echo "Configure/Build failed"
  exit -1
fi

FAILED_TEST=$(find _build -name "Test.xml" -type f | xargs grep "<Test Status=\"failed\">" -c)

if [ "$FAILED_TEST" -gt "0" ]; then
  echo "$FAILED_TEST tests failed"
  exit -1
fi

# seems like everything worked
exit 0
