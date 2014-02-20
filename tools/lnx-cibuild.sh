#!/bin/bash

# helper function to convert build name into a proper, cdash compatible format
function cdashify()
{
  cdash_version=$(echo $1 | sed 's/\//_/g')
  echo ${cdash_version}
}

# get some infos from git to embed it in the build name
export SOURCE_DIRECTORY=`pwd`
mkdir _build

# define the build name
if [ "$COVERITY_SCAN_BRANCH" == "1" ]; then
  exit 0
elif [ "$TRAVIS_PULL_REQUEST" != "false" ]; then
  export BUILD_NAME=$TRAVIS_PULL_REQUEST
elif [[ -n "$TRAVIS_COMMIT_RANGE" ]]; then
  export BUILD_NAME=$TRAVIS_COMMIT_RANGE
else
  _build_name=${_build_name}"-"$(cdashify ${TRAVIS_COMMIT})
fi

# append compiler info
_build_name=${_build_name}"-"${CXX}

# add style to build name if requested
if [ "${ENABLE_STYLE_TESTING}" = "On" ]; then
  _build_name=${_build_name}"-coding-style"
fi

# add unity to build name if requested
if [ "${ENABLE_UNITYBUILD}" = "On" ]; then
  _build_name=${_build_name}"-unity-build"
fi

# add unity to build name if requested
if [ "${WITH_GUI}" = "Off" ]; then
  _build_name=${_build_name}"-no-gui"
fi

# we will use this in the cmake script
export BUILD_NAME=${_build_name}

# we need an X-server for building the documentation and some tests
# so we start xvfb
export DISPLAY=:99.0
sh -e /etc/init.d/xvfb start

# add search engines to PATH
export PATH=${SOURCE_DIRECTORY}/_searchengines/MyriMatch:$PATH
export PATH=${SOURCE_DIRECTORY}/_searchengines/OMSSA:$PATH
export PATH=${SOURCE_DIRECTORY}/_searchengines/XTandem:$PATH
export PATH=${SOURCE_DIRECTORY}/_searchengines/MSGFPlus:$PATH

# if we perform style tests, add cppcheck to path
if [ $ENABLE_STYLE_TESTING == "On" ]; then
  export PATH=${SOURCE_DIRECTORY}/cppcheck:$PATH
fi

ctest -VV -S tools/lnx-cibuild.cmake

# tell the user where he can find the results
echo "Please check the build results at: http://cdash.openms.de/index.php?project=OpenMS&date="$(date +"%y-%m-%d")"#Continuous"
echo "This build has the name: ${BUILD_NAME}"

# we indicate build failures if ctest experienced any errors
if [ -f ${SOURCE_DIRECTORY}/failed ]; then
  echo "Configure/Build failed"
  exit -1
fi

FAILED_TEST=$(find _build -name "Test.xml" -type f | xargs grep "<Test Status=\"failed\">" -c)

if [ "$FAILED_TEST" -gt "0" ]; then
  echo "$FAILED_TEST tests failed"
  exit -1
fi

# it seems like everything worked
exit 0
