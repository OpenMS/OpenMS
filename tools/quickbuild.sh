#!/bin/sh
set -e

###################################
# 1. download and build contrib (if not present already)
if [ ! -d "contrib" ] ; then
  git clone https://github.com/OpenMS/contrib.git
fi

# set contrib absolute path for configure
CONTRIB_PATH=`pwd`/contrib-build

if [ ! -d "contrib-build" ] ; then
  mkdir -p contrib-build
  cd contrib-build
  cmake -DBUILD_TYPE=ALL ../contrib
  cd ..
fi

###################################
# 2. build OpenMS
#  - default builds all but first argument is passed to make
mkdir -p openms-build
cd openms-build

# contrib path needs to be an absolute path!
cmake -DOPENMS_CONTRIB_LIBS="$CONTRIB_PATH" -DBOOST_USE_STATIC=On ../
make $@

cd ..

