#!/bin/sh
set -e

echo `pwd`

###################################
# 1. download and build contrib (if not present already)
git submodule update --init contrib

# set contrib absolute path for configure
CONTRIB_PATH=`pwd`/contrib-build

if [ ! -d "contrib-build" ] ; then
  mkdir -p contrib-build
  cd contrib-build
  cmake -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_C_COMPILER=clang -DBUILD_TYPE=ALL ../contrib
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

