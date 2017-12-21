#!/bin/sh
set -e

# first argument used for number of jobs
if [ ! -z "$1" ]
then 
  numberOfJobs=$1
else 
  numberOfJobs=1
fi

###################################
# 1. download and build contrib (if not present already)
git submodule update --init contrib

# set contrib absolute path for configure
CONTRIB_PATH=`pwd`/contrib-build

if [ ! -d "contrib-build" ] ; then
  mkdir -p contrib-build
  cd contrib-build
  cmake -DBUILD_TYPE=ALL -DNUMBER_OF_JOBS=$numberOfJobs ../contrib
  cd ..
fi

###################################
# 2. build OpenMS
#  - default builds all but first argument is passed to make
mkdir -p openms-build
cd openms-build

# contrib path needs to be an absolute path!
cmake -DOPENMS_CONTRIB_LIBS="$CONTRIB_PATH" -DNUMBER_OF_JOBS=$numberOfJobs -DBOOST_USE_STATIC=On ../
make -j $numberOfJobs

cd ..

