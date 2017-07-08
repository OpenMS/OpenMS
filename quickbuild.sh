#!/bin/sh
set -e

###################################
# 1. download and build contrib (if not present already)
if [ ! -d "contrib" ] ; then
  git clone https://github.com/OpenMS/contrib.git
fi
mkdir -p contrib-build
cd contrib-build

CONTRIB_PATH=`pwd`
cmake -DBUILD_TYPE=ALL ../contrib

###################################
# 2. build OpenMS
#  - default builds all but first argument is passed to make
cd ..
mkdir -p openms-build
cd openms-build

# needs absolute path!
cmake -DOPENMS_CONTRIB_LIBS="$CONTRIB_PATH" ../
if [ -z "$1" ]; then
  make $1
else
  make
fi

