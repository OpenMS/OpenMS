#!/bin/bash
# Time ThermoRawFile::transform (counting consumer) for three versions of the three Thermo source files:
# develop, PR v1 (7bebaff053), PR v2 (HEAD). Rebuilds libOpenMS in the PR build dir for each, restores HEAD at the end.
set -e
R=/home/singjc/Documents/github/OpenMS-thermo-readahead
B=/home/singjc/Documents/github/OpenMS-thermo-readahead-build
T=/home/singjc/.jcode/scratch/thermo_time
RAW=/home/singjc/Documents/datasets/proteobench-astral-entrapment/raw/LFQ_Astral_DIA_15min_50ng_Human_01.raw
export DOTNET_ROOT=/home/singjc/Documents/github/OpenMS-wp00-build/dependencies/dotnet
FILES="src/openms/include/OpenMS/FORMAT/HANDLERS/ThermoRawScanReader.h src/openms/source/FORMAT/HANDLERS/ThermoRawScanReader.cpp src/openms/source/FORMAT/ThermoRawFile.cpp"
cd $R
cat $RAW > /dev/null
for V in upstream/develop 7bebaff053 HEAD; do
  git checkout -q $V -- $FILES
  (cd $B && nice ninja OpenMS > /dev/null 2>&1)
  for rep in 1 2 3; do
    for args in "2 4" "1 4" "2 1"; do
      echo "$V rep$rep $($T/thermo_time $RAW $args)"
    done
  done
done
git checkout -q HEAD -- $FILES
(cd $B && nice ninja OpenMS > /dev/null 2>&1)
git status --short
echo done
