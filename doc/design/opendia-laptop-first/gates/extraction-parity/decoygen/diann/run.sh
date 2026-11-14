#!/bin/bash
# DIA-NN 2.2 with library-provided decoys: the PAIRLIB library with its old decoys (pseudo-reverse + K/R swap, copied
# target intensities and RT) and the arm D library (same sequences, own predicted intensities and RT), raws 01 and 02.
# Settings as diann-pairlib-20261113 (and diann-ourlib-20261112). --report-decoys made DIA-NN 2.2 segfault (exit 139)
# while writing the report on all five first runs (empty report.parquet), so it is not used; library decoy use is checked
# from the log (precursor count of the loaded library). --gen-decoys there was not a DIA-NN 2.2 option (ignored).
# Also DIA-NN's own decoys on raw02 (targets only), which the earlier run lacked.
# Usage: run.sh <tag> <lib dir> <raw 01|02> <threads>
cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/diann-decoys-20261114
tag=$1; lib=$2; raw=$3; th=$4
out=$tag/raw$raw; mkdir -p $out/tmp
/usr/bin/time -v docker run --rm --cpus $th -u $(id -u):$(id -g) \
  -v /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/dd-plots-20261111:/in:ro -v $PWD:/w \
  -v /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107:/b:ro \
  --entrypoint /usr/diann-2.2.0/diann-linux diann:2.2.0 \
  --f /in/raw$raw.mzML --lib /w/$lib/lib.tsv --threads $th --verbose 1 \
  --out /w/$out/report.parquet --qvalue 1 \
  --unimod4 --min-pr-charge 1 --max-pr-charge 5 --temp /w/$out/tmp > $out/diann.log 2> $out/time.log
echo $? > $out/EXIT
rm -rf $out/tmp
