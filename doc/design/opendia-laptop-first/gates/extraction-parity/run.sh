#!/bin/bash
# Same DIA-NN 2.2.0 raw01 search as diann-ourlib-20261112/run.sh, plus --xic 60 (fragment + MS1 traces +-60 s around
# each identified apex) and --export-quant (fragment heights). Reference for OpenDIA extraction parity.
cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/extraction-parity-20261112/diann
trap 'echo $? > EXIT' EXIT
echo running > STATUS
mkdir -p tmp
docker run --rm --cpus 14 -u $(id -u):$(id -g) \
  -v /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/dd-plots-20261111:/in:ro \
  -v /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/diann-ourlib-20261112:/speclib:ro -v $PWD:/w \
  --entrypoint /usr/diann-2.2.0/diann-linux diann:2.2.0 \
  --f /in/raw01.mzML --lib /speclib/lib_targets.tsv --threads 14 --verbose 1 \
  --out /w/report_raw01.parquet --qvalue 1 --gen-decoys \
  --unimod4 --min-pr-charge 1 --max-pr-charge 5 --temp /w/tmp --xic 60 --export-quant > diann_raw01.log 2>&1
echo done > STATUS
