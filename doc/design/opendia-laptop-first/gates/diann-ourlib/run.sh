#!/bin/bash
cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/diann-ourlib-20261112
trap 'echo $? > EXIT' EXIT
echo running > STATUS
docker run --rm --cpus 14 -u $(id -u):$(id -g) \
  -v /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/dd-plots-20261111:/in:ro -v $PWD:/w \
  --entrypoint /usr/diann-2.2.0/diann-linux diann:2.2.0 \
  --f /in/raw01.mzML --lib /w/lib_targets.tsv --threads 14 --verbose 1 \
  --out /w/report_raw01.parquet --qvalue 1 --gen-decoys \
  --unimod4 --min-pr-charge 1 --max-pr-charge 5 --temp /w/tmp > diann_raw01.log 2>&1
echo done > STATUS
