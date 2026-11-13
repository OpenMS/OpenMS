#!/bin/bash
# PLAN_DECOYGEN library step: per-arm library TSV (build_lib.py) and PQP (TargetedFileConverter, as PAIRLIB lib/convert_predict.sh).
cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/decoygen-20261113
PY=/home/singjc/.jcode/scratch/torchenv/.venv/bin/python
BIN=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
for a in ${ARMS:-A B C D}; do
  L=arm$a/home/lib
  [ -s $L/library.tsv ] || $PY scripts/build_lib.py $a > /dev/null || { echo "build $a failed"; exit 1; }
  /usr/bin/time -v $BIN/TargetedFileConverter -in $L/library.tsv -in_type tsv -out $L/library.pqp -out_type pqp > $L/convert.stdout.log 2> $L/convert.stderr.log
  echo $? > $L/CONVERT_EXIT
  rm -f $L/library.tsv
  cp ../pairlib-20261113/lib/pred12.parquet $L/
  (cd $L && sha256sum library.pqp pred12.parquet decoy12.tsv > lib.sha256)
  echo "arm $a: convert exit $(cat $L/CONVERT_EXIT)"
done
