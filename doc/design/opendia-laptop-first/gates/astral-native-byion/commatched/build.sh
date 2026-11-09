#!/bin/bash
cd "$(dirname "$0")"; B=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
set -o pipefail
$JCODE_PY build_lib.py > build.log 2>&1 && \
$B/TargetedFileConverter -in library.tsv -in_type tsv -out library.pqp -out_type pqp > convert.log 2>&1
echo $? > build.EXIT
