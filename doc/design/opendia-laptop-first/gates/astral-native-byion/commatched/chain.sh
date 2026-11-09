#!/bin/bash
cd "$(dirname "$0")"
./run.sh "$1" > run.log 2>&1; echo $? > run.EXIT
"$1" evaluate.py raw01 > EVALUATION.txt 2>&1; "$1" evaluate.py raw02 >> EVALUATION.txt 2>&1; echo $? > eval.EXIT
