#!/bin/bash
# Raw 02 replication of the OSW-like teacher-free scorer (design diagnostic; user approved raw 02 reuse).
cd "$(dirname "$0")"
PY=$JCODE_SCRATCH_DIR/torchenv/.venv/bin/python
export OPENDIA_RUN=02 NUMBA_NUM_THREADS=12
$PY build_xic.py > raw02_build.log 2>&1 && \
$PY oswlike_feats.py > raw02_feats.txt 2>&1 && \
$PY oswlike_score.py > raw02_score.txt 2>&1
echo $? > raw02.EXIT
