#!/bin/bash
cd "$(dirname "$0")"
PY=$JCODE_SCRATCH_DIR/torchenv/.venv/bin/python
( for R in 01 02; do OPENDIA_RUN=$R $PY calib_mitig.py 2>/dev/null; done ) > calib_mitig.txt
echo $? > calib_mitig.EXIT
