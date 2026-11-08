#!/bin/bash
cd "$(dirname "$0")"
( for R in 01 02; do OPENDIA_RUN=$R $JCODE_SCRATCH_DIR/torchenv/.venv/bin/python calib_libreal.py 2>/dev/null; done ) > calib_libreal.txt
echo $? > calib_libreal.EXIT
