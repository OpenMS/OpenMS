#!/bin/bash
cd "$(dirname "$0")"
( for R in 01 02; do OPENDIA_RUN=$R $JCODE_SCRATCH_DIR/torchenv/.venv/bin/python calib_scan.py 2>&1 | grep -v Warn; done ) > calib_scan.txt
echo $? > calib_scan.EXIT
