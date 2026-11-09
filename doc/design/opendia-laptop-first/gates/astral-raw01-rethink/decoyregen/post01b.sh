#!/bin/bash
cd "$(dirname "$0")"; while [ ! -f post01.EXIT ]; do sleep 60; done
OPENDIA_RUN=R01 $JCODE_SCRATCH_DIR/torchenv/.venv/bin/python score4.py > score4_01.txt 2>&1; echo $? > post01b.EXIT
