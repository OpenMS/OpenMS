#!/bin/bash
# After chain01: the post-hoc variants on raw 01 (seed fallback, clean-pair set, shared-fragment tail).
cd "$(dirname "$0")"; P=$JCODE_SCRATCH_DIR/torchenv/.venv/bin/python
while [ ! -f chain01.EXIT ]; do sleep 60; done
OPENDIA_RUN=R01 $P score2_seed.py > score01_seed.txt 2>&1
OPENDIA_RUN=R01 $P s_tail.py > s_tail01.txt 2>&1
OPENDIA_RUN=R01 $P score3.py > score3_01.txt 2>&1
echo $? > post01.EXIT
