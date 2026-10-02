#!/bin/sh
cd "$(dirname "$0")"
until grep -q '^rc=' run_fchr.log 2>/dev/null; do sleep 30; done
python3 run.py --arms a35_local_hr a35_local_E a35_raw_local_E --workers 2 > run_local.log 2>&1; echo "rc=$?" >> run_local.log
