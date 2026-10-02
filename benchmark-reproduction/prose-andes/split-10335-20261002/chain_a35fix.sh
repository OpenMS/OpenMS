#!/bin/sh
cd "$(dirname "$0")"
until grep -q '^rc=\|failed' run_dd.log 2>/dev/null; do sleep 30; done
python3 run.py --resume --arms a35_hs_single a35_hs_multi --workers 2 velos_5000_R2 > run_a35c.log 2>&1; echo "rc=$?" >> run_a35c.log
