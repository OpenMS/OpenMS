#!/bin/sh
cd "$(dirname "$0")"
until grep -q '^rc=' run_opt.log 2>/dev/null; do sleep 30; done
python3 run.py --arms d_fc_hr --workers 2 > run_fchr.log 2>&1; echo "rc=$?" >> run_fchr.log
