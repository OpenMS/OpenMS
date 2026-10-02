#!/bin/sh
cd "$(dirname "$0")"
timeout 2400 python3 build.py fc_val --jobs 2 > build_fc_val.log 2>&1; echo "rc=$?" >> build_fc_val.log
until grep -q '^rc=' run_a35c.log 2>/dev/null; do sleep 30; done
if grep -q "'tests_passed': True" build_fc_val.log; then
  python3 run.py --arms d_fc --workers 2 > run_fc.log 2>&1; echo "rc=$?" >> run_fc.log
else
  echo "fc_val build or tests failed" > run_fc.log
fi
python3 run.py --arms a35_raw a35_raw_local a35_local --workers 2 > run_opt.log 2>&1; echo "rc=$?" >> run_opt.log
