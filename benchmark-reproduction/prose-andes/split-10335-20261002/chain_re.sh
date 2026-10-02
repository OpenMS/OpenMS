#!/bin/sh
cd "$(dirname "$0")"
until grep -q '^rc=' build_re_val.log 2>/dev/null; do sleep 30; done
until grep -q '^rc=' run_local.log 2>/dev/null; do sleep 30; done
if grep -q "'tests_passed': True" build_re_val.log; then
  python3 run.py --arms re_off re_rl re_r --workers 2 > run_re.log 2>&1; echo "rc=$?" >> run_re.log
else
  echo "re_val build or tests failed" > run_re.log; echo "rc=1" >> run_re.log
fi
