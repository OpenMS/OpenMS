#!/bin/sh
cd "$(dirname "$0")"
timeout 3600 python3 build.py re_on --jobs 2 > build_re_on3.log 2>&1; echo "rc=$?" >> build_re_on3.log
until grep -q '^rc=' run_re.log 2>/dev/null; do sleep 30; done
if grep -q "'tests_passed': True" build_re_on3.log; then
  python3 run.py --arms re_auto --workers 2 > run_auto.log 2>&1; echo "rc=$?" >> run_auto.log
else
  echo "re_on build or tests failed" > run_auto.log; echo "rc=1" >> run_auto.log
fi
