#!/bin/sh
# 2026-10-02: open PRs on develop with the #10391 deisotoping fix.
cd "$(dirname "$0")"
python3 run.py --arms d_base d_priors d_79 d_79_mass --workers 2 > run_d.log 2>&1; echo "rc=$?" >> run_d.log
until grep -q '^rc=' build_p35_val.log; do sleep 30; done
if grep -q '^rc=0' build_p35_val.log; then
  python3 run.py --arms d_35 --workers 2 > run_d35.log 2>&1; echo "rc=$?" >> run_d35.log
else
  echo "p35_val build failed" > run_d35.log
fi
