#!/bin/sh
# #10335 low-resolution scoring split, then #10335's peptide deduplication alone on develop.
cd "$(dirname "$0")"
python3 run.py --resume --arms a35_hs_single a35_hs_multi a35_cal_single --workers 2 > run_a35b.log 2>&1; echo "rc=$?" >> run_a35b.log
until grep -q '^rc=' build_dd_val.log 2>/dev/null; do sleep 30; done
if grep -q "'tests_passed': True" build_dd_val.log; then
  python3 run.py --arms d_dedup d_dedup_off d_dedup_E --workers 2 > run_dd.log 2>&1; echo "rc=$?" >> run_dd.log
else
  echo "dd_val build or tests failed" > run_dd.log
fi
