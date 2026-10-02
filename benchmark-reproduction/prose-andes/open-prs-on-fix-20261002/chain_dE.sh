#!/bin/sh
cd "$(dirname "$0")"
until grep -q '^rc=\|build failed' run_d35.log 2>/dev/null; do sleep 30; done
python3 run.py --arms d_priors_E d_79_E d_base_E --workers 2 > run_dE.log 2>&1; echo "rc=$?" >> run_dE.log
