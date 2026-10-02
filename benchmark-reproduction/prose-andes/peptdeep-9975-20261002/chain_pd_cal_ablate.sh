#!/bin/bash
# after chain_pd_cal.sh: Percolator-only feature ablation of the calibration-fix arms
cd "$(dirname "$0")"
until grep -q CHAIN_CAL_DONE logs_chain_pd_cal.txt; do sleep 30; done
set -x
python3 ablate_pd.py --arms pd_cal pd_cal_E --workers 3
echo CHAIN_ABLATE_DONE
