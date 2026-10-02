#!/bin/bash
# after chain_pd.sh: the calibration-set fix on all files and the shuffled entrapment databases
cd "$(dirname "$0")"
until grep -q CHAIN_DONE logs_chain_pd.txt; do sleep 30; done
until python3 -c "import json,sys; sys.exit(0 if json.load(open('build_records.json')).get('pdr_cal',{}).get('tests_passed') else 1)"; do sleep 30; done
set -x
python3 run.py --arms pd_cal --resume --workers 2
python3 run.py --arms pd_cal_E --resume --workers 2
echo CHAIN_CAL_DONE
