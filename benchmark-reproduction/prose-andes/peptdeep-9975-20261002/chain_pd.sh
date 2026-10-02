#!/bin/bash
# #9975 benchmark: identity first (cheap), then the main arm, entrapment and the instrument variant.
cd "$(dirname "$0")"
set -x
python3 run.py --arms pd_off --resume --workers 2
python3 run.py --arms pd_on --resume --workers 2
python3 run.py --arms re_auto_E pd_on_E --resume --workers 2
python3 run.py --arms pd_inst --resume --workers 2
echo CHAIN_DONE
