#!/bin/bash
# searches each run once its download is complete (file ends with </indexedmzML>)
source $(dirname $0)/env.sh
for r in $RUNS; do
  until tail -c 64 $DATA/$r.mzML 2>/dev/null | grep -q "</indexedmzML>"; do sleep 20; done
  s=$(date +%s); $(dirname $0)/search_one.sh $r > $WORK/search/$r.log 2>&1; echo "$r exit=$? $(( $(date +%s)-s )) s"
done
echo SEARCH_DONE
