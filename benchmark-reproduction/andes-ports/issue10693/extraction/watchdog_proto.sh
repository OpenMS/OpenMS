#!/bin/bash
# stop master_proto.sh and the combine arms if a prototype detection fails (they would re-detect with 4 threads)
L=/home/user/logs/master_proto.log; n=$(grep -c . $L)
while true; do
  if tail -n +$((n+1)) $L | grep -q "FAILED .*det_proto"; then
    for p in $(pgrep -f "bash ./master_proto.sh"); do kill $p; done
    sleep 2
    for p in $(pgrep -x ProteomicsLFQ); do kill $p; done
    echo "$(date -u +%FT%TZ) watchdog: stopped chain after failed detection" >> $L; exit 0
  fi
  tail -n +$((n+1)) $L | grep -q -E "done .*det_proto" && exit 0
  sleep 10
done
