#!/bin/bash
# kill any ProSE whose RSS exceeds LIMIT_KB (default 12.5 GiB) so the container never swaps/OOMs
LIMIT_KB=${LIMIT_KB:-13107200}
while true; do
  for p in $(pgrep -x "ProSE|andes"); do
    r=$(awk '/VmRSS/{print $2}' /proc/$p/status 2>/dev/null)
    if [ -n "$r" ] && [ "$r" -gt "$LIMIT_KB" ]; then kill -9 $p; echo "$(date +%T) killed $p rss=$r" >> /home/user/bench/watchdog.log; fi
  done; sleep 0.5
done
