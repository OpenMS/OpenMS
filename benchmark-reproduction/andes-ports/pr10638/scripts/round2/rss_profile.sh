#!/bin/bash
# Samples the RSS of the running ProSE once per second with the last progress/log line: rss_profile.sh <prose.log> <out>
LOG=$1; OUT=$2; : > "$OUT"
until pgrep -x ProSE >/dev/null; do sleep 0.2; done
t=0
while pid=$(pgrep -x ProSE); do
  r=$(ps -o rss= -p $pid 2>/dev/null | tr -d ' ')
  l=$(grep -a -v "unknown parameter" "$LOG" | tr '\r' '\n' | grep -a -v '^\s*$' | tail -1 | cut -c1-90)
  echo "$t $r $l" >> "$OUT"; t=$((t+1)); sleep 1
done
