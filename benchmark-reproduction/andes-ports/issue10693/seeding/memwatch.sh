#!/bin/bash
# log peak RSS of ProteomicsLFQ processes every 20 s
while true; do ps -C ProteomicsLFQ -o rss=,etime=,args= 2>/dev/null | awk -v t="$(date -u +%FT%TZ)" '{print t, $1/1048576 " GB", $2, $NF}' ; sleep 20; done
