#!/bin/bash
set -ex
cd /home/singjc/.jcode/scratch/compact-20261124
search=/home/singjc/.jcode/scratch/perf/lib/search.opendia.parquet
scoring=/home/singjc/.jcode/scratch/perf/lib/scoring.opendia.parquet
# Concurrent micro timings: identity raw01 with 16 threads may compete for CPU and memory bandwidth.
date -Is
ps -eo pid,etimes,pcpu,rss,args | grep -E 'OpenDIA|identity' || true
sha256sum compact_heap compact_probe /home/singjc/Documents/github/OpenMS-r1-full-ninja-build/lib/libOpenDIA_impl.a
/usr/bin/time -v ./compact_heap "$search" "$scoring" > compact-heap.log 2>&1
for trial in 1 2 3; do
  /usr/bin/time -v ./old/probe analysis "$search" "$scoring" > "old-probe-analysis-$trial.log" 2>&1
  /usr/bin/time -v ./compact_probe analysis "$search" "$scoring" > "compact-probe-analysis-$trial.log" 2>&1
done
date -Is
sha256sum old/* compact_heap compact_probe *.cpp *.sh *.log > final-sha256.txt
