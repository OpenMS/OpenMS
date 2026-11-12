#!/bin/bash
# Identical DIA-NN 2.2.0 raw01 search, timed. GNU time runs inside the container on the host binary's libc-free path is not
# guaranteed, so time is measured two ways: /usr/bin/time -v around the docker client (wall clock only is meaningful) and
# the container cgroup's memory.peak / cpu.stat, read by a sidecar loop while the container runs.
cd /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/diann-ourlib-20261112/run2
trap 'echo $? > EXIT' EXIT
echo running > STATUS
date -u +%FT%TZ > start.txt
CID=$(docker run -d --cpus 14 -u $(id -u):$(id -g) \
  -v /home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/dd-plots-20261111:/in:ro -v $PWD/..:/speclib:ro -v $PWD:/w \
  --entrypoint /usr/diann-2.2.0/diann-linux diann:2.2.0 \
  --f /in/raw01.mzML --lib /speclib/lib_targets.tsv --threads 14 --verbose 1 \
  --out /w/report_raw01.parquet --qvalue 1 --gen-decoys \
  --unimod4 --min-pr-charge 1 --max-pr-charge 5 --temp /w/tmp)
echo $CID > cid.txt
CG=$(find /sys/fs/cgroup -maxdepth 4 -type d -name "*${CID}*" | head -1); echo "$CG" > cgroup_path.txt
peak=0
while [ "$(docker inspect -f '{{.State.Running}}' $CID)" = true ]; do
  if [ -n "$CG" ] && [ -r $CG/memory.current ]; then c=$(cat $CG/memory.current); [ $c -gt $peak ] && peak=$c; fi
  [ -n "$CG" ] && [ -r $CG/memory.peak ] && cat $CG/memory.peak > memory_peak_bytes.txt
  [ -n "$CG" ] && [ -r $CG/cpu.stat ] && cat $CG/cpu.stat > cpu.stat.txt
  sleep 1
done
echo $peak > memory_current_max_bytes.txt
docker inspect -f '{{.State.StartedAt}} {{.State.FinishedAt}} {{.State.ExitCode}}' $CID > container_times.txt
docker logs $CID > diann_raw01.log 2>&1
docker rm $CID > /dev/null
date -u +%FT%TZ > end.txt
echo done > STATUS
