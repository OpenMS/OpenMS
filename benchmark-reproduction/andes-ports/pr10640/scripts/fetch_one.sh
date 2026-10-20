#!/bin/bash
# download one raw file (resumable, 5 tries); wait while >= 4 downloaded raw files wait for conversion
name=$1; size=$2; RAW=/home/user/data/pxd022787/raw
[ -e /home/user/data/pxd022787/mzml/${name%.raw}.mzML ] && exit 0
[ -e $RAW/$name.done ] && exit 0
waiting() { local n=0; for d in $RAW/*.raw.done; do [ -e "${d%.done}" ] && n=$((n+1)); done; echo $n; }
while [ $(waiting) -ge 4 ]; do sleep 15; done
cd $RAW
for i in 1 2 3 4 5; do
  curl -sS -C - -o $name https://ftp.pride.ebi.ac.uk/pride/data/archive/2024/05/PXD022787/$name
  [ "$(stat -c %s $name)" -eq "$size" ] && { touch $name.done; exit 0; }
  sleep $((2**i))
done
echo "FAILED $name"
