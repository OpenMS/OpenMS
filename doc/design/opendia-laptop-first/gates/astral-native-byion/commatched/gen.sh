#!/bin/bash
# k11 candidate decoys: shuffles with switchKR (precursor mass shift kept, y-series moved), 4 seeds, from the same targets.
cd "$(dirname "$0")"; B=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
T=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/astral-design-raw01-20261107/rethink/decoyregen/targets.pqp
set -o pipefail
for s in 1 2 3 4; do
  $B/OpenSwathDecoyGenerator -in $T -out shuf$s.tsv -method shuffle -shuffle_seed $s -decoy_tag DECOY_ -separate -threads 8 > gen$s.log 2>&1 || { echo 1 > gen.EXIT; exit 1; }
done
echo 0 > gen.EXIT
