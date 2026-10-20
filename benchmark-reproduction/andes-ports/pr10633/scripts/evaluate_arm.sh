#!/bin/bash
# Evaluate one arm: transfer counts / entrapment FDP from the consensusXML, ProteoBench metrics (min 3 runs) on the
# unnormalized (ProteomicsLFQ -out_msstats) and rt_local-normalized MSstats tables.
# usage: evaluate_arm.sh <db> <arm>   -> lfq/<db>/<arm>/eval_*.tsv
set -euo pipefail
B=/home/user/bench; d=$B/lfq/$1/$2
fasta=()
[ "$1" = entrap ] && fasta=(--fasta $B/db/HYE_entrap_td.fasta)
python3 -I $B/scripts/transfer_stats.py <(zcat -f $d/lfq.consensusXML*) "${fasta[@]}" --label "$1/$2" --header > $d/eval_transfers.tsv
for t in msstats rtlocal.msstats; do
  $B/pbenv/bin/python -I $B/scripts/score_proteobench.py $d/$t.csv --label "$1/$2/$t" --header 2>/dev/null > $d/eval_pb_$t.tsv
done
cat $d/eval_*.tsv
