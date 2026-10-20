#!/bin/bash
# usage: id_all.sh <db label> <fasta> [threads]
for r in $(source /home/user/bench/scripts/env.sh; echo $RUNS); do /home/user/bench/scripts/id_one.sh $1 $2 $r ${3:-4} || echo "FAILED $r"; done
echo "ID_DONE $1"
