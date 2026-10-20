#!/bin/bash
B=$(dirname $(readlink -f $0)); export WORK2=${WORK2:-/home/user/work2} DATA2=${DATA2:-/home/user/data2}; L=$WORK2/lfq; mkdir -p $L
for r in 202505616 202505617 202505618 202505622 202505623 202505624; do
  $B/search_one_agilent.sh $r > $WORK2/$r.search.log 2>&1; echo "$(date +%T) search $r exit=$?"; done
run() { name=$1; shift; echo "$(date +%T) start $name"; BIN=/opt/openms-pr/bin $B/plfq_agilent.sh $L/$name "$@" > $L/$name.log 2>&1; echo "$(date +%T) $name exit=$?"; }
run T_off $L/feat_T -targeted_only true
run T_on  $L/feat_T -targeted_only true -consensus_elution_width true
run P_off $L/feat_P -targeted_only false -pip_echo true
run P_on  $L/feat_P -targeted_only false -pip_echo true -consensus_elution_width true
echo ALL_DONE
