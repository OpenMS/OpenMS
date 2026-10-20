#!/bin/bash
# All arms, one at a time (PIP-ECHO detection peaks at 11-13 GB). Option on/off share IDs and -feat_dir checkpoints.
B=$(dirname $(readlink -f $0)); L=${WORK:-/home/user/work}/lfq; mkdir -p $L
run() { name=$1; bin=$2; shift 2; echo "$(date +%T) start $name"; BIN=$bin $B/plfq.sh $L/$name "$@" > $L/$name.log 2>&1; echo "$(date +%T) $name exit=$?"; }
PR=/opt/openms-pr/bin; BASE=/opt/openms-base/bin
run T_off  $PR   $L/feat_T_pr  -targeted_only true
run T_on   $PR   $L/feat_T_pr  -targeted_only true  -consensus_elution_width true
run P_off  $PR   $L/feat_P_pr  -targeted_only false -pip_echo true
run P_on   $PR   $L/feat_P_pr  -targeted_only false -pip_echo true -consensus_elution_width true
run T_base $BASE $L/feat_T_base -targeted_only true
run P_base $BASE $L/feat_P_base -targeted_only false -pip_echo true
echo ALL_DONE
