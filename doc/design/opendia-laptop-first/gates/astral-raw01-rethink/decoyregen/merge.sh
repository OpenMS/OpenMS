#!/bin/bash
# Merge targets + three decoy families (DECOY_ = original pseudo-reverse with switchKR, control; DECOYP_ = pseudo-reverse
# without switchKR; DECOYS_ = shuffle without switchKR) into one PQP for a single extraction.
cd "$(dirname "$0")"; B=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
set -o pipefail
$B/TargetedFileConverter -in targets.pqp -out targets.tsv -out_type tsv > merge_t.log 2>&1 && \
( cat targets.tsv; tail -n +2 pr_orig.tsv; tail -n +2 pr_noKR.tsv; tail -n +2 shuf.tsv ) > merged.tsv && \
$B/TargetedFileConverter -in merged.tsv -in_type tsv -out merged.pqp -out_type pqp > merge_c.log 2>&1
echo $? > merge.EXIT
