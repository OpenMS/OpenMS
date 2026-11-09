#!/bin/bash
# Decoy-family regeneration (design diagnostic, k3). Same targets, three decoy families, merged for one extraction.
cd "$(dirname "$0")"; B=/home/singjc/Documents/github/OpenMS-r1-full-ninja-build/bin
set -o pipefail
$B/OpenSwathDecoyGenerator -in targets.pqp -out pr_noKR.tsv -method pseudo-reverse -switchKR false -decoy_tag DECOYP_ -separate -threads 8 > gen_pr.log 2>&1 && \
$B/OpenSwathDecoyGenerator -in targets.pqp -out shuf.tsv -method shuffle -switchKR false -decoy_tag DECOYS_ -separate -threads 8 > gen_sh.log 2>&1 && \
$B/OpenSwathDecoyGenerator -in targets.pqp -out pr_orig.tsv -method pseudo-reverse -decoy_tag DECOY_ -separate -threads 8 > gen_or.log 2>&1
echo $? > gen.EXIT
