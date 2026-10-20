#!/bin/bash
# HYE target-decoy database as in #10585
set -e; source $(dirname $0)/env.sh
$BIN/DecoyDatabase -in $DATA/ProteoBenchFASTA_MixedSpecies_HYE.fasta -out $WORK/db/HYE_td.fasta \
  -decoy_string DECOY_ -decoy_string_position prefix -method reverse -enzyme Trypsin -no_progress
