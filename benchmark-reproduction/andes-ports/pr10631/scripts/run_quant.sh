#!/bin/bash
# Isobaric quantification of one run with one arm; default parameters, same in both arms.
# usage: run_quant.sh <arm: base|pr> <mzML> <idXML> <type> <outdir>
arm=$1; mzml=$2; idxml=$3; type=$4; out=$5
source $(dirname "$0")/env.sh $arm
mkdir -p $out; stem=$(basename $mzml .mzML)
IsobaricAnalyzer -in $mzml -out $out/$stem.IA.consensusXML -type $type -threads 1 > $out/$stem.IA.log 2>&1
IsobaricWorkflow -in $mzml -in_id $idxml -out $out/$stem.IW.consensusXML -type $type -calculate_id_purity -threads 1 > $out/$stem.IW.log 2>&1
echo "$arm $stem IA=$(grep -c 'took' $out/$stem.IA.log) IW=$(grep -c 'took' $out/$stem.IW.log)"
