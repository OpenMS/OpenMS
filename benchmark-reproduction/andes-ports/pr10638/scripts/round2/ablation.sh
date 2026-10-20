#!/bin/bash
# HF-X ablations of the andes-matching changes (local build with PROSE_CHIMERIC_ABLATION switches).
S=/tmp/claude-0/-home-user-OpenMS/1f1cec10-611d-54cf-aaaf-cdae17ef92a8/scratchpad
export PERCOLATOR=$S/tp/Linux/x86_64/Percolator/percolator
NAME=${NAME:-hfx}; M=${M:-/home/user/data/LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01.mzML}; DB=/home/user/data/db/hye_entrap_td.fasta
cd /home/user/runs2
for arm in "$@"; do
  O=$NAME/abl_$arm
  PROSE_CHIMERIC_ABLATION=$arm SKIP_PERC=1 $S/bench/run_search.sh /opt/openms-abl abl $M $DB $O -Search:chimeric:enabled
  rm -f $O/search.idXML
  $S/bench/perc.sh $O
  python3 -I $S/bench/evaluate.py $O | grep MEAN >> $NAME/ablation.jsonl
  echo "$arm $(grep -E 'Elapsed' $O/time.txt | awk '{print $NF}') $(grep 'Maximum resident' $O/time.txt | awk '{print $NF}') $(grep -h 'co-isolated precursors searched' $O/prose.log)" >> $NAME/ablation_runs.txt
  rm -f $O/input.pin
done
echo ABL_DONE >> $NAME/ablation_runs.txt
