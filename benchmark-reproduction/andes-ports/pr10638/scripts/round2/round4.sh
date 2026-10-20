#!/bin/bash
# chimeric:top_hits arms (/opt/openms-pr4): whole Astral run with 3, HF-X with 1, 2, 3.
S=/tmp/claude-0/-home-user-OpenMS/1f1cec10-611d-54cf-aaaf-cdae17ef92a8/scratchpad
export PERCOLATOR=$S/tp/Linux/x86_64/Percolator/percolator
DB=/home/user/data/db/hye_entrap_td.fasta
A=/home/user/data/LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.mzML
H=/home/user/data/LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01.mzML
cd /home/user/runs2
run() { # mzML outdir top_hits
  local MZ=$1 O=$2 N=$3
  mkdir -p $O; ($S/bench/rss_profile.sh $O/prose.log $O.rss &)
  SKIP_PERC=1 $S/bench/run_search.sh /opt/openms-pr4 pr4 $MZ $DB $O -Search:chimeric:enabled -Search:chimeric:top_hits $N
  rm -f $O/search.idXML
  [ -s $O/search.pin.tsv ] && $S/bench/perc.sh $O && python3 -I $S/bench/evaluate.py $O | grep MEAN >> round4.jsonl
  echo "$O $(grep -E 'Elapsed' $O/time.txt | awk '{print $NF}') $(grep 'Maximum resident' $O/time.txt | awk '{print $NF}') $(grep -c 'signal 9' $O/time.txt)" >> round4_runs.txt
  rm -f $O/input.pin
}
run $A astral/pr4_t3 3
run $H hfx/pr4_t1 1
run $H hfx/pr4_t2 2
run $H hfx/pr4_t3 3
echo ROUND4_DONE >> round4_runs.txt
