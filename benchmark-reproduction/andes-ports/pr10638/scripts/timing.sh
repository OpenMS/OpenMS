#!/bin/bash
# Interleaved timing: 3 rounds of (pr_off, pr_on) per file, search only, nothing else running; outputs deleted after.
S=/tmp/claude-0/-home-user-OpenMS/1f1cec10-611d-54cf-aaaf-cdae17ef92a8/scratchpad
DB=/home/user/data/db/hye_entrap_td.fasta
cd /home/user/runs
for rep in 1 2 3; do
  for f in hfx:/home/user/data/LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01.mzML astral:/home/user/data/astral_A1_rt623-1048.mzML; do
    name=${f%%:*}; M=${f#*:}
    for arm in pr_off pr_on; do
      x=""; [ $arm = pr_on ] && x=-Search:chimeric:enabled
      O=timing/$name/$arm/t$rep
      SKIP_PERC=1 $S/bench/run_search.sh /opt/openms-pr $arm $M $DB $O $x
      echo "$name $arm $rep $(grep -E 'Elapsed' $O/time.txt | awk '{print $NF}') $(grep 'Maximum resident' $O/time.txt | awk '{print $NF}') $(grep -o 'search [0-9.]*s' $O/prose.log)" >> timing/timing.txt
      rm -f $O/search.idXML $O/search.pin.tsv
    done
  done
done
echo TIMING_DONE >> timing/timing.txt
