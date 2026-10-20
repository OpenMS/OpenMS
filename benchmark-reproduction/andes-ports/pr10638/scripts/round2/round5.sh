#!/bin/bash
# Final build (f2800c61, /opt/openms-pr4): FAIMS, Astral slice (also andes-faithful pr2 there), HF-X in-process
# Percolator, then interleaved timing (3 rounds of off/on on HF-X and the whole Astral run).
S=/tmp/claude-0/-home-user-OpenMS/1f1cec10-611d-54cf-aaaf-cdae17ef92a8/scratchpad
export PERCOLATOR=$S/tp/Linux/x86_64/Percolator/percolator
DB=/home/user/data/db/hye_entrap_td.fasta
A=/home/user/data/LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.mzML
SL=/home/user/data/astral_A1_rt623-1048.mzML
F=/home/user/data/ExplorisFAIMS-DDA-LFQ.mzML
H=/home/user/data/LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01.mzML
cd /home/user/runs2
run() { # prefix mzML outdir args
  local P=$1 MZ=$2 O=$3; shift 3
  SKIP_PERC=1 $S/bench/run_search.sh $P x $MZ $DB $O "$@"
  rm -f $O/search.idXML
  $S/bench/perc.sh $O && python3 -I $S/bench/evaluate.py $O >> round5.jsonl
  rm -f $O/input.pin
}
run /opt/openms-pr4 $F faims/pr4_on -Search:chimeric:enabled
run /opt/openms-pr4 $SL slice/pr4_on -Search:chimeric:enabled
run /opt/openms-pr2 $SL slice/pr2_on -Search:chimeric:enabled
SKIP_PERC=1 $S/bench/run_inproc.sh /opt/openms-pr4 x $H $DB hfx_inproc/pr4_on -Search:chimeric:enabled
python3 -I $S/bench/count_idxml.py hfx_inproc/pr4_on/search.idXML > hfx_inproc/pr4_on.json; rm -f hfx_inproc/pr4_on/search.idXML
echo SEARCHES_DONE
for rep in 1 2 3; do
  for f in hfx:$H astral:$A; do
    name=${f%%:*}; M=${f#*:}
    for arm in off on; do
      x=""; [ $arm = on ] && x=-Search:chimeric:enabled
      O=timing/$name/$arm/t$rep
      SKIP_PERC=1 $S/bench/run_search.sh /opt/openms-pr4 $arm $M $DB $O $x
      echo -e "$name\t$arm\t$rep\t$(grep -E 'Elapsed' $O/time.txt | awk '{print $NF}')\t$(grep 'Maximum resident' $O/time.txt | awk '{print $NF}')" >> timing/timing.tsv
      rm -f $O/search.idXML $O/search.pin.tsv
    done
  done
done
echo ROUND5_DONE
