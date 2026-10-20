#!/bin/bash
# Validation round: whole Astral run and FAIMS, base2 / pr2_off / pr2_on (andes-faithful, f28b682d) / pr3_on (0e0a9640).
S=/tmp/claude-0/-home-user-OpenMS/1f1cec10-611d-54cf-aaaf-cdae17ef92a8/scratchpad
export PERCOLATOR=$S/tp/Linux/x86_64/Percolator/percolator
DB=/home/user/data/db/hye_entrap_td.fasta
A=/home/user/data/LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.mzML
F=/home/user/data/ExplorisFAIMS-DDA-LFQ.mzML
H=/home/user/data/LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01.mzML
cd /home/user/runs2
search() { # prefix arm mzML outdir [args]
  local P=$1 ARM=$2 MZ=$3 O=$4; shift 4
  mkdir -p $O; ($S/bench/rss_profile.sh $O/prose.log $O.rss &)
  SKIP_PERC=1 $S/bench/run_search.sh $P $ARM $MZ $DB $O "$@"
  rm -f $O/search.idXML
}
search /opt/openms-pr3 pr3_on $A astral/pr3_on -Search:chimeric:enabled
search /opt/openms-pr2 pr2_off $A astral/pr2_off
search /opt/openms-base2 base2 $A astral/base2
search /opt/openms-pr3 pr3_on $H hfx/pr3_on -Search:chimeric:enabled
search /opt/openms-base2 base2 $F faims/base2
search /opt/openms-pr2 pr2_off $F faims/pr2_off
search /opt/openms-pr2 pr2_on $F faims/pr2_on -Search:chimeric:enabled
search /opt/openms-pr3 pr3_on $F faims/pr3_on -Search:chimeric:enabled
sha256sum astral/base2/search.pin.tsv astral/pr2_off/search.pin.tsv faims/base2/search.pin.tsv faims/pr2_off/search.pin.tsv hfx/pr3_on/search.pin.tsv hfx/abl_topn,iso,allpeaks,tol10/search.pin.tsv > identity.sha256
rm -f astral/base2/search.pin.tsv faims/base2/search.pin.tsv
echo SEARCHES_DONE
for O in astral/pr3_on astral/pr2_off faims/pr2_off faims/pr2_on faims/pr3_on hfx/pr3_on; do
  $S/bench/perc.sh $O &
  while [ $(pgrep -c -f "Percolator/percolator") -ge 2 ]; do sleep 10; done
done
wait
echo ROUND3_DONE
