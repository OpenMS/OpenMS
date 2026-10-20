#!/bin/bash
# Round 2 for one file: base2 (842adef9), pr2_off and pr2_on (f28b682d) searches, RSS profile of pr2_on, Percolator.
# usage: round2_file.sh <name> <mzML>
S=/tmp/claude-0/-home-user-OpenMS/1f1cec10-611d-54cf-aaaf-cdae17ef92a8/scratchpad
export PERCOLATOR=$S/tp/Linux/x86_64/Percolator/percolator
NAME=$1; M=$2; DB=/home/user/data/db/hye_entrap_td.fasta
cd /home/user/runs2
SKIP_PERC=1 $S/bench/run_search.sh /opt/openms-base2 base2 $M $DB $NAME/base2
SKIP_PERC=1 $S/bench/run_search.sh /opt/openms-pr2 pr2_off $M $DB $NAME/pr2_off
mkdir -p $NAME/pr2_on; ($S/bench/rss_profile.sh $NAME/pr2_on/prose.log $NAME/pr2_on.rss &)
SKIP_PERC=1 $S/bench/run_search.sh /opt/openms-pr2 pr2_on $M $DB $NAME/pr2_on -Search:chimeric:enabled
sha256sum $NAME/base2/search.pin.tsv $NAME/pr2_off/search.pin.tsv
for a in base2 pr2_off; do grep -v 'date\|creation' $NAME/$a/search.idXML | sha256sum; done
rm -f $NAME/base2/search.idXML $NAME/pr2_off/search.idXML $NAME/pr2_on/search.idXML
$S/bench/perc.sh $NAME/pr2_on
[ -s $NAME/../$NAME.off_done ] || $S/bench/perc.sh $NAME/pr2_off
rm -f $NAME/base2/search.pin.tsv
echo ROUND2_DONE $NAME
