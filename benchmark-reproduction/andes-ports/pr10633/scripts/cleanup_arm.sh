#!/bin/bash
# After evaluation: drop the normalized consensusXML and the mzTab, gzip the consensusXML (kept for re-analysis).
d=/home/user/bench/lfq/$1/$2
rm -f $d/rtlocal.consensusXML $d/lfq.mzTab
[ -f $d/lfq.consensusXML ] && gzip -f $d/lfq.consensusXML
