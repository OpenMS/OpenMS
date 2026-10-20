#!/bin/bash
# Usage: search_all.sh <bin_dir> <td_fasta> <out_dir>   (runs search_one.sh for the 6 runs, sequentially)
B=$1; DB=$2; OUT=$3; S=$(dirname "$0"); export LD_LIBRARY_PATH=$(dirname $B)/lib
TP=/home/user/ext/THIRDPARTY/Linux/x86_64
for c in A B; do for i in 01 02 03; do
  r=LFQ_Orbitrap_DDA_Condition_${c}_Sample_Alpha_${i}
  [ -s $OUT/$r.idXML ] && continue
  bash -c 'time "$@"' _ $S/search_one.sh $B /home/user/data/mzml_pb/$r.mzML $DB $OUT $TP/Comet/comet.exe $TP/Percolator/percolator > $OUT/$r.log 2>&1 || echo "FAILED $r"
  echo "done $r $(date)"
done; done
echo SEARCH_ALL_DONE
