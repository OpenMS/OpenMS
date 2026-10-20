#!/bin/bash
# ProteoBench's own mzML (as used by #10585), not the archive's re-conversion
cd "$1"
for f in A_Sample_Alpha_01 A_Sample_Alpha_02 A_Sample_Alpha_03 B_Sample_Alpha_01 B_Sample_Alpha_02 B_Sample_Alpha_03; do
  echo "https://proteobench.cubimed.rub.de/raws/DDA/mzml/LFQ_Orbitrap_DDA_Condition_$f.mzML"
done | xargs -P 3 -n 1 curl -s -O --retry 5 --retry-delay 5
echo DL_DONE
