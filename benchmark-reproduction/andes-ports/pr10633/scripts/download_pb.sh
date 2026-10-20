#!/bin/bash
# ProteoBench DDA quant module (QExactive HF-X) mzML, as used by #10585, and the ProteoBench HYE FASTA.
set -euo pipefail
cd ${1:-.}
for f in A_Sample_Alpha_01 A_Sample_Alpha_02 A_Sample_Alpha_03 B_Sample_Alpha_01 B_Sample_Alpha_02 B_Sample_Alpha_03; do
  echo "https://proteobench.cubimed.rub.de/raws/DDA/mzml/LFQ_Orbitrap_DDA_Condition_$f.mzML"
done > urls.txt
xargs -P 3 -n 1 curl -sS -O --retry 5 --retry-delay 5 < urls.txt
curl -sS -O https://proteobench.cubimed.rub.de/fasta/ProteoBenchFASTA_MixedSpecies_HYE.zip
