#!/bin/bash
# Calibration search (dev): keep scored rows and the columns cal_robust.py needs, gzip. Usage: reduce.sh <features.tsv> <out.tsv.gz>
set -euo pipefail
awk -F'\t' -v OFS='\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; n=split("precursor_id is_decoy score d_lib d_max i_nclaim",want," ");
  for(j=1;j<=n;j++){if(!(want[j] in c)){print "missing "want[j] > "/dev/stderr"; exit 3}; k[j]=c[want[j]]}; s=c["score"]}
  {if(NR>1 && $s+0 <= -1e29) next; line=$k[1]; for(j=2;j<=n;j++) line=line OFS $k[j]; print line}' "$1" | gzip -1 > "$2"
