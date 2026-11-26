#!/bin/bash
# Download two never-used PXD027359 HeLa diaPASEF raws (new held-out set for PRECHECK_EVSHARE).
set -euo pipefail
cd /nfs/research/juan/singj/opendia/evshare-20261126/data
for f in $(awk '{print $2}' pride.sha1); do
  test -e "$f" || wget -q "https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/07/PXD027359/$f"
done
sha1sum -c pride.sha1
for f in *.d.zip; do test -e "${f%.zip}" || unzip -q "$f"; done
du -sh *.d
