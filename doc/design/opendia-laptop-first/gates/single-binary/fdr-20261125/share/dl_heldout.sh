#!/bin/bash
# Download PXD017703 HeLa_Evosep_diaPASEF_RAW.zip (never-used 60SPD and 100SPD HeLa diaPASEF runs) for the share-rule precheck.
set -euo pipefail
cd /nfs/research/juan/singj/opendia/share-20261126/data
test -e HeLa_Evosep_diaPASEF_RAW.zip || wget -q https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/12/PXD017703/HeLa_Evosep_diaPASEF_RAW.zip
sha1sum -c pride.sha1
unzip -q -o HeLa_Evosep_diaPASEF_RAW.zip '60SPD_py8/*' '100SPD_py8/*'
du -sh 60SPD_py8/*.d 100SPD_py8/*.d
