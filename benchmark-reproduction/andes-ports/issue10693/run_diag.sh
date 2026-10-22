#!/bin/bash
export PIPECHO_DIAG=/home/user/runs/diag/link.tsv OPENMS_DATA_PATH=/home/user/OpenMS/share/OpenMS
D=/home/user/data; R=/home/user/runs
/home/user/bench/scripts/plfq.sh /opt/openms-diag $D/db/HYE_td.fasta $D/ids_hye $R/hye/feat_pe $R/diag/hye_off_q05 -targeted_only false -pip_echo true -PipEcho:shape_features false && echo DIAG_DONE || echo DIAG_FAILED
