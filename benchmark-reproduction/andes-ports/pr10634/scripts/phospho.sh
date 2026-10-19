#!/bin/bash
# PXD007653 phospho file (andes's phospho benchmark), mouse + 1:1 shuffled entrapment + reversed decoys
set -e
R=/home/user/bench
for arm in "$@"; do
  case $arm in
    off)  B=/opt/openms-pr;   X=() ;;
    on)   B=/opt/openms-pr;   X=(--set Search:report:collapse_positional_isomers=true) ;;
    base) B=/opt/openms-base; X=() ;;
  esac
  [ -f $R/runs/phospho/$arm/metrics.json ] && continue
  python3 -I $R/scripts/run_arm.py --bin $B --params $R/cfg/phospho.json --mzml $R/phospho/control2_phospho.mzML \
    --fasta $R/phospho_db/mouse_entrap_td.fasta --out $R/runs/phospho/$arm "${X[@]}"
done
