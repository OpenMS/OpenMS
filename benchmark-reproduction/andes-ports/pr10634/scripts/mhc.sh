#!/bin/bash
# MHCquant spike-in (PXD011628, 100 fmol, msms52/53/54): PR build, option off vs on, interleaved per file
R=/home/user/bench
for n in 52 53 54; do for arm in off on; do
  X=(); [ $arm = on ] && X=(--set Search:report:collapse_positional_isomers=true)
  [ -f $R/runs/mhc/${arm}_$n/metrics.json ] && continue
  python3 -I $R/scripts/run_arm.py --bin /opt/openms-pr --params $R/cfg/mhc.json --mzml $R/mhc/mhc100fm_$n.mzML \
    --fasta $R/mhc_db/human_sp_td.fasta --out $R/runs/mhc/${arm}_$n "${X[@]}"
done; done
