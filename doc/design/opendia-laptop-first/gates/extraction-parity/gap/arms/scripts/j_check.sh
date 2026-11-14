#!/bin/bash
#SBATCH -t 02:00:00
#SBATCH --mem=48G
#SBATCH -c 12
# Post-hoc parent-artefact control: featJ (block I without same-composition claimants) as block I of the frozen trainer.
G=/nfs/research/juan/singj/opendia/gap-20261114
W=/nfs/research/juan/singj/opendia/pairlib-20261113
LOC=/home/singjc/Documents/github/OpenMS-r1-benchmark-20261107/pairlib-20261113
PY=/home/singjc/.jcode/scratch/torchenv/.venv/bin/python
S="singularity exec --cleanenv --bind $G/home:$LOC $W/container/opendia-pairlib-a8be5bed3c.sif"
J=$G/home/score/jcheck; mkdir -p $J
for r in raw01 raw02; do
  $S $PY $LOC/scripts_gap/features_J.py $r $LOC $LOC/lib12/library12_pop.pqp $LOC/iso_windows.tsv $LOC/$r/pipeline/search.features.tsv > $J/featJ_$r.log 2>&1 || exit 1
done
$S $PY $LOC/scripts_gap/j_train.py > $J/train.log 2>&1 || exit 1
$S $PY $LOC/scripts_gap/compare_gap.py D11FDWT D11FDIWT D11FDJWT > $J/COMPARE_J.txt 2> $J/compare.err || exit 1
cd $G/home/scripts 2>/dev/null
for a in D11FDJWT; do $S bash -c "cd $LOC/scripts && $PY diagnose.py $a $LOC xic" > $J/DIAG_$a.txt 2>&1; done
$S $PY $LOC/scripts_gap/pooled_rho.py $LOC D11FDJWT > $J/RHO_D11FDJWT.txt 2>&1
