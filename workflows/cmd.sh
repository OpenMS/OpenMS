#!/bin/usr/sh

# activate python that has the snakemake module
# source /project/bin/pyenv38/bin/activate 

# Load singularity for docker containers
## You may need to load singularity if you're using a computing cluster
# module load singularity/3.8

# Configuration files to run snakemake on SLURM clusters
snakemake_config=config/

## To Run Snakemake WFs, uncomment the ones you want to run below.

# Run msconvert WF
# NOTE: env variable TMPMYWINEPREFIX and singularity arg --writable-tmpfs is needed to run the wine version of MSConvert, which allows for vendor format converisons such as Thermo.
#
#TMPMYWINEPREFIX=`mkdir -p /dev/shm/mywineprefix`
# snakemake --profile $snakemake_config --snakefile Snakefile.msconvert --use-singularity -j 16 --singularity-args "-B /localscratch:/localscratch --writable-tmpfs" 

# Run MOBI-DIK diapysef conversion WF
# snakemake --profile $snakemake_config --snakefile Snakefile.diapasef_convert --use-singularity -j 16 --singularity-args "-B /localscratch:/localscratch"

# Create small datasets with OpenMS FileFilter
# snakemake --profile $snakemake_config --snakefile Snakefile.filefilter --use-singularity -j 16 --singularity-args "-B /localscratch:/localscratch" 

# Run DIA-Umpire WF
# snakemake --profile $snakemake_config --snakefile Snakefile.diau_workflow --use-singularity -j 16 --singularity-args "-B /localscratch:/localscratch"

# Run DDA - Library Generation WF
# snakemake --profile $snakemake_config --snakefile Snakefile.dda_workflow --use-singularity -j 16 --singularity-args "-B /localscratch:/localscratch"

# Run OpenSWATH WF
# snakemake --profile $snakemake_config --snakefile Snakefile.dia_workflow --use-singularity -j 16 --singularity-args "-B /localscratch:/localscratch" 

