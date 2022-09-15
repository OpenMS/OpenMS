# Snakemake Mass Spectrometry Data Analysis Workflows

This repository is meant to house various snakemake workflows for performing mass spectrometry data analyses.

[Snakemake](https://snakemake.readthedocs.io/en/stable/) is a workflow manager implemented and distributed via python, that allows the ease of scaling up workflows from single-core machines to large compute clusters. Snakemake follows the philosophy of human reable workflows that are easy to implement and allow for interoperatability when re-analyzing data, or moving to a different computing system. [1] Furthermore, workflows are set up as a list of rules to run, which can generate directed acyclic graphs to identify which rules have dependencies of prior rules exectuing successfully before self execution. This makes snakemake ideal for running pipelines for analyzing mass spectrometry data, because it allows for mutiple analyses to be run in parallel for many MS experiment files across several computing nodes.

# Using a cluster

If you plan on submitting jobs using a computing cluster, take note of the config folder and it's contents. The current config is an example profile for running workflows on a cluster that usings a SLURM scheduler. The **config/config.yaml** file contains some some instructions of what files to look for to use for job submissions, and how many times to restart a job if it fails (`restart-times`), as well as other additional arguments that control how long to wait for an expected output file to be generated before after a job has finished. The **config/.envs/res.json** file contains the resource request allocation for each job submission. You can tailor the requested amount of resources for specific rules that get submitted for jobs, as well as, have a default amount of resources requested if none is defined for a specfific cluster.

**NOTE**: Snakemake workflow rules allow for singularity containers that can be downloaded and build from images stored on dockerhub. However, if you are working on a cluster, depending on network connectivity and the size of the image, it may take a while to download and stored these as singularity images. In these cases, you can build a singularity image locally, and reference that singularity image in the snakemake workflow rule instead.

# Data Conversion

## ProteoWizard MSConvert

`Snakefile.msconvert` is a snakemake workflow for converting between various mass spec data formats using [ProteoWizards MSConvert](https://proteowizard.sourceforge.io/download.html) tool contained in the `chambm/pwiz-skyline-i-agree-to-the-vendor-licenses` dokcer image.

### Example

<details>

Make sure you have a python enviornment activated that has the snakemake module installed. 

**Note:** you can specify the number of files to be processed in parallel using the `-j` flag. In the example below, we supply 16, meaning 16 raw MS data files will be converted in parallel.

**NOTE:** We need to initialize `TMPMYWINEPREFIX` and add the `--writable-tmpfs` to the singularity arguments, in order to properly use the wine version of MSConvert, which allows for the proprietary vendor file conversions.

```
TMPMYWINEPREFIX=`mkdir -p /dev/shm/mywineprefix`
snakemake --snakefile Snakefile.msconvert --use-singularity -j 16 --singularity-args "-B /localscratch:/localscratch --writable-tmpfs"
```

</details>

## Mobi-DIK diapysef

If you have timsTOF ion mobility diaPASEF data, you also have the option of using the conversion tool in [diapysef](https://github.com/Roestlab/dia-pasef). This can be performed using the `Snakefile.diapasef_convert` snakemake workflow.

### Example

<details>

```
snakemake --snakefile Snakefile.diapasef_convert --use-singularity -j 16 --singularity-args "-B /localscratch:/localscratch --writable-tmpfs"
```

</details>

# Pseudo DDA Spectra Generation from DIA Data using DIA-Umpire

To generate pseudo spectra from DIA data, you can use DIA-Umpire's signal extraction module. This may be useful if you don't have any DDA runs to generate a spectral library, or if you want to generate run specific irt libraries that can be used for linear and non-linear alignment in OpenSwathWorkflow. You can pass a parameter file to DIA-Umpire's signal extraction module, an example for Thermo and SCIEX data can be found in the params folder. We can use the `Snakefile.diau_workflow` snakemake file to genera pseudo spectra for DIA data.

# Spectral Library Generation from DDA Data

To analyse DDA data and for spectral library generation for corresponding DIA data, we can use the `Snakefile.dda_workflow`.

**NOTE:** This workflow uses [MSFragger](https://github.com/Nesvilab/MSFragger/wiki/Preparing-MSFragger#Downloading-MSFragger) for database peptide searching. MSFragger is freely available for academic research, non-commercial and education purposes, which means docker images for this tool needs to be private. You will need to build your own version of an MSFragger docker image.

# Spectral Library Generation from DIA Data

To analyse DIA data, we can use the `Snakefile.dia_workflow`, that will perform targeted OpenSwathWorkflow data extraction, perform stastical validation using pyprophet, and perform RT alignment for recapturing missed signals and improved quantification using DIAlignR.

# Contributors

* Justin Sing
* George Rosenberger

# References

1, Mölder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., Forster, J., Lee, S., Twardziok, S.O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., Köster, J., 2021. Sustainable data analysis with Snakemake. F1000Res 10, 33.