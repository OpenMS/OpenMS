#!/bin/sh

#######################################################################
##	Rename sqMass files
## 
## The purpose of this script is to rename OpenSwathWorkflow output
## extracted ion chromatogram files to add a .chrom.sqMass 
## naming convetion, which is expected by DIAlignR

data_dir=${1:-results/xics/*.sqMass}
file_regex=${2:-'[[:alnum:]][.]chrom[.]sqMass'}
target_name=${3:-'${f%.sqMass}.chrom.sqMass'}

echo -e "-------------- Parameters --------------"
echo -e "data directory: $data_dir"
echo -e "file regex pattern: $file_regex"
echo -e "target rename: $target_name"
echo -e "----------------------------------------"

for f in $data_dir; do if ! grep -q $file_regex "$f" ; then eval new_name="$target_name"; mv -v -- "$f" $new_name; fi ; done