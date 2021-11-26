#!/bin/bash

# find the location of FLASHDeconv file
export FLASHDeconvFile=$(find . -path "*/bin/FLASHDeconv")

# usage function
usage() {
	echo "  Usage:
	./$(basename "$0") <options> 
  
  Options (mandatory options marked with '*'):		
	-in <file>* 		Input file(*.mzML) or directory
	-out <file>*		Output file(*.tsv) or directory
	-out_spec <0/1>		true:1, false:0 - to write spectrum level deconvoluted masses per ms level (default: 1)
	-out_mzml <0/1>		true:1, false:0 - to write spectral deconvolution in mzml output format (default: 0)
	-out_topFD <0/1>	true:1, false:0 - to write output file in topFD format (default: 0)
	-out_promex <0/1>	true:1, false:0 - to write output file in promex format (default: 0)
	-max_MS_level		maximum MS level (inclusive) for deconvolution (default: 2)
	-h 			Shows script options
	--help			Shows FLASHDeconv options
	--helphelp		Shows all FLASHDeconv options (including advanced)

  Except for these options, other options corresponds to FLASHDeconv.
  "

	exit 0
}

# process output file options
run_FLASH_with_options() {

	# add output options based on max_MS_level
	[ -z "$max_MS_level" ] && max_MS_level=2

	# if out_mzml is true
	[ ! -z $out_mzml ] && [ $out_mzml = "1" ] && out_mzml_file=" -out_mzml \"${out_prefix}_deconv.mzML\""

	# if out_promex is true
	[ ! -z $out_promex ] && [ $out_promex = "1" ] && out_promex_file=" -out_promex \"${out_prefix}.ms1ft\""

	# for spectral deconvolution outputs
	if [ -z $out_spec ] || [ $out_spec = "1" ]; then
		out_spec_files=" -out_spec"
		for num_ms in $(seq 1 $max_MS_level); do
			out_spec_files="${out_spec_files} \"${out_prefix}_ms${num_ms}.tsv\""
		done
	fi

	# if out_topFD is true
	if [ ! -z $out_topFD ] && [ $out_topFD = "1" ]; then
		out_topFD_files=" -out_topFD"
		for num_ms in $(seq 1 $max_MS_level); do
			out_topFD_files="${out_topFD_files} \"${out_prefix}_ms${num_ms}.msalign\""
		done
	fi

	# hidden argument - train
	[ ! -z $hidden_train ] && [ $hidden_train = "1" ] && train_files=" -in_train \"${in_file%.*}_ms2_toppic_prsm.tsv\" -out_train \"${out_prefix}_att.csv\""

	# hidden argument - read_flashida_log
	if [ ! -z $hidden_flashida ] && [ $hidden_flashida = "1" ]; then
		in_file_name=$(basename "$in_file")
		in_file_name=${in_file_name%.*}
		in_flashida_file=" -in_log \"$(dirname "$in_file")/IDALog_${in_file_name}.log\""
	fi

	file_options="-in \"${in_file}\" -out \"${out_prefix}.tsv\"${out_spec_files}${out_mzml_file}${out_topFD_files}${out_promex_file} -max_MS_level ${max_MS_level}${train_files}${in_flashida_file}"

	# running FLASHDeconv
	echo " " # for blank line, indicating new file will be processed.
	echo "$FLASHDeconvFile ${file_options}$1"
	eval "$FLASHDeconvFile ${file_options}$1"
}

# # parsing arguments
while [ "${1:-}" != "" ]; do
	case "$1" in
		"-h")
			usage
			;;
		"--help")
			$FLASHDeconvFile --help
			exit 1
			;;
		"--helphelp")
			$FLASHDeconvFile --helphelp
			exit 1
			;;
		"-in")
			shift # move pointer to argument
			in="$1"
			;;
		"-out")
			shift
			out="$1"
			;;
		"-out_spec")
			shift
			out_spec=$1
			;;
		"-out_mzml")
			shift
			out_mzml=$1
			;;
		"-out_topFD")
			shift
			out_topFD=$1
			;;
		"-out_promex")
			shift
			out_promex=$1
			;;	
		"-max_MS_level")
			shift
			max_MS_level=$1
			;;
		"-train")
			shift
			hidden_train=$1
			;;
		"-read_flashida_log")
			shift
			hidden_flashida=$1
			;;
		*)
			other_args="$other_args $1"
			;;
	esac
	shift # move pointer to next option
done


# if $in or $out were given
if [ -z "$in" ] || [ -z "$out" ]; then
	usage
fi

# if in or out does not exist
if [ ! -f "$in" ] && [ ! -d "$in" ]; then
	echo "ERROR : $in does not exist."
	exit 0
fi

# check if out is file format or directory
is_out_file=""  # initialize
if [ -f "$out" ]; then
	[ "${out##*.}" = "tsv" ] && is_out_file=true || usage
elif [ -d "$out" ]; then
	is_out_file=false;
else
	# out doesn't exist. create one.
	if [[ $out != *.* ]]; then
		# if out doesn't contain period (directory)
		is_out_file=false;
	else
		# if out is file
		out_format=$(echo "${out##*.}" | tr '[:upper:]' '[:lower:]')
		[ $out_format = "tsv" ] && is_out_file=true || usage
	fi
fi


# if in is directory but out is file
if [ -d "$in" ] && [ $is_out_file = true ]; then
	echo "ERROR: if [-in] is directory, [-out] should be directory."
	exit 0

# if in is file
elif [ -f "$in" ]; then

	# check format of in file
	in_format=$(echo "${in##*.}" | tr '[:upper:]' '[:lower:]')
	[ $in_format != "mzml" ] && usage
	in_file=$in # variable name used in function 'get_output_file_options'

	# 1. if out is file
	if [ $is_out_file = true ]; then
		# get output prefix
		out_prefix="${out%.*}"
	# 2. if out is directory
	else 
		# use in as file name, at out directory
		out_prefix=$(basename $in)
		out_prefix="${out%%/}/${out_prefix%.*}"
	fi
	
	run_FLASH_with_options "$other_args"

# if both in and out are directories 
elif [ -d "$in" ] && [ $is_out_file = false ]; then
	[ -d "$out" ] || mkdir "$out" # if out directory does not exist, create one. 

	# for all mzml files in "in directory"
	for in_file in "${in%%/}"/*.mzML; do
		out_prefix=$(basename "$in_file")
		out_prefix="${out%%/}/${out_prefix%.*}"

		run_FLASH_with_options "$other_args"
	done

# cannot reach here. something's wrong.
else
	usage
fi