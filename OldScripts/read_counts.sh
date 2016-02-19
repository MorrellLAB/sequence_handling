#!/bin/bash

set -e
set -u
set -o pipefail

#   This is a script to count the depths
#   of reads defined by a list of samples

#   This script uses bioawk to count the
#   read depth, please make sure bioawk
#   is installed before running this script

usage() {
    echo -e "\
Usage: ./read_counts.sh sample_info outdirectory \n\
where:  sample_info is a list of samples to be processed \n\
\n\
        outdirectory is the directory where the file should be placed \n\
" >&2
    exit 1
}


if [ "$#" -lt 2 ]; then
    usage;
fi

#   List of samples to be processed
sample_info=$1

#   Specify path to outdirectory
OUTDIR=$2

#   Check to see if bioawk is installed
if `command -v bioawk > /dev/null 2> /dev/null`
then
    echo "Bioawk is installed"
else
    echo "Please install Bioawk and add it to your PATH"
    echo
    echo "Running 'installer.sh bioawk' will do this"
    exit 1 
fi

#   Truncate sample info file into output file name
outfile=$(basename $sample_info .txt)

echo "$outfile"

echo "$sample_info"

#   Create a bash array of sample names for files
sample_names=($(cut -f 1 "$sample_info"))

#   Create output file
touch ${OUTDIR}/${outfile}_out.txt

echo ${sample_names[*]}

#   Iterate over each of the sample names and calculate
#   Read depth 
for sample in ${sample_names[*]}
    do
    		count="$(bioawk -cfastx 'END{print NR}' $sample)"
    		printf %s"$sample \t $count \n" >> ./${outfile}_out.txt
    done

echo Results can be found at "${OUTDIR}"/"$outfile"_out.txt
