#!/bin/bash

set -e
set -u
set -o pipefail

#   This is a simple script to create a list
#   of samples in a directory and all of its
#   subdirectories.

usage() {
    echo -e "\
Usage: ./sample_list_generator.sh file_ext reads_dir out_dir out_name \n\
where:  file_ext is the extenstion of samples being found \n\
            example: \n\
                .txt.gz \n\
                .fq.gz \n\
                .sam \n\
\n\
        reads_dir is the directory in which all samples are found or \n\
            the root directory for subdirectories containing samples \n\
\n\
        out_dir is the desired out directory for the list. \n\
            NOTE: this should NOT be the same as reads_dir \n\
\n\
        out_name is the desired name for the sample list \n\
" >&2
    exit 1
}

if [ "$#" -lt 4 ]; then
    usage;
fi

#   File extension of samples
#       Example:
#           .fastq.gz
#           .txt.gz
#           .fq.gz
FILE_EXT="$1"

#   Full path to directory in
#   which ALL samples are stored
#   NOTE: subdirectories for indivdidual
#   samples may be within this directory.
#   They will still be searched through
READS_DIR="$2"

#   Desired path and name of outfile
#   Should not be same as $READS_DIR
OUT_DIR="$3"
OUT_NAME="$4"

find "$READS_DIR" -name "*$FILE_EXT" | sort > ${OUT_DIR}/${OUT_NAME}

if [[ -f "check_sample_list.sh" ]]
then
    bash check_sample_list.sh "${OUT_DIR}"/"${OUT_NAME}"
else
    echo "To ensure proper sample list generation, use 'check_sample_list.sh', bundled with 'sequence_handling'"
    echo "Find 'sequence_handling' at https://github.com/MorrellLab/sequence_handling"
fi
