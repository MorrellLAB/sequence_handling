#!/bin/bash

set -e
set -u
set -o pipefail

#   This is a simple script to create a list
#   of samples in a directory.

usage() {
    echo -e "\
Usage: ./sample_list_generator.sh file_extension directory out_name \n\
where:  file_extension is the suffix of the samples you want included in the list \n\
            examples: \n\
                .fastq.gz \n\
                .sam \n\
\n\
        directory is the directory in which all samples are found \n\
\n\
        out_name is the desired name for the sample list \n\
" >&2
    exit 1
}

if [ "$#" -lt 3 ]; then
    usage;
fi

#   File extension of samples
FILE_EXT="$1"

#   Full path to directory in which ALL samples are stored
INPUT_DIR="$2"
cd "${INPUT_DIR}"
READS_DIR=$(pwd -P)

#   Desired path and name of outfile
OUT_NAME="$3"

find "$READS_DIR" -maxdepth 1 -name "*$FILE_EXT" | sort > ${READS_DIR}/${OUT_NAME}

NUMBER=$(cat ${READS_DIR}/${OUT_NAME} | wc -l)

echo "A list of ${NUMBER} ${FILE_EXT} samples was generated at ${READS_DIR}/${OUT_NAME}"
echo "Please double check that all of your samples were included properly"