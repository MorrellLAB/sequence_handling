#!/bin/bash

#   This script calculates coverage from FastQC output
#   Please run the 'Assess_Quality.sh' script before running

#PBS -l mem=4000mb,nodes=1:ppn=4,walltime=6:00:00
#PBS -m abe
#PBS -M user@example.com
#PBS -q lab

set -e
set -u
set -o pipefail

module load parallel

#   This script is a QSub submission for calculating the read depths of a batch of samples
#   To use, on line 5, change the 'user@example.com' to your own email address
#       to get notifications on start and completion for this script
#   Add the full file path to list of samples on the 'SAMPLE_INFO' field on line 38
#       This should look like:
#           SAMPLE_INFO=${HOME}/Directory/list.txt
#       Use ${HOME}, as it is a link that the shell understands as your home directory
#           and the rest is the full path to the actual list of samples
#   Put the full directory path for the output in the 'SCRATCH' field on line 41
#       This should look like:
#           SCRATCH="${HOME}/Out_Directory
#       Adjust for your own out directory.
#   Name the project in the 'PROJECT' field on line 44
#       This should look lke:
#           PROJECT=Barley
#       The full output path will be ${SCRATCH}/${PROJECT}/Read_Depths
#   Define the target size for your samples in the 'TARGET' field on line 47
#       This should look like:
#           TARGET=60000000
#   Run the script using the qsub command
#       qsub
#   This script outputs a text file with all of the read depths for each sample within it.

#   List of samples to be processed
SAMPLE_INFO=

#   Full path to out directory
SCRATCH=

#   Project name
PROJECT=

#   Target size for samples
TARGET=

#   Make the out directory
OUT=${SCRATCH}/${PROJECT}/Read_Depths
mkdir -p ${OUT}

#   Define a function to unzip the FastQC report files, extract
#       total number of sequences and sequence length, calculate
#       read depth, and get rid of unzipped FastQC report files
function read_depths() {
    #   Define variables in relation to function
    #   In order: zipped FastQC report, target size for read depths,
    #       out directory, and project name
    ZIPFILE="$1"
    TARG="$2"
    OUTDIR="$3"
    PROJ="$4"
    #   Find name of uzipped FastQC report directory, sample,
    #       and directory where zipfiles are contained
    ZIP_DIR=`basename "${ZIPFILE}" .zip`
    SAMPLE_NAME=`basename "${ZIPFILE}" _fastqc.zip`
    ROOT=`dirname "${ZIPFILE}"`
    #   Change to directory where zipfiles are contained
    cd "${ROOT}"
    #   Unzip and change into FastQC report directory
    unzip "${ZIPFILE}"
    cd "${ZIP_DIR}"
    #   Extract total number of sequences and sequence length
    TOTAL_SEQUENCE=`grep 'Total Sequences' fastqc_data.txt | cut -f 2`
    SEQUENCE_LENGTH=`grep 'Sequence length' fastqc_data.txt | cut -f 2 | rev | cut -d '-' -f 1 | rev`
    #   Do math and write to output file
    echo -e "${SAMPLE_NAME} \t $[${TOTAL_SEQUENCE} * ${SEQUENCE_LENGTH} / ${TARG}]" >> "${OUTDIR}"/"${PROJ}"_read_depths.txt
    #   Get rid of unzipped FastQC report files
    cd "${ROOT}"
    rm -rf "${ZIP_DIR}"
}

#   Export the function to be used by GNU parallel
export -f read_depths

#   Pass variables to and run thye function in parallel
cat ${SAMPLE_INFO} | parallel "read_depths {} ${TARGET} ${OUT} ${PROJECT}"

#   Where can the final output be found?
echo "Final output file can be found at"
echo "${OUT}/${PROJECT}_read_depths.txt"
