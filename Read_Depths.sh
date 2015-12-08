#!/bin/bash

#   This script calculates coverage from FastQC output
#   Please run the 'Assess_Quality.sh' script before running

#PBS -l mem=4000mb,nodes=1:ppn=4,walltime=6:00:00
#PBS -m abe
#PBS -M user@example.com
#PBS -q lab

set -e
set -o pipefail

# module load parallel

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
# SAMPLE_INFO=
#
# #   Full path to out directory
# SCRATCH=
#
# #   Project name
# PROJECT=
#
# #   Target size for samples
# TARGET=
#
# #   Make the out directory
# OUT=${SCRATCH}/${PROJECT}/Read_Depths
# mkdir -p ${OUT}

#   Define a function to unzip the FastQC report files, extract
#       total number of sequences and sequence length, calculate
#       read depth, and get rid of unzipped FastQC report files
function getCoverage() {
    local zipFile="$1" # What is our zipped file?
    local zipDirectory=`basename "${zipFile}" .zip` # Get the basename of the file/name of unzipped directory
    local sampleName=`basename "${zipFile}" _fastqc.zip` # Get the name of the sample
    local ROOT=`dirname "${zipFile}"` # Where is the zipfile stored?
    cd "${ROOT}" # Change to directory where zipfiles are contained
    unzip "${zipFile}" # Unzip the FastQC report directory
    cd "${zipDirectory}" # Change into the report directory
    local totalSequence=`grep 'Total Sequences' fastqc_data.txt | cut -f 2` # Extract total number of sequences
    local sequenceLength=`grep 'Sequence length' fastqc_data.txt | cut -f 2 | rev | cut -d '-' -f 1 | rev` # Extract sequence length
    echo -e "${sampleName} \t $[${totalSequence} * ${sequenceLength} / ${target}]" >> "${outDirectory}"/"${project}"_depths.txt # Do math and write to output file
    #   Get rid of unzipped FastQC report files
    cd "${ROOT}"
    rm -rf "${zipDirectory}"
}

#   Export the coverage function
export -f getCoverage

#   Define a function to get the read counts
function getCounts() {
    local sample="$1" # What is our sample?
    local name="$( basename ${sample} )" # Get the name of the sample
    local count="$( bioawk -cfastx 'END{print NR}' $sample )" # Find the read counts
    printf %s"$name \t $count \n" # Write the results all pretty like
}

#   Export counts function
export -f getCounts

#   Define a function to run all of this
function Read_Depths() {
    local rawSamples="$1" # FastQ samples for getCounts
    local zipSamples="$2" # Zipped files for getCoverage
    local outDirectory="$3"/Read_Depths # The out directory to store our files
    local project="$4" # Name for the output files
    local target="$5" # Target size for getCoverage
    mkdir -p "${outDirectory}" # Make our output directory
    parallel getCoverage {} :::: "${zipSamples}" # Run getCoverage in parallel
    parallel "getCounts {} :::: ${rawSamples}" > "${outDirectory}"/"${project}"_counts.txt # Run getCounts in parallel and write results to an output file
}

#   Export the function
export -f Read_Depths

#   Pass variables to and run thye function in parallel
# cat ${SAMPLE_INFO} | parallel "read_depths {} ${TARGET} ${OUT} ${PROJECT}"
#
# #   Where can the final output be found?
# echo "Final output file can be found at"
# echo "${OUT}/${PROJECT}_read_depths.txt"
