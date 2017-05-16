#!/bin/bash

set -o pipefail

#   What are the dependencies for Read_Depths?
declare -a Read_Depths_Dependencies=(parallel)

#   Define a function to unzip the FastQC report files, extract
#       total number of sequences and sequence length, calculate
#       read depths, and get rid of unzipped FastQC report files
function getCoverage() {
    local zipFile="$1" # What is our zipped file?
    local -i target="$2" # What is our target size?
    local outDirectory="$3"
    local project="$4"
    local zipDirectory=`basename "${zipFile}" .zip` # Get the basename of the file/name of unzipped directory
    local sampleName=`basename "${zipFile}" _fastqc.zip` # Get the name of the sample
    local ROOT=`dirname "${zipFile}"` # Where is the zipfile stored?
    cd "${ROOT}" # Change to directory where zipfiles are contained
    unzip "${zipFile}" # Unzip the FastQC report directory
    cd "${zipDirectory}" # Change into the report directory
    local -i totalSequence=`grep 'Total Sequences' fastqc_data.txt | cut -f 2` # Extract total number of sequences
    local -i sequenceLength=`grep 'Sequence length' fastqc_data.txt | cut -f 2 | rev | cut -d '-' -f 1 | rev` # Extract sequence length
    echo -e "${sampleName}\t$(( totalSequence * sequenceLength / target ))" >> "${outDirectory}"/"${project}"_depths.txt # Do math and write to output file
    #   Get rid of unzipped FastQC report files
    cd "${ROOT}"
    rm -rf "${zipDirectory}"
}

#   Export the coverage function
export -f getCoverage

#   Define a function to run all of this
function Read_Depths() {
    local rawSamples="$1" # FastQ samples for getDepths
    local zipSamples="$2" # Zipped files for getCoverage
    local outDirectory="$3"/Read_Depths # The out directory to store our files
    local project="$4" # Name for the output files
    local -i target="$5" # Target size for getCoverage
    mkdir -p "${outDirectory}" # Make our output directory
    parallel getCoverage {} ${target} ${outDirectory} ${project} :::: "${zipSamples}" # Run getCoverage in parallel
}

#   Export the function
export -f Read_Depths
