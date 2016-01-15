#!/bin/bash

#   This script calculates coverage from FastQC output
#   Please run the 'Assess_Quality.sh' script before running

set -e
set -o pipefail

#   Define a function to unzip the FastQC report files, extract
#       total number of sequences and sequence length, calculate
#       read depth, and get rid of unzipped FastQC report files
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
    local -i target="$5" # Target size for getCoverage
    mkdir -p "${outDirectory}" # Make our output directory
    parallel getCoverage {} ${target} ${outDirectory} ${project} :::: "${zipSamples}" # Run getCoverage in parallel
    #parallel "getCounts {} :::: ${rawSamples}" #> "${outDirectory}"/"${project}"_counts.txt # Run getCounts in parallel and write results to an output file
}

#   Export the function
export -f Read_Depths
