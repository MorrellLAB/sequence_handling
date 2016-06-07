#!/bin/sh

#   This script runs FastQC on a series of samples
#   Please have FastQC installed before using this script

set -e
set -o pipefail

#   What are the dependencies Assess Quality?
declare -a Quality_Assessment_Dependencies=(fastqc parallel)

#   A function to run the quality assesment
function Quality_Assessment() {
    local sampleList="$1" # What is our list of samples?
    local outDir="$2" # Where are we storing our results?
    local proj="$3" # What do we call our results?
    local out="${outDir}"/Quality_Assessment # Full path to output directory
    mkdir -p "${out}" # Make our output directory
    cat "${sampleList}" | parallel "fastqc --outdir ${out} {}" # Run FastQC in parallel
    find "${out}" -name "*.zip" | sort > "${out}"/"${proj}"_FastQC_ZipFiles.txt # Write a list of zipped files
}

export -f Quality_Assessment
