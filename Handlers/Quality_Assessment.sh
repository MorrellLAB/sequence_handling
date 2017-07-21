#!/bin/sh

#   This script runs FastQC on a series of samples

set -o pipefail

#   What are the dependencies for Quality_Assessment?
declare -a Quality_Assessment_Dependencies=(fastqc parallel)

#   A function to run the quality assessment
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
