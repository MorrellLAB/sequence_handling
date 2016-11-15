#!/bin/bash

#   This script demultiplexes FastQ files based on barcode sequences
#   using FastX_Barcode_Splitter. Inspired from Jeff Neyhart's GBarleyS
#   pipeline, available on GitHub: https://github.com/neyhartj/GBarleyS

set -o pipefail

#   What are the dependencies for GBS_Demultiplexer?
declare -a GBS_Demultiplex_Dependencies=(Rscript fastx_barcode_splitter.pl parallel)

#   A function to determine if our barcodes are at the beginning or end of the sequence
function checkLineEnding() {
    local lineEnding="$1" # Get the line ending
    if [[ "${lineEnding}" == 'beginning' ]] # If the barcodes are at the beginning?
    then
        local endArg='--bol' # Set the argument as such
    elif [[ "${lineEnding}" == 'end' ]] # If the barcodes are at the end?
    then
        local endArg='--eol' # Set the argument as such
    else # Otherwise...
        echo "Incorrect ending specified!" >&2
        exit 86 # Exit with error
    fi
    echo "${endArg}" # Return the argument
}

#   Export the function
export -f checkLineEnding

#   A function to check if a file is empty
function testExistence() {
    local sample="$1"
    if [[ ! -s "${sample}" ]]
    then
        rm "${sample}"
    fi
}

#   Export the function
export -f testExistence

#   A function to handle the demultiplexing
function demultiplexFastQ() {
    local sample="$1" # What sample are we working with?
    local extension="$2" # What's the extension of our sample?
    local barcodeFile="$3" # Where is our barcode file?
    local outDirectory="$4" # Where are we storing our output files?
    local endArg="$5" # What is our argument for line endings?
    local mismatches="$6" # How many mismatches are allowed?
    local partial="$7" # Allow partial overlap?
    local sampleName=$(basename "${sample}" "${extension}")
    local baseBarcode=$(basename "${barcodeFile}" .barcode) # Get the name of the barcode file
    local barcodeLength="${baseBarcode: -1}" # The last letter of the barcode file is the length of the barcodes
    local prefix="${outDirectory}"/"${sampleName}"_"${barcodeLength}"_ # Make a prefix for fastx_barcode_spliter.pl
    #   Check the compression level of our sample
    if [[ $( echo "${sample}" | rev | cut -f 1 -d '.' | rev) == 'gz' ]] # If gzipped...
    then
        local toDemultiplex="${outDirectory}/${sampleName}_${barcodeLength}_sample_PIPE" # Create a name for the pipe
        rm -f "${toDemultiplex}" # Remove any existing pipes
        mkfifo "${toDemultiplex}" # Make the pipe
        gzip -cd "${sample}" > "${toDemultiplex}" # Uncompress the file and write to pipe
    elif [[ $( echo "${sample}" | rev | cut -f 1 -d '.' | rev) == 'bz2' ]] # If bzipped...
    then
        local toDemultiplex="${outDirectory}/${sampleName}_${barcodeLength}_sample_PIPE" # Create a name for the pipe
        rm -f "${toDemultiplex}" # Remove any existing pipes
        mkfifo "${toDemultiplex}" # Make the pipe
        bzip2 -cd "${sample}" > "${toDemultiplex}" & # Uncompress the file and write to pipe
    else # Otherwise
        local toDemultiplex="${sampleName}" # Use the name of the sample as 'toDemultiplex'
    fi
    #   Send the multiplexed sample to fastx_barcode_spliter.pl
    cat "${toDemultiplex}" | fastx_barcode_splitter.pl \
        --bcfile "${barcodeFile}" \
        --prefix "${prefix}" \
        --suffix .fastq \
        "${endArg}" \
        --mismatches "${mismatches}" \
        --partial "${partial}" \
        --quiet
    #   Find all demultiplexed samples output from fastx_barcode_spliter.pl
    local -a demultiplexed=($(find "${outDirectory}" -name "${sampleName}_${barcodeLength}_*.fastq"))
    #   Remove the empty files
    parallel testExistence {} ::: "${demultiplexed[@]}"
    #   Find all samples that didn't get removed
    local -a stillThere=($(find "${outDirectory}" -name "${sampleName}_${barcodeLength}_*.fastq"))
    #   Gzip our files
    parallel gzip {} ::: "${stillThere[@]}"
}

#   Export the function
export -f demultiplexFastQ

#   A function to split the key file based on barcode length
function barcodeGenerator() {
    local helperScripts="$1"/HelperScripts
    local keyfile="$2"
    local outDirectory="$3"/GBS_Demultiplex/barcodes
    local project="$4"
    mkdir -p "${outDirectory}" # Make our output directory
    if [[ ! -d "${helperScripts}" ]]; then echo "Cannot find directory with helper scripts!" >&2; exit 1; fi
    #   Create barcode files from the keyfile based on barcode length
    Rscript "${helperScripts}"/collect_barcodes.R "${keyfile}" "${outDirectory}" "${project}"
    #   Create a list with all the barcode files
    find "${outDirectory}" -name "*.barcode" > "${outDirectory}"/"${project}"_barcode_list.txt
}

#   Export the function
export -f barcodeGenerator

#   A function to run the demultiplexing
function GBS_Demultiplex() {
    local barcode="$1" # What is our barcode file?
    local sampleList="$2" # What is our list of samples?
    local outDirectory="$3"/GBS_Demultiplex # Where are we storing the demultiplexed files?
    local lineEnding="$4" # Which end of the sequence are the barcodes at?
    local mismatches="${5:-1}" # How many mismatches are allowed?
    local partial="${6:-0}" # Allow partial overlap?
    local extension="$7" # What's the extension for the samples?
    local project="$8" # What's the name of our project?
    mkdir -p "${outDirectory}" # Make our output directory
    #   Create an array of sample names
    #local -a sampleNames=($(parallel basename {} "${extension}" :::: "${sampleList}"))
    #   Get the argument for the line endings
    local endArg=$(checkLineEnding "${lineEnding}")
    #   Do the demultiplexing
    parallel demultiplexFastQ {} "${extension}" "${barcode}" "${outDirectory}" "${endArg}" "${mismatches}" "${partial}" :::: "${sampleList}"
}

#   Export the function
export -f GBS_Demultiplex
