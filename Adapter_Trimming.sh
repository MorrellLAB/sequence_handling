#!/bin/bash

#   This script performs adapter trimming
#   on a series of FASTQ samples using Scythe
#   Please install these before use.

set -e
set -o pipefail

#   What are the dependencies for Adapter_Trimming?
declare -a Adapter_Trimming_Dependencies=(scythe Rscript)

#   A function to perform the trimming
#       Adapted from Tom Kono
function trimAdapters() {
    #   Set the arguments for the trimming
    local sample="$1" # What sample are we working with?
    local outPrefix="$2" #  Outdirectory
    local adapters="$3" # Adapter file
    local prior="$4" # Prior
    local platform="$5" # Quality encoding platform
    local helperScripts="$6" # Directory with helper scripts
    local forwardNaming="$7" # What is the forward naming scheme?
    local reverseNaming="$8" # What is the reverse naming scheme?
    #   Check to see if we have a forward or reverse sample
    if [[ $(echo "${sample}" | grep "${forwardNaming}") ]] # Is this a forward sample?
    then # If yes
        local name="$(basename ${sample} ${forwardNaming})"_Forward # Make the name say forward
        local out="${outPrefix}/$(basename ${sample} ${forwardNaming})" # Make a name for the outdirectory
    elif [[ $(echo "${sample}"| grep "${reverseNaming}") ]] # Is this the reverse sample?
    then # If yes
        local name="$(basename ${sample} ${reverseNaming})"_Reverse # Make the name say reverse
        local out="${outPrefix}/$(basename ${sample} ${reverseNaming})" # Make a name for the outdirectory
    else # If this is neither
        local name="$(basename ${sample} | cut -f 1 -d '.')"_Single # Make the name say single
        local out="${outPrefix}/$(basename ${sample} | cut -f 1 -d '.')" # Make a name for the outdirectory
    fi
    #   Make the outdirectory
    mkdir -p "${out}"
    #   Is our sample compressed?
    if [[ "$(echo ${sample} | rev | cut -f 1 -d '.' | rev)" == "gz" ]] # Is this gzipped?
    then # If so
        local toTrim="${out}/${name}_PIPE" # Make a name for the pipe which will be passed to the trimmer
        rm -f "${toTrim}" # Remove any pipes with the same name
        mkfifo "${toTrim}" # Make the pipe
        gzip -cd "${sample}" > "${toTrim}" & # Uncompress our sample to the pipe
    elif [[ "$(echo ${sample} | rev | cut -f 1 -d '.' | rev)" == "bz2" ]] # Is this bzipped?
    then # If so
        local toTrim="${out}/${name}_PIPE" # Make a name for the pipe which will be passed to the trimmer
        rm -f "${toTrim}" # Remove any pipes with the same name
        mkfifo "${toTrim}" # Make the pipe
        bzip2 -cd "${sample}" > "${toTrim}" & # Uncompress our sample to the pipe
    else # Otherwise
        local toTrim="${sample}" # Our name will be the sample itself
    fi
    #   Make a named for the trimmed sample
    local trimmed="${out}/${name}_ScytheTrimmed.fastq.gz"
    #   Trim the sample
    scythe -a "${adapters}" -p "${prior}" --quiet "${toTrim}" | gzip -c > "${trimmed}"
}

#   Export the function
export -f trimAdapters

#   A function to run the quality trimming
function Adapter_Trimming() {
    local rawSamples="$1" # What is our list of samples?
    local outDirectory="$2"/Adapter_Trimming # Where are we storing our results?
    local project="$3" # What do we call our results?
    local forwardNaming="$4" # What is the extension indicating a forward read?
    local reverseNaming="$5" # What is the extension indicating a reverse read?
    local adapters="$6" # What is our adapter file?
    local prior="$7" # What is Scythe's prior?
    local platform="$8" # What platform did we sequence on?
    local helperScripts="$9"/Helper_Scripts # Where are the helper scripts stored?
    # checkCounts "${rawSamples}" "${forwardNaming}" "${reverseNaming}" # Check the counts of forward and reverse reads
    if [[ "$?" -ne 0 ]]; then echo "Unbalanced forward and reverse reads" >&2; exit 1; fi # If not an equal amount, exit out with error
    mkdir -p outDirectory # Make our out directory
    parallel trimAutoplot {} "${outDirectory}" "${adapters}" "${prior}" "${platform}" "${helperScripts}" "${forwardNaming}" "${reverseNaming}" :::: "${rawSamples}" # Perform the trim
    find "${outDirectory}" -type p -exec rm {} \; # Clean up all pipes
    find "${outDirectory}" -name "*_ScytheTrimmed.fastq.gz" | sort > "${outDirectory}"/"${project}"_trimmed_adapters.txt # Create our list of trimmmed files
}