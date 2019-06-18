#!/bin/bash

#   This script performs adapter trimming on a series of FASTQ samples using scythe.
#   Please install scythe before use.

set -o pipefail

#   What are the dependencies for Adapter_Trimming?
declare -a Adapter_Trimming_Dependencies=(scythe parallel)

function Adapter_Trimming_One() {
    local sample="$1"
    local out="$2" # Where are we storing our results?    
    local project="$3" # What do we call our results?
    local forwardNaming="$4" # What is the extension indicating a forward read?
    local reverseNaming="$5" # What is the extension indicating a reverse read?
    local adapters="$6" # What is our adapter file?
    local prior="$7" # What is Scythe's prior?
    local platform="$8" # What platform did we sequence on?
    #   Check to see if we have a forward or reverse sample
    if [[ ! -z "${forwardNaming}" && $(echo "${sample}" | grep "${forwardNaming}") ]]
    then # If forward
        local name="$(basename ${sample} ${forwardNaming})"_Forward # Make the name say forward
    elif [[ ! -z "${reverseNaming}" && $(echo "${sample}"| grep "${reverseNaming}") ]]
    then # If reverse
        local name="$(basename ${sample} ${reverseNaming})"_Reverse # Make the name say reverse
    else # If this is neither
        local name="$(basename ${sample} | cut -f 1 -d '.')"_Single # Make the name say single
    fi
    #   Trim the sample
    (set -x; scythe -a "${adapters}" -p "${prior}" -q "${platform}" "${sample}" > "${out}/${name}_ScytheTrimmed.fastq"
    #   Gzip the output fastq file
    gzip "${out}/${name}_ScytheTrimmed.fastq" )
}

#   Export the function
export -f Adapter_Trimming_One

#   A function to run the adapter trimming, adapted from Tom Kono
function Adapter_Trimming() {
    local sample_list="$1" # What is our list of samples?
    local out="$2/Adapter_Trimming" # Where are we storing our results?
    local project="$3" # What do we call our results?
    local forwardNaming="$4" # What is the extension indicating a forward read?
    local reverseNaming="$5" # What is the extension indicating a reverse read?
    local adapters="$6" # What is our adapter file?
    local prior="$7" # What is Scythe's prior?
    local platform="$8" # What platform did we sequence on?
    #   Make our out directory
    mkdir -p "${out}"
    if [[ "${USE_PBS}" == "true" ]]
    then
        #   Make an array of samples from the sample list
        declare -a sample_array=($(grep -E ".fastq|.fastq.gz" "${sample_list}"))
        #   Get which sample in the list we are working on	
        local sample="${sample_array[${PBS_ARRAYID}]}"
        Adapter_Trimming_One ${sample} ${out} ${project} ${forwardNaming} ${reverseNaming} ${adapters} ${prior} ${platform}
    else # Not with PBS, use parallel
        grep -E ".fastq|.fastq.gz" "${sample_list}" | parallel "Adapter_Trimming_One {} ${out} ${project} ${forwardNaming} ${reverseNaming} ${adapters} ${prior} ${platform}"
        sort "${out}/${project}_trimmed_adapters.txt" > "${out}/${project}_trimmed_adapters.sorted.txt"
        rm "${out}/${project}_trimmed_adapters.txt"
        mv "${out}/${project}_trimmed_adapters.sorted.txt" "${out}/${project}_trimmed_adapters.txt"
    fi
}

#   Export the function
export -f Adapter_Trimming
