#!/bin/env bash

#   This script performs quality trimming
#   on a series of FASTQ samples using sickle and seqqs.
#   Please install these before use.

set -e
set -o pipefail

#   What are the dependencies for Adapter_Trimming?
declare -a Quality_Trimming_Dependencies=(sickle seqqs Rscript)

#   A function to perform the trimming and plotting
#       Adapted from Tom Kono and Peter Morrell
function trimAutoplot() {
    #   Set the arguments for trimming
    local sampleName="$1" #   Name of the sample
    local forward="$2" #  Forward file
    local reverse="$3" #  Reverse file
    local out="$4"/"${sampleName}" #  Outdirectory
    local threshold="$5" #    Threshold Value
    local encoding="$6" # Platform for sequencing
    local seqHand="$7" #  The sequence_handling directory
    if [[ -d "${seqHand}"/HelperScripts ]] #   Check to see if helper scripts directory exists
    then
        local helper="${seqHand}"/HelperScripts #   The directory for the helper scripts
    else
        echo "Cannot find directory with helper scripts!"
        exit 1
    fi
    #   Make the out directories
    stats="${out}"/stats
    plots="${stats}"/plots
    mkdir -p "${plots}"
    #   Run seqqs on the raw samples
    seqqs -q "${encoding}" -p "${stats}"/raw_"${sampleName}"_R1 "${forward}"
    seqqs -q "${encoding}" -p "${stats}"/raw_"${sampleName}"_R2 "${reverse}"
    #   Trim the sequences based on quality
    sickle pe -t "${encoding}" -q "${threshold}" \
        -f "${forward}" \
        -r "${reverse}" \
        -o "${out}"/"${sampleName}"_R1_trimmed.fastq \
        -p "${out}"/"${sampleName}"_R2_trimmed.fastq \
        -s "${out}"/"${sampleName}"_singles_trimmed.fastq
    #   Run seqqs on the trimmed samples
    seqqs -q "${encoding}" -p "${stats}"/trimmed_"${sampleName}"_R1 "${out}"/"${sampleName}"_R1_trimmed.fastq
    seqqs -q "${encoding}" -p "${stats}"/trimmed_"${sampleName}"_R2 "${out}"/"${sampleName}"_R2_trimmed.fastq
    #   Gzip the trimmed files
    gzip "${out}"/"${sampleName}"_R1_trimmed.fastq
    gzip "${out}"/"${sampleName}"_R2_trimmed.fastq
    gzip "${out}"/"${sampleName}"_singles_trimmed.fastq
    #   Fix the quality scores
    "${helper}"/fix_quality.sh "${stats}"/raw_"${sampleName}"_R1_qual.txt
    "${helper}"/fix_quality.sh "${stats}"/raw_"${sampleName}"_R2_qual.txt
    "${helper}"/fix_quality.sh "${stats}"/trimmed_"${sampleName}"_R1_qual.txt
    "${helper}"/fix_quality.sh "${stats}"/trimmed_"${sampleName}"_R2_qual.txt
    #   Make the plots
    Rscript "${helper}"/plot_seqqs.R "${stats}" "${sampleName}"
}

#   Export the function
export -f trimAutoplot

#   A function to run the quality trimming
function Quality_Trimming() {
    local sampleName="$1" #   Name of the sample
    local forward="$2" #  Forward file
    local reverse="$3" #  Reverse file
    local out="$4"/"Quality_Trimming" #  Outdirectory
    local threshold="$5" #    Threshold Value
    local encoding="$6" # Platform for sequencing
    local seqHand="$7" #  The sequence_handling directory
    #   Test to see if there are equal numbers of forward and reverse reads
    local FORWARD_COUNT="`grep -cE "${forward}" ${sampleName}`"
    local REVERSE_COUNT="`grep -cE "${reverse}" ${sampleName}`"
    if [[ "${FORWARD_COUNT}" -eq "${REVERSE_COUNT}" ]]
    then
        echo "Equal numbers of forward and reverse samples."
    else
        exit 1
    fi
    #   Make the out directory
    mkdir -p ${out}
    #   Create arrays of forward and reverse samples
    local -a FORWARD=(`grep -E "${forward}" "${sampleName}"`)
    local -a REVERSE=(`grep -E "${reverse}" "${sampleName}"`)
    #   Create an array of sample names
    local -a SAMPLE_NAMES=()
    local counter=0
    for s in "${FORWARD[@]}"
    do
        SAMPLE_NAMES[`echo "$counter"`]="`basename $s ${forward}`"
        let "counter += 1"
    done
    parallel --xapply trimAutoplot {1} {2} {3} ${out} ${threshold} ${encoding} ${seqHand} ::: ${SAMPLE_NAMES[@]} ::: ${FORWARD[@]} ::: ${REVERSE[@]}
}

#   Export the function
export -f Quality_Trimming